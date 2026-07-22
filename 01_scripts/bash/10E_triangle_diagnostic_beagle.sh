#!/bin/bash
#SBATCH --job-name=tri_diag
#SBATCH --partition=medium
#SBATCH --time=12:00:00
#SBATCH --mem=60G
#SBATCH --cpus-per-task=8
#SBATCH --output=99_logs/10E_%j.out
#SBATCH --error=99_logs/10E_%j.err

# =============================================================================
# 10E_triangle_diagnostic_beagle.sh -- Identify diagnostic SNPs (delta-AF > 0.5 between LLS and JAY)
#       AND extract Beagle GLs at those sites for ALL individuals
#
# Inputs : 26_triangle/maf_anchors/{LLS,JAY}.mafs.gz  (per-anchor MAFs)
# Outputs: 26_triangle/diagnostic_snps.tsv            (chr, pos, major, minor,
#                                                       freq_LLS, freq_JAY, dAF)
#          26_triangle/diagnostic_sites.{sites,chrs}  (ANGSD sites format)
#          26_triangle/triangle.beagle.gz             (GLs all ind, diag SNPs)
# =============================================================================

set -euo pipefail
source config/00_config.sh

OUTDIR="26_triangle"
THRESHOLD=0.5

LLS_MAF="$OUTDIR/maf_anchors/LLS.mafs.gz"
JAY_MAF="$OUTDIR/maf_anchors/JAY.mafs.gz"

[[ -s "$LLS_MAF" ]] || { echo "Missing $LLS_MAF" >&2; exit 1; }
[[ -s "$JAY_MAF" ]] || { echo "Missing $JAY_MAF" >&2; exit 1; }

echo "=== Identifying diagnostic SNPs (|MAF_LLS - MAF_JAY| > $THRESHOLD) ==="

# --- Join LLS and JAY MAFs on chr+pos+major+minor ----------------------------
# .mafs.gz columns: chromo, position, major, minor, ref/anc?, knownEM, nInd
# We use 'knownEM' as the major-allele frequency since -doMajorMinor 3 was used
# in 10D_triangle_anchor_maf.sh (consistent major/minor across pops).
python3 - <<'PYEOF'
import gzip, os, sys

def read_maf(path):
    d = {}
    with gzip.open(path, 'rt') as f:
        header = f.readline().rstrip().split('\t')
        # find columns
        try:
            i_chr   = header.index('chromo')
            i_pos   = header.index('position')
            i_maj   = header.index('major')
            i_min   = header.index('minor')
            i_freq  = header.index('knownEM')
        except ValueError:
            i_chr, i_pos, i_maj, i_min = 0, 1, 2, 3
            # fallback: pick first frequency-like column
            for j, h in enumerate(header):
                if h.lower() in ('knownem', 'phat', 'freq'):
                    i_freq = j
                    break
        for line in f:
            t = line.rstrip().split('\t')
            key = (t[i_chr], t[i_pos], t[i_maj], t[i_min])
            d[key] = float(t[i_freq])
    return d

print("Reading LLS MAF...", flush=True)
lls = read_maf("26_triangle/maf_anchors/LLS.mafs.gz")
print(f"  LLS sites: {len(lls):,}", flush=True)

print("Reading JAY MAF...", flush=True)
jay = read_maf("26_triangle/maf_anchors/JAY.mafs.gz")
print(f"  JAY sites: {len(jay):,}", flush=True)

shared = set(lls) & set(jay)
print(f"  Shared sites: {len(shared):,}", flush=True)

THRESHOLD = 0.5
n_diag = 0
with open("26_triangle/diagnostic_snps.tsv", "w") as out:
    out.write("chromo\tposition\tmajor\tminor\tfreq_LLS\tfreq_JAY\tdeltaAF\n")
    for k in shared:
        f_lls = lls[k]; f_jay = jay[k]
        d = abs(f_lls - f_jay)
        if d > THRESHOLD:
            out.write(f"{k[0]}\t{k[1]}\t{k[2]}\t{k[3]}\t{f_lls:.4f}\t{f_jay:.4f}\t{d:.4f}\n")
            n_diag += 1

print(f"  Diagnostic SNPs (dAF > {THRESHOLD}): {n_diag:,}", flush=True)
PYEOF

# --- Build ANGSD sites file & chr list for diagnostic SNPs -------------------
tail -n +2 "$OUTDIR/diagnostic_snps.tsv" | \
    awk 'BEGIN{OFS="\t"} {print $1, $2, $3, $4}' \
    > "$OUTDIR/diagnostic_sites.sites"
cut -f1 "$OUTDIR/diagnostic_sites.sites" | sort -u > "$OUTDIR/diagnostic_sites.chrs"

${ANGSD:-angsd} sites index "$OUTDIR/diagnostic_sites.sites"

N_DIAG=$(wc -l < "$OUTDIR/diagnostic_sites.sites")
echo ""
echo "Diagnostic sites file: $OUTDIR/diagnostic_sites.sites"
echo "  N SNPs: $N_DIAG"
echo "  N chromosomes: $(wc -l < "$OUTDIR/diagnostic_sites.chrs")"

# --- Extract Beagle GLs for ALL individuals at diagnostic SNPs ---------------
echo ""
echo "=== Extracting Beagle GLs at diagnostic SNPs (Kivalliq + LLS + JAY) ==="

BAMLIST="02_info/bamlists/triangle.bamlist"
N_IND=$(wc -l < "$BAMLIST")
MIN_IND=$(awk -v n="$N_IND" -v p="$PERCENT_IND" 'BEGIN{printf "%d", n*p+0.999}')

${ANGSD:-angsd} \
    -bam "$BAMLIST" \
    -ref "$GENOME" \
    -out "$OUTDIR/triangle" \
    -sites "$OUTDIR/diagnostic_sites.sites" \
    -rf    "$OUTDIR/diagnostic_sites.chrs" \
    -GL 1 \
    -doMajorMinor 3 \
    -doMaf 1 \
    -doGlf 2 \
    -doCounts 1 \
    -minMapQ 30 -minQ 20 \
    -minInd "$MIN_IND" \
    -nThreads "$SLURM_CPUS_PER_TASK"

echo ""
echo "=== Outputs ==="
echo "  Diagnostic SNP table: $OUTDIR/diagnostic_snps.tsv"
echo "  Beagle file:          $OUTDIR/triangle.beagle.gz"
echo "  N SNPs in Beagle:     $(zcat "$OUTDIR/triangle.beagle.gz" | tail -n +2 | wc -l)"
