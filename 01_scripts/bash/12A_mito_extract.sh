#!/bin/bash
# =============================================================================
# 12A_mito_extract.sh — Mitochondrial consensus sequence per individual
#
#   Mitogenome in the reference: NC_000861.1 (S. alpinus, 16,659 bp)
#   Coverage is very high (~750x for a 2x nuclear sample), so a simple
#   majority consensus per fish is reliable.
#   Insertions are not written and deletions are written as "-", so every
#   sequence is 16,659 bp and the sequences are already aligned (no MAFFT).
#
# BAM list: 02_info/bam_mito.filelist = bam.filelist (the 405 Kivalliq fish,
#   from 00_setup.sh) + JAY (Arctic reference) + LLS (Atlantic reference)
#   + DV (outgroup, used only to root the tree in 12B). Built by "submit".
#
# From angsd_pipeline2/:
#   1) bash 01_scripts/bash/12A_mito_extract.sh submit   (builds the list, submits the array)
#   2) when all jobs are done:
#      bash 01_scripts/bash/12A_mito_extract.sh merge
# =============================================================================
#SBATCH -J "12A_mito"
#SBATCH -o logs/12A_mito_%A_%a.out
#SBATCH -c 1
#SBATCH -p small
#SBATCH --time=0-01:00
#SBATCH --mem=4G

set -uo pipefail
cd "${SLURM_SUBMIT_DIR:-$(pwd)}"
source config/00_config.sh
set +e
module load samtools 2>/dev/null

BAMLIST="${INFO_DIR}/bam_mito.filelist"
MITO_CHR="NC_000861.1"
MITO_LEN=16659
OUT_DIR="${BASE_DIR}/24_mito"
CONS_DIR="${OUT_DIR}/consensus"
mkdir -p "$CONS_DIR" logs
printf "%s\t%s\n" "$MITO_CHR" "$MITO_LEN" > "${OUT_DIR}/mito_scaffold.txt"

# Sample ID from the BAM name. Drop the processing suffixes (.trimmed, .sorted,
# .bam...), then keep the last dot-separated field, so that
#   NS.1979.003.N701---S520.saalJAYs_013-17_1.trimmed.sorted.bam...  -> saalJAYs_013-17_1
#   Ba-A18-21_1.trimmed.sorted.bam...                                -> Ba-A18-21_1
sample_id() { basename "$1" | sed -E 's/\.(trimmed|sorted|dedup|no_overlap|merged|markdup|rmdup|realigned|filtered|bam).*$//; s/^.*\.//'; }

# ---- submit mode: build the BAM list and submit the array --------------------
if [[ "${1:-}" == "submit" ]]; then
    cat "${INFO_DIR}/bam.filelist" \
        "${INFO_DIR}/bamlists/JAY.bamlist" \
        "${INFO_DIR}/bamlists/LLS.bamlist" \
        "${INFO_DIR}/bamlists/DV.bamlist" | sed '/^\s*$/d' | sort -u > "$BAMLIST"
    N=$(wc -l < "$BAMLIST")
    NDUP=$(while read -r b; do sample_id "$b"; done < "$BAMLIST" | sort | uniq -d | wc -l)
    log_msg "BAM list: $BAMLIST ($N BAMs, $NDUP duplicated sample IDs)"
    if [[ "$NDUP" -gt 0 ]]; then
        log_msg "ERROR: duplicated IDs. BAMs concerned (first 6):"
        DUPS=$(while read -r b; do sample_id "$b"; done < "$BAMLIST" | sort | uniq -d)
        while read -r b; do id=$(sample_id "$b"); grep -qxF "$id" <<< "$DUPS" && echo "  $id <- $b"; done < "$BAMLIST" | head -6
        exit 1
    fi
    sbatch --array=1-${N}%50 "$0"
    exit 0
fi

# ---- merge mode: one multi-FASTA + a QC table --------------------------------
if [[ "${1:-}" == "merge" ]]; then
    ALN="${OUT_DIR}/mito_all.fasta"
    cat "${CONS_DIR}"/*.fa > "$ALN"
    echo -e "sample\tlength\tN_count\tN_prop\tmean_depth" > "${OUT_DIR}/mito_qc.tsv"
    for f in "${CONS_DIR}"/*.fa; do
        id=$(basename "$f" .fa)
        seq=$(grep -v "^>" "$f" | tr -d '\n')
        len=${#seq}; n=$(echo -n "$seq" | tr -cd 'Nn' | wc -c)
        dp=$(cat "${CONS_DIR}/${id}.depth" 2>/dev/null || echo NA)
        echo -e "${id}\t${len}\t${n}\t$(echo "scale=4; $n/$len" | bc)\t${dp}"
    done >> "${OUT_DIR}/mito_qc.tsv"
    log_msg "Merged $(grep -c '^>' "$ALN") sequences into $ALN"
    log_msg "Lengths found: $(awk 'NR>1{print $2}' "${OUT_DIR}/mito_qc.tsv" | sort -u | tr '\n' ' ')"
    exit 0
fi

# ---- one individual per array task -------------------------------------------
BAM=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$BAMLIST")
[[ -z "$BAM" ]] && { log_msg "No BAM at line ${SLURM_ARRAY_TASK_ID}"; exit 0; }
[[ -f "$BAM" ]] || { log_msg "Missing BAM: $BAM"; exit 1; }
ID=$(sample_id "$BAM")
OUT_FA="${CONS_DIR}/${ID}.fa"
[[ -s "$OUT_FA" ]] && { log_msg "$ID already done"; exit 0; }

# mean depth on the mitogenome (QC)
samtools depth -a -r "$MITO_CHR" -Q 20 "$BAM" | awk '{s+=$3} END {printf "%.0f\n", s/NR}' > "${CONS_DIR}/${ID}.depth"

# simple majority consensus: a base is called if >= 70% of reads agree and
# depth >= 10, otherwise N. No insertions; deletions kept as gaps ("*" -> "-")
# so that every position of the reference is kept -> fixed length.
samtools consensus -r "$MITO_CHR" -a --show-ins no --show-del yes \
    --mode simple --call-fract 0.7 --min-depth 10 \
    --min-MQ 20 --min-BQ 20 -f fasta "$BAM" |
    awk -v id="$ID" '/^>/ {print ">" id; next} {gsub(/\*/, "-"); print}' > "${OUT_FA}.tmp" && mv "${OUT_FA}.tmp" "$OUT_FA"

LEN=$(grep -v "^>" "$OUT_FA" | tr -d '\n' | wc -c)
NN=$(grep -v "^>" "$OUT_FA" | tr -cd 'Nn' | wc -c)
log_msg "${ID}: ${LEN} bp, ${NN} N, mean depth $(cat "${CONS_DIR}/${ID}.depth")x"
[[ "$LEN" -ne "$MITO_LEN" ]] && log_msg "WARNING: ${ID} length ${LEN} != ${MITO_LEN}"
exit 0
