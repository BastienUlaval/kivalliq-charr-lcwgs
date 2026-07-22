#!/bin/bash
#SBATCH --job-name=tri_plot
#SBATCH --partition=small
#SBATCH --time=1-00:00
#SBATCH --mem=24G
#SBATCH --cpus-per-task=2
#SBATCH --output=99_logs/10F_%j.out
#SBATCH --error=99_logs/10F_%j.err

# =============================================================================
# 10F_triangle_plot.sh -- Run triangle_plot_v2.R: hybrid index + interspecific
# heterozygosity (Fig_Triangle_combined)
# =============================================================================

set -euo pipefail
source config/00_config.sh

OUTDIR="26_triangle"
FIG_DIR="04_figures/triangle"

mkdir -p "$FIG_DIR"

[[ -s "$OUTDIR/triangle.beagle.gz"     ]] || { echo "Missing Beagle. Run 10E_triangle_diagnostic_beagle.sh." >&2; exit 1; }
[[ -s "$OUTDIR/diagnostic_snps.tsv"   ]] || { echo "Missing diagnostic_snps.tsv. Run 10E_triangle_diagnostic_beagle.sh." >&2; exit 1; }
[[ -s "02_info/bamlists/triangle.popmap" ]] || { echo "Missing popmap. Run 10A_triangle_bamlist.sh." >&2; exit 1; }

# Restrict diagnostic SNPs to the ngsParalog-canonical (non-deviant) set --
# see triangle_plot_v2.R header. Falls back to no masking (with a warning)
# if this file doesn't exist yet.
export CANON_SITES="${INFO_DIR}/sites_all_${SUFFIX}_canonical"

"$RSCRIPT" 01_scripts/R/triangle_plot_v2.R \
    "$OUTDIR" \
    "02_info/bamlists/triangle.popmap" \
    "$FIG_DIR"

echo ""
echo "Figures in: $FIG_DIR"
ls -la "$FIG_DIR"
