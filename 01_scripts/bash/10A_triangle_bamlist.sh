#!/bin/bash
#SBATCH --job-name=tri_bamlist
#SBATCH --partition=small
#SBATCH --time=1-00:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --output=99_logs/10A_%j.out
#SBATCH --error=99_logs/10A_%j.err

# =============================================================================
# 10A_triangle_bamlist.sh — Build combined bamlist: Kivalliq (13 pops) + LLS + JAY
# Output: 02_info/bamlists/triangle.bamlist (1 BAM path per line)
#         02_info/bamlists/triangle.popmap   (BAM \t POP \t REGION)
# =============================================================================

set -euo pipefail

source config/00_config.sh
mkdir -p 99_logs 02_info/bamlists

OUT_BAM="02_info/bamlists/triangle.bamlist"
OUT_POP="02_info/bamlists/triangle.popmap"

KIVALLIQ_POPS=(AKL AUL CRB DIA MEL HOR ITI KGJ NOP PAM SUP TIN WHI)
ANCHOR_POPS=(LLS JAY)

> "$OUT_BAM"
> "$OUT_POP"

# Helper: tag region from pop code
region_for () {
    case "$1" in
        AKL|AUL|CRB|DIA|MEL) echo "Rankin" ;;
        HOR)                  echo "Baker" ;;
        ITI|KGJ|NOP|PAM|SUP|TIN|WHI) echo "Naujaat" ;;
        LLS) echo "Atlantic_anchor" ;;
        JAY) echo "Arctic_anchor" ;;
        *)   echo "Unknown" ;;
    esac
}

n_total=0
for POP in "${KIVALLIQ_POPS[@]}" "${ANCHOR_POPS[@]}"; do
    BAMLIST="02_info/bamlists/${POP}.bamlist"
    if [[ ! -s "$BAMLIST" ]]; then
        echo "ERROR: Missing or empty bamlist: $BAMLIST" >&2
        exit 1
    fi
    REG=$(region_for "$POP")
    n=0
    while IFS= read -r BAM; do
        [[ -z "$BAM" ]] && continue
        echo "$BAM"             >> "$OUT_BAM"
        echo -e "${BAM}\t${POP}\t${REG}" >> "$OUT_POP"
        n=$((n+1))
    done < "$BAMLIST"
    printf "  %-4s  %-15s  %3d ind\n" "$POP" "$REG" "$n"
    n_total=$((n_total + n))
done

echo ""
echo "Total individuals: $n_total"
echo "Bamlist: $OUT_BAM"
echo "Popmap : $OUT_POP"
