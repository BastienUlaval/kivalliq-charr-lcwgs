#!/bin/bash
# =============================================================================
# 10_RUN_triangle_pipeline.sh -- Master orchestration (v3 -- array-based + merge)
#
# Pipeline:
#   10A -- Build combined bamlist (Kivalliq + LLS + JAY)
#   10B -- Call SNPs per chromosome (ARRAY)
#   10C -- Consolidate per-chr MAFs into site list
#   10D -- Compute MAF for LLS and JAY anchors (ARRAY)
#   10E -- Filter diagnostic SNPs + extract Beagle
#   10F -- Run triangle_plot_v2.R (Fig_Triangle_combined)
#
# Usage: bash 01_scripts/bash/10_RUN_triangle_pipeline.sh
#
# SLURM array dependency note:
#   afterok:JOBID waits for ALL tasks of an array to succeed when JOBID
#   refers to the whole array (no _N suffix). This is correct SLURM behavior.
#   However sbatch --parsable on an array returns "JOBID" only (no underscore).
#   So afterok:$JOB_B1 correctly waits for the full array.
# =============================================================================

set -euo pipefail
source config/00_config.sh

SCRIPTS="01_scripts/bash"
mkdir -p 99_logs 26_triangle 04_figures/triangle

die() { echo "ERROR: $*" >&2; exit 1; }

echo "=== Submitting triangle plot pipeline ==="
echo ""

# --- Step A: bamlist (instant, runs locally) ---------------------------------
echo "Step 10A: Building combined bamlist..."
bash "${SCRIPTS}/10A_triangle_bamlist.sh" \
    || die "10A failed"
[[ -s "02_info/bamlists/triangle.bamlist" ]] \
    || die "triangle.bamlist not created by 10A"
echo "  OK -- triangle.bamlist ready"

# --- Step B step 1: ARRAY -- call SNPs per chromosome -------------------------
echo ""
echo "Step 10B: Calling SNPs per chromosome (array)..."
JOB_B1=$(sbatch --parsable "${SCRIPTS}/10B_triangle_call_sites.sh") \
    || die "sbatch 10B failed"
[[ -n "$JOB_B1" ]] || die "sbatch 10B returned empty job ID"
echo "  Submitted array job: $JOB_B1"

# --- Step B merge: consolidate per-chr files ---------------------------------
# afterok:JOBID on a full array ID waits for ALL tasks to finish successfully.
echo ""
echo "Step 10C: Merging per-chr files..."
JOB_BM=$(sbatch --parsable \
    --dependency=afterok:${JOB_B1} \
    "${SCRIPTS}/10C_triangle_merge_sites.sh") \
    || die "sbatch 10C failed"
[[ -n "$JOB_BM" ]] || die "sbatch 10C returned empty job ID"
echo "  Submitted: $JOB_BM  (after array $JOB_B1)"

# --- Step B anchor: per-anchor MAF (LLS, JAY) --------------------------------
echo ""
echo "Step 10D: Computing per-anchor MAF (array LLS + JAY)..."
JOB_B=$(sbatch --parsable \
    --dependency=afterok:${JOB_BM} \
    "${SCRIPTS}/10D_triangle_anchor_maf.sh") \
    || die "sbatch 10D failed"
[[ -n "$JOB_B" ]] || die "sbatch 10D returned empty job ID"
echo "  Submitted: $JOB_B  (after $JOB_BM)"

# --- Step C: diagnostic SNPs + Beagle extraction -----------------------------
echo ""
echo "Step 10E: Diagnostic SNPs + Beagle..."
JOB_C=$(sbatch --parsable \
    --dependency=afterok:${JOB_B} \
    "${SCRIPTS}/10E_triangle_diagnostic_beagle.sh") \
    || die "sbatch 10E failed"
[[ -n "$JOB_C" ]] || die "sbatch 10E returned empty job ID"
echo "  Submitted: $JOB_C  (after $JOB_B)"

# --- Step D: R triangle plot --------------------------------------------------
echo ""
echo "Step 10F: R triangle plot..."
JOB_D=$(sbatch --parsable \
    --dependency=afterok:${JOB_C} \
    "${SCRIPTS}/10F_triangle_plot.sh") \
    || die "sbatch 10F failed"
[[ -n "$JOB_D" ]] || die "sbatch 10F returned empty job ID"
echo "  Submitted: $JOB_D  (after $JOB_C)"

echo ""
echo "=== Pipeline submitted successfully ==="
printf "  %-35s %s\n" "10B (per-chr SNP call array)" "$JOB_B1"
printf "  %-35s %s\n" "10C (merge per-chr sites)"                    "$JOB_BM"
printf "  %-35s %s\n" "10D (anchor MAF, LLS+JAY array)"     "$JOB_B"
printf "  %-35s %s\n" "10E (diagnostic + Beagle)"          "$JOB_C"
printf "  %-35s %s\n" "10F (R triangle plot)"              "$JOB_D"
echo ""
echo "Monitor with: squeue -u \$USER"
echo "Check logs:   ls 99_logs/10*"