#!/bin/bash
# =============================================================================
# 05_admixture.sh — NGSadmix (global, Rankin, Naujaat) with replicate runs
# Submits 3 independent jobs (one per dataset) to respect 7-day time limits.
#
# Optimal K (global=7, Rankin=2, Naujaat=5) was selected via the
# StructureSelector web tool from the resulting .qopt/.log files, not by a
# script in this repo (see manuscript Methods). Once K is chosen, the
# admixture barplots are built directly by plot_panel_global_structure.R
# and plot_panel_regional_structure.R (Fig. 2c, 3e-f).
# Usage: bash 01_scripts/bash/05_admixture.sh
# =============================================================================

cd "${SLURM_SUBMIT_DIR:-$(pwd)}"
source config/00_config.sh

GLOBAL_BEAGLE="${LD_DIR}/all_${SUFFIX}.prunednosex.beagle.gz"
RANKIN_BEAGLE="${LD_DIR}/rankin_${SUFFIX}.prunednosex.beagle.gz"
NAUJAAT_BEAGLE="${LD_DIR}/naujaat_${SUFFIX}.prunednosex.beagle.gz"

require_file "$NGSADMIX_PATH" "NGSadmix binary"
require_file "$GLOBAL_BEAGLE" "Global Beagle"
require_file "$RANKIN_BEAGLE" "Rankin Beagle"
require_file "$NAUJAAT_BEAGLE" "Naujaat Beagle"

mkdir -p "${ADMIX_DIR}/global" "${ADMIX_DIR}/rankin" "${ADMIX_DIR}/naujaat"

log_msg "=== SUBMITTING 3 ADMIXTURE JOBS ==="

# ─── Global (K=2-15, 10 reps) ───────────────────────────────────────────────
JID_GLOBAL=$(sbatch --parsable \
    -J "admix_global" \
    -o logs/admix_global_%j.out \
    -e logs/admix_global_%j.err \
    -c 5 -p medium --time=7-00:00 --mem=100G \
    -D "$BASE_DIR" \
    --wrap "
source config/00_config.sh
for k in \$(seq ${K_MIN} ${K_MAX_GLOBAL}); do
    for rep in \$(seq 1 ${NGSADMIX_NREP}); do
        seed=\$((k * 1000 + rep))
        echo \"Global K=\${k} rep=\${rep} seed=\${seed}\"
        ${NGSADMIX_PATH} -P 5 \
            -likes ${GLOBAL_BEAGLE} \
            -minMaf ${MIN_MAF} \
            -K \${k} \
            -seed \${seed} \
            -o ${ADMIX_DIR}/global/global_K\${k}_rep\${rep}
    done
done
echo 'Global DONE'
")

log_msg "  Global: job ${JID_GLOBAL} (K=${K_MIN}-${K_MAX_GLOBAL}, ${NGSADMIX_NREP} reps)"

# ─── Rankin (K=2-7, 10 reps) ────────────────────────────────────────────────
JID_RANKIN=$(sbatch --parsable \
    -J "admix_rankin" \
    -o logs/admix_rankin_%j.out \
    -e logs/admix_rankin_%j.err \
    -c 5 -p medium --time=7-00:00 --mem=50G \
    -D "$BASE_DIR" \
    --wrap "
source config/00_config.sh
for k in \$(seq ${K_MIN} ${K_MAX_RANKIN}); do
    for rep in \$(seq 1 ${NGSADMIX_NREP}); do
        seed=\$((k * 1000 + rep))
        echo \"Rankin K=\${k} rep=\${rep} seed=\${seed}\"
        ${NGSADMIX_PATH} -P 5 \
            -likes ${RANKIN_BEAGLE} \
            -minMaf ${MIN_MAF} \
            -K \${k} \
            -seed \${seed} \
            -o ${ADMIX_DIR}/rankin/rankin_K\${k}_rep\${rep}
    done
done
echo 'Rankin DONE'
")

log_msg "  Rankin: job ${JID_RANKIN} (K=${K_MIN}-${K_MAX_RANKIN}, ${NGSADMIX_NREP} reps)"

# ─── Naujaat (K=2-9, 10 reps) ───────────────────────────────────────────────
JID_NAUJAAT=$(sbatch --parsable \
    -J "admix_naujaat" \
    -o logs/admix_naujaat_%j.out \
    -e logs/admix_naujaat_%j.err \
    -c 5 -p medium --time=7-00:00 --mem=50G \
    -D "$BASE_DIR" \
    --wrap "
source config/00_config.sh
for k in \$(seq ${K_MIN} ${K_MAX_NAUJAAT}); do
    for rep in \$(seq 1 ${NGSADMIX_NREP}); do
        seed=\$((k * 1000 + rep))
        echo \"Naujaat K=\${k} rep=\${rep} seed=\${seed}\"
        ${NGSADMIX_PATH} -P 5 \
            -likes ${NAUJAAT_BEAGLE} \
            -minMaf ${MIN_MAF} \
            -K \${k} \
            -seed \${seed} \
            -o ${ADMIX_DIR}/naujaat/naujaat_K\${k}_rep\${rep}
    done
done
echo 'Naujaat DONE'
")

log_msg "  Naujaat: job ${JID_NAUJAAT} (K=${K_MIN}-${K_MAX_NAUJAAT}, ${NGSADMIX_NREP} reps)"

log_msg ""
log_msg "=== ALL 3 JOBS SUBMITTED ==="
log_msg "When all 3 finish: upload the .qopt/.log files to StructureSelector"
log_msg "to choose K, then run plot_panel_global_structure.R / plot_panel_regional_structure.R"