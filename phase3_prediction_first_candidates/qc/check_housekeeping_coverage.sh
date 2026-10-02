#!/bin/bash
# ============================================================================
# Phase 3 - Prediction-first candidate selection
#
# Is a suspect donor low everywhere (RNA quality) or only in the gene of interest (real biology)?
#
# Original location in the analysis project:
#     0716/code/10_qc_smht029_brain.sh
#
# NOTE: paths below may still point at the original environment
#       (local Mac or Broad cluster). See the phase README.
# ============================================================================
# 10_qc_smht029_brain.sh
# ======================
# Sanity check: is SMHT029 brain generally low-coverage (RIN/quality issue)?
# Or is RBFOX1 specifically down in this donor?
#
# Compares mean depth in 3 housekeeping genes (GAPDH/ACTB/GUSB) vs RBFOX1 zoom
# across all 12 samples in the 3AL cohort.
#
# If SMHT029 is LOW everywhere -> QC issue, discard RBFOX1 finding
# If SMHT029 is NORMAL in housekeeping but LOW in RBFOX1 -> real biology

BAM_DIR=/medpop/esp2/yxliu/SMAHT/specific_BAM_for_Yilin
CARRIER_BAM=$BAM_DIR/SMHT029-3AL-006A2-M72-A101-nygc-SMAFI4T62QHB-star_2.7.10b_GRCh38_gencode_v47.aligned.sorted.bam

echo "=== Housekeeping vs RBFOX1 coverage: SMHT029 (carrier) vs cohort ==="
for GENE in "chr12:6533927-6538371:GAPDH" \
            "chr7:5527151-5563902:ACTB" \
            "chr7:65960687-65982150:GUSB" \
            "chr16:7600000-7720000:RBFOX1_zoom"; do
    REGION=$(echo "$GENE" | cut -d':' -f1-2)
    NAME=$(echo "$GENE" | cut -d':' -f4)
    printf "\n--- %s (%s) ---\n" "$NAME" "$REGION"
    for BAM in "$BAM_DIR"/*-3AL-*.bam; do
        SAMPLE=$(basename "$BAM" | awk -F'-' '{print $1"-"$2"-"$4"-"$6}')
        DEPTH=$(samtools depth -a -r "$REGION" "$BAM" 2>/dev/null | \
                awk '{s+=$3;n++} END{if(n>0) printf "%.2f",s/n; else print "0"}')
        MARKER=""
        if [[ "$BAM" == "$CARRIER_BAM" ]]; then MARKER=" ⭐ CARRIER"; fi
        printf "  %-40s mean_depth=%s%s\n" "$SAMPLE" "$DEPTH" "$MARKER"
    done
done

echo ""
echo "Done. Interpretation:"
echo "  If SMHT029 is low in RBFOX1 but NORMAL in housekeeping -> real RBFOX1 down"
echo "  If SMHT029 is low everywhere -> RNA quality issue, discard"
