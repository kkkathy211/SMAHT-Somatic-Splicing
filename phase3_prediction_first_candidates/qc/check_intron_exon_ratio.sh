#!/bin/bash
# ============================================================================
# Phase 3 - Prediction-first candidate selection
#
# Intron / exon coverage ratio: distinguishes intron retention from simple over-expression.
#
# Original location in the analysis project:
#     0716/code/12_intron_exon_ratio_SLC10A6.sh
#
# NOTE: paths below may still point at the original environment
#       (local Mac or Broad cluster). See the phase README.
# ============================================================================
# 12_intron_exon_ratio_SLC10A6.sh
# ================================
# Compare intron 1 coverage vs exon 1 coverage per donor.
#
# If SMHT029 has abnormally HIGH intron/exon ratio -> real intron retention (Story A)
# If SMHT029's ratio is similar to cohort -> just high expresser (Story B)
#
# Regions:
#   Intron 1: chr4:86,834,000-86,847,000  (interior of intron 1, avoiding
#             splice site edges — 13 kb window)
#   Exon 1:   chr4:86,848,739-86,849,384  (canonical exon 1)

BAM_DIR=/medpop/esp2/yxliu/SMAHT/specific_BAM_for_Yilin
INTRON_R="chr4:86834000-86847000"
EXON_R="chr4:86848739-86849384"

printf "\n%-15s  %-12s  %-12s  %-12s  %s\n" "Donor" "Intron_mean" "Exon_mean" "Int/Exo" "Marker"
printf "%-15s  %-12s  %-12s  %-12s  %s\n" "-----" "-----------" "---------" "-------" "------"

for BAM in $BAM_DIR/*-3AD-*.bam; do
    SAMPLE=$(basename "$BAM" | cut -d'-' -f1)
    INT=$(samtools depth -a -r "$INTRON_R" "$BAM" 2>/dev/null | \
          awk '{s+=$3;n++} END{if(n>0) printf "%.3f",s/n; else print "0"}')
    EX=$(samtools depth -a -r "$EXON_R" "$BAM" 2>/dev/null | \
         awk '{s+=$3;n++} END{if(n>0) printf "%.3f",s/n; else print "0"}')
    RATIO=$(awk -v i="$INT" -v e="$EX" 'BEGIN{if(e>0) printf "%.5f",i/e; else print "NA"}')
    MARKER=""
    [[ "$BAM" == *SMHT029* ]] && MARKER=" ⭐ CARRIER1 (cryptic ψ=0.032)"
    [[ "$BAM" == *SMHT039* ]] && MARKER=" ⭐ CARRIER2 (cryptic ψ=0.018)"
    printf "  %-13s  %-12s  %-12s  %-12s  %s\n" "$SAMPLE" "$INT" "$EX" "$RATIO" "$MARKER"
done

echo ""
echo "Interpretation:"
echo "  Story A (TR-driven intron retention): SMHT029 has HIGHER Int/Exo ratio than cohort median"
echo "  Story B (just high expresser):         SMHT029's Int/Exo ratio is similar to cohort"
