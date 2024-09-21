## To intersect SVs with genebody, TSS, TTS for rloop and nonrloop regions
# Start with intersection score 200 regions
genebody="input/Rloop_intersected_score200/rloop_genebody.bed"
nogenebody="input/Rloop_intersected_score200/norloop_genebody.bed"
tss="input/Rloop_intersected_score200/rloop_tss.bed"
notss="input/Rloop_intersected_score200/norloop_tss.bed"
tts="input/Rloop_intersected_score200/rloop_tts.bed"
notts="input/Rloop_intersected_score200/norloop_tts.bed"

# The bedpe files for
erpos="input/BRCA-EU_SV/serena_erpos_SV.bedpe"
trpneg="input/BRCA-EU_SV/serena_tripneg_SV.bedpe"

## Section to intersect for ERpos tumors
mkdir -p output/intersect_serena_erpos_rloops/
# Run the intersection for the breakpoints
bedtools pairtobed -a $erpos -b $tss > output/intersect_serena_erpos_rloops/rloop_tss-sc200.bedpe
bedtools pairtobed -a $erpos -b $notss > output/intersect_serena_erpos_rloops/norloop_tss-sc200.bedpe
bedtools pairtobed -a $erpos -b $tts > output/intersect_serena_erpos_rloops/rloop_tts-sc200.bedpe
bedtools pairtobed -a $erpos -b $notts > output/intersect_serena_erpos_rloops/norloop_tts-sc200.bedpe
bedtools pairtobed -a $erpos -b $genebody > output/intersect_serena_erpos_rloops/rloop_genebody-sc200.bedpe
bedtools pairtobed -a $erpos -b $nogenebody > output/intersect_serena_erpos_rloops/norloop_genebody-sc200.bedpe

## To intersect SVs with genebody, TSS, TTS for rloop and nonrloop regions
genebody="input/output_isaac_score100/rloop_genebody.bed"
nogenebody="input/output_isaac_score100/norloop_genebody.bed"
tss="input/output_isaac_score100/rloop_tss.bed"
notss="input/output_isaac_score100/norloop_tss.bed"
tts="input/output_isaac_score100/rloop_tts.bed"
notts="input/output_isaac_score100/norloop_tts.bed"

# Run the intersection sc100 for the breakpoints
bedtools pairtobed -a $erpos -b $tss > output/intersect_serena_erpos_rloops/rloop_tss-sc100.bedpe
bedtools pairtobed -a $erpos -b $notss > output/intersect_serena_erpos_rloops/norloop_tss-sc100.bedpe
bedtools pairtobed -a $erpos -b $tts > output/intersect_serena_erpos_rloops/rloop_tts-sc100.bedpe
bedtools pairtobed -a $erpos -b $notts > output/intersect_serena_erpos_rloops/norloop_tts-sc100.bedpe
bedtools pairtobed -a $erpos -b $genebody > output/intersect_serena_erpos_rloops/rloop_genebody-sc100.bedpe
bedtools pairtobed -a $erpos -b $nogenebody > output/intersect_serena_erpos_rloops/norloop_genebody-sc100.bedpe

## Section to intersect with consensus 100 and 200 region
# Intersect the sv bedpe file with the rloop consensus score 100 regions
bedtools pairtobed -a $erpos -b input/consensus_sc_gt_100.bed > output/intersect_serena_erpos_rloops/consensus_sc100.bedpe
bedtools pairtobed -a $erpos -b input/consensus_sc_gt_200.bed > output/intersect_serena_erpos_rloops/consensus_sc200.bedpe
bedtools pairtobed -a $erpos -b input/consensus_no_cutoff.bed > output/intersect_serena_erpos_rloops/consensus_no-cutoff.bedpe
# The following is to obtain the negative bedpe set which are not in the
# consensus region
bedtools pairtobed -a $erpos -b input/consensus_no_cutoff.bed -type neither > output/intersect_serena_erpos_rloops/no_consensus.bedpe

## Section to intersect for trpneg tumors
mkdir -p output/intersect_serena_trpneg_rloops/
# Run the intersection for the breakpoints
bedtools pairtobed -a $trpneg -b $tss > output/intersect_serena_trpneg_rloops/rloop_tss-sc200.bedpe
bedtools pairtobed -a $trpneg -b $notss > output/intersect_serena_trpneg_rloops/norloop_tss-sc200.bedpe
bedtools pairtobed -a $trpneg -b $tts > output/intersect_serena_trpneg_rloops/rloop_tts-sc200.bedpe
bedtools pairtobed -a $trpneg -b $notts > output/intersect_serena_trpneg_rloops/norloop_tts-sc200.bedpe
bedtools pairtobed -a $trpneg -b $genebody > output/intersect_serena_trpneg_rloops/rloop_genebody-sc200.bedpe
bedtools pairtobed -a $trpneg -b $nogenebody > output/intersect_serena_trpneg_rloops/norloop_genebody-sc200.bedpe

## To intersect SVs with genebody, TSS, TTS for rloop and nonrloop regions
genebody="input/output_isaac_score100/rloop_genebody.bed"
nogenebody="input/output_isaac_score100/norloop_genebody.bed"
tss="input/output_isaac_score100/rloop_tss.bed"
notss="input/output_isaac_score100/norloop_tss.bed"
tts="input/output_isaac_score100/rloop_tts.bed"
notts="input/output_isaac_score100/norloop_tts.bed"

# Run the intersection sc100 for the breakpoints
bedtools pairtobed -a $trpneg -b $tss > output/intersect_serena_trpneg_rloops/rloop_tss-sc100.bedpe
bedtools pairtobed -a $trpneg -b $notss > output/intersect_serena_trpneg_rloops/norloop_tss-sc100.bedpe
bedtools pairtobed -a $trpneg -b $tts > output/intersect_serena_trpneg_rloops/rloop_tts-sc100.bedpe
bedtools pairtobed -a $trpneg -b $notts > output/intersect_serena_trpneg_rloops/norloop_tts-sc100.bedpe
bedtools pairtobed -a $trpneg -b $genebody > output/intersect_serena_trpneg_rloops/rloop_genebody-sc100.bedpe
bedtools pairtobed -a $trpneg -b $nogenebody > output/intersect_serena_trpneg_rloops/norloop_genebody-sc100.bedpe

## Section to intersect with consensus 100 and 200 region
# Intersect the sv bedpe file with the rloop consensus score 100 regions
bedtools pairtobed -a $trpneg -b input/consensus_sc_gt_100.bed > output/intersect_serena_trpneg_rloops/consensus_sc100.bedpe
bedtools pairtobed -a $trpneg -b input/consensus_sc_gt_200.bed > output/intersect_serena_trpneg_rloops/consensus_sc200.bedpe
bedtools pairtobed -a $trpneg -b input/consensus_no_cutoff.bed > output/intersect_serena_trpneg_rloops/consensus_no-cutoff.bedpe
# The following is to obtain the negative bedpe set which are not in the
# consensus region
bedtools pairtobed -a $trpneg -b input/consensus_no_cutoff.bed -type neither > output/intersect_serena_trpneg_rloops/no_consensus.bedpe

