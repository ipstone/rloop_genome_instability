# Clean previous calculation result 
rm -rf tmp/final_output_beds
rm tmp/*.bed 

# Intersect for the positive and negative R-loop regions

neg_bed=input/Consensus_no_cutoffs.txt.bed_SORTED.bed
pos_bed=input/Consensus_sc_gt_200.txt_length_lt_5_num_gt_5.bed_SORTED.bed
exome=input/exome_chr-sorted.bed

## TSS
# intersecting the negative reference set
bedtools intersect -a input/tss-1kb-window-sorted.bed \
    -b $neg_bed -v > \
    tmp/norloop_tss.bed
# intersecting the pos reference set
bedtools intersect -a input/tss-1kb-window-sorted.bed \
    -b $pos_bed -wa > \
    tmp/rloop_tss.bed

## TTS
bedtools intersect -a input/tts-1kb-window-sorted.bed \
    -b $neg_bed -v > \
    tmp/norloop_tts.bed
bedtools intersect -a input/tts-1kb-window-sorted.bed \
    -b $pos_bed -wa > \
    tmp/rloop_tts.bed

# TSS transcribed
bedtools intersect -a input/tss-transcribed-1kb-window-sorted.bed \
    -b $neg_bed -v > \
    tmp/norloop_tss_transcribed.bed
# intersecting the pos reference set
bedtools intersect -a input/tss-transcribed-1kb-window-sorted.bed \
    -b $pos_bed -wa > \
    tmp/rloop_tss_transcribed.bed

## TTS transcribed
bedtools intersect -a input/tts-transcribed-1kb-window-sorted.bed \
    -b $neg_bed -v > \
    tmp/norloop_tts_transcribed.bed
bedtools intersect -a input/tts-transcribed-1kb-window-sorted.bed \
    -b $pos_bed -wa > \
    tmp/rloop_tts_transcribed.bed

#########################
## intersecting final tss,tts,transcribed with exome bed
for r in rloop norloop ; do
    for region in tss tts tss_transcribed tts_transcribed ; do 

    echo working on tmp/"$r"_"$region".bed 
    bedtools intersect -a tmp/"$r"_"$region".bed -b $exome > \
        tmp/"$r"_"$region"_exome.bed

    done
done
