README for input data files:

* Rloop region bed files: 
    Generated: in R-loop_consensus_regions/output

These files are linked here from the above output folder:

    Consensus_no_cutoffs.txt.bed_SORTED.bed.gz (compressed to save space)

    Consensus_sc_gt_100.txt_length_lt_5_num_gt_5.bed_SORTED.bed
    Consensus_sc_gt_200.txt_length_lt_5_num_gt_5.bed_SORTED.bed

* Annotation bed file for different regions
    Around protein coding genes
        genebody-sorted.bed (protein coding region)
        lincRNA-sorted.bed
        pseudogene-sorted.bed
        tss-1kb-window-sorted.bed
        tts-1kb-window-sorted.bed

* Downloaded exom target file from: https://gdc.cancer.gov/about-data/publications/mc3-2017
gencode.v19.basic.exome.bed:
    this is the standard exome target file

gaf_20111020Plusbroad_wex_1.1_hg19.bed:
    this isn't very clear: the regions covered here is only about half of the
    size of the exome regions.

* Sorted TSS/TES(TTS) bed regions (sorted cleaned for intersection with Rloop regions)
    tss-1kb-window-sorted.bed
    tss-transcribed-1kb-window-sorted.bed
    tts-1kb-window-sorted.bed
    tts-transcribed-1kb-window-sorted.bed
