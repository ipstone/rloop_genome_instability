The code and data in this folder are used to calculate the TSS (Transcription start site) and TES (Transcription end site, also refered as tts in the code/file-naming conventions) Rloop vs control regions.

* The bed files generated in the 'output' folder are used to intersect mutation sites (SNVs, Indels, SVs) to generate the stat tests for the figures in the paper.


* Below are the steps to prepare the bed files (in the output folder) for further analysis:

    step1_sort_input_transcribed_tss_tts.sh
    step2_score200_intersect_tss-tts-transcribed-exome.sh
    step3_process_bedfiles_clean_sort_merge.sh

* Note:

This subfolder in the input folder contains the code and notes to generate the bed files for the regions of different gene related annotations. Please refer to the README in the folder for more details.

    input
    ├── prepare_annotation_bed_files
    │   ├── output
    │   │   ├── genebody.bed
    │   │   ├── lincRNA.bed
    │   │   ├── pseudogene.bed
    │   │   ├── transcript_protein-coding_gencode-v19.bed
    │   │   ├── transcript_protein-coding_gencode-v19.tsv
    │   │   ├── tss-1kb-window.bed
    │   │   ├── tss-transcribed-1kb-window.bed
    │   │   ├── tts-1kb-window.bed
    │   │   └── tts-transcribed-1kb-window.bed
