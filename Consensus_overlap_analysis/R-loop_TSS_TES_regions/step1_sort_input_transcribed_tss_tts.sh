bedtools sort -i input/prepare_annotation_bed_files/output/tss-1kb-window.bed > \
    input/tss-1kb-window-sorted.bed

bedtools sort -i input/prepare_annotation_bed_files/output/tts-1kb-window.bed > \
    input/tts-1kb-window-sorted.bed

bedtools sort -i input/prepare_annotation_bed_files/output/tss-transcribed-1kb-window.bed > \
    input/tss-transcribed-1kb-window-sorted.bed

bedtools sort -i input/prepare_annotation_bed_files/output/tts-transcribed-1kb-window.bed > \
    input/tts-transcribed-1kb-window-sorted.bed

bedtools sort -i input/exome_chr.bed > \
    input/exome_chr-sorted.bed

