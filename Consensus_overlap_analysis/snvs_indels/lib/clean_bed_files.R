# Convert the provided rloop related bedtools to compatable bedfiles
rm(list = ls())
library(data.table)
library(tidyverse)
library(ipfun)

convert_types <- expand.grid(
    c("rloop_", "norloop_"),
    c("genebody", "intergenic", "tss", "tts")
)
convert_types_list = (paste0(convert_types$Var1, convert_types$Var2))

clean_bed <- function(x) {
    # convert the input bedfile x, to cleaned bed file (return data.table)
    # rloop_exon <- fread("input/rloop_bed_files/positive_set/HIGH_CONFIDENCE_CONSENSUS_PEAK_SCORE_GT_100_NUM_GT_5_ovlp_EXONS.txt_converted_to_bed.bed")
    ## bed files also need conversion from the chrX... to X format
    normal_chroms <- c(paste0("chr", seq(1, 22)), "chrX", "chrY")
    rloop_exon <- fread(x)
    rloop_exon <- rloop_exon[V1 %in% normal_chroms]
    rloop_exon$chrom <- gsub("chr", "", rloop_exon$V1, fixed = T)
    rloop_sorted <- rloop_exon[order(chrom, V2), .(chrom, V2, V3)]
    return(rloop_sorted)
}

r = list()
r$norloop_genebody<-clean_bed( "MAY_17TH_2021_FILES_FOR_ISAAC/negative_set/GENE_BODY_REFERENCE_SET_SUBTRACT.bed")
r$norloop_intergenic<-clean_bed( "MAY_17TH_2021_FILES_FOR_ISAAC/negative_set/INTERGENIC_REFERENCE_SET.bed")
r$norloop_tss<-clean_bed( "MAY_17TH_2021_FILES_FOR_ISAAC/negative_set/TSS_REFERENCE_SET.bed")
r$norloop_tts<-clean_bed( "MAY_17TH_2021_FILES_FOR_ISAAC/negative_set/TTS_REFERENCE_SET.bed")
r$rloop_genebody<-clean_bed( "MAY_17TH_2021_FILES_FOR_ISAAC/positive_set_score_gt_200/GENEBODY_R_LOOP_POSITIVE_SCORE_GT_200.bed")
r$rloop_intergenic<-clean_bed( "MAY_17TH_2021_FILES_FOR_ISAAC/positive_set_score_gt_200/INTERGENIC_R_LOOP_POSITIVE_SCORE_GT_200.bed")
r$rloop_tss<-clean_bed( "MAY_17TH_2021_FILES_FOR_ISAAC/positive_set_score_gt_200/TSS_R_LOOP_POSITIVE_SCORE_GT_200.bed")
r$rloop_tts<-clean_bed( "MAY_17TH_2021_FILES_FOR_ISAAC/positive_set_score_gt_200/TTS_R_LOOP_POSITIVE_SCORE_GT_200.bed")


save_clean_rloop_bedfiles = function(x){
    print(paste("Working on: ", x))
    fwrite(r[[x]],paste0(x, ".bed"), col.names = FALSE, sep = "\t")
} 

sapply(convert_types_list, save_clean_rloop_bedfiles)
