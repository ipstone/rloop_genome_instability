# Convert the provided rloop related bedtools to compatable bedfiles
rm(list = ls())
library(data.table)
library(tidyverse)
library(ipfun)

# 2 arguments are needed for the script: 1 - input file, 2- output folder
args <- commandArgs(trailingOnly = TRUE)
input_file = args[1]
output_folder = args[2]
mkdirp(output_folder)

output_file = paste0(output_folder, input_file)

# cat(args, sep = "\n") # returned list of given argument
# -- except script name and other common optons
    
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

fwrite(clean_bed(input_file), output_file , col.names = FALSE, sep = "\t")


