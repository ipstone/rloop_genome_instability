# Here are some misc codes previously used to clean bed files format

## bed files also need conversion from the chrX... to X format
normal_chroms <- c(paste0("chr", seq(1, 22)), "chrX", "chrY")

clean_bed <- function(x) {
    # convert the input bedfile x, to cleaned bed file (return data.table)
    # rloop_exon <- fread("input/rloop_bed_files/positive_set/HIGH_CONFIDENCE_CONSENSUS_PEAK_SCORE_GT_100_NUM_GT_5_ovlp_EXONS.txt_converted_to_bed.bed")
    rloop_exon <- fread(x)
    rloop_exon <- rloop_exon[V1 %in% normal_chroms]
    rloop_exon$chrom <- gsub("chr", "", rloop_exon$V1, fixed = T)
    rloop_sorted <- rloop_exon[order(chrom, V2), .(chrom, V2, V3, V4)]
    return(rloop_sorted)
}

save_cleaned_bed_files <- function() {
    # Save all the rloop, norloop bedfiles
    # Rloop positive regions
    clean_bed("input/rloop_bed_files/positive_set/HIGH_CONFIDENCE_CONSENSUS_PEAK_SCORE_GT_100_NUM_GT_5_ovlp_EXONS.txt_converted_to_bed.bed") %>%
        fwrite("input/rloop_bed_files/rloop_exon.bed", col.names = FALSE, sep = "\t")

    clean_bed("input/rloop_bed_files/positive_set/HIGH_CONFIDENCE_CONSENSUS_PEAK_SCORE_GT_100_NUM_GT_5_ovlp_GENE_BODY.txt_converted_to_bed.bed") %>%
        fwrite("input/rloop_bed_files/rloop_genebody.bed", col.names = FALSE, sep = "\t")

    clean_bed("input/rloop_bed_files/positive_set/HIGH_CONFIDENCE_CONSENSUS_PEAK_SCORE_GT_100_NUM_GT_5_ovlp_TSS.txt_converted_to_bed.bed") %>%
        fwrite("input/rloop_bed_files/rloop_tss.bed", col.names = FALSE, sep = "\t")

    clean_bed("input/rloop_bed_files/positive_set/HIGH_CONFIDENCE_CONSENSUS_PEAK_SCORE_GT_100_NUM_GT_5_ovlp_TTS.txt_converted_to_bed.bed") %>%
        fwrite("input/rloop_bed_files/rloop_tts.bed", col.names = FALSE, sep = "\t")

    # Rloop negative regions
    clean_bed("input/rloop_bed_files/negative_set/EXONS_REFERENCE_SET.bed") %>%
        fwrite("input/rloop_bed_files/norloop_exon.bed", col.names = FALSE, sep = "\t")

    clean_bed("input/rloop_bed_files/negative_set/GENE_BODY_REFERENCE_SET.bed") %>%
        fwrite("input/rloop_bed_files/norloop_genebody.bed", col.names = FALSE, sep = "\t")

    clean_bed("input/rloop_bed_files/negative_set/TSS_REFERENCE_SET.bed") %>%
        fwrite("input/rloop_bed_files/norloop_tss.bed", col.names = FALSE, sep = "\t")

    clean_bed("input/rloop_bed_files/negative_set/TTS_REFERENCE_SET.bed") %>%
        fwrite("input/rloop_bed_files/norloop_tts.bed", col.names = FALSE, sep = "\t")
}

# Get the region bp length for different rloop bed file region and snv counts on these
# -- TODO: probably this is not needed as we are moving these code to the lib
# load dataprep R code
get_rloop_type <- function() {
    # Get the different rloop type and genomic region stat info

    # Check genome regions covered by each bed files
    convert_types <- expand.grid(
        c("rloop_", "norloop_"),
        c("exon", "genebody", "tss", "tts")
    )
    convert_types$all <- paste0(convert_types$Var1, convert_types$Var2)
    convert_types <- data.table(convert_types)

    # Calculate different rloop types genome region length
    count_bed_gnome_region <- function(rloop_type) {
        fp <- paste0("input/rloop_bed_files/", rloop_type, ".bed")
        b <- fread(fp)
        b$len <- abs(b$V3 - b$V2)
        gsum <- sum(b$len)
        return(gsum)
    }

    convert_types$genome_covered <- sapply(convert_types$all, count_bed_gnome_region)

    # Calculate the snv density along rloop regions
    count_intersected_snv <- function(rloop_type) {
        fp <- paste0("analysis/ICGC_OV/", rloop_type, "/intersected.txt")
        b <- fread(fp)
        bsvn <- b[mut_type == "SNP"]
        return(nrow(bsvn))
    }
    convert_types$snv_count <- sapply(convert_types$all, count_intersected_snv)
    convert_types$snv_density <- convert_types$snv_count / convert_types$genome_covered
    return(convert_types)
}
