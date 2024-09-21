# these are library functions for loading differnt prep data types for
# different purposes.
library(data.table)

# Loading differnt rloop type and some stats numbers for the regions
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
        fp <- paste0("input/rloop_bed_files/_archive/bedfiles_2021-05-08/", rloop_type, ".bed")
        b <- fread(fp)
        b$len <- abs(b$V3 - b$V2)
        gsum <- sum(b$len)
        return(gsum)
    }

    convert_types$genome_covered <- sapply(convert_types$all, count_bed_gnome_region)

    # Calculate the snv density along rloop regions
    count_intersected_snv <- function(rloop_type) {
        # fp <- paste0("analysis/ICGC_OV/", rloop_type, "/intersected.txt")
        fp <- paste0("_archive/analysis/analysis_2021-05-14_using_bedtools_intersection/ICGC_OV/", rloop_type, "/input/intersected.txt")
        # -- Using the bedtools intersected mutations data to count diferent values
        # _archive/analysis_2021-05-14_using_bedtools_intersection/ICGC_OV/norloop_exon/input/intersected.txt
        b <- fread(fp)
        bsvn <- b[mut_type == "SNP"]
        return(nrow(bsvn))
    }
    convert_types$snv_count <- sapply(convert_types$all, count_intersected_snv)
    # Multiple with 1 million - so the density will be how many changes per
    # million nucleotides
    convert_types$snv_density <- convert_types$snv_count / convert_types$genome_covered * 1e6

    # Obtain the gc content for the different rloop regions
    calc_region_gcpct <- function(region) {
        d <- fread(paste0("input/rloop_bed_files/_archive/bedfiles_2021-05-08/test_bedtools_nuc/bedtools_nuc_results/nuc_", region, ".bed_result.tsv"))
        # d <- fread("input/rloop_bed_files/test_bedtools/bedtools_nuc_results/nuc_rloop_exon.bed_result.tsv")
        d$length <- d[[3]] - d[[2]]
        d$gc_count <- d[[6]] * d$length # gc pct is at column 6
        gcpct <- sum(d$gc_count) / sum(d$length)
        return(gcpct)
    }
    convert_types$gcpct <- sapply(convert_types$all, calc_region_gcpct)

    return(convert_types)
}

# Loading differnt rloop type and some stats numbers for the regions
get_rloop_type2 <- function() {
    # This is a modified version of rloop types, with only tss, and tts
    # Get the different rloop type and genomic region stat info

    # Check genome regions covered by each bed files
    convert_types <- expand.grid(
        c("rloop_", "norloop_"),
        c("tss", "tts", "genebody", "pseudogene", "lincRNA")
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

    # Obtain the gc content for the different rloop regions
    calc_region_gcpct <- function(region) {
        d <- fread(paste0("input/rloop_bed_files/bedtools_nuc_results/nuc_", region, ".bed_result.tsv"))
        # d <- fread("input/rloop_bed_files/test_bedtools/bedtools_nuc_results/nuc_rloop_exon.bed_result.tsv")
        d$length <- d[[3]] - d[[2]]
        # d$gc_count <- d[[6]] * d$length # gc pct is at column 6
        d$gc_count <- d[[5]] * d$length # Notice that different version of bedtools have different nuc gc content output. The current output is at column 5 for the bedfile version  2021-05-20
        gcpct <- sum(d$gc_count) / sum(d$length)
        return(gcpct)
    }
    convert_types$gcpct <- sapply(convert_types$all, calc_region_gcpct)

    return(convert_types)
}

calc_region_gcpct_field <- function(input_file, gc_field_num = 5) {
    # Different befiles have different numbers of columns, which affect the
    # position of the gc pct calculated fields
    # d <- fread(paste0("input/rloop_bed_files/bedtools_nuc_results/nuc_", region, ".bed_result.tsv"))
    d <- fread(input_file)
    # d <- fread("input/rloop_bed_files/test_bedtools/bedtools_nuc_results/nuc_rloop_exon.bed_result.tsv")
    d$length <- d[[3]] - d[[2]]
    # d$gc_count <- d[[6]] * d$length # gc pct is at column 6
    d$gc_count <- d[[gc_field_num]] * d$length # Notice that different version of bedtools have different nuc gc content output. The current output is at column 5 for the bedfile version  2021-05-20
    gcpct <- sum(d$gc_count) / sum(d$length)
    return(gcpct)
}

# Loading differnt rloop type and some stats numbers for the regions
get_rloop_type3 <- function() {
    # This is a modified version of rloop types, with only tss, and tts
    # Get the different rloop type and genomic region stat info

    # Check genome regions covered by each bed files
    convert_types <- expand.grid(
        c("rloop_", "norloop_"),
        c("intergenic", "genebody", "tss", "tts")
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

    # Obtain the gc content for the different rloop regions
    calc_region_gcpct <- function(region) {
        d <- fread(paste0("input/rloop_bed_files/bedtools_nuc_results/nuc_", region, ".bed_result.tsv"))
        # d <- fread("input/rloop_bed_files/test_bedtools/bedtools_nuc_results/nuc_rloop_exon.bed_result.tsv")
        d$length <- d[[3]] - d[[2]]
        # d$gc_count <- d[[6]] * d$length # gc pct is at column 6
        d$gc_count <- d[[5]] * d$length # Notice that different version of bedtools have different nuc gc content output. The current output is at column 5 for the bedfile version  2021-05-20
        gcpct <- sum(d$gc_count) / sum(d$length)
        return(gcpct)
    }
    convert_types$gcpct <- sapply(convert_types$all, calc_region_gcpct)

    return(convert_types)
}

# Loading differnt rloop type and some stats numbers for the regions
get_rloop_transcribed <- function() {
    # This is a modified version of rloop types, with only tss, and tts, with
    # transcribed intersected, vs. no transcribed region. tx or -notx

    # Check genome regions covered by each bed files
    convert_types <- expand.grid(
        c("rloop", "norloop"),
        c(
            "_tss", "_tts", "_tss_transcribed", "_tts_transcribed",
            "_tss_exome", "_tts_exome", "_tss_transcribed_exome", "_tts_transcribed_exome",
            "_consensus_sc200"
        )
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

    # Obtain the gc content for the different rloop regions
    calc_region_gcpct <- function(region) {
        d <- fread(paste0("input/rloop_bed_files/bedtools_nuc_results/nuc_", region, ".bed_result.tsv"))
        # d <- fread("input/rloop_bed_files/test_bedtools/bedtools_nuc_results/nuc_rloop_exon.bed_result.tsv")
        d$length <- d[[3]] - d[[2]]
        # d$gc_count <- d[[6]] * d$length # gc pct is at column 6
        d$gc_count <- d[[5]] * d$length # Notice that different version of bedtools have different nuc gc content output. The current output is at column 5 for the bedfile version  2021-05-20
        gcpct <- sum(d$gc_count) / sum(d$length)
        return(gcpct)
    }

    convert_types$gcpct <- sapply(convert_types$all, calc_region_gcpct)

    return(convert_types)
}

get_rloop_pseudo_tss_tts <- function() {
    # This is a modified version of rloop types, with only tss, and tts, with
    # transcribed intersected, vs. no transcribed region. tx or -notx

    # Check genome regions covered by each bed files
    convert_types <- expand.grid(
        c("rloop", "norloop"),
        c("_tss", "_tts"),
        c("_gene", "_pseudogene", "_lincRNA")
    )
    convert_types$all <- paste0(
        convert_types$Var1,
        convert_types$Var2,
        convert_types$Var3
    )

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

    # Obtain the gc content for the different rloop regions
    calc_region_gcpct <- function(region) {
        d <- fread(paste0("input/rloop_bed_files/bedtools_nuc_results/nuc_", region, ".bed_result.tsv"))
        # d <- fread("input/rloop_bed_files/test_bedtools/bedtools_nuc_results/nuc_rloop_exon.bed_result.tsv")
        d$length <- d[[3]] - d[[2]]
        # d$gc_count <- d[[6]] * d$length # gc pct is at column 6
        d$gc_count <- d[[5]] * d$length # Notice that different version of bedtools have different nuc gc content output. The current output is at column 5 for the bedfile version  2021-05-20
        gcpct <- sum(d$gc_count) / sum(d$length)
        return(gcpct)
    }

    convert_types$gcpct <- sapply(convert_types$all, calc_region_gcpct)

    return(convert_types)
}

# Obtain the gc content for the different rloop regions
calc_region_gcpct <- function(path, region) {
    # d <- fread(paste0("input/rloop_bed_files/bedtools_nuc_results/nuc_", region, ".bed_result.tsv"))
    d <- fread(paste0(path, "/gcpct_bedtools/nuc_", region, ".bed_result.tsv"))
    # d <- fread("input/rloop_bed_files/test_bedtools/bedtools_nuc_results/nuc_rloop_exon.bed_result.tsv")
    d$length <- d[[3]] - d[[2]]
    # d$gc_count <- d[[6]] * d$length # gc pct is at column 6
    d$gc_count <- d[[5]] * d$length # Notice that different version of bedtools have different nuc gc content output. The current output is at column 5 for the bedfile version  2021-05-20
    gcpct <- sum(d$gc_count) / sum(d$length)
    return(gcpct)
}

# Calculate different rloop types genome region length
count_bed_gnome_region <- function(path, rloop_type) {
    # fp <- paste0("input/rloop_bed_files/", rloop_type, ".bed")
    fp <- paste0(path, "/", rloop_type, ".bed")
    b <- fread(fp)
    b$len <- abs(b$V3 - b$V2)
    gsum <- sum(b$len)
    return(gsum)
}

get_rloop_wholeGenome <- function(path = "analysis/structural_variants/input/whole-genome-rloops-bed/final_output_beds") {
    # This is a modified version of rloop types for whole genome regions
    # - Including consensus no cutoff region, score 100/200, and negative regions
    convert_types <- data.table(all = c(
        "consensus_no_cutoff",
        "consensus_sc_gt_100",
        "consensus_sc_gt_200",
        "negative_regions_without_consensus_rloop"
    ))

    local_count_bed_gnome_region <- function(x) count_bed_gnome_region(path, x)
    local_calc_region_gcpct <- function(x) calc_region_gcpct(path, x)

    convert_types$genome_covered <- sapply(
        convert_types$all,
        local_count_bed_gnome_region
    )
    convert_types$gcpct <- sapply(
        convert_types$all,
        local_calc_region_gcpct
    )

    return(convert_types)
}
