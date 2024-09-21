# With the given input folder, will convert the intersected bed in the folder to simple snv txt file

rm(list = ls())
library(data.table)
library(magrittr)
library(ipfun)


args <- commandArgs(trailingOnly = TRUE)
# cat(args, sep = "\n") # returned list of given argument
# -- except script name and other common optons
print(args)
project <- args[1]
region_type = args[2] # input for the rloop region type


# project <- "ICGC_OV"
# project <- "Serena_ER"
# project <- "Serena-Davies_ER_Biallelic"

# convert the bed file back to simple snv txt
convert_bed_intersected <- function(x) {
    #d <- fread("input/ICGC_OV/icgc_ov_simple_snv.txt")
    d <- fread(x)
    d <- d[, 4:14]
    names(d) <- c("Project", "Sample", "ID", "Genome", "mut_type", "chrom", "pos_start", "pos_end", "ref", "alt", "Type")
    # Remove any duplications
    d1 <- d[!duplicated(d)]
    return(d1)
}

save_converted <- function(intertype = "norloop_exon") {
    # save converted
    input <- paste0("input/_intersected_snv_beds/", project, "/", intertype, "/intersected.bed")
    output_folder <- paste0("analysis/mutation_matrix/", project, "/", intertype)
    mkdirp(output_folder)
    output <- paste0("analysis/mutation_matrix/", project, "/", intertype, "/intersected.txt")
    convert_bed_intersected(input) %>%
        fwrite(output, sep = "\t")
}

# Constructing combination of rloop, norloop types
#convert_types <- expand.grid(
    ## c("norloop_"),
    ## c("intergenic", "genebody")
    #c("rloop_", "norloop_"),
    #c("intergenic", "genebody", "tss", "tts")
#)
#convert_types$all <- paste0(convert_types$Var1, convert_types$Var2)

# Convert all the intersected into simple snv text file
#lapply(convert_types$all, save_converted)
save_converted(region_type)

