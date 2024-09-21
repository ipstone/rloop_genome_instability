# Prepare mutation dataset for biallelic/control samples for intersted tumor types

rm(list = ls())
library(data.table)
library(magrittr)
library(ipfun)

# Loading ER+ Serena dataset
save_Serena_biallelic <- function() {
    davies <- readRDS("input/biallelic_sample_lists/Davies_Serena_biallelic_samples.rds")

    # Save biallelic and control samples
    er <- fread("input/serena_erpos_simple_snv.txt")
    er.biallelic <- er[Sample %in% davies$biallelic]
    er.control <- er[Sample %in% davies$control]
    fwrite(er.biallelic, "input/Serena-Davies_ER_Biallelic/davies_ER_Biallelic_simple_mutation.txt", sep = "\t")
    fwrite(er.control, "input/Serena-Davies_ER_Control/davies_ER_Control_simple_mutation.txt", sep = "\t")


    tri <- fread("input/serena_tripneg_simple_snv.txt")
    tri.biallelic <- tri[Sample %in% davies$biallelic]
    tri.control <- tri[Sample %in% davies$control]
    fwrite(tri.biallelic, "input/Serena-Davies_TRP_Biallelic/davies_TRP_Biallelic_simple_mutation.txt", sep = "\t")
    fwrite(tri.control, "input/Serena-Davies_TRP_Control/davies_TRP_Control_simple_mutation.txt", sep = "\t")
}

save_icgc_ov_biallelic <- function() {
    d <- fread("input/ICGC_OV/icgc_ov_simple_snv.txt")
    ov <- fread("input/biallelic_sample_lists/icgc_ov_biallelic_cat.csv")
    ov.biallelic <- ov[cat == "Biallelic"]
    ov.control <- ov[cat == "control"]
    d.biallelic <- d[Sample %in% ov.biallelic$icgc_donor_id]
    d.control <- d[Sample %in% ov.control$icgc_donor_id]
    fwrite(d.biallelic, "input/ICGC_OV_Biallelic/ICGC_OV_Biallelic_simple_mutation.tsv", sep = "\t")
    fwrite(d.control, "input/ICGC_OV_Control/ICGC_OV_Control_simple_mutation.tsv", sep = "\t")
}