# dataprep for SV to be converted to bedpe format for intersection
rm(list = ls())
library(data.table)
library(magrittr)
library(tidyverse)
library(ipfun)

get_serena_samples = function() {
    # Loading Serena breast tumor type info
    sample <- fread("../../input/bed-version-2021-05-08/BRCA-EU/serena_nature_2016_560_breast_tumor_table1_clindata.csv")
    sample$sampleid <- gsub("(a|b).*$", "", sample$sample_name)

    # This is the cohort we are using for final testing 320 ER+/HER2-
    er_pos_sample <- sample[final.ER == "positive" & final.HER2 == "negative"]
    erPosSamples <- unique(er_pos_sample$sampleid)
    trp_neg_sample <- sample[final.ER == "negative" &
        final.HER2 == "negative" &
        final.PR == "negative"]
    trpNegSamples <- unique(trp_neg_sample$sampleid)

    rl = list()
    rl$erPosSamples = erPosSamples 
    rl$trpNegSamples = trpNegSamples
    return(rl)
}

save_serena_sv_bedpe_format <- function() {
    # Prepare Serena's SNV dataset for analysis
    d <- fread("input/BRCA-EU_SV/structural_somatic_mutation.BRCA-EU.tsv")
    d$sampleid <- gsub("(a|b).*$", "", d$submitted_sample_id) # There are 569 unique sample_id, after removing a|b, it's still 569 unique sample id
    table(d$chr_from)
    table(d$chr_to)

    serena= get_serena_samples()

    # Recode the SV file to the bedpe format styles - try to fulfill all bedpe
    # definitions
    d$Sample <- d$sampleid # For Serena's dataset, better use sampleid as we will use it to differentiate between biallelic and control cases
    d$svclass <- recode(d$variant_type,
        "deletion" = "DEL",
        "inversion" = "INV",
        "tandem duplication" = "DUP",
        "interchromosomal rearrangement - unknown type" = "TRA"
    )

    d$chrom1 <- d$chr_from
    d$start1 <- d$chr_from_bkpt -1
    d$end1 <- d$chr_from_bkpt 

    d$chrom2 <- d$chr_to
    d$start2 <- d$chr_to_bkpt -1
    d$end2 <- d$chr_to_bkpt 

    #d$sv_id   # sv_id is already given
    d$pe_support = 0 # there's no supporting SV breakpoints counts in this data 
    d$strand1 = ifelse(d$chr_from_strand == 1, "+", "-")
    d$strand2 = ifelse(d$chr_from_strand == 2, "+", "-")
    d$svmethod = "BRASS"

    #table(duplicated(d$sv_id))
    # -- no duplicated SV calls are in this dataset


    # Start to subset mutation data and save ER, tripNeg simple mutation file.
    d.er <- d[sampleid %in% serena$erPosSamples, .(chrom1, start1, end1, chrom2, start2, end2, sv_id, pe_support, strand1, strand2, svclass,svmethod, Sample)]
    fwrite(d.er, "input/BRCA-EU_SV/serena_erpos_SV.bedpe", sep = "\t", col.names=F)

    d.tri <- d[sampleid %in% serena$trpNegSamples, .(chrom1, start1, end1, chrom2, start2, end2, sv_id, pe_support, strand1, strand2, svclass,svmethod, Sample)]
    fwrite(d.tri, "input/BRCA-EU_SV/serena_tripneg_SV.bedpe", sep = "\t", col.names=F)

}
#save_serena_sv_bedpe_format()

test_serena_BRCA_EU_samples <- function () {
    # Test whether all Serena's samples are within the BRCA-EU cases 
    # -- tested, all the 320 samples are within the BRCA-EU cases
    # -- also tesssted, all 163 triple negative breast tumors are within the
    # list
    serena = get_serena_samples()
    eu = fread("input/BRCA-EU_SV/sample.BRCA-EU.tsv")
    eu$sampleid <- gsub("(a|b).*$", "", eu$submitted_sample_id)
    #return(table(serena$erPosSamples %in% eu$sampleid))
    return(table(serena$trpNegSamples %in% eu$sampleid))
    
}
