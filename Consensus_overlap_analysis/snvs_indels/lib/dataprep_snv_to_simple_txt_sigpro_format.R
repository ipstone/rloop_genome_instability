# Convert the SNV files from ICGC input to sigPorfileMatrixGenerator input text format
rm(list = ls())
library(data.table)
library(magrittr)
library(ipfun)
library(dplyr)

# Subset/save ovarian simple SNV data
save_ovarian_icgc_simple_maf <- function() {
    # Loading from the source of ICGC OV tumors' WGS SNV file
    d <- fread("/data/share/icgc/pcwag/consensus_snv_indel/final_consensus_passonly.snv_mnv_indel.icgc.public.maf")
    od <- d[Project_Code == "Ovary-AdenoCA"]
    fwrite(od, "/data/projects/peix/rloop_project/input/icgc_ovary_final_consensus_passonly_public.maf", sep = "\t")
    # -- it turns out the format is not correct
}

save_ovarian_snv_simple_text_format <- function() {
    # Convert the snv to the correct simple text for]mat
    d <- fread("input/icgc_ovary_final_consensus_passonly_public.maf")
    d$Project <- d$Project_Code
    d$Sample <- d$Donor_ID
    d$ID <- d$Hugo_Symbol # We don't have an ID in the input, just use Hugo_Symbol in this place
    d$Genome <- "GRCh37"
    d$mut_type <- d$Variant_Type
    d$chrom <- d$Chromosome
    d$pos_start <- d$Start_position
    d$pos_end <- d$End_position
    d$ref <- d$Reference_Allele
    d$alt <- d$Tumor_Seq_Allele2
    d$Type <- "SOMATIC"
    d$vaf <- d$t_alt_count / (d$t_alt_count + d$t_ref_count)

    # Save the renamed d as d1 with selected fields
    # d1 <- d[, .(Project, Sample, ID, Genome, mut_type, chrom, pos_start, pos_end, ref, alt, Type)]
    # fwrite(d1, "analysis/icgc_ov_analysis/icgc_ov_simple_snv.txt", sep = "\t")
    # d.t = fread("analysis/icgc_ov_analysis/icgc_ov_simple_snv.txt") # -- testing data by loading it

    d2 <- d[vaf >= 0.26] # Looking at the high VAF SNVs
    d3 <- d2[, .(Project, Sample, ID, Genome, mut_type, chrom, pos_start, pos_end, ref, alt, Type)]
    fwrite(d3, "input/ICGC_OV_highVAF_26/icgc_ov_simple_snv_highVAF_26.txt", sep = "\t")
}

convert_icgc_maf_simple_mut_textfile <- function(maf, simpf) {
    # Convert the icgc open maf file to the correct simple text for]mat
    # d <- fread("input/icgc_ovary_final_consensus_passonly_public.maf")
    d <- fread(maf)
    d$Project <- d$Project_Code
    d$Sample <- d$Donor_ID
    d$ID <- d$Hugo_Symbol # We don't have an ID in the input, just use Hugo_Symbol in this place
    d$Genome <- "GRCh37"
    d$mut_type <- d$Variant_Type
    d$chrom <- d$Chromosome
    d$pos_start <- d$Start_position
    d$pos_end <- d$End_position
    d$ref <- d$Reference_Allele
    d$alt <- d$Tumor_Seq_Allele2
    d$Type <- "SOMATIC"

    # Save the renamed d as d1 with selected fields
    d1 <- d[, .(Project, Sample, ID, Genome, mut_type, chrom, pos_start, pos_end, ref, alt, Type)]
    fwrite(d1, simpf, sep = "\t")
    # d2 = fread("analysis/icgc_ov_analysis/icgc_ov_simple_snv.txt") # -- testing data by loading it
}

save_other_tumor_simple_mutation_file <- function() {
    # Save the prosate , pancreatic, liver and icgc-breast
    #mkdirp("input/ICGC_Prostate")
    #mkdirp("input/ICGC_Pancreatic")
    #mkdirp("input/ICGC_Liver")
    #mkdirp("input/ICGC_Breast")
    mkdirp("input/ICGC_CNS-Medullo")

    #convert_icgc_maf_simple_mut_textfile(
        #"input/ICGC_Prostate_consensus_passonly.maf",
        #"input/ICGC_Prostate/ICGC_Prostate_simple_snv.txt"
    #)

    #convert_icgc_maf_simple_mut_textfile(
        #"input/ICGC_Pancreatic_adeno_consensus_passonly.maf",
        #"input/ICGC_Pancreatic/ICGC_Pancreatic_simple_snv.txt"
    #)

    #convert_icgc_maf_simple_mut_textfile(
        #"input/ICGC_Liver_consensus_passonly.maf",
        #"input/ICGC_Liver/ICGC_Liver_simple_snv.txt"
    #)

    #convert_icgc_maf_simple_mut_textfile(
        #"input/ICGC_Breast_consensus_passonly.maf",
        #"input/ICGC_Breast/ICGC_Breast_simple_snv.txt"
    #)

    convert_icgc_maf_simple_mut_textfile(
        "input/ICGC_CNS-Medullo_consensus_passonly.maf",
        "input/ICGC_CNS-Medullo/ICGC_CNS-Medullo_simple_snv.txt"
    )
}

# Function to convert the snv files to bed files for intersection
convert_snv_bedfile_for_intersection <- function(x) {
    # x is input file
    # y is output bedfile
    d <- fread(x)
    # d <- fread("input/ICGC_OV/icgc_ov_simple_snv.txt")
    d.names <- names(d)
    d$V1 <- d$chrom
    d$V2 <- d$pos_start
    d$V3 <- d$pos_end
    new.names <- c("V1", "V2", "V3", d.names)
    d1 <- d[, ..new.names]
    return(d1)
}

convert_icgc_tumors_beds = function(){
    ## Here is a collection of previously converting different ICGC
    ## tumors code

    # convert_snv_bedfile_for_intersection("input/ICGC_OV/icgc_ov_simple_snv.txt") %>%
    #     fwrite("input/ICGC_OV/icgc_ov_simple_snv.txt.bed", col.names = FALSE, sep = "\t")


    ## Convert all Serena's SNVs to bed format for intersection
    # convert_snv_bedfile_for_intersection("input/BRCA-EU/serena_erpos_simple_snv_no_duplication.txt") %>%
    #     fwrite("input/Serena_ER/icgc_ov_simple_snv.txt.bed", col.names = FALSE, sep = "\t")

    ## Convert all Serena's SNVs to bed format for intersection
    # convert_snv_bedfile_for_intersection("input/BRCA-EU/serena_erpos_simple_snv_multiple_transcripts.txt") %>%
    #     fwrite("input/Serena_ER/icgc_ov_simple_snv.txt.bed", col.names = FALSE, sep = "\t")

    ## Testing biallelic samples in Serena's ER biallelic samples
    # convert_snv_bedfile_for_intersection("input/bed-version-2021-05-08/Serena-Davies_ER_Biallelic/davies_ER_Biallelic_simple_mutation.txt") %>%
    #     fwrite("input/Serena-Davies_ER_Biallelic/icgc_ov_simple_snv.txt.bed", col.names = FALSE, sep = "\t")

    #     ## Testing ICGC-Breast
    #     mkdirp("input/ICGC_Breast")
    #     convert_snv_bedfile_for_intersection("input/bed-version-2021-05-08/ICGC_Breast/ICGC_Breast_simple_snv.txt") %>%
    #         fwrite("input/ICGC_Breast/icgc_ov_simple_snv.txt.bed", col.names = FALSE, sep = "\t")

    #     ## Testing ICGC-Liver
    #     mkdirp("input/ICGC_Liver")
    #     convert_snv_bedfile_for_intersection("input/bed-version-2021-05-08/ICGC_Liver/ICGC_Liver_simple_snv.txt") %>%
    #         fwrite("input/ICGC_Liver/icgc_ov_simple_snv.txt.bed", col.names = FALSE, sep = "\t")

    #     ## Testing ICGC-Pancreatic
    #     mkdirp("input/ICGC_Pancreatic")
    #     convert_snv_bedfile_for_intersection("input/bed-version-2021-05-08/ICGC_Pancreatic/ICGC_Pancreatic_simple_snv.txt") %>%
    #         fwrite("input/ICGC_Pancreatic/icgc_ov_simple_snv.txt.bed", col.names = FALSE, sep = "\t")

    #     ## Testing ICGC-Prostate
    #     mkdirp("input/ICGC_Prostate")
    #     convert_snv_bedfile_for_intersection("input/bed-version-2021-05-08/ICGC_Prostate/ICGC_Prostate_simple_snv.txt") %>%
    #         fwrite("input/ICGC_Prostate/icgc_ov_simple_snv.txt.bed", col.names = FALSE, sep = "\t")

    ## Testing Serena_ER
    # mkdirp("input/Serena_ER")
    # convert_snv_bedfile_for_intersection("input/bed-version-2021-05-08/serena_erpos_simple_snv.txt") %>%
    #     fwrite("input/Serena_ER/icgc_ov_simple_snv.txt.bed", col.names = FALSE, sep = "\t")

    ## Save high vaf icgc OV snvs
    #save_ovarian_snv_simple_text_format()
    #convert_snv_bedfile_for_intersection("input/ICGC_OV_highVAF_26/icgc_ov_simple_snv_highVAF_26.txt") %>%
    #fwrite("input/ICGC_OV_highVAF_26/icgc_ov_simple_snv.txt.bed", col.names = FALSE, sep = "\t")

    ## Save CNS-Medullo
    #save_other_tumor_simple_mutation_file()
    #convert_snv_bedfile_for_intersection("input/ICGC_CNS-Medullo/ICGC_CNS-Medullo_simple_snv.txt") %>%
    #fwrite("input/ICGC_CNS-Medullo/icgc_ov_simple_snv.txt.bed", col.names = FALSE, sep = "\t")

}

simple_icgc_tumor_save = function(tumor.name, maf.file) {
    # Combine several function together to convert icgc maf data to snv text
    # format and bed format for intersection and further calculation down the
    # pipeline
    folder = paste0("input/ICGC_", tumor.name, "/")
    snv.file = paste0(folder, "icgc_", tumor.name, "_simple_snv.txt")
    bed.file = paste0(folder, "icgc_ov_simple_snv.txt.bed")
    mkdirp(folder)
    convert_icgc_maf_simple_mut_textfile(maf.file, snv.file)
    convert_snv_bedfile_for_intersection(snv.file) %>%
        fwrite(bed.file, col.names= FALSE, sep="\t")
}


save_serena_snv_simple_text_format <- function() {
    # Prepare Serena's SNV dataset for analysis
    d <- fread("input/bed-version-2021-05-08/BRCA-EU/simple_somatic_mutation.open.BRCA-EU.tsv")
    d$sampleid <- gsub("(a|b).*$", "", d$submitted_sample_id) # There are 569 unique sample_id, after removing a|b, it's still 569 unique sample id

    # Loading Serena breast tumor type info
    sample <- fread("input/bed-version-2021-05-08/BRCA-EU/serena_nature_2016_560_breast_tumor_table1_clindata.csv")
    sample$sampleid <- gsub("(a|b).*$", "", sample$sample_name)
    # This is the cohort we are using for final testing 320 ER+/HER2-
    er_pos_sample <- sample[final.ER == "positive" & final.HER2 == "negative"]
    erPosSamples <- unique(er_pos_sample$sampleid)
    trp_neg_sample <- sample[final.ER == "negative" &
        final.HER2 == "negative" &
        final.PR == "negative"]
    trpNegSamples <- unique(trp_neg_sample$sampleid)

    # Recode the mutation file to the simple text format
    d$Project <- d$project_code
    # d$Sample <- d$icgc_donor_id
    d$Sample <- d$sampleid # For Serena's dataset, better use sampleid as we will use it to differentiate between biallelic and control cases
    d$ID <- d$gene_affected # We don't have an ID in the input, just use Hugo_Symbol in this place
    d$Genome <- "GRCh37"
    d$mut_type <- recode(d$mutation_type,
        "deletion of <=200bp" = "DEL",
        "insertion of <=200bp" = "INS",
        "multiple base substitution (>=2bp and <=200bp)" = "DNP",
        "single base substitution" = "SNP"
    )

    d$chrom <- d$chromosome
    d$pos_start <- d$chromosome_start
    d$pos_end <- d$chromosome_end
    d$ref <- d$reference_genome_allele
    d$alt <- d$mutated_to_allele
    d$Type <- "SOMATIC"

    # Processing duplicated samples from
    d$label <- paste(d$sampleid, d$mut_type, d$chrom, d$pos_start, d$pos_end, d$ref, d$alt, sep = "-")
    d$label2 <- paste(d$label, d$transcript_affected) # Try to identify SNVs with mulitple transcript

    # One way is to only include verified sequencing result
    d2 <- d[verification_status == "tested and verified"]
    d2 <- d2[!duplicated(d2$label)]

    # Another way is to include variants with multiple tanscripts
    d3 <- d[!duplicated(d$label2)] # all unique tarnscript
    d4 <- d3[duplicated(d3$label)] # Only include those SNVs with duplicated position but with mutliple tanscripts
    d4 <- d4[!duplicated(label)] # This is needed to remove cases with multiple entries at unique positions

    # Start to subset mutation data and save ER, tripNeg simple mutation file.
    d.er <- d2[sampleid %in% erPosSamples, .(Project, Sample, ID, Genome, mut_type, chrom, pos_start, pos_end, ref, alt, Type)]
    fwrite(d.er, "input/BRCA-EU/serena_erpos_simple_snv_no_duplication.txt", sep = "\t")
    d.tri <- d2[sampleid %in% trpNegSamples, .(Project, Sample, ID, Genome, mut_type, chrom, pos_start, pos_end, ref, alt, Type)]
    fwrite(d.tri, "input/BRCA-EU/serena_tripneg_simple_snv_no_duplication.txt", sep = "\t")

    # Start to subset mutation data and save ER, tripNeg simple mutation file.
    # d.er <- d4[sampleid %in% erPosSamples, .(Project, Sample, ID, Genome, mut_type, chrom, pos_start, pos_end, ref, alt, Type)]
    # fwrite(d.er, "input/BRCA-EU/serena_erpos_simple_snv_multiple_transcripts.txt", sep = "\t")
    # d.tri <- d4[sampleid %in% trpNegSamples, .(Project, Sample, ID, Genome, mut_type, chrom, pos_start, pos_end, ref, alt, Type)]
    # fwrite(d.tri, "input/BRCA-EU/serena_tripneg_simple_snv_multiple_transcripts.txt", sep = "\t")
}


#######################################

save_TCGA_mutation_input <- function () {
    ## Coding mc3 tumor types
    #setwd("..")
    mc3 = fread("/data/share/tcga/mc3.v0.2.8.PUBLIC.maf")
    mc3$sample=substr(mc3$Tumor_Sample_Barcode, 1, 12)

    # Loading ER sample list
    sl = readRDS("input/biallelic_sample_lists/biallelic_tcga_sample_lists_biallelic-cn-project.rds")

    # Loading TCGA tumor type sample list
    tcga = fread("input/TCGA-pan-cancer-survival_2018_Cell_data.csv")
    tcga$sample = tcga$bcr_patient_barcode
    table(tcga$type)


    # When coding the fields into SigProfilerMutation matrix gen input:
    # -- we can first focus on making the simple text format for the input 
    # -- then using the convert_snv_bedfile_for_intersection function to save the
    # format to the bed file format needed for the current pipeline

    convert_TCGA_maf_simple_mut_textfile <- function(maf.dt, 
                                                     project_code="TCGA",
                                                     simpf) {
        # Using the input maf data.table object
        d <- maf.dt
        d$Project <- project_code
        d$Sample <- d$sample
        d$ID <- d$Hugo_Symbol # We don't have an ID in the input, just use Hugo_Symbol in this place
        d$Genome <- "GRCh37"
        d$mut_type <- d$Variant_Type
        d$chrom <- d$Chromosome
        d$pos_start <- d$Start_Position
        d$pos_end <- d$End_Position
        d$ref <- d$Reference_Allele
        d$alt <- d$Tumor_Seq_Allele2
        d$Type <- "SOMATIC"

        # Save the renamed d as d1 with selected fields
        d1 <- d[, .(Project, Sample, ID, Genome, mut_type, chrom, pos_start, pos_end, ref, alt, Type)]
        fwrite(d1, simpf, sep = "\t")
    }

    convert_mc3_tumorData_snvTxt_bed = function(dt, project.code="TCGA", snv_file, bed_file){
        # With the input MC3 data.table subset, save the snv text format and
        # bed format for the pipeline to run intersection and signature matrix
        # generation and extractions
        er.dt =dt
        convert_TCGA_maf_simple_mut_textfile(er.dt, project_code=project.code, snv_file)
        convert_snv_bedfile_for_intersection(snv_file) %>%
            fwrite(bed_file, col.names = FALSE, sep = "\t")

        return(dt)
    }

    # Specify / narrow down to interested samples
    ov = tcga[ type == "OV" ]
    liver = tcga[ type =="LIHC"]
    er = tcga[ sample %in% sl$er_positive ]

    ov.dt = mc3[ sample %in% ov$sample ]
    liver.dt = mc3[ sample %in% liver$sample ]
    er.dt = mc3[ sample %in% er$sample ]

    #mkdirp("input/TCGA_OV/")
    #mkdirp("input/TCGA_Liver/")
    #mkdirp("input/TCGA_ER/")

    #convert_mc3_tumorData_snvTxt_bed(ov.dt, project.code="TCGA_OV",
                                    #"input/TCGA_OV/icgc_ov_simple_snv.txt",
                                    #"input/TCGA_OV/icgc_ov_simple_snv.txt.bed")
    #convert_mc3_tumorData_snvTxt_bed(liver.dt, project.code="TCGA_Liver",
                                    #"input/TCGA_Liver/icgc_ov_simple_snv.txt",
                                    #"input/TCGA_Liver/icgc_ov_simple_snv.txt.bed")
    #convert_mc3_tumorData_snvTxt_bed(er.dt, project.code="TCGA_ER",
                                    #"input/TCGA_ER/icgc_ov_simple_snv.txt",
                                    #"input/TCGA_ER/icgc_ov_simple_snv.txt.bed")

    # Additional function to make call each tumor type easier
    simple_save = function(tumor) {
        s = tcga[ type == tumor ]
        f.dt = mc3[ sample %in% s$sample ]
        project = paste0("TCGA_", tumor)
        folder = paste0("input/TCGA_", tumor, "/")
        mkdirp(folder)
        snv = paste0(folder, "icgc_ov_simple_snv.txt")
        bed = paste0(folder, "icgc_ov_simple_snv.txt.bed")
        convert_mc3_tumorData_snvTxt_bed(f.dt, project.code=project, snv, bed)
    }

    # Save additional tumor types data
    #simple_save("PRAD")
    #simple_save("LUAD")
    #simple_save("LUSC")
    #simple_save("SKCM")
    simple_save("PAAD")
}

## main execution section
main <- function() {
    # save_cleaned_bed_files()
    # save_serena_snv_simple_text_format()
    # save_other_tumor_simple_mutation_file()

    ## Save Serena's SNV using verified and remove duplications
    # save_serena_snv_simple_text_format()

    save_TCGA_mutation_input()
    #simple_icgc_tumor_save("Melanoma",
                           #"input/ICGC_Melanoma_consensus_passonly.maf")

}

main()
