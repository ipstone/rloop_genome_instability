# Extract all signatures through SigProfilerExtractor
# -- using 
tumors = [
    "ICGC_Breast", "ICGC_Liver", 
    "ICGC_OV", "ICGC_CNS-Medullo", 
    "TCGA_OV", "TCGA_ER", "TCGA_Liver",
    "ICGC_Pancreatic", "ICGC_Prostate", 
    "Serena_ER", "Serena-Davies_ER_Biallelic",
    "TCGA_PRAD", "TCGA_LUAD", "TCGA_LUSC", "TCGA_SKCM",
    "ICGC_Melanoma", "TCGA_PAAD"
    # "Serena-Davies_ER_Control", 
    # "Serena-Davies_TRP_Biallelic",
    # "Serena-Davies_TRP_Control", 
    # "Serena_ER",
    # "Serena_TRP"
    # "ICGC_OV_Biallelic", "ICGC_OV_Control", 
]

# tumors = [ "ICGC_OV"]
# tumors = [ "Serena_ER"]

# analysis_folder = "analysis/Serena_TRP/"
# analysis_folder = "analysis/Serena_ER/"
# result_file = "/extracted_signatures/SBS96/SBS96_selection_plot.pdf"

all_results = expand("output/SBS_compare_rloop_regions/{tumor}/", tumor=tumors)

rule all:
    input:
        directory(all_results)

rule gen_sigs:
    output:
        directory("output/SBS_compare_rloop_regions/{tumor}/")

    shell:
        # "mkdir -p output/SBS_compare_rloop_regions/{wildcards.tumor}"
        "Rscript lib/analyze_SBS_freq_diff_rloops.R {wildcards.tumor}"

