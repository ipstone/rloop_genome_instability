# Generate all mutation matrix through sigProfiler Matrix Gen for all type of
# tumors
# -- using 
#tumors = [
    #"ICGC_Breast", "ICGC_Liver", 
    #"ICGC_OV", "ICGC_CNS-Medullo", 
    #"ICGC_Pancreatic", "ICGC_Prostate", 
    #"Serena_ER", "Serena-Davies_ER_Biallelic"
#]

#rloop_type = ["rloop_", "norloop_"]

#region_type =['tss', 'tts', 
              #'genebody', 'pseudogene', 'lincRNA']

# -- using 
tumors = [
    "ICGC_Breast", "ICGC_Liver", 
    "ICGC_OV", "ICGC_CNS-Medullo", 
    "TCGA_OV", "TCGA_ER", "TCGA_Liver",
    "ICGC_Pancreatic", "ICGC_Prostate", 
    "Serena_ER", "Serena-Davies_ER_Biallelic"
    # "Serena-Davies_ER_Control", 
    # "Serena-Davies_TRP_Biallelic",
    # "Serena-Davies_TRP_Control", 
    # "Serena_ER",
    # "Serena_TRP"
    # "ICGC_OV_Biallelic", "ICGC_OV_Control", 
]
rloop_type = ["rloop_", "norloop_"]
region_type =['tss', 'tts', 
              'genebody', 'pseudogene', 'lincRNA']

# Expand to obtain all rloop region types
regions = expand("{rloop}{rtype}", rloop=rloop_type, rtype=region_type)

rule all:
    input:
        expand("input/_archive/intersected/{tumor}/{region}/", tumor=tumors, region=regions)

rule clean_intersected_snv:
    output:
        directory("input/_archive/intersected/{tumor}/{region}/")
        
    run:
        shell("mkdir -p input/_archive/intersected/{wildcards.tumor}/{wildcards.region}")
        shell("mv input/{wildcards.tumor}/{wildcards.region}/ input/_archive/intersected/{wildcards.tumor}/{wildcards.region}")
