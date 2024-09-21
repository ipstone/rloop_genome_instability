# Generate all mutation matrix through sigProfiler Matrix Gen for all type of
# tumors
# -- using 
tumors = [
    "ICGC_Breast", "ICGC_Liver", 
    "ICGC_OV", "ICGC_CNS-Medullo", 
    "TCGA_OV", "TCGA_ER", "TCGA_Liver",
    "ICGC_Pancreatic", "ICGC_Prostate", 
    "Serena_ER", "Serena-Davies_ER_Biallelic",
    "TCGA_PRAD", "TCGA_LUAD", "TCGA_LUSC", "TCGA_SKCM",
    "ICGC_Melanoma", "TCGA_PAAD"
]

rloop_type = ["rloop_", "norloop_"]
region_type =['tss', 'tts', 'tss_transcribed', 'tts_transcribed',
              'tss_exome', 'tts_exome', 'tss_transcribed_exome', 'tts_transcribed_exome',
              'consensus_sc200']
#region_type =['tss', 'tts', 'tss_transcribed', 'tts_transcribed',
              #'tss_exome', 'tts_exome', 'tss_transcribed_exome', 'tts_transcribed_exome']
#region_type =['tss', 'tts', 
              #'genebody', 'pseudogene', 'lincRNA']

# Expand to obtain all rloop region types
regions = expand("{rloop}{rtype}", rloop=rloop_type, rtype=region_type)

#region_type =['tss', 'tts'] 
#gene_types =['_gene', '_pseudogene', '_lincRNA']
#regions = expand("{rloop}{rtype}{gene_type}", 
                #rloop=rloop_type, 
                #rtype=region_type, 
                #gene_type=gene_types)


# Rules to generate all different SBS/indel mutation count
rule all:
    input:
        expand("analysis/mutation_matrix/{tumor}/{region}/output/SBS/icgc_ov.SBS18.all", tumor=tumors, region=regions)
        #expand("analysis/mutation_matrix/{tumor}/{region}/output/ID/icgc_ov.ID28.all", tumor=tumors, region=regions)

rule gen_mutmatrix:
    input: 
        "input/_intersected_snv_beds/{tumor}/{region}/intersected.bed"
    output: 
        "analysis/mutation_matrix/{tumor}/{region}/output/SBS/icgc_ov.SBS18.all"
        #"analysis/mutation_matrix/{tumor}/{region}/output/ID/icgc_ov.ID28.all"
    run:
        shell("mkdir -p analysis/mutation_matrix/{wildcards.tumor}/{wildcards.region}")
        shell("Rscript lib/convert_bed_back_simple_snv_txt.R {wildcards.tumor} {wildcards.region}")
        shell("python lib/run_sigProfileMatrixGen_datafolder.py analysis/mutation_matrix/{wildcards.tumor}/{wildcards.region}")

# Generate all mutation matrix through sigProfiler Matrix Gen for all type of
# tumors

#rule intersect_bed:
    #input:
        #expand("input/{tumor}/{region}/intersected.bed", tumor=tumors, region=regions)

rule intersect_snv:
    input: 
        "input/{tumor}/icgc_ov_simple_snv.txt.bed"
    output: 
        "input/_intersected_snv_beds/{tumor}/{region}/intersected.bed"
    run:
        shell("mkdir -p input/_intersected_snv_beds/{wildcards.tumor}/{wildcards.region}")
        shell("bedtools intersect -wa -a {input} -b input/rloop_bed_files/{wildcards.region}.bed > {output} ")

