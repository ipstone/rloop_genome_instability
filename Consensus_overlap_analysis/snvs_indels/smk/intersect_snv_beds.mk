# project = "ICGC_OV"
# project = Serena_ER
# project = Serena-Davies_ER_Biallelic
# project = ICGC_Breast
# project = ICGC_Liver
# project = ICGC_Pancreatic
# project = ICGC_Prostate
#project = "ICGC_OV_highVAF_26"
project = "ICGC_CNS-Medullo"

## To use bedtools to intersect snv files with different bed files.
snv = "input/$(project)/icgc_ov_simple_snv.txt.bed"

.PHONY: all clean create_folders

# all: input/ICGC_OV/norloop_intergenic/intersected.bed \
# 	input/ICGC_OV/norloop_genebody/intersected.bed 

all: create_folders input/$(project)/rloop_intergenic/intersected.bed \
	input/$(project)/rloop_genebody/intersected.bed \
	input/$(project)/rloop_tss/intersected.bed \
	input/$(project)/rloop_tts/intersected.bed \
	input/$(project)/norloop_intergenic/intersected.bed \
	input/$(project)/norloop_genebody/intersected.bed \
	input/$(project)/norloop_tss/intersected.bed \
	input/$(project)/norloop_tts/intersected.bed

create_folders: 
	mkdir -p input/$(project)/rloop_intergenic/
	mkdir -p input/$(project)/rloop_genebody/
	mkdir -p input/$(project)/rloop_tss/
	mkdir -p input/$(project)/rloop_tts/
	mkdir -p input/$(project)/norloop_intergenic/
	mkdir -p input/$(project)/norloop_genebody/
	mkdir -p input/$(project)/norloop_tss/
	mkdir -p input/$(project)/norloop_tts/



input/$(project)/rloop_intergenic/intersected.bed:
	bedtools intersect -wa -a $(snv) -b input/rloop_bed_files/rloop_intergenic.bed > $@

input/$(project)/rloop_genebody/intersected.bed:
	bedtools intersect -wa -a $(snv) -b input/rloop_bed_files/rloop_genebody.bed > $@

input/$(project)/rloop_tss/intersected.bed:
	bedtools intersect -wa -a $(snv) -b input/rloop_bed_files/rloop_tss.bed > $@

input/$(project)/rloop_tts/intersected.bed:
	bedtools intersect -wa -a $(snv) -b input/rloop_bed_files/rloop_tts.bed > $@

input/$(project)/norloop_intergenic/intersected.bed:
	bedtools intersect -wa -a $(snv) -b input/rloop_bed_files/norloop_intergenic.bed > $@

input/$(project)/norloop_genebody/intersected.bed:
	bedtools intersect -wa -a $(snv) -b input/rloop_bed_files/norloop_genebody.bed > $@

input/$(project)/norloop_tss/intersected.bed:
	bedtools intersect -wa -a $(snv) -b input/rloop_bed_files/norloop_tss.bed > $@

input/$(project)/norloop_tts/intersected.bed:
	bedtools intersect -wa -a $(snv) -b input/rloop_bed_files/norloop_tts.bed > $@

clean: 
	rm input/$(project)/norloop_intergenic/intersected.bed \
	input/$(project)/norloop_genebody/intersected.bed 
	rm input/$(project)/rloop_intergenic/intersected.bed \
	input/$(project)/rloop_genebody/intersected.bed \
	input/$(project)/rloop_tss/intersected.bed \
	input/$(project)/rloop_tts/intersected.bed \
	input/$(project)/norloop_tss/intersected.bed \
	input/$(project)/norloop_tts/intersected.bed 
