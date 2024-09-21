# Extract all signatures through SigProfilerExtractor
gen_sigs = lib/run_sigProfExtact_SBS96.py
archive_place = _archive/analysis/extracted_signatures_2021-05-10

.PHONY: list gen_sigs clean

list:
	@echo "This makefile is to run signature extraction using sigProfilerExtractor"
	@$(MAKE) -pRrq -f $(lastword $(MAKEFILE_LIST)) : 2>/dev/null | awk -v RS= -F: '/^# File/,/^# Finished Make data base/ {if ($$1 !~ "^[#.]") {print $$1}}' | sort | egrep -v -e '^[^[:alnum:]]' -e '^$@$$'

gen_sigs:
	python $(gen_sigs) analysis/icgc_ov_analysis
	python $(gen_sigs) analysis/ICGC_OV/norloop_exon
	python $(gen_sigs) analysis/ICGC_OV/rloop_exon
	python $(gen_sigs) analysis/ICGC_OV/norloop_genebody
	python $(gen_sigs) analysis/ICGC_OV/norloop_tss
	python $(gen_sigs) analysis/ICGC_OV/norloop_tts
	python $(gen_sigs) analysis/ICGC_OV/rloop_genebody
	python $(gen_sigs) analysis/ICGC_OV/rloop_tss
	python $(gen_sigs) analysis/ICGC_OV/rloop_tts

clean:
	mkdir -p $(archive_place)
	###################################################################
	## This is not the best approach : need better cleaning method
	##################################################################
	mv -b analysis/icgc_ov_analysis/extracted_signatures $(archive_place)/icgc_ov_analysis
	mv -b analysis/ICGC_OV/norloop_exon/extracted_signatures $(archive_place)/norloop_exon
	mv -b analysis/ICGC_OV/norloop_genebody/extracted_signatures $(archive_place)/norloop_genebody
	mv -b analysis/ICGC_OV/norloop_tss/extracted_signatures  $(archive_place)/norloop_tss
	mv -b analysis/ICGC_OV/norloop_tts/extracted_signatures $(archive_place)/norloop_tts
	mv -b analysis/ICGC_OV/rloop_exon/extracted_signatures $(archive_place)/rloop_exon
	mv -b analysis/ICGC_OV/rloop_genebody/extracted_signatures $(archive_place)/rloop_genebody
	mv -b analysis/ICGC_OV/rloop_tss/extracted_signatures $(archive_place)/rloop_tss
	mv -b analysis/ICGC_OV/rloop_tts/extracted_signatures $(archive_place)/rloop_tts
