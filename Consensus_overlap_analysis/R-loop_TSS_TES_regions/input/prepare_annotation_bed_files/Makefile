# Make file to parse out the downloaded genecode.v19 annotation gtf file.
#
gencode-v19-extracted.txt: gencode.v19.annotation.gtf-no-header
	python extract_gencode_fields.py gencode.v19.annotation.gtf-no-header > $@ 

gencode.v19.annotation.gtf-no-header: gencode.v19.annotation.gtf
	grep -v '###' $< > $@


