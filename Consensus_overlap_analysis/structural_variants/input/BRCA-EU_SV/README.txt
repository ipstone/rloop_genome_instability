The input structural variants are from PCAWG release 28 BRCA-EU data.

- There are 544 patients included in this dataset (Serena's ER and triple
  negative breast tumors data): compared with the 1040~ patients included in
  the folder, perhaps half? of the patients don't have SV called ...?
      - that could be a safe assumption for now.

- All these variants are labeled as "not tested", are all called by BRASS
  algorithm

- the variant types are labeled as:
      deletion
      interchromosomal rearramgent (TRA)
      inversion
      tandem duplication(duplication)


* Analysis 
      We analyze with similar steps as for the ICGC SV, counting the
      overlapping breakpoints with the Rloop respective regions.

      Step 1. Oragnize the SV data to similar format/BEDPE format as the ICGC
              data.
              - There could be cases with 0 SV events, which are not included
                in the sv file.


* Data files not included in the code repo:
  These data were downloaded from ICGC data portal for BRCA-EU projects.
    sample.BRCA-EU.tsv
    structural_somatic_mutation.BRCA-EU.tsv

  The following data were the converted bedpe format from BRCA-EU SVs data, separated to ER-pos breast tumors and Triple-negative breast tumors.
  (converted by the ../../dataprep_serena_sv_bedpe.R script)

    serena_erpos_SV.bedpe
    serena_tripneg_SV.bedpe
