Analyzing structural variants / breakpoints with the Rloop regions vs. non-Rloop regions


* scripts to intersect the SVs with the different bed files 
    intersect_consensus.sh
    intersect_no_consensus.sh
    intersect_rloop_regions_sc100.sh
    intersect_rloop_regions_sc200.sh
    intersect_sc100.sh
    intersect_sc200.sh
    intersect_sereana_rloops_all_regions.sh

    Some other helper scripts to manage the data:
        combine_all_raw_sv_result.py
        count_intersected_region_sv.sh

* Intersection approach : use bedtools pairtobed -a bedpefile -b bed 
      default option is to report the bedpe file intersect with bed file at
      either ends.
             bedtools pairtobed -a example.bedpe -b consensus_no_cutoff.bed > result_pairtobed_all_consensus.txt


      We can also use -type=neither to find that the bedpe file which didn't
      intersect with the bedpe files at all.
             sed 's/^...//' Consensus_sc_gt_100.txt_length_lt_5_num_gt_5.bed_SORTED.bed > consensus_sc_gt_100.bed
             bedtools pairtobed -a example.bedpe -b consensus_no_cutoff.bed -type neither > result-pairtobed-neither-consensus.txt

* input
  This folder holds the different type of inputs used in this analysis. Some inputs are omitted here due to size limitation of the code repo. 

* R
  This folder contains R code to analyze the intersected SV breakpoints with genomics regions with or without R-loops.
  These analysis codes are written in R, should be run at this folder level.

