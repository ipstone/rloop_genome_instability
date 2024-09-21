* Overview

** Highlights:

   The 3 bed files in the output folder are the consensus region generated, which is used for downstream analysis:
   output/
        Consensus_no_cutoffs.txt.bed_SORTED.bed                            # any R-loop region, which is used to be substracted, to generate R-loop negative reigons
        Consensus_sc_gt_100.txt_length_lt_5_num_gt_5.bed_SORTED.bed        # positive R-loop region with score greater than 100
        Consensus_sc_gt_200.txt_length_lt_5_num_gt_5.bed_SORTED.bed        # positive R-loop region with score greater than 200


** Steps to calculate consensus regions:

(1)	Downloaded bigwig files from ucsc
(2)	Manually mapped the ucsc track name to the bigwig file name
(3)	Converted bigwig files into bed format
(4)	Took the max intensity score between the posiive and negatuvce
(5)	Merged replicates
(6)	Make correlation plot
(7)	Remove dodgy files
(8)	Called peaks on the bed files
(9)	Got consensus regions
(10)	Make reference dataset
(11)	Make QC plots


* Download bigwig files outlines in Steps_to_download_bigwig_files.pdf

(1) Navigate to http://genome.ucsc.edu/s/fredericchedinlab/hg19_DRIP_Correlation
(2) Navigate to R-loop Cross Check Data
(3) click on blue track name 
(4) Click on schema
(5) until you get to the page "Binary file of type bigWig stored at https://s3-us-west-1.amazonaws.com/muhucsc/HeLa_DRIPc_WT_LS61A_rep1_pos.bw"
(6) copy and paste the link in a new browser window - this will download the bigwig files

I have mapped the bigwig files to the appropriate track names - this is in the file : 

./peaks_called_from_bedgraph_files/Map_bigwig_files_to_UCSC_track_names_paper_list_for_consensus.txt_only_WT_files.txt


* Scripts to submit jobs/processing data in the ./scripts folder:

 get_max.bsub
 merge_replicates.bsub
 merge_strands.bsub
 convert_to_bed.bsub

convert_to_bed.bsub	(Run_script_to_convert_to_bed.sh)
merge_strands.bsub      (Run_script_to_merge_strands.sh)
get_max.bsub		(get_max_intensity_score_across_strands.pl)
merge_replicates.bsub	(Run_script_to_merge_replicates_after_merging_strands.sh)
merge_replicates_no_strand.bsub	(Run_script_to_merge_replicates_no_strand.sh) 
run_calculate_mean_intensity_score.bsub	(Run_calculate_mean_intensity_score.sh, calculate_mean_intensity_score.pl)

** list of how the replicate and strand files were combined is in the excel file : figure_out_replicates_only_WT_files.xls

FINAL FILES : final_23_files/ 

** call peaks

call_peaks.sub  (Run_macs2_bdgpeakcall_on_all_bed_files_after_merging_across_replicates_and_strands.sh)

###
cd peaks_files/
run the rest of thee commands within the directory peaks_files/
###


BC_DRIP_WT_GSE117671.bw_converted_to_bed.bed_bdgpeakcall_peaks.txt
CHLA10_DRIP_WT_GSE68847.bw_converted_to_bed.bed_bdgpeakcall_peaks.txt
GSM4474217_NC.bw_converted_to_bed.bed_bdgpeakcall_peaks.txt
HEK293_RDIP_WT_GSE68953_merged_strands.bed_after_getting_max_across_strands.bed_bdgpeakcall_peaks.txt
HeLa_DRIP_WT_GSE93368_mean_of_replicates.bed_after_computing_average.bed_bdgpeakcall_peaks.txt
HeLa_DRIP_WT_LS2A_mean_of_replicates.bed_after_computing_average.bed_bdgpeakcall_peaks.txt
HeLa_DRIP_WT_Pasero0520_mean_of_replicates.bed_after_computing_average.bed_bdgpeakcall_peaks.txt
HeLa_DRIPc_WT_LS61ACH_mean_of_replicates.bed_after_computing_average.bed_bdgpeakcall_peaks.txt
HeLa_RDIP_WT_GSE120371_mean_of_replicates.bed_after_computing_average.bed_bdgpeakcall_peaks.txt
HeLa_qDRIP_WT_mean_of_replicates.bed_after_computing_average.bed_bdgpeakcall_peaks.txt
IMR90_RDIP_WT_GSE68953_merged_strands.bed_after_getting_max_across_strands.bed_bdgpeakcall_peaks.txt
K562_DRIP_WT_GSE127979.bw_converted_to_bed.bed_bdgpeakcall_peaks.txt
K562_DRIP_WT_LS42C.bw_converted_to_bed.bed_bdgpeakcall_peaks.txt
K562_DRIPc_WT_GSE127979_mean_of_replicates.bed_after_computing_average.bed_bdgpeakcall_peaks.txt
K562_DRIPc_WT_LS42B_merged_strands.bed_after_getting_max_across_strands.bed_bdgpeakcall_peaks.txt
K562_DRIPc_WT_LS60AB_mean_of_replicates.bed_after_computing_average.bed_bdgpeakcall_peaks.txt
NT2_DRIP_WT_LS15A_mean_of_replicates.bed_after_computing_average.bed_bdgpeakcall_peaks.txt
NT2_DRIPc_WT_SX011C_mean_of_replicates.bed_after_computing_average.bed_bdgpeakcall_peaks.txt
TC32_DRIP_WT_GSE68847.bw_converted_to_bed.bed_bdgpeakcall_peaks.txt
U2OS_DRIP_WT_EMTAB6318.bw_converted_to_bed.bed_bdgpeakcall_peaks.txt
U2OS_DRIP_WT_GSE115957_mean_of_replicates.bed_after_computing_average.bed_bdgpeakcall_peaks.txt
U2OS_DRIP_WT_GSE129907.bw_converted_to_bed.bed_bdgpeakcall_peaks.txt
hPSC_DRIP_WT_PRJNA474076_mean_of_replicates.bed_after_computing_average.bed_bdgpeakcall_peaks.txt


** Then I use only the 19 file after rembving the 4 dodgy files to get consensus region swith no cut-off

./Run_use_bedtools_to_get_consensus_no_score_cut_off_FINAL_19_FILES.sh

output file : Consensus_no_cutoffs.txt

Now I convert this to bed file and remove regions that are more than 5kb long

perl convert_consensus_to_bed_remove_greater_5kb.pl 

output file : Consensus_no_cutoffs.txt.bed

Now sort this file so that I can run bedtoolss on it

sort -k1,1 -k2,2n Consensus_no_cutoffs.txt.bed > Consensus_no_cutoffs.txt.bed_SORTED.bed

Now find the part of the genome that contains no R loops at all

#perl make_bed_file_whole_genome.pl 

output file : hg19_chromosome_sizes.txt_converted_to_bed.bed

sort -k1,1 -k2,2n hg19_chromosome_sizes.txt_converted_to_bed.bed > hg19_chromosome_sizes.txt_converted_to_bed.bed_SORTED.bed

~/share/usr/bin/bedtools subtract -a hg19_chromosome_sizes.txt_converted_to_bed.bed_SORTED.bed -b Consensus_no_cutoffs.txt.bed_SORTED.bed | cut -f 1-3 | sort | uniq > WHOLE_GENOME_without_R_LOOPS.bed

# Sanity check ######

#Consensus_no_cutoffs.txt.bed_SORTED.bed + INTERGENIC_REFERENCE_SET_SUBTRACT.bed + GENE_BODY_REFERENCE_SET_SUBTRACT.bed + TSS_REFERENCE_SET.bed + TTS_REFERENCE_SET.bed + GENCODE_V19_TSS_1kb.bed_SORTED.bed + GENCODE_V19_TTS_1kb.bed_SORTED.bed = should equal whole genome

cut -f 1-3 Consensus_no_cutoffs.txt.bed_SORTED.bed | sort | uniq > sanity_check.bed
cut -f 1-3 INTERGENIC_REFERENCE_SET_SUBTRACT.bed  | sort | uniq >> sanity_check.bed 
cut -f 1-3 GENE_BODY_REFERENCE_SET_SUBTRACT.bed | sort | uniq >> sanity_check.bed 
cut -f 1-3 TSS_REFERENCE_SET.bed | sort | uniq >> sanity_check.bed 
cut -f 1-3 GENCODE_V19_TSS_1kb.bed_SORTED.bed | sort | uniq >> sanity_check.bed 
cut -f 1-3 TTS_REFERENCE_SET.bed  | sort | uniq >> sanity_check.bed 
cut -f 1-3 GENCODE_V19_TTS_1kb.bed_SORTED.bed  | sort | uniq >> sanity_check.bed 

sort -k1,1 -k2,2n sanity_check.bed  > sanity_check.bed_SORTED.bed
bedtools merge -i sanity_check.bed_SORTED.bed > sanity_check.bed_SORTED.bed_COLLAPSED.bed
perl calculate_total_distance_covered.pl 
2,982,704,466

bedtools subtract -a hg19_chromosome_sizes.txt_converted_to_bed.bed_SORTED.bed -b sanity_check.bed_SORTED.bed_COLLAPSED.bed | cut -f 1-3 | sort | uniq > missing_regions.bed
sort -k1,1 -k2,2n missing_regions.bed > missing_regions.bed_SORTED.bed

No missing regions:
wc missing_regions.bed*
       0       0       0 missing_regions.bed
       0       0       0 missing_regions.bed_SORTED.bed

#################################
** Now I need to make the high confidence peaks - with score cut-off of 100 and present in more than 5 files

perl filter_peaks_with_score_gte_100.pl

** Now get consensus regions for peaks > 100

./Run_use_bedtools_to_get_consensus_no_score_cut_off_FINAL_19_FILES_SCORE_GT_100.sh

output file : Consensus_sc_gt_100.txt

Now I need to convert to bed and only keep regions that overlap > 5 file and length < 5000

perl convert_consensus_to_bed_remove_greater_5kb_score_gt_100.pl 

output file : Consensus_sc_gt_100.txt_length_lt_5_num_gt_5.bed

sort -k1,1 -k2,2n Consensus_sc_gt_100.txt_length_lt_5_num_gt_5.bed > Consensus_sc_gt_100.txt_length_lt_5_num_gt_5.bed_SORTED.bed

** Now get consensus regions for peaks > 200

perl filter_peaks_with_score_gte_200.pl 

./Run_use_bedtools_to_get_consensus_no_score_cut_off_FINAL_19_FILES_SCORE_GT_200.sh

output file : Consensus_sc_gt_200.txt

perl convert_consensus_to_bed_remove_greater_5kb_score_gt_200.pl

output file: Consensus_sc_gt_200.txt_length_lt_5_num_gt_5.bed

sort -k1,1 -k2,2n Consensus_sc_gt_200.txt_length_lt_5_num_gt_5.bed > Consensus_sc_gt_200.txt_length_lt_5_num_gt_5.bed_SORTED.bed


###########################################################
* Final r-loop positive sets

Consensus_sc_gt_100.txt_length_lt_5_num_gt_5.bed_SORTED.bed
Consensus_sc_gt_200.txt_length_lt_5_num_gt_5.bed_SORTED.bed

###########################################################

