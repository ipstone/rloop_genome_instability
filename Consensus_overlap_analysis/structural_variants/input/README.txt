 README.org

* chromosome bed file
 wget https://gist.githubusercontent.com/ianyfchang/dd9a5c93533ca65d73ac/raw/07269e6b1e7c8440f68d88cdf8d4ddc4d0f4c5ec/hg19.chrlen.bed
 cp hg19.chrlen.bed hg19_chromsomes.bed
 edit the file to make the bed  format: hg19_chromsomes.bed

* Consensus Rloop region bed files
  There bedfiles used in the script are the same output as:
    R-loop_consensus_regions/output/

  These bedfiles are modified to remove 'chr' from the chromosome
    consensus_no_cutoff.bed (not included in the code repo to save space)
    consensus_sc_gt_100.bed (not included)
    consensus_sc_gt_200.bed (included)

    negative_regions_without_consensus_rloop.bed

* nature2020_panCancerWGS_suppl_table1.txt
  This table list information for sample's tumor types etc. which is used in our analysis and plots.

* Rloop_intersected_score200
  Contains the bed files for TSS, TTS, gene-body regions that are intersected with consensus Rloops.(the specific version is the score 200 version)
