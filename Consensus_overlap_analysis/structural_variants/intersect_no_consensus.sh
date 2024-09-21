# Use parallel to intersect the consensus Rloop regions for non-matched
# regions: which are the structural variants not in the Rloop region

find input/sv_bedpe_files/ -name "*.bedpe" | xargs -I{} basename '{}' | \
    parallel "bedtools pairtobed -a input/sv_bedpe_files/{} -b input/consensus_no_cutoff.bed -type neither > output/intersect_no_consensus/{} "


