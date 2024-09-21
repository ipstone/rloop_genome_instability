# Intersect the sv bedpe file with the rloop consensus score 100 regions

find input/sv_bedpe_files/ -name "*.bedpe" | xargs -I{} basename '{}' | \
    parallel "bedtools pairtobed -a input/sv_bedpe_files/{} -b input/consensus_sc_gt_200.bed > output/intersect_sc200/{} "


