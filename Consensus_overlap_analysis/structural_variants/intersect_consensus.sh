testfor() {
    for i in input/sv_bedpe_files/*.bedpe; do 
        fname="$(basename -- $i)"
        echo working on $fname
        bedtools pairtobed -a "$i" -b input/consensus_no_cutoff.bed > \
            output/intersect_consensus/"$fname"
        done
}

find input/sv_bedpe_files/ -name "*.bedpe" | xargs -I{} basename '{}' | \
    parallel "bedtools pairtobed -a input/sv_bedpe_files/{} -b input/consensus_no_cutoff.bed > output/intersect_consensus/{} "


