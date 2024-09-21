## To intersect SVs with genebody, TSS, TTS for rloop and nonrloop regions
genebody="input/Rloop_intersected_score200/rloop_genebody.bed"
nogenebody="input/Rloop_intersected_score200/norloop_genebody.bed"
tss="input/Rloop_intersected_score200/rloop_tss.bed"
notss="input/Rloop_intersected_score200/norloop_tss.bed"
tts="input/Rloop_intersected_score200/rloop_tts.bed"
notts="input/Rloop_intersected_score200/norloop_tts.bed"

#mkdir -p output/intersect_regions_sc200/rloop_genebody/
#mkdir -p output/intersect_regions_sc200/norloop_genebody/
mkdir -p output/intersect_regions_sc200/rloop_tss/
mkdir -p output/intersect_regions_sc200/norloop_tss/
mkdir -p output/intersect_regions_sc200/rloop_tts/
mkdir -p output/intersect_regions_sc200/norloop_tts/

# Use parallel to intersect the consensus Rloop regions for non-matched
# regions: which are the structural variants not in the Rloop region

#find input/sv_bedpe_files/ -name "*.bedpe" | xargs -I{} basename '{}' | \
    #parallel "bedtools pairtobed -a input/sv_bedpe_files/{} -b $genebody > output/intersect_regions_sc200/rloop_genebody/{} "
#find input/sv_bedpe_files/ -name "*.bedpe" | xargs -I{} basename '{}' | \
    #parallel "bedtools pairtobed -a input/sv_bedpe_files/{} -b $nogenebody > output/intersect_regions_sc200/norloop_genebody/{} "

find input/sv_bedpe_files/ -name "*.bedpe" | xargs -I{} basename '{}' | \
    parallel "bedtools pairtobed -a input/sv_bedpe_files/{} -b $tss > output/intersect_regions_sc200/rloop_tss/{} "
find input/sv_bedpe_files/ -name "*.bedpe" | xargs -I{} basename '{}' | \
    parallel "bedtools pairtobed -a input/sv_bedpe_files/{} -b $notss > output/intersect_regions_sc200/norloop_tss/{} "

find input/sv_bedpe_files/ -name "*.bedpe" | xargs -I{} basename '{}' | \
    parallel "bedtools pairtobed -a input/sv_bedpe_files/{} -b $tts > output/intersect_regions_sc200/rloop_tts/{} "
find input/sv_bedpe_files/ -name "*.bedpe" | xargs -I{} basename '{}' | \
    parallel "bedtools pairtobed -a input/sv_bedpe_files/{} -b $notts > output/intersect_regions_sc200/norloop_tts/{} "
