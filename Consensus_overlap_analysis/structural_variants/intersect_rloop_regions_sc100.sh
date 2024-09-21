## To intersect SVs with genebody, TSS, TTS for rloop and nonrloop regions
genebody="input/output_isaac_score100/rloop_genebody.bed"
nogenebody="input/output_isaac_score100/norloop_genebody.bed"
tss="input/output_isaac_score100/rloop_tss.bed"
notss="input/output_isaac_score100/norloop_tss.bed"
tts="input/output_isaac_score100/rloop_tts.bed"
notts="input/output_isaac_score100/norloop_tts.bed"

mkdir -p output/intersect_regions_sc100/rloop_genebody/
mkdir -p output/intersect_regions_sc100/norloop_genebody/
mkdir -p output/intersect_regions_sc100/rloop_tss/
mkdir -p output/intersect_regions_sc100/norloop_tss/
mkdir -p output/intersect_regions_sc100/rloop_tts/
mkdir -p output/intersect_regions_sc100/norloop_tts/

# Use parallel to intersect the consensus Rloop intersect_regions_sc100 for non-matched
# intersect_regions_sc100: which are the structural variants not in the Rloop region

find input/sv_bedpe_files/ -name "*.bedpe" | xargs -I{} basename '{}' | \
    parallel "bedtools pairtobed -a input/sv_bedpe_files/{} -b $genebody > output/intersect_regions_sc100/rloop_genebody/{} "
find input/sv_bedpe_files/ -name "*.bedpe" | xargs -I{} basename '{}' | \
    parallel "bedtools pairtobed -a input/sv_bedpe_files/{} -b $nogenebody > output/intersect_regions_sc100/norloop_genebody/{} "

find input/sv_bedpe_files/ -name "*.bedpe" | xargs -I{} basename '{}' | \
    parallel "bedtools pairtobed -a input/sv_bedpe_files/{} -b $tss > output/intersect_regions_sc100/rloop_tss/{} "
find input/sv_bedpe_files/ -name "*.bedpe" | xargs -I{} basename '{}' | \
    parallel "bedtools pairtobed -a input/sv_bedpe_files/{} -b $notss > output/intersect_regions_sc100/norloop_tss/{} "

find input/sv_bedpe_files/ -name "*.bedpe" | xargs -I{} basename '{}' | \
    parallel "bedtools pairtobed -a input/sv_bedpe_files/{} -b $tts > output/intersect_regions_sc100/rloop_tts/{} "
find input/sv_bedpe_files/ -name "*.bedpe" | xargs -I{} basename '{}' | \
    parallel "bedtools pairtobed -a input/sv_bedpe_files/{} -b $notts > output/intersect_regions_sc100/norloop_tts/{} "
