# For each regions in a subfolder of the argv $1, generate the line counts of
# each file (line-counts_) and combined_intersected_regions.tsv

#echo $1
#subfolders=$(find output/intersect_regions_sc200 -maxdepth 1 -mindepth 1 -type d  -printf '%f\n')
subfolders=$(find "$1" -maxdepth 1 -mindepth 1 -type d)
echo $subfolders

for i in $subfolders; do
    base=$(basename $i)
    echo working on $base
    python combine_intersected_sv_result.py $i > $1/combined_intersected_"$base".tsv
    cd $i 
    wc -l *.bedpe | head -n -1 > ../line-counts_"$base".txt 
    cd -
done

