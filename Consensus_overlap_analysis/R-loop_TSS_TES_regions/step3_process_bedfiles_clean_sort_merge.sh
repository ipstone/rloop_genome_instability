## The followings are steps to clean and prepare bedfiles for usage/input for
## my snv calculations etc.
## -- Start with the specific bedfiles version folder

# Clean the previous final_output 
mv -b output/final_output_beds _archive/
rm -rf output/cleaned_bedfiles
mkdir -p output/cleaned_bedfiles


## Start to clean bedfiles and then calculate gc content
cd tmp/

# Step 1: clean/locate only chromosome bed file entries
for i in *.bed; do Rscript ../lib/clean_input_bedfile.R $i ./cleaned_bedfiles/ ; done

# Step 2: sort and merge bedfiles 
cd ./cleaned_bedfiles   # move to the cleaned_bedfiles folder
mkdir final_output_beds
for i in *.bed; do bedtools sort -i $i | bedtools merge > final_output_beds/$i; done

cp -a final_output_beds ../output

