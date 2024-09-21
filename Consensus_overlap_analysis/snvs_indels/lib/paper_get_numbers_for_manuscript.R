# To get the number of cases used in the analysis for the manuscript
rm(list = ls())
library(data.table)
library(magrittr)
library(ipfun)
library(dplyr)

sv <- fread("output/fig_plot_data/SV_counts_all-tumors.tsv")
table(sv$tumor_type)
sbs <- fread("output/fig_plot_data/SBS_counts_short-tumor-list.tsv")
lu(sbs$variable)
table(sbs$tumor)
for (i in unique(sbs$tumor)) {
    print(i)
    tumor <- sbs[sbs$tumor == i, ]
    print(lu(tumor$variable))
}
# [1] "ICGC_Liver"
# [1] 264
# [1] "ICGC_OV"
# [1] 71
# [1] "ICGC_Pancreatic"
# [1] 239
# [1] "ICGC_Prostate"
# [1] 189
# [1] "Serena_ER"
# [1] 320
# [1] "ICGC_Melanoma"
# [1] 70

combined <- fread("output/combined_sv_indel_density_rloop_sc200.tsv")
table(combined$Project_Code)
