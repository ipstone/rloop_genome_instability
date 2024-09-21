# Add chr to the chromsome in gencode bed files
rm(list = ls())
library(data.table)
library(magrittr)
#library(tidyverse)
library(ipfun) 

d = fread("gencode.v19.basic.exome.bed")
d$V1 = paste0('chr', d$V1)
fwrite(d,"exome_chr.bed", sep="\t", col.names =FALSE )

