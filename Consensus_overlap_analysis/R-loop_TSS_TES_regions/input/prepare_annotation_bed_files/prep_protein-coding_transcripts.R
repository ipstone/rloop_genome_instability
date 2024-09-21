# Prepare data for protein-coding transcript bed files for use / also TSS/TTS
# files.

rm(list = ls())
library(data.table)
library(magrittr)
#library(tidyverse)
library(ipfun)

# Filter and prepare the right type of protein coding transcript.
d = fread("./extract_gencode_annotation/gencode-v19-extracted.txt")

save_transcript_tss_tts <- function () {
    
    tp0 = d[ field_type=="transcript" ]
    data.table(table(tp0$gene_type))[order(N),]

    tp = d[ field_type=="transcript" & transcript_type=="protein_coding" & gene_type =="protein_coding" ]
    tp = tp[ chrom != "chrM" ]
#fwrite(tp, 'output/transcript_protein-coding_gencode-v19.tsv', sep="\t")
    fwrite(tp[, .(chrom, start-1, end)], 'output/transcript_protein-coding_gencode-v19.bed', sep="\t", col.names=F )
# -- the start position is adjusted to 0 based. End positin in bed is 1
# based-not-included, which would translate into end - 1 + 1 (the first
# not-included : so adding 1)



# Further coding TSS, and TTS with +/- 1000bp window
# -- for gencodeV19, the index is 1 based.
# -- the export for bedtools, start is 0 based, end is 1-based-not-included.
# TSS
    tp$tss_minus_1kb = ifelse(tp$strand =="+", tp$start - 1 - 1000 , tp$end - 1 - 1000 )
    tp$tss_plus_1kb =  ifelse(tp$strand =="+", tp$start + 1000 -1 , tp$end + 1000 -1 )
    tp$tss_minus_1kb = ifelse(tp$tss_minus_1kb <0, 0, tp$tss_minus_1kb)
#table(tp$tss_minus_1kb < tp$tss_plus_1kb)
    fwrite(tp[, .(chrom, tss_minus_1kb, tss_plus_1kb, 
                  strand, gene_name, gene_id, transcript_id)],
           "output/tss-1kb-window.bed", sep="\t", col.names = F)

# TTS
    tp$tts_minus_1kb = ifelse(tp$strand =="+", tp$end - 1- 1000 , tp$start - 1 -1000 )
    tp$tts_plus_1kb =  ifelse(tp$strand =="+", tp$end + 1000 -1 , tp$start + 1000 -1 )
    tp$tts_minus_1kb = ifelse(tp$tts_minus_1kb <0, 0, tp$tts_minus_1kb)
#table(tp$tts_minus_1kb < tp$tts_plus_1kb)
    fwrite(tp[, .(chrom, tts_minus_1kb, tts_plus_1kb, 
                  strand, gene_name, gene_id, transcript_id)],
           "output/tts-1kb-window.bed", sep="\t", col.names = F)


}

save_transcript_tss_tts_transcribed <- function () {
   # Save the TSS/TTS region but only covers the transcribed region 
    #tp0 = d[ field_type=="transcript" ]
    #data.table(table(tp0$gene_type))[order(N),]

    tp = d[ field_type=="transcript" & transcript_type=="protein_coding" & gene_type =="protein_coding" ]
    tp = tp[ chrom != "chrM" ]

    # -- the start position is adjusted to 0 based. End positin in bed is 1
    # based-not-included, which would translate into end - 1 + 1 (the first
    # not-included : so adding 1)

# Further coding TSS, and TTS with +/- 1000bp window
# -- for gencodeV19, the index is 1 based.
# -- the export for bedtools, start is 0 based, end is 1-based-not-included.
# TSS
    tp$tss_start = as.integer(ifelse(tp$strand =="+", tp$start - 1 , tp$end - 1 - 1000 ))
    tp$tss_end =  as.integer(ifelse(tp$strand =="+", tp$start + 1000 -1 , tp$end -1 ))
    #table(tp$tss_start - tp$tss_end)
    fwrite(tp[, .(chrom, tss_start, tss_end, 
                  strand, gene_name, gene_id, transcript_id)],
           "output/tss-transcribed-1kb-window.bed", sep="\t", col.names = F)

# TTS
    tp$tts_start = as.integer(ifelse(tp$strand =="+", tp$end - 1- 1000 , tp$start - 1 ))
    tp$tts_end =  as.integer(ifelse(tp$strand =="+", tp$end - 1 , tp$start + 1000 -1 ))
    #table(tp$tts_end - tp$tts_start)
    fwrite(tp[, .(chrom, tts_start, tts_end, 
                  strand, gene_name, gene_id, transcript_id)],
           "output/tts-transcribed-1kb-window.bed", sep="\t", col.names = F)

}

save_pseudogene_tss_tts <- function () {
    ## Save pseudogene tss tts start/end region
    
    #tp0 = d[ field_type=="transcript" ]
    #data.table(table(tp0$gene_type))[order(N),]
    #data.table(table(tp0$transcript_type))[order(N),]

    tp = d[ field_type=="transcript" & gene_type =="pseudogene" ]
    tp = tp[ chrom != "chrM" ]
#fwrite(tp, 'output/transcript_protein-coding_gencode-v19.tsv', sep="\t")
    #fwrite(tp[, .(chrom, start-1, end)], 'output/transcript_pseudogene_gencode-v19.bed', sep="\t", col.names=F )
    # -- the start position is adjusted to 0 based. End positin in bed is 1
    # based-not-included, which would translate into end - 1 + 1 (the first
    # not-included : so adding 1)



# Further coding TSS, and TTS with +/- 1000bp window
# -- for gencodeV19, the index is 1 based.
# -- the export for bedtools, start is 0 based, end is 1-based-not-included.
# TSS
    tp$tss_minus_1kb = ifelse(tp$strand =="+", tp$start - 1 - 1000 , tp$end - 1 - 1000 )
    tp$tss_plus_1kb =  ifelse(tp$strand =="+", tp$start + 1000 -1 , tp$end + 1000 -1 )
    tp$tss_minus_1kb = ifelse(tp$tss_minus_1kb <0, 0, tp$tss_minus_1kb)
#table(tp$tss_minus_1kb < tp$tss_plus_1kb)
    fwrite(tp[, .(chrom, tss_minus_1kb, tss_plus_1kb, 
                  strand, gene_name, gene_id, transcript_id)],
           "output/pseudogene_tss-1kb-window.bed", sep="\t", col.names = F)

# TTS
    tp$tts_minus_1kb = ifelse(tp$strand =="+", tp$end - 1- 1000 , tp$start - 1 -1000 )
    tp$tts_plus_1kb =  ifelse(tp$strand =="+", tp$end + 1000 -1 , tp$start + 1000 -1 )
    tp$tts_minus_1kb = ifelse(tp$tts_minus_1kb <0, 0, tp$tts_minus_1kb)
#table(tp$tts_minus_1kb < tp$tts_plus_1kb)
    fwrite(tp[, .(chrom, tts_minus_1kb, tts_plus_1kb, 
                  strand, gene_name, gene_id, transcript_id)],
           "output/pseudogene_tts-1kb-window.bed", sep="\t", col.names = F)
}

save_lincRNA_tss_tts <- function () {
    ## Save lincRNA tss tts start/end region
    
    #tp0 = d[ field_type=="transcript" ]
    #data.table(table(tp0$gene_type))[order(N),]
    #data.table(table(tp0$transcript_type))[order(N),]

    tp = d[ field_type=="transcript" & transcript_type =="lincRNA" ]
    tp = tp[ chrom != "chrM" ]
    # -- the start position is adjusted to 0 based. End positin in bed is 1
    # based-not-included, which would translate into end - 1 + 1 (the first
    # not-included : so adding 1)



# Further coding TSS, and TTS with +/- 1000bp window
# -- for gencodeV19, the index is 1 based.
# -- the export for bedtools, start is 0 based, end is 1-based-not-included.
# TSS
    tp$tss_minus_1kb = ifelse(tp$strand =="+", tp$start - 1 - 1000 , tp$end - 1 - 1000 )
    tp$tss_plus_1kb =  ifelse(tp$strand =="+", tp$start + 1000 -1 , tp$end + 1000 -1 )
    tp$tss_minus_1kb = ifelse(tp$tss_minus_1kb <0, 0, tp$tss_minus_1kb)
#table(tp$tss_minus_1kb < tp$tss_plus_1kb)
    fwrite(tp[, .(chrom, tss_minus_1kb, tss_plus_1kb, 
                  strand, gene_name, gene_id, transcript_id)],
           "output/lincRNA_tss-1kb-window.bed", sep="\t", col.names = F)

# TTS
    tp$tts_minus_1kb = ifelse(tp$strand =="+", tp$end - 1- 1000 , tp$start - 1 -1000 )
    tp$tts_plus_1kb =  ifelse(tp$strand =="+", tp$end + 1000 -1 , tp$start + 1000 -1 )
    tp$tts_minus_1kb = ifelse(tp$tts_minus_1kb <0, 0, tp$tts_minus_1kb)
#table(tp$tts_minus_1kb < tp$tts_plus_1kb)
    fwrite(tp[, .(chrom, tts_minus_1kb, tts_plus_1kb, 
                  strand, gene_name, gene_id, transcript_id)],
           "output/lincRNA_tts-1kb-window.bed", sep="\t", col.names = F)
}

save_genebody = function(){
    # Use field_type == gene, to select protein coding gene start/end pos on
    # chromosomes

    #gp = d[ field_type=="gene"]
    #data.table(table(gp$gene_type))[order(N),]

    gp = d[ field_type=="gene" & gene_type=="protein_coding"]
    data.table(table(gp$transcript_type))[order(N),]

    gp$gene_start = gp$start - 1
    gp$gene_end =  gp$end - 1

    fwrite(gp[, .(chrom, gene_start, gene_end, 
                  strand, gene_name, gene_id, transcript_id)],
           "output/genebody.bed", sep="\t", col.names = F)

}

save_lincRNA = function(){
    # Use the transcript table, to select the lincRNA and peusdo genes to save
    # respective bed files

    tp = d[ field_type=="transcript" & transcript_type =="lincRNA" ]
    #table(tp$transcript_type)
    #data.table(table(tp$transcript_type))[order(N),]

    tp$rna_start = tp$start - 1
    tp$rna_end =  tp$end - 1

    fwrite(tp[, .(chrom, rna_start, rna_end, 
                  strand, gene_name, gene_id, transcript_id)],
           "output/lincRNA.bed", sep="\t", col.names = F)

}

save_pseudogene = function(){
    # Use field_type == gene, to select peudogene coding gene start/end pos on
    # chromosomes

    gp = d[ field_type=="gene" & gene_type=="pseudogene"]
    gp$gene_start = gp$start - 1
    gp$gene_end =  gp$end - 1

    fwrite(gp[, .(chrom, gene_start, gene_end, 
                  strand, gene_name, gene_id, transcript_id)],
           "output/pseudogene.bed", sep="\t", col.names = F)

}

main <- function () {
    #save_genebody()
    #save_lincRNA()
    #save_pseudogene()
    #save_pseudogene_tss_tts()
    #save_lincRNA_tss_tts()
    #save_transcript_tss_tts()

    save_transcript_tss_tts_transcribed() #There is an interesting bug that double 1e7 was read by bedtools as 1
    # So to fix, we are using integer format here
}

main()
