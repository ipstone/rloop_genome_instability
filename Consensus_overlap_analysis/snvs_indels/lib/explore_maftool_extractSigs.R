# Explore extractSignatures function from maftools to extract signatures
rm(list = ls())
library(data.table)
library(magrittr)
library(ipfun)
library(maftools)
library(NMF)
library(pheatmap)

explore_overall_ov_sigs <- function() {
    # Exploring overall signatures with all OV SBS sigantures
    # Loading 96SBS matrix generated from sigProfilerMatrixGenerator
    d <- fread("analysis/icgc_ov_analysis/output/SBS/icgc_ov.SBS96.all")
    d[[2]] <- NULL
    dm <- as.matrix(d, rownames = "MutationType")
    laml.tnm <- list(nmf_matrix = t(dm)) # Provide data in similar format as maftools

    # Run extractSignatures function from maftools
    laml.sign <- estimateSignatures(mat = laml.tnm, nTry = 6)
    png("output/maftools_signatures/icgc_ov_cophenetic.png")
    plotCophenetic(res = laml.sign)
    dev.off()
    # from viewing teh cophenetic plots, n=4 drops most significantly

    laml.sig <- extractSignatures(mat = laml.tnm, n = 4)

    # Compate against original 30 signatures
    laml.og30.cosm <- compareSignatures(nmfRes = laml.sig, sig_db = "legacy")
    # Save the similarity to orginal cosmic 30 sigantures
    library("pheatmap")
    pheatmap(
        mat = laml.og30.cosm$cosine_similarities, cluster_rows = FALSE, main = "cosine similarity against validated signatures",
        filename = "output/maftools_signatures/icgc_ov_cosm_original30.pdf"
    )

    # Save the signatures itself
    png("output/maftools_signatures/icgc_ov_denovo_sigs.png")
    plotSignatures(nmfRes = laml.sig, title_size = 1.2, sig_db = "SBS")
    dev.off()
}

saveplot_cophenetic <- function(matrix_path, trynumber = 9) {
    # Loading SBS96 from given matrix path, and generate cophenetic plots
    f <- paste0("analysis/ICGC_OV/", matrix_path, "/output/SBS/icgc_ov.SBS96.all")

    # now loading and coding data to run
    d <- fread(f)
    d[[2]] <- NULL
    dm <- as.matrix(d, rownames = "MutationType")
    laml.tnm <- list(nmf_matrix = t(dm)) # Provide data in similar format as maftools

    # Run extractSignatures function from maftools
    laml.sign <- estimateSignatures(mat = laml.tnm, nTry = trynumber)
    outf <- paste0("output/maftools_signatures/icgc_ov_cophenetic_", matrix_path, ".png")
    png(outf)
    plotCophenetic(res = laml.sign)
    dev.off()
}

# Constructing combination of rloop, norloop types
convert_types <- expand.grid(
    c("rloop_", "norloop_"),
    c("exon", "genebody", "tss", "tts")
)
convert_types$all <- paste0(convert_types$Var1, convert_types$Var2)
# NOTE: norloop_genebody has little counts, need to remove
convert_types <- convert_types[convert_types$all != "norloop_genebody", ]

# Run saveplot for each rloop region combintations
run_save_all_groups_cophetic <- function() lapply(convert_types$all, saveplot_cophenetic)

## Function to run each condition for denovo signatures (set n=4)
saveplot_extractSigs <- function(matrix_path, sig.n = 4) {
    writeLines(paste("Working on:", matrix_path))
    # Loading SBS96 from given matrix path, and generate cophenetic plots
    f <- paste0("analysis/ICGC_OV/", matrix_path, "/output/SBS/icgc_ov.SBS96.all")
    # now loading and coding data to run
    d <- fread(f)
    d[[2]] <- NULL
    dm <- as.matrix(d, rownames = "MutationType")
    laml.tnm <- list(nmf_matrix = t(dm)) # Provide data in similar format as maftools

    # Save the extracted signature plot
    laml.sig <- extractSignatures(mat = laml.tnm, n = sig.n)

    # Compate against original 30 signatures
    heatmapfile <- paste0("output/maftools_signatures/icgc_ov_pheatmap_cosm30_", matrix_path, "_nsig", sig.n, ".pdf")
    laml.og30.cosm <- compareSignatures(nmfRes = laml.sig, sig_db = "legacy")
    # Save the similarity to orginal cosmic 30 sigantures
    pheatmap(
        mat = laml.og30.cosm$cosine_similarities, cluster_rows = FALSE, main = "cosine similarity against validated signatures",
        filename = heatmapfile
    )

    # Save the signatures plot itself
    denovo_sigfile <- paste0("output/maftools_signatures/icgc_ov_denovoSigs_", matrix_path, "_nsig", sig.n, ".png")
    png(denovo_sigfile)
    plotSignatures(nmfRes = laml.sig, title_size = 1.2, sig_db = "SBS")
    dev.off()

    # Save the denovo signature extracted and the contribution of eadh signatures to each samples
    table_signature <- paste0("output/maftools_signatures/tables/denovoSigs_", matrix_path, "_nsig", sig.n, ".csv")
    fwrite(data.table(laml.sig$signatures, keep.rownames = T), table_signature)
    # Save the contributions
    table_contribution <- paste0("output/maftools_signatures/tables/denovoSigs_", matrix_path, "_nsig", sig.n, "_contribution.csv")
    fwrite(data.table(laml.sig$contributions, keep.rownames = T), table_contribution)
}

# Run saveplot for each rloop region combintations
run_save_all_groups_extractSigs <- function(n.sig) {
    n_saveplot <- function(x) saveplot_extractSigs(x, sig.n = n.sig)
    lapply(convert_types$all, n_saveplot)
}

###############################################################
# main excution part of the script
main <- function() {
    # run_save_all_groups_cophetic()
    run_save_all_groups_extractSigs(4)
    run_save_all_groups_extractSigs(3)
}
main()