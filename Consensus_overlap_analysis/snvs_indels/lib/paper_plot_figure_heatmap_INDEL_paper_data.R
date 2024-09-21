# Plot the heatmap for different rloop region compare results

rm(list = ls())
library(data.table)
library(magrittr)
library(tidyverse)
library(ipfun)
library(ggpubr)

tumors_list <- c(
    # "ICGC_Breast",
    "ICGC_Pancreatic",
    "TCGA_PAAD",
    "Serena_ER",
    "TCGA_ER",
    "ICGC_Prostate",
    "TCGA_PRAD",
    "ICGC_OV",
    "TCGA_OV",
    "ICGC_Liver",
    "TCGA_Liver",
    # "Serena-Davies_ER_Biallelic",
    "TCGA_LUAD",
    "TCGA_LUSC",
    "ICGC_Melanoma",
    "TCGA_SKCM"
    # "ICGC_CNS-Medullo"
    # , "Serena-Davies_ER_Biallelic"
)

# The short tumor list is for the current verison of the figures
short_tumors_list <- c(
    # "ICGC_Breast",
    "ICGC_Pancreatic",
    "Serena_ER",
    "ICGC_Prostate",
    "ICGC_OV",
    "ICGC_Liver",
    # "Serena-Davies_ER_Biallelic",
    "ICGC_Melanoma"
    # "ICGC_CNS-Medullo"
    # , "Serena-Davies_ER_Biallelic"
)

# tumors <- c(
# "ICGC_Breast", "ICGC_Liver",
# "ICGC_OV", "ICGC_CNS-Medullo",
# "ICGC_Pancreatic", "ICGC_Prostate",
# "Serena_ER")


collect_all_sbs_results <- function(testtype = "SBS6") {
    # Collect all different tumors SBS test results
    collect_tumor <- function(tumor) {
        input_file <- paste0(
            "SBS_compare_rloop_regions/",
            tumor, "/", testtype, "_test_results.tsv"
        )
        d <- fread(input_file)
        d$tumor_type <- tumor
        return(d)
    }
    alld <- rbindlist(lapply(tumors_list, collect_tumor))
    alld$tumor_type <- factor(alld$tumor_type, levels = tumors_list)
    return(alld)
}
# d = collect_all_sbs_results("SBS18")
# a = collect_all_sbs_results("SBS6")
# b = a[ tumor_type=="ICGC_CNS-Medullo"]
# bs= b[ MutationType=="T>G" & region == "tss_transcribed_exome"]
# bs= b[ region == "tss_transcribe_exome"]
# table(b$region)

# table(d$tumor_type)

collect_all_indel_results <- function() {
    # Collect all different tumors SBS test results
    collect_tumor <- function(tumor) {
        input_file <- paste0(
            "indel_compare_rloop_regions/",
            tumor, "/indel_test_results.tsv"
        )
        d <- fread(input_file)
        d$tumor_type <- tumor
        return(d)
    }
    alld <- rbindlist(lapply(tumors_list, collect_tumor))
    alld$tumor_type <- factor(alld$tumor_type, levels = tumors_list)
    return(alld)
}
# d = collect_all_indel_results()
# table(d$tumor_type)
plot_pval <- function(x) {
    # Return the heatmap plot object for pval with the range
    d <- x
    d$plot_pval <- ifelse(d$pval < 0.00001, 0.00001, d$pval)
    p1 <- ggplot(d, aes(x = MutationType, y = region)) +
        geom_tile(aes(fill = log10(plot_pval)), color = "white") +
        scale_fill_distiller(palette = "YlGnBu", limits = c(-5, 0)) +
        theme(axis.text.x = element_text(angle = 90))
    return(p1)
}

plot_foldChange <- function(x, mask.by.pval = F, fold.change.median = F) {
    # Return the heatmap plot object for pval with the range
    # If the fold.change.median is given by T, then return the
    # foldChangeMedian

    d <- x
    if (fold.change.median == T) {
        d$foldChange <- d$foldChangeMedian
    }

    d$foldChange <- ifelse(d$foldChange <= 0.125, 0.125, d$foldChange)
    d$foldChange <- ifelse(d$foldChange >= 8, 8, d$foldChange)
    # zCuts <- seq(3, -3, -0.5)

    # If given mask.by.pval ==T, mask these fold change plot
    if (mask.by.pval == T) {
        d$foldChange <- ifelse(d$pval > 0.07, NA, d$foldChange)
    }

    p2 <- ggplot(d, aes(
        x = MutationType, y = region,
        # fill = cut(log2(foldChange), zCuts))) +
        fill = log2(foldChange)
    )) +
        geom_tile(color = "white") +
        # scale_fill_brewer(palette="RdBu") +
        scale_fill_distiller(palette = "RdBu", limits = c(-3, 3)) +
        # scale_fill_distiller(palette = "RdBu", limits = c(-3, 3), na.value = NA) +
        # scale_fill_gradient2(midpoint=0) +
        # scale_fill_distiller(palette="RdBu") +
        # scale_fill_distiller(palette="Spectral") +
        theme(axis.text.x = element_text(angle = 90))
    return(p2)
}

plot_foldChange_byTumor <- function(x, mask.by.pval = F, fold.change.median = F) {
    # This version put tumor type on the X axis, then using facet to wrap for
    # different mutation types
    # Return the heatmap plot object for pval with the range
    # If the fold.change.median is given by T, then return the
    # foldChangeMedian

    d <- x
    if (fold.change.median == T) {
        d$foldChange <- d$foldChangeMedian
    }

    d$foldChange <- ifelse(d$foldChange <= 0.125, 0.125, d$foldChange)
    d$foldChange <- ifelse(d$foldChange >= 8, 8, d$foldChange)
    # zCuts <- seq(3, -3, -0.5)

    # If given mask.by.pval ==T, mask these fold change plot
    if (mask.by.pval == T) {
        d$foldChange <- ifelse(d$pval > 0.07, NA, d$foldChange)
    }

    p2 <- ggplot(d, aes(
        x = tumor_type, y = region,
        # fill = cut(log2(foldChange), zCuts))) +
        fill = log2(foldChange)
    )) +
        geom_tile(color = "white") +
        # scale_fill_brewer(palette="RdBu") +
        scale_fill_distiller(palette = "RdBu", limits = c(-3, 3)) +
        # scale_fill_distiller(palette = "RdBu", limits = c(-3, 3), na.value = NA) +
        # scale_fill_gradient2(midpoint=0) +
        # scale_fill_distiller(palette="RdBu") +
        # scale_fill_distiller(palette="Spectral") +
        theme(axis.text.x = element_text(angle = 90))
    return(p2)
}

saveplot_alltumors_heatmap <- function(testtype = "SBS6", column = 2) {
    # Save all tumors heatmap plots based on the testtype
    if (testtype == "indel") {
        get_data_fun <- collect_all_indel_results
        folderpath <- "indel_compare_rloop_regions/"
        fd <- get_data_fun()
    } else {
        get_data_fun <- collect_all_sbs_results
        folderpath <- "SBS_compare_rloop_regions/"
        fd <- get_data_fun(testtype)
    }

    # Filter out the region type to be plotted
    # fd = fd[ region =="tss" | region == "tts"]
    fd <- fd[region %in% c("tss", "tts", "consensus_sc200")]
    fd <- fd[tumor_type %in% short_tumors_list]
    fd$tumor_type <- factor(fd$tumor_type, levels = short_tumors_list)


    # Specifiy plot file
    pval_file <- paste0(folderpath, "all-tumors_", testtype, "_selected_pval_heatmap.png")
    fchg_file <- paste0(folderpath, "all-tumors_", testtype, "_selected_fold-change_heatmap.png")
    fchg_file2 <- paste0(folderpath, "all-tumors_", testtype, "_selected_fold-change_significant-only_heatmap.png")
    # -- this is the fold change plot with only significant wilcox test
    both_file <- paste0(folderpath, "all-tumors_", testtype, "_selected_both_heatmap.png")
    fchg_file_median <- paste0(folderpath, "all-tumors_", testtype, "_selected_fold-change-median_heatmap.png")
    # -- this is the fold change plot withthe median value fold change

    print(pval_file)
    p1 <- plot_pval(fd) +
        facet_wrap(~tumor_type, ncol = column)
    ggsave(pval_file, p1, width = 6, height = 8)

    print(fchg_file)
    p2 <- plot_foldChange(fd) +
        facet_wrap(~tumor_type, ncol = column)
    ggsave(fchg_file, p2, width = 6, height = 8)

    print(fchg_file2)
    p.foldchange.2 <- plot_foldChange(fd, mask.by.pval = T) +
        facet_wrap(~tumor_type, ncol = column)
    ggsave(fchg_file2, p.foldchange.2, width = 6, height = 8)


    print(both_file)
    pboth <- ggarrange(p1, p2, nrow = 1)
    # ggsave(both_file, pboth, width=16, height=12)
    ggsave(both_file, pboth, width = 12, height = 12)

    print(fchg_file_median)
    p.foldchange.median <- plot_foldChange(fd, fold.change.median = T) +
        facet_wrap(~tumor_type, ncol = column)
    ggsave(fchg_file_median, p.foldchange.median, width = 6, height = 8)
}

saveplot_alltumors_heatmap_byMutationTypes <- function(testtype = "SBS6", column = 2) {
    # Save all tumors heatmap plots based on the testtype
    # get_data_fun <- collect_all_sbs_results
    # folderpath <- "SBS_compare_rloop_regions/"
    # fd <- get_data_fun("SBS6")
    ## table(fd$region)

    # Save all tumors heatmap plots based on the testtype
    if (testtype == "indel") {
        get_data_fun <- collect_all_indel_results
        folderpath <- "indel_compare_rloop_regions/"
        fd <- get_data_fun()
    } else {
        get_data_fun <- collect_all_sbs_results
        folderpath <- "SBS_compare_rloop_regions/"
        fd <- get_data_fun(testtype)
    }

    # Filter out the region type to be plotted
    fd <- fd[region %in% c("tss", "tts", "consensus_sc200")]
    fd <- fd[tumor_type %in% short_tumors_list]
    fd$tumor_type <- factor(fd$tumor_type, levels = short_tumors_list)

    # Filter out the region type to be plotted
    # fd = fd[ region =="tss" | region == "tts"]

    # Specifiy plot file
    fchg_file <- paste0(folderpath, "all-tumors_", testtype, "_fold-change_selected_heatmap_byMutationType.png")
    fchg_file2 <- paste0(folderpath, "all-tumors_", testtype, "_fold-change_selected_heatmap_byMutationType_significant-only.png")

    print(fchg_file)
    p2 <- plot_foldChange_byTumor(fd) +
        facet_wrap(~MutationType, ncol = column)
    ggsave(fchg_file, p2, width = 6, height = 8)

    print(fchg_file2)
    p.foldchange.2 <- plot_foldChange_byTumor(fd, mask.by.pval = T) +
        facet_wrap(~MutationType, ncol = column)
    ggsave(fchg_file2, p.foldchange.2, width = 6, height = 8)
}

main <- function() {
    # saveplot_alltumors_heatmap("SBS6", column = 2)
    # saveplot_alltumors_heatmap_byMutationTypes("SBS6", column = 2)
    saveplot_alltumors_heatmap("indel", column = 2)
    saveplot_alltumors_heatmap_byMutationTypes("indel", column = 2)
}

main()

save_sbs_indel_result_data_for_paper<- function() {
    # a <- collect_all_sbs_results()
    # fwrite(a, "/data/projects/peix/rloop_project/output/fig_plot_data/SBS6_fold-change_test-rsults_all-tumors.tsv", sep = "\t")
    b <- collect_all_indel_results()
    fwrite(b, "/data/projects/peix/rloop_project/output/fig_plot_data/indel_fold-change_test-rsults_all-tumors.tsv", sep = "\t")
}
# save_sbs_indel_result_data_for_paper()
