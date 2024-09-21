# Plot the heatmap for different rloop region compare results
# From original icgc tumor plot to adapt for Serena's tumors specifically

rm(list = ls())
library(data.table)
library(magrittr)
library(tidyverse)
library(ipfun)
library(ggpubr)

# SV results are organized differently with SBS/indel results
# We have the results saved in the output/sv_test_results folder.

collect_all_sv_results <- function() {
    # Collect all different results in the output/sv_test_results folder,
    # adding tumor type and region tags
    # Note: in order to have this function working, I renamed tss_100 to
    # tss-100, so that the regions are all together when splitted
    files <- list.files("output/serena_sv_test_results/", pattern = "csv$")
    tumor <- sapply(strsplit(files, "_", fixed = T), `[`, 1)
    tags <- sapply(strsplit(files, "_", fixed = T), `[`, 3)
    regions <- gsub(".csv", "", tags, fixed = T)

    collect_tumor <- function(tumor_file) {
        input_file <- paste0(
            "output/serena_sv_test_results/",
            tumor_file
        )

        d <- fread(input_file)
        files <- tumor_file
        tumor <- sapply(strsplit(files, "_", fixed = T), `[`, 1)
        tags <- sapply(strsplit(files, "_", fixed = T), `[`, 3)
        regions <- gsub(".csv", "", tags, fixed = T)

        d$tumor_type <- tumor
        d$region <- regions
        d$file <- tumor_file
        return(d)
    }

    alld <- rbindlist(lapply(files, collect_tumor))
    return(alld)
}
# a = collect_all_sv_results()

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

saveplot_alltumors_sv_heatmap <- function(column = 2) {
    # Save all tumors SV heatmap plots based on the testtype
    fd <- collect_all_sv_results()
    folderpath <- "output/serena_heatmap_plots/"
    mkdirp(folderpath)
    testtype <- "SV"

    # Specifiy plot file
    pval_file <- paste0(folderpath, "all-tumors_", testtype, "_pval_heatmap.png")
    fchg_file <- paste0(folderpath, "all-tumors_", testtype, "_fold-change_heatmap.png")
    fchg_file2 <- paste0(folderpath, "all-tumors_", testtype, "_fold-change_significant-only_heatmap.png")
    # -- this is the fold change plot with only significant wilcox test
    both_file <- paste0(folderpath, "all-tumors_", testtype, "_both_heatmap.png")
    fchg_file_median <- paste0(folderpath, "all-tumors_", testtype, "_fold-change-median_heatmap.png")
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

saveplot_alltumors_sv_heatmap_byMutationTypes <- function(column = 2) {
    # Save all tumors heatmap plots based on the testtype

    fd <- collect_all_sv_results()
    folderpath <- "output/serena_heatmap_plots/"
    mkdirp(folderpath)
    testtype <- "SV"

    # Filter out the region type to be plotted
    # fd = fd[ region =="tss" | region == "tts"]

    # Specifiy plot file
    fchg_file <- paste0(folderpath, "all-tumors_", testtype, "_fold-change_heatmap_byMutationType.png")
    fchg_file2 <- paste0(folderpath, "all-tumors_", testtype, "_fold-change_heatmap_byMutationType_significant-only.png")

    print(fchg_file)
    p2 <- plot_foldChange_byTumor(fd) +
        facet_wrap(~MutationType, ncol = column)
    ggsave(fchg_file, p2, width = 6, height = 8)

    print(fchg_file2)
    p.foldchange.2 <- plot_foldChange_byTumor(fd, mask.by.pval = T) +
        facet_wrap(~MutationType, ncol = column)
    ggsave(fchg_file2, p.foldchange.2, width = 6, height = 8)
}

saveplot_selected_tumors_sv_heatmap_byMutationTypes <- function(column = 2) {
    # Save selected tumors heatmap plots based on the testtype

    fd <- collect_all_sv_results()
    fd <- fd[region != "consensus"]

    folderpath <- "output/serena_heatmap_plots/"
    mkdirp(folderpath)
    testtype <- "SV"

    # Filter out the region type to be plotted
    # fd = fd[ region =="tss" | region == "tts"]

    # Specifiy plot file
    fchg_file <- paste0(folderpath, "selected-tumors_", testtype, "_fold-change_heatmap_byMutationType.png")
    fchg_file2 <- paste0(folderpath, "selected-tumors_", testtype, "_fold-change_heatmap_byMutationType_significant-only.png")

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
    # Save the SV summary results heatmaps
    # saveplot_alltumors_sv_heatmap()
    # saveplot_alltumors_sv_heatmap_byMutationTypes()
    saveplot_selected_tumors_sv_heatmap_byMutationTypes()
}

# main()


save_allelic_foldChange_heatmap <- function(test.type = "rloop-vs-norloop") {
    a <- fread("output/serena_biallelic_tests/all_serena_biallelic_test_results.csv")
    # a$region <- paste0(a$test_group, "_", a$test_type, "_", a$group)
    a$region <- paste0(a$test_group, "_", a$group)
    a$tumor_type <- a$tumor
    a <- a[test_type == test.type]

    folderpath <- "output/serena_heatmap_plots/"
    mkdirp(folderpath)
    testtype <- test.type

    # Filter out the region type to be plotted
    # fd = fd[ region =="tss" | region == "tts"]

    # Specifiy plot file
    fchg_file <- paste0(folderpath, "allelic_compare-rloop_", testtype, "_fold-change_heatmap_byMutationType.png")
    fchg_file2 <- paste0(folderpath, "allelic_compare-rloop_", testtype, "_fold-change_heatmap_byMutationType_significant-only.png")

    print(fchg_file)
    p2 <- plot_foldChange_byTumor(a) +
        facet_wrap(~MutationType, ncol = 2)
    ggsave(fchg_file, p2, width = 6, height = 8)

    print(fchg_file2)
    p.foldchange.2 <- plot_foldChange_byTumor(a, mask.by.pval = T) +
        facet_wrap(~MutationType, ncol = 2)
    ggsave(fchg_file2, p.foldchange.2, width = 6, height = 8)
}

save_allelic_foldChange_heatmap("rloop-vs-norloop")
save_allelic_foldChange_heatmap("biallelic-vs-control")