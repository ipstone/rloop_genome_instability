# Combine the different tumor SBS data, and plot violin plots for the TTS/TSS
# comparing Rloops vs. Non-Rloops regions
setwd("/data/projects/peix/rloop_project")

rm(list = ls())
library(data.table)
library(magrittr)
library(ipfun)
library(ggplot2)
library(ggpubr)
source("lib/lib_load_dataprep.R")

###############################################################################
## Global variable and setting sections
## Actually for totalSNV, it is not adjusted with GC content, as the result would be complicated to interpretate.
global_gcadjust <- FALSE # To choose whether we use gc adjustment for the SNV calculations
## -- this use of global gcadjustment would be logical ok, though might not needed at all when we look at totalSNV

# Combine all the different tumors SNV 6 categry data into one long data.table
tumors <- c(
    "ICGC_Breast", "ICGC_Liver", "ICGC_OV", "ICGC_CNS-Medullo",
    "TCGA_OV", "TCGA_ER", "TCGA_Liver", "ICGC_Pancreatic", "ICGC_Prostate",
    "Serena_ER", "Serena-Davies_ER_Biallelic", "TCGA_PRAD", "TCGA_LUAD",
    "TCGA_LUSC", "TCGA_SKCM", "ICGC_Melanoma", "TCGA_PAAD"
)

# Use grch37 non-n genome size to normal whole genome snv
human_genome_size <- 2897293955
# rloop_type <- get_rloop_transcribed()
# rloop_type <- get_rloop_type2() # tts, tss, genebody, pseudogene, lincRNA
# region_types =  c("tss", "tts", "genebody", "pseudogene", "lincRNA")

rloop_type <- get_rloop_transcribed() # tts, tss, genebody, pseudogene, lincRNA
region_types <- c(
    "tss", "tts", "tss_exome", "tts_exome",
    "tss_transcribed", "tts_transcribed",
    "tss_transcribed_exome", "tts_transcribed_exome",
    "consensus_sc200"
)

################################################################################
# Loading the SBS muttype counts from sbspath, return data.table with adjusted
# freq by dividing the adjsize
get_sbs_count_muttype_persample <- function(sbspath,
                                            sbstype = "SBS6",
                                            adjsize = human_genome_size,
                                            label = "wgs") {
    # Loading SBS data, from 6 types of SBS first
    dfile <- paste0(
        "analysis/mutation_matrix/", sbspath, "/output/SBS/icgc_ov.",
        sbstype, ".all"
        # sbstype, ".region" # The new intersected counts ends with extension .region
    )
    # d = fread("analysis/mutation_matrix/ICGC_OV/norloop_tss/output/SBS/icgc_ov.SBS6.all")
    d <- fread(dfile)
    dm <- melt(d, id.vars = "MutationType")

    dm <- dm[variable != "="]
    dm$adjusted_count <- dm$value / adjsize * 1e6 # Adjust count to per million nucleotides

    # Adding total SBS count here as a MutationType
    total_dm <- dm[, .(totalSNV = sum(value)), by = "variable"]
    total_dm$MutationType <- "totalSNV"
    total_dm$adjusted_count <- total_dm$totalSNV / adjsize * 1e6
    total_dm$value <- total_dm$totalSNV
    total_dm_sel <- total_dm[, .(MutationType, variable, value, adjusted_count)]

    # Combinding total SNV counts with the SBS types counts
    dm <- rbind(dm, total_dm_sel)
    dm$type <- label
    return(dm)
}
# Testing the function
# td = get_sbs_count_muttype_persample("icgc_ov_analysis")

create_test_dataset <- function(sbs.type = "SBS6", a.file, b.file,
                                a.label, b.label,
                                a.adjsize, b.adjsize) {
    # Loading two different SBS dataset and then combine before plotting and
    # testing differences.
    a <- get_sbs_count_muttype_persample(a.file,
        sbstype = sbs.type,
        adjsize = a.adjsize,
        label = a.label
    )
    b <- get_sbs_count_muttype_persample(b.file,
        sbstype = sbs.type,
        adjsize = b.adjsize,
        label = b.label
    )
    both <- rbind(a, b)
    return(both)
}


create_testdata_region <- function(sbs.type = "SBS6", label, tumor, gc.adjust = global_gcadjust) {
    # When given the region type, call the create_test_dataset to creat eh
    # regions
    # Return testdataset for rloop/norloop type: label, sbs.type
    # The label parameter is for the region types
    a.label <- paste0("rloop_", label)
    b.label <- paste0("norloop_", label)

    a.file <- paste0(tumor, "/", a.label)
    print(paste0("-- ", "Working on ", a.file))
    b.file <- paste0(tumor, "/", b.label)
    print(paste0("-- ", "Working on ", b.file))

    a.adjsize <- rloop_type$genome_covered[rloop_type$all == a.label]
    b.adjsize <- rloop_type$genome_covered[rloop_type$all == b.label]

    # Create the merged region of rloop, norloop
    m <- create_test_dataset(
        sbs.type,
        a.file, b.file,
        a.label, b.label,
        a.adjsize, b.adjsize
    )

    # Adjusting with rloop region gc pct (C-> adjust with gc, T-> adjust with 1-gcpct)
    if (gc.adjust == T) {
        # -- When adjusting gc content, adjust it to the assumed 50% by dividing 2
        # -- We could also *0.409 for the whole human genome gc content
        # -- that will give the ratioed but essentially eventual same result when comparing
        a.adj.gc <- rloop_type$gcpct[rloop_type$all == a.label]
        b.adj.gc <- rloop_type$gcpct[rloop_type$all == b.label]
        m[
            type == a.label & grepl("C>", MutationType),
            adjusted_count := .(adjusted_count / a.adj.gc / 2)
        ] # Adjust the gc content
        m[
            type == a.label & grepl("T>", MutationType),
            adjusted_count := .(adjusted_count / (1 - a.adj.gc) / 2)
        ]
        m[
            type == b.label & grepl("C>", MutationType),
            adjusted_count := .(adjusted_count / b.adj.gc / 2)
        ] # Adjust the gc content
        m[
            type == b.label & grepl("T>", MutationType),
            adjusted_count := .(adjusted_count / (1 - b.adj.gc) / 2)
        ]
    }

    # Adding the tumor type information to the SBS type count data.table
    m$tumor <- tumor
    return(m)
}
td <- create_testdata_region("SBS6", "tts", "ICGC_OV")
# table(td$type)


obtain_all_tumor_rloop_counts <- function(sbs.type = "SBS6", paired = FALSE) {
    # Compare between rloop vs non-rloop for different mutationType
    # Return all test results of all the different rlooptypes
    # ---------------------------------------
    # fr <- list()
    # for (i in region_types) {
    # tdata <- create_testdata_region(sbs.type, i)
    # fr[[i]] <- test_muttype_region(tdata, label = i, paired.test = paired)
    # }
    # return(fr)
    # ---------------------------------------

    tumor_data_list <- list()
    fr <- list()

    # Two levels of loop: First tumors, Second regions
    for (j in tumors) {
        for (i in region_types) {
            fr[[i]] <- create_testdata_region(sbs.type, i, j)
            fr[[i]]$region <- i
        }
        tumor_data_list[[j]] <- rbindlist(fr)
    }

    all_tumor_region_data <- rbindlist(tumor_data_list)
    return(all_tumor_region_data)
}
# a = obtain_all_tumor_rloop_counts()
# table(a$MutationType)
# table(a$tumor)
# table(a$type)


save_all_tumor_rloop_count_data <- function() {
    # save_all_tumor_rloop_count_data()
    a <- obtain_all_tumor_rloop_counts()
    mkdirp("output/SBS_violin_plots")
    fwrite(a, "output/SBS_violin_plots/all_tumors_SBS6_rloop_counts.csv")
}


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

long_tumors_list <- c(
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

plot_figure_violin <- function() {
    # 2021-08-10 update: According to paper's figure outlines to plot all
    # tumors' SNV counts together
    d <- fread("output/SBS_violin_plots/all_tumors_SBS6_rloop_counts.csv")
    d$rloop <- sapply(strsplit(d$type, "_", fixed = T), `[`, 1)
    # table(d$rloop)

    # table(d$tumor)
    # table(d$type)
    # Only keep the tss and tts region, for now only plot totalSNV
    d <- d[region %in% c("tss", "tts") & MutationType == "totalSNV"]
    # d <- d[region %in% c("consensus_sc200") & MutationType == "totalSNV"]
    d <- d[tumor %in% long_tumors_list]
    d$tumor <- factor(d$tumor, levels = long_tumors_list)
    d$rloop <- factor(d$rloop, levels = c("rloop", "norloop"))

    # Select the short tumor list data
    d.s <- d[tumor %in% short_tumors_list]

    my_violin <- function(dset, fname, plot.fun = ggviolin, logscale = T) {
        # with the given dataset, save the plot in fname
        p <- plot.fun(dset,
            x = "rloop", y = "adjusted_count",
            # x = "tumor", y = "adjusted_count",
            color = "rloop", palette = "npg",
            size = 0.5,
            alpha = 0.3,
            ylab = "log10 - Adjusted SNV per million",
            title = "ICGC tumors Rloops vs Non-Rloops regions",
            add = c("median_iqr", "jitter"),
        ) + font("x", size = 8) +
            # yscale("log10", .format = T) +
            # ylim(y.lim[1], y.lim[2]) +
            font("legend.text", size = 6) +
            theme(axis.text.x = element_text(angle = 90)) +
            stat_compare_means(aes(group = rloop),
                method = "wilcox.test",
                label = "..p.format..",
                # label = "..p.signif..",
                # paired = T,
                size = 3
            ) + facet_wrap(~region, ncol = 2)

        if (logscale == T) {
            p <- p + yscale("log10", .format = T)
        }
        ggsave(fname, p, width = 10, height = 5)
    }

    my_violin(d.s, "output/SBS_violin_plots/ggviolin_SBS6_ICGC_tumors_short-list.png")
    my_violin(d, "output/SBS_violin_plots/ggviolin_SBS6_ICGC_tumor_long-lists.png")
    my_violin(d.s, "output/SBS_violin_plots/ggbox_SBS6_ICGC_tumors_short-list.png", ggboxplot)
    my_violin(d, "output/SBS_violin_plots/ggbox_SBS6_ICGC_tumor_long-lists.png", ggboxplot)
    my_violin(d.s, "output/SBS_violin_plots/ggbox_nolog_SBS6_ICGC_tumors_short-list.png",
        ggboxplot,
        logscale = F
    )
    my_violin(d, "output/SBS_violin_plots/ggbox_nolog_SBS6_ICGC_tumor_long-lists.png",
        ggboxplot,
        logscale = F
    )
}

save_sbs_plot_data <- function() {
    # this is to provide plot data to Manisa
    d <- fread("output/SBS_violin_plots/all_tumors_SBS6_rloop_counts.csv")
    d$rloop <- sapply(strsplit(d$type, "_", fixed = T), `[`, 1)
    # table(d$rloop)

    # table(d$tumor)
    # table(d$type)
    # Only keep the tss and tts region, for now only plot totalSNV
    d <- d[region %in% c("tss", "tts") & MutationType == "totalSNV"]
    # d <- d[region %in% c("consensus_sc200") & MutationType == "totalSNV"]
    d <- d[tumor %in% long_tumors_list]
    d$tumor <- factor(d$tumor, levels = long_tumors_list)
    d$rloop <- factor(d$rloop, levels = c("rloop", "norloop"))

    # Select the short tumor list data
    d.s <- d[tumor %in% short_tumors_list]
    mkdirp("output/fig_plot_data")
    fwrite(d.s, "output/fig_plot_data/SBS_counts_short-tumor-list.tsv", sep = "\t")
    return(d.s)
}


main <- function() {
    # save_all_tumor_rloop_count_data()
    # plot_figure_violin()
    # save_sbs_plot_data()
}
main()
a <- save_sbs_plot_data()
