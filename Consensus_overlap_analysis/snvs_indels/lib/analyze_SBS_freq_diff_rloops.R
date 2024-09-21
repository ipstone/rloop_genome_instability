# To analyze the SBS mutaion frequency in different context of rloops regions
# - Compare rloop genomic regons with non-rloop regions
# - Compare rloop genomic regions with the overall genome regions

rm(list = ls())
library(data.table)
library(magrittr)
library(ipfun)
library(ggplot2)
library(ggpubr)
source("lib/lib_load_dataprep.R")

args <- commandArgs(trailingOnly = TRUE)
# cat(args, sep = "\n") # returned list of given argument
# -- except script name and other common optons
print(args)
project_folder <- paste0(args[1], "/") # using the given parameter 1 as the project folder
#project_folder = "ICGC_Liver/"

# Specify which project folder to process the data
# project_folder <- "ICGC_OV/"
# project_folder <- "Serena_ER/"
# project_folder <- "Serena_TRP/"
# project_folder <- "Serena-Davies_ER_Biallelic/"
# project_folder <- "Serena-Davies_ER_Control/"
# project_folder <- "ICGC_Breast/"
# project_folder <- "ICGC_Liver/"
# project_folder <- "ICGC_Pancreatic/"
# project_folder <- "ICGC_Prostate/"

# project_folder <- "Serena-Davies_TRP_Biallelic/"
# project_folder <- "Serena-Davies_TRP_Control/"
# project_folder <- "ICGC_OV_Biallelic/"
# project_folder <- "ICGC_OV_Control/"
# project_folder <- "ICGC_OV_highVAF_26/"

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

## Section to compare with whole genome -----------------------------
old_testing_whole_genome_rloop <- function() {
    # Testing functions with the whole genome results
    create_testdata_all_region <- function(sbs.type, label) {
        a.file <- "icgc_ov_analysis"
        a.adjsize <- human_genome_size
        a.label <- "wgs"
        b.file <- paste0("ICGC_OV/", label)
        b.adjsize <- rloop_type$genome_covered[rloop_type$all == label]
        b.label <- label
        return(create_test_dataset(
            sbs.type, a.file, b.file,
            a.label, b.label,
            a.adjsize, b.adjsize
        ))
    }
    # td = create_testdata_all_region("SBS6","rloop_exon")

    # Wilcox test each category of MutationType in the given dataset
    test_muttype_wholegenome <- function(testdata, label = "rloop_exon") {
        # Return a list/result object from the testdata MutationType by comparing
        # each type in the test data
        fr <- list() # functions result
        for (i in unique(testdata$MutationType)) {
            subdata <- testdata[MutationType == i]
            fr[[paste0(label, "_", i)]] <-
                wilcox.test(adjusted_count ~ type, data = subdata)$p.value
        }
        return(fr)
    }
    # tr = test_muttype_wholegenome(td, label="rloop_exon")

    # Obtain all rloop regions type tests results with whole genome
    obtain_wholegenome_all_rloop_types_tests <- function(sbs.type = "SBS6") {
        # Return all test results of all the different rlooptypes
        rtypes <- rloop_type$all
        fr <- list()
        for (i in rtypes) {
            tdata <- create_testdata_all_region(sbs.type, i)
            fr[[i]] <- test_muttype_wholegenome(tdata, label = i)
        }
        return(fr)
    }
    # tr <- obtain_wholegenome_all_rloop_types_tests("SBS6")
}

#---------------------------------------- 
## Section of function to compare rloop vs non-rloop for each genome type
create_testdata_region <- function(sbs.type, label, gc.adjust = TRUE) {
    # Return testdataset for rloop/norloop type: label, sbs.type
    a.label <- paste0("rloop_", label)
    b.label <- paste0("norloop_", label)

    a.file <- paste0(project_folder, a.label)
    print(paste0("-- ", "Working on ", a.file))
    b.file <- paste0(project_folder, b.label)
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
        a.adj.gc <- rloop_type$gcpct[rloop_type$all == a.label]
        b.adj.gc <- rloop_type$gcpct[rloop_type$all == b.label]
        m[
            type == a.label & grepl("C>", MutationType),
            adjusted_count := .(adjusted_count / a.adj.gc / 2)
        ] # Adjust the gc content
        # -- When adjusting gc content, adjust it to the assumed 50% by dividing 2
        # -- We could also *0.409 for the whole human genome gc content
        # -- that will give the ratioed but essentially eventual same result when comparing
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
    return(m)
}

# td = create_testdata_region("SBS6","exon")
# table(td$type)

test_muttype_region <- function(testdata, label = "genebody", paired.test = FALSE) {
    # Compare rloop vs norloop relevant regions
    # Return a list/result object from the testdata MutationType by comparing
    # each type in the test data
    # ---------------------------------------
    # -- previous code to return a list
    # fr <- list() # functions result
    # for (i in unique(testdata$MutationType)) {
    # subdata <- testdata[MutationType == i]
    # fr[[paste0(label, "_", i)]] <-
    # wilcox.test(adjusted_count ~ type, data = subdata)$p.value
    # }
    # return(fr)
    # ---------------------------------------

    meano <- function(x) {
        y <- ifelse(is.na(x), 0, x)
        return(mean(y) + 5e-324)
    }
    mediano <- function(x) {
        y <- ifelse(is.na(x), 0, x)
        return(median(y) + 5e-324)
    }

    mediano <- function(x) median(na.omit(x)) + 5e-324

    allrows <- data.table() # keep all rows as value to return

    for (i in unique(testdata$MutationType)) {
        # Obtain each mutationType data for testing/calculations
        subdata <- testdata[MutationType == i]
        rloop <- subdata[type == paste0("rloop_", label)]
        norloop <- subdata[type == paste0("norloop_", label)]
        rowd.pvalue <- wilcox.test(adjusted_count ~ type,
            data = subdata
        )$p.value
        # alternative = "less")$p.value
        mean_rloop <- meano(rloop$adjusted_count)
        mean_norloop <- meano(norloop$adjusted_count)
        rowd.fold <- mean_rloop / mean_norloop

        # Obtain the median count/ratio
        median_rloop <- mediano(rloop$adjusted_count)
        median_norloop <- mediano(norloop$adjusted_count)

        rowd.median.fold <- median_rloop / median_norloop
        rowdt <- data.table(
            MutationType = i, pval = rowd.pvalue, foldChange = rowd.fold,
            rloop_mean = mean_rloop, norloop_mean = mean_norloop,
            rloop_median = median_rloop, norloop_median = median_norloop,
            foldChangeMedian = rowd.median.fold
        )
        allrows <- rbind(allrows, rowdt)
    }
    return(allrows)
}

obtain_all_kinds_rloop_tests <- function(sbs.type = "SBS6", paired = FALSE) {
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

    fr <- list()
    for (i in region_types) {
        tdata <- create_testdata_region(sbs.type, i)
        fr[[i]] <- test_muttype_region(tdata, label = i, paired.test = paired)
        fr[[i]]$region <- i
    }
    allregions <- rbindlist(fr)
    return(allregions)
}

save_all_test_result <- function(sbs.type = "SBS6") {
    tr <- obtain_all_kinds_rloop_tests(sbs.type = sbs.type)
    project <- substr(project_folder, 1, nchar(project_folder) - 1)

    # Adding the following lines to make sure the output folder exist
    output_folder <- paste0("output/SBS_compare_rloop_regions/", project_folder)
    mkdirp(output_folder)

    outputfile <- paste0(
        "output/SBS_compare_rloop_regions/", project,
        "/", sbs.type, "_test_results.tsv"
    )
    fwrite(tr, outputfile, sep = "\t")
}

# tr <- obtain_all_kinds_rloop_tests("SBS18")
# tr.paired <- obtain_all_kinds_rloop_tests("SBS18", paired = T)

#------------------------------------------------------------
## Section to create plots for each of the SBS categories (ggviolin or
## ggboxplot
saveplot_rloop_types_SBS_plot <- function(sbs.type = "SBS6", n.row = 2, adj.gc = T, y.lim = c(0, 10)) {
    # Will create the sbs.type plots for each mutationType category comparison
    # of the type (rloop vs norloop) pannels, there will be 4 pannels created
    output_folder <- paste0("output/SBS_compare_rloop_regions/", project_folder)
    mkdirp(output_folder)
    gc_label <- ""
    if (adj.gc == T) gc_label <- "_gc_adjusted"

    pdf(paste0(output_folder, sbs.type, gc_label, "_compare_rloop_regons.pdf"))

    for (i in region_types) {
        tdata <- create_testdata_region(sbs.type, i, gc.adjust = adj.gc)
        print(i) # -- debug/track run progress
        print(adj.gc) # -- debug/track run progress
        p <- ggboxplot(tdata,
            x = "type", y = "adjusted_count",
            color = "type", palette = "npg",
            # size = 0.1
            alpha = 0.7,
            ylab = "SNV per million",
            title = i,
            add = c("jitter", "median_iqr"),
        ) + rremove("x.text") +
            font("x", size = 8) +
            ylim(y.lim[1], y.lim[2]) +
            font("legend.text", size = 8) +
            stat_compare_means(aes(group = type),
                method = "wilcox.test",
                label = "..p.format..",
                # label = "..p.signif..",
                # paired = T,
                size = 3
            ) + facet_wrap(~MutationType, nrow = n.row)
        print(p)
    }
    dev.off()
}

plot_rloop_type_compare <- function(sbs.type,
                                    rloop_region,
                                    adj.gc = T,
                                    y.lim = c(0, 10),
                                    n.row = 2) {
    # With the input rloop_region: return ggplot object of comparing
    # rloop_region is one of the following 4 types: c("exon", "genebody", "tss", "tts")
    i <- rloop_region

    tdata <- create_testdata_region(sbs.type, i, gc.adjust = adj.gc)
    p <- ggboxplot(tdata,
        x = "type", y = "adjusted_count",
        color = "type", palette = "npg",
        # size = 0.1
        alpha = 0.7,
        ylab = "SNV per million",
        title = i,
        add = c("jitter", "median_iqr"),
    ) + rremove("x.text") +
        font("x", size = 8) +
        ylim(y.lim[1], y.lim[2]) +
        font("legend.text", size = 8) +
        stat_compare_means(aes(group = type),
            method = "wilcox.test",
            label = "..p.format..",
            # label = "..p.signif..",
            # paired = T, # Added to test with paired wilcox test
            size = 3
        ) + facet_wrap(~MutationType, nrow = n.row)
    return(p)
}

saveplot_3regionType_compare <- function(sbs.type, ylim = c(0, 3), nrow = 2) {
    p.1 <- plot_rloop_type_compare(sbs.type, "intergenic", y.lim = ylim, n.row = nrow)
    # p.1 <- plot_rloop_type_compare(sbs.type, "exon", y.lim = ylim, n.row = nrow)
    p.2 <- plot_rloop_type_compare(sbs.type, "tss", y.lim = ylim, n.row = nrow)
    p.3 <- plot_rloop_type_compare(sbs.type, "tts", y.lim = ylim, n.row = nrow)
    p <- ggarrange(p.1, p.2, p.3, nrow = 1)

    output_folder <- paste0("output/SBS_compare_rloop_regions/", project_folder)
    mkdirp(output_folder)
    gc_label <- "_gc_adjusted"

    ofile <- paste0(output_folder, sbs.type, gc_label, "_compare_all3_regons.pdf")
    ggsave(ofile, p, width = 18, height = 12)
}

saveplot_2regionType_compare <- function(sbs.type, ylim = c(0, 3), nrow = 2) {
    # p.1 <- plot_rloop_type_compare(sbs.type, "exon", y.lim = ylim, n.row = nrow)
    p.2 <- plot_rloop_type_compare(sbs.type, "_tss", y.lim = ylim, n.row = nrow)
    p.3 <- plot_rloop_type_compare(sbs.type, "_tts", y.lim = ylim, n.row = nrow)
    p <- ggarrange(p.2, p.3, nrow = 1)

    output_folder <- paste0("output/SBS_compare_rloop_regions/", project_folder)
    mkdirp(output_folder)
    gc_label <- "_gc_adjusted"

    ofile <- paste0(output_folder, sbs.type, gc_label, "_compare_all2_regons.pdf")
    ggsave(ofile, p, width = 18, height = 12)
}

saveplot_4regionType_compare <- function(sbs.type, ylim = c(0, 3), nrow = 2) {
    p.1 <- plot_rloop_type_compare(sbs.type, "intergenic", y.lim = ylim, n.row = nrow)
    p.2 <- plot_rloop_type_compare(sbs.type, "genebody", y.lim = ylim, n.row = nrow)
    p.3 <- plot_rloop_type_compare(sbs.type, "tss", y.lim = ylim, n.row = nrow)
    p.4 <- plot_rloop_type_compare(sbs.type, "tts", y.lim = ylim, n.row = nrow)
    p <- ggarrange(p.1, p.2, p.3, p.4, nrow = 1)

    mkdirp("output/SBS_compare_rloop_regions/")
    output_folder <- paste0("output/SBS_compare_rloop_regions/", project_folder)
    mkdirp(output_folder)
    gc_label <- "_gc_adjusted"

    ofile <- paste0(output_folder, sbs.type, gc_label, "_compare_all4_regons.pdf")
    ggsave(ofile, p, width = 24, height = 12)
}

main <- function() {
    # Run different type of SBS
    # saveplot_rloop_types_SBS_plot("SBS6", y.lim = c(0, 5))
    # saveplot_rloop_types_SBS_plot("SBS18", n.row = 3, y.lim = c(0, 3))
    # saveplot_rloop_types_SBS_plot("SBS24", n.row = 4, y.lim = c(0, 3))

    # Plot with gc_adjusted
    # saveplot_rloop_types_SBS_plot("SBS6", adj.gc = T, y.lim = c(0, 5))
    # saveplot_rloop_types_SBS_plot("SBS18",
    #     n.row = 3, adj.gc = T,
    #     y.lim = c(0, 3)
    # )
    # saveplot_rloop_types_SBS_plot("SBS24",
    #     n.row = 4, adj.gc = T,
    #     y.lim = c(0, 3)
    # )

    # Save the 3 regions together types; these are all with GC adjustment
    # saveplot_3regionType_compare("SBS6", ylim = c(0, 5), nrow = 2)
    # saveplot_3regionType_compare("SBS18", ylim = c(0, 5), nrow = 3)
    # saveplot_2regionType_compare("SBS6", ylim = c(0, 5), nrow = 2)
    # saveplot_2regionType_compare("SBS18", ylim = c(0, 5), nrow = 3)
    # saveplot_2regionType_compare("SBS6", ylim = c(0, 5), nrow = 2)
    # saveplot_2regionType_compare("SBS18", ylim = c(0, 5), nrow = 3)

    saveplot_rloop_types_SBS_plot("SBS6", y.lim = c(0, 5))
    saveplot_rloop_types_SBS_plot("SBS18", n.row = 3, y.lim = c(0, 3))
    save_all_test_result("SBS6")
    save_all_test_result("SBS18")
}
main()
