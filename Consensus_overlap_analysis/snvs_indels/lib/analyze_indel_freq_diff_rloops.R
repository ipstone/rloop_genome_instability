# To analyze the Indel mutaion frequency in different context of rloops regions
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

# Specify which project folder to process the data
# project_folder <- "ICGC_OV/"

# Use grch37 non-n genome size to normal whole genome snv
human_genome_size <- 2897293955
# rloop_type <- get_rloop_transcribed()
# rloop_type <- get_rloop_type2() # tts, tss, genebody, pseudogene, lincRNA
# region_types =  c("tss", "tts", "genebody", "pseudogene", "lincRNA")
# region_types =  c("tss", "tts", "genebody", "pseudogene")
# rloop_type <- get_rloop_pseudo_tss_tts() # tts, tss, genebody, pseudogene, lincRNA
# region_types =  c("tss_gene", "tts_gene",
# "tss_pseudogene", "tts_pseudogene",
# "tss_lincRNA", "tts_lincRNA")

rloop_type <- get_rloop_transcribed() # tts, tss, genebody, pseudogene, lincRNA
region_types <- c(
    "tss", "tts", "tss_exome", "tts_exome",
    "tss_transcribed", "tts_transcribed",
    "tss_transcribed_exome", "tts_transcribed_exome",
    "consensus_sc200"
)

# The type of indel fields to test
indel_category <- c(
    "adj_del_count", "adj_ins_count", "adj_long_del",
    "adj_long_ins", "adj_mh", "adj_complex", "adj_total",
    "adj_long_indel", "adj_long_simple_indel"
)

# Loading the indel muttype counts from sbspath, return data.table with adjusted
# freq by dividing the adjsize
get_indel_count_muttype_persample <- function(sbspath,
                                              sbstype = "ID28",
                                              adjsize = human_genome_size,
                                              label = "wgs") {
    # Loading SBS data, from 6 types of SBS first
    dfile <- paste0(
        "analysis/mutation_matrix/", sbspath, "/output/ID/icgc_ov.",
        sbstype, ".all"
        # sbstype, ".region" # The new intersected counts ends with extension .region
    )
    # -- When called for specific rloop type, the rloop type subfolder will be
    # included in the sbspath when given to this function

    d <- fread(dfile)
    # d <- fread("analysis/mutation_matrix/ICGC_Breast/norloop_tts/output/ID/icgc_ov.ID28.all") # examplary loading data for exploration
    dm <- melt(d, id.vars = "MutationType")
    dm <- dm[variable != "="]

    # Start to count different categories of indels
    dm_del_count <- dm[grepl("1:Del:", MutationType),
        .(del_count = sum(value)),
        by = "variable"
    ]
    dm_ins_count <- dm[grepl("1:Ins:", MutationType),
        .(ins_count = sum(value)),
        by = "variable"
    ]
    dm_long_del <- dm[grepl("long_Del", MutationType),
        .(long_del = sum(value)),
        by = "variable"
    ]
    dm_long_ins <- dm[grepl("long_Ins", MutationType),
        .(long_ins = sum(value)),
        by = "variable"
    ]
    dm_mh <- dm[grepl("MH", MutationType),
        .(mh = sum(value)),
        by = "variable"
    ]
    dm_complex <- dm[grepl("complex", MutationType),
        .(complex_count = sum(value)),
        by = "variable"
    ]
    dm_count <- merge(dm_del_count, dm_ins_count, by = "variable")
    dm_count <- merge(dm_count, dm_long_del, by = "variable")
    dm_count <- merge(dm_count, dm_long_ins, by = "variable")
    dm_count <- merge(dm_count, dm_mh, by = "variable")
    dm_count <- merge(dm_count, dm_complex, by = "variable")

    # Obtain the total count, note that: sum would sum different columns
    # together
    dm_count[, total_indel := .(del_count + ins_count +
        long_del + long_ins + mh + complex_count)]

    dm_count$adj_del_count <- dm_count$del_count / adjsize * 1e6
    dm_count$adj_ins_count <- dm_count$ins_count / adjsize * 1e6
    dm_count$adj_long_del <- dm_count$long_del / adjsize * 1e6
    dm_count$adj_long_ins <- dm_count$long_ins / adjsize * 1e6
    dm_count$adj_mh <- dm_count$mh / adjsize * 1e6
    dm_count$adj_complex <- dm_count$complex_count / adjsize * 1e6
    dm_count$adj_total <- dm_count$total_indel / adjsize * 1e6
    # dm$adjusted_count <- dm$value / adjsize * 1e6 # Adjust count to per million nucleotides
    # dm_count$type <- label

    ## Adding additional category for testing: long-indel (anything but 1bp indles), long indels no-MH or complex-indel
    dm_count$adj_long_simple_indel <- (dm_count$long_del +
        dm_count$long_ins) / adjsize * 1e6
    dm_count$adj_long_indel <- (dm_count$long_del +
        dm_count$long_ins +
        dm_count$mh) / adjsize * 1e6

    # Select dm fields for casting
    selected_fields <- c("variable", indel_category)
    sdm <- dm_count[, ..selected_fields]
    sdm$sample <- sdm$variable
    sdm$variable <- NULL

    dmt <- melt(sdm, id.vars = "sample")
    names(dmt) <- c("sample", "MutationType", "adjusted_count")
    dmt$type <- label
    return(dmt)
}
# d = get_indel_count_muttype_persample("ICGC_OV/rloop_tss") # Test code

create_test_dataset <- function(sbs.type = "ID28", a.file, b.file,
                                a.label, b.label,
                                a.adjsize, b.adjsize) {
    # Loading two different SBS dataset and then combine before plotting and
    # testing differences.
    a <- get_indel_count_muttype_persample(a.file,
        sbstype = sbs.type,
        adjsize = a.adjsize,
        label = a.label
    )
    b <- get_indel_count_muttype_persample(b.file,
        sbstype = sbs.type,
        adjsize = b.adjsize,
        label = b.label
    )
    both <- rbind(a, b)
    return(both)
}


#---------------------------------------- 
## Section of function to compare rloop vs non-rloop for each genome type
create_indel_testdata_region <- function(label) {
    # Return testdataset for rloop/norloop type: label, sbs.type
    a.label <- paste0("rloop_", label)
    b.label <- paste0("norloop_", label)

    # a.file <- paste0("ICGC_OV/", a.label)
    # b.file <- paste0("ICGC_OV/", b.label)

    a.file <- paste0(project_folder, a.label)
    print(paste0("-- ", "Working on ", a.file))
    b.file <- paste0(project_folder, b.label)
    print(paste0("-- ", "Working on ", b.file))

    a.adjsize <- rloop_type$genome_covered[rloop_type$all == a.label]
    b.adjsize <- rloop_type$genome_covered[rloop_type$all == b.label]

    # Create the merged region of rloop, norloop
    m <- create_test_dataset(
        "ID28",
        a.file, b.file,
        a.label, b.label,
        a.adjsize, b.adjsize
    )

    return(m)
}

# td = create_testdata_region("SBS6","exon")
# table(td$type)

test_muttype_region <- function(testdata,
                                label = "genebody", paired.test = FALSE) {
    # Compare rloop vs norloop relevant regions
    # Return a list/result object from the testdata MutationType by comparing
    # each type in the test data

    # meano <- function(x) mean(na.omit(x)) + 5e-324
    meano <- function(x) {
        y <- ifelse(is.na(x), 0, x)
        return(mean(y) + 5e-324)
    }

    mediano <- function(x) {
        y <- ifelse(is.na(x), 0, x)
        return(median(y) + 5e-324)
    }

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

        rowdt <- data.table(
            MutationType = i, pval = rowd.pvalue, foldChange = rowd.fold,
            rloop_mean = mean_rloop, norloop_mean = mean_norloop
        )
        allrows <- rbind(allrows, rowdt)
    }
    return(allrows)
}

obtain_all_kinds_rloop_tests <- function(paired = FALSE) {
    # Compare between rloop vs non-rloop for different mutationType
    # Return all test results of all the different rlooptypes
    # tumor.
    fr <- list()
    for (i in region_types) {
        tdata <- create_indel_testdata_region(i)
        fr[[i]] <- test_muttype_region(tdata, label = i, paired.test = paired)
        fr[[i]]$region <- i
    }
    # return(fr)
    allregions <- rbindlist(fr)
    return(allregions)
}

save_all_test_result <- function() {
    tr <- obtain_all_kinds_rloop_tests()
    project <- substr(project_folder, 1, nchar(project_folder) - 1)
    outputfile <- paste0(
        "output/indel_compare_rloop_regions/", project,
        "/indel_test_results.tsv"
    )
    fwrite(tr, outputfile, sep = "\t")
}
# tr <- obtain_all_kinds_rloop_tests()
# tr.paired <- obtain_all_kinds_rloop_tests("SBS18", paired = T)

#------------------------------------------------------------
## Section to create plots for simple ID28 indel plots (ggviolin or
## ggboxplot

experiment_plot_single_indel_count <- function(tdata, ptitle) {
    # Return a ggplot object for box plot for type~adjusted_count

    p <- ggboxplot(tdata,
        x = "type", y = "adjusted_count",
        color = "type", palette = "npg",
        # size = 0.1
        alpha = 0.7,
        ylab = "SNV per million",
        title = ptitle,
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
        )
    return(p)
}

saveplot_rloop_types_indel_plot <- function(n.row = 1, y.lim = c(0, 10), ylog10 = F) {
    # Will create the sbs.type plots for each mutationType category comparison
    # of the type (rloop vs norloop) pannels, there will be 4 pannels created
    output_folder <- paste0("output/indel_compare_rloop_regions/", project_folder)
    mkdirp(output_folder)
    filetag <- ""
    smallest_pos_value <- 5e-324

    if (ylog10 == T) {
        filetag <- "log10_"
    }
    pdf(paste0(output_folder, filetag, "indel_compare_rloop_regons.pdf"))

    for (i in region_types) {
        # region_types is a global variable for this script
        tdata <- create_indel_testdata_region(i)
        if (ylog10 == T) {
            tdata$adjusted_count <- tdata$adjusted_count + smallest_pos_value
        }

        print(i) # -- debug/track run progress
        # p <- ggviolin(tdata,
        p <- ggboxplot(tdata,
            x = "type", y = "adjusted_count",
            color = "type", palette = "npg",
            # size = 0.1
            alpha = 0.7,
            ylab = "event per MB",
            title = i,
            add = c("jitter", "median_iqr"),
            # add = c("jitter", "mean_sd"),
            # error.plot="crossbar",
        ) + rremove("x.text") +
            font("x", size = 8) +
            # coord_trans(y = "log10") +
            # ylim(y.lim[1], y.lim[2]) +
            font("legend.text", size = 8) +
            facet_wrap(~MutationType, nrow = n.row) +
            stat_compare_means(aes(group = type),
                method = "wilcox.test",
                label = "..p.format..",
                # label = "..p.signif..",
                # paired = T,
                size = 3
            )

        if (ylog10 == T) {
            p <- p + yscale("log10")
        }
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

    tdata <- create_indel_testdata_region(i)
    p <- ggboxplot(tdata,
        x = "type", y = "adjusted_count",
        color = "type", palette = "npg",
        # size = 0.1
        alpha = 0.7,
        ylab = "events per MB",
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
    p.2 <- plot_rloop_type_compare(sbs.type, "tss", y.lim = ylim, n.row = nrow)
    p.3 <- plot_rloop_type_compare(sbs.type, "tts", y.lim = ylim, n.row = nrow)
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

##############################################################################
explore_wilcox_test <- function() {
    # Some testing code to explore the wilcox test result to find odds ratio
    # etc.
    project_folder <- "ICGC_OV/"
    d <- create_indel_testdata_region("tss")
}

main <- function() {
    # Run different type of indel
    ## -- Note!: when plotting indels, if there are a lot empty indel cases
    ## ---- the plot function currently can not handel it well.
    ## ---- In these situations, I only run the test to generate the test results for the heatmaps.

    #  saveplot_rloop_types_indel_plot()
    #  saveplot_rloop_types_indel_plot(ylog10=T)
    save_all_test_result()
}
main()