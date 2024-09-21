# Analyze the structural variants for rloop regions
rm(list = ls())
library(data.table)
library(magrittr)
library(ipfun)
library(ggpubr)

short_tumor_list <- c(
    "Panc-AdenoCA", "Breast-AdenoCA", "Prost-AdenoCA",
    "Ovary-AdenoCA", "Liver-HCC", "Skin-Melanoma"
)

load_assemble_sv_count <- function(sv.line.file,
                                   sv.count.file,
                                   label = "no_consensus",
                                   region_length) {
    # Assemble the total sv counts per sample and different category of sv
    # count
    a <- fread(sv.line.file)
    # a=fread('output/line-counts_intersected_no_consensus.txt')
    a$sample <- sapply(strsplit(a$V2, ".", fixed = T), `[`, 1)
    names(a) <- c("total_sv", "file", "sample")

    # Cast the different category of SV counts
    d <- fread(sv.count.file)
    # d=fread("output/combined_intersected_no_consensus.tsv")
    d2 <- d[, .N, by = .(V1, V12)]
    d3 <- dcast(d2, V1 ~ V12, value.var = "N", fill = 0)

    # Combine total count and count of each category
    b <- merge(a, d3, by.x = "sample", by.y = "V1", all.x = T)
    b$DEL <- ifelse(is.na(b$DEL), 0, b$DEL)
    b$DUP <- ifelse(is.na(b$DUP), 0, b$DUP)
    b$TRA <- ifelse(is.na(b$TRA), 0, b$TRA)
    b$h2hINV <- ifelse(is.na(b$h2hINV), 0, b$h2hINV)
    b$t2tINV <- ifelse(is.na(b$t2tINV), 0, b$t2tINV)
    b$INV <- b$h2hINV + b$t2tINV

    # Further transformation and melt
    b$total_sv <- as.numeric(b$total_sv)
    m <- melt(b, id.vars = c(1), measure.vars = c(2, 4, 5, 6, 7, 8, 9))
    m$adjusted_value <- m$value / region_length * 1e6
    m$label <- label

    return(m)
}

get_region_length <- function(region_file) {
    # Return the genomic region length
    # r = fread("input/negative_regions_without_consensus_rloop.bed")
    r <- fread(region_file)
    return(sum(r$V3 - r$V2))
}

# Get the different region total length
get_all_region_length <- function() {
    region_length <- list()
    region_length$no_consensus_region <- get_region_length("input/negative_regions_without_consensus_rloop.bed")
    region_length$consensus_region <- get_region_length("input/consensus_no_cutoff.bed")
    region_length$sc100_region <- get_region_length("input/consensus_sc_gt_100.bed")
    region_length$sc200_region <- get_region_length("input/consensus_sc_gt_200.bed")

    # Additional diffferent types of regions: rloop vs no-rloops
    region_length$genebody <- get_region_length("input/output_isaac_score100/rloop_genebody.bed")
    region_length$no_genebody <- get_region_length("input/output_isaac_score100/norloop_genebody.bed")
    region_length$tss <- get_region_length("input/output_isaac_score100/rloop_tss.bed")
    region_length$no_tss <- get_region_length("input/output_isaac_score100/norloop_tss.bed")
    region_length$tts <- get_region_length("input/output_isaac_score100/rloop_tts.bed")
    region_length$no_tts <- get_region_length("input/output_isaac_score100/norloop_tts.bed")

    region_length$tss_sc200 <- get_region_length("input/Rloop_intersected_score200/rloop_tss.bed")
    region_length$no_tss_sc200 <- get_region_length("input/Rloop_intersected_score200/norloop_tss.bed")
    region_length$tts_sc200 <- get_region_length("input/Rloop_intersected_score200/rloop_tts.bed")
    region_length$no_tts_sc200 <- get_region_length("input/Rloop_intersected_score200/norloop_tts.bed")
    return(region_length)
}
all_region_length <- get_all_region_length()


get_all_test_data <- function() {
    # Return a list of differeent groups of testdata
    testdata <- list()

    # Get the sv intersected data from different  regions
    no_consensus <- load_assemble_sv_count(
        "output/line-counts_intersected_no_consensus.txt",
        "output/combined_intersected_no_consensus.tsv",
        "norloop_consensus",
        all_region_length$no_consensus_region
    )
    consensus <- load_assemble_sv_count(
        "output/line-counts_intersected_consensus.txt",
        "output/combined_intersected_consensus.tsv",
        "rloop_consensus",
        all_region_length$consensus_region
    )
    sc100 <- load_assemble_sv_count(
        "output/intersect_all_regions_scores/line-counts_intersect_sc100.txt",
        "output/intersect_all_regions_scores/combined_intersected_intersect_sc100.tsv",
        "rloop_sc100",
        all_region_length$sc100_region
    )

    sc200 <- load_assemble_sv_count(
        "output/intersect_all_regions_scores/line-counts_intersect_sc200.txt",
        "output/intersect_all_regions_scores/combined_intersected_intersect_sc200.tsv",
        "rloop_sc200",
        all_region_length$sc200_region
    )
    # Combine positive sets and negative sets for testing
    testdata$consensus_data <- rbind(consensus, no_consensus)

    # Adapt the sc100 dataset with the right label for norloop_sc100
    no_consensus2 <- no_consensus
    no_consensus2$label <- "norloop_sc100" # The negative control for sc100 is also the same no_consensus regions
    testdata$sc100_data <- rbind(sc100, no_consensus2)

    no_consensus3 <- no_consensus
    no_consensus3$label <- "norloop_sc200" # The negative control for sc100 is also the same no_consensus regions
    testdata$sc200_data <- rbind(sc200, no_consensus3)


    # Assemble other TSS/TTS related region data
    testdata$genebody_data <- rbind(
        load_assemble_sv_count(
            "output/intersect_regions_sc100/line-counts_rloop_genebody.txt",
            "output/intersect_regions_sc100/combined_intersected_rloop_genebody.tsv",
            "rloop_genebody",
            all_region_length$genebody
        ),
        load_assemble_sv_count(
            "output/intersect_regions_sc100/line-counts_norloop_genebody.txt",
            "output/intersect_regions_sc100/combined_intersected_norloop_genebody.tsv",
            "norloop_genebody",
            all_region_length$no_genebody
        )
    )

    testdata$tss_data_score100 <- rbind(
        load_assemble_sv_count(
            "output/intersect_regions_sc100/line-counts_rloop_tss.txt",
            "output/intersect_regions_sc100/combined_intersected_rloop_tss.tsv",
            "rloop_tss",
            all_region_length$tss
        ),
        load_assemble_sv_count(
            "output/intersect_regions_sc100/line-counts_norloop_tss.txt",
            "output/intersect_regions_sc100/combined_intersected_norloop_tss.tsv",
            "norloop_tss",
            all_region_length$no_tss
        )
    )

    testdata$tts_data_score100 <- rbind(
        load_assemble_sv_count(
            "output/intersect_regions_sc100/line-counts_rloop_tts.txt",
            "output/intersect_regions_sc100/combined_intersected_rloop_tts.tsv",
            "rloop_tts",
            all_region_length$tts
        ),
        load_assemble_sv_count(
            "output/intersect_regions_sc100/line-counts_norloop_tts.txt",
            "output/intersect_regions_sc100/combined_intersected_norloop_tts.tsv",
            "norloop_tts",
            all_region_length$no_tts
        )
    )

    testdata$tts_data_score200 <- rbind(
        load_assemble_sv_count(
            "output/intersect_regions_sc200/line-counts_rloop_tts.txt",
            "output/intersect_regions_sc200/combined_intersected_rloop_tts.tsv",
            "rloop_tts",
            all_region_length$tts_sc200
        ),
        load_assemble_sv_count(
            "output/intersect_regions_sc200/line-counts_norloop_tts.txt",
            "output/intersect_regions_sc200/combined_intersected_norloop_tts.tsv",
            "norloop_tts",
            all_region_length$no_tts_sc200
        )
    )

    testdata$tss_data_score200 <- rbind(
        load_assemble_sv_count(
            "output/intersect_regions_sc200/line-counts_rloop_tss.txt",
            "output/intersect_regions_sc200/combined_intersected_rloop_tss.tsv",
            "rloop_tss",
            all_region_length$tss_sc200
        ),
        load_assemble_sv_count(
            "output/intersect_regions_sc200/line-counts_norloop_tss.txt",
            "output/intersect_regions_sc200/combined_intersected_norloop_tss.tsv",
            "norloop_tss",
            all_region_length$no_tss_sc200
        )
    )

    # Return all the assembled data

    return(testdata)
}
svdata <- get_all_test_data()


# p <- ggdensity(sc100_data[variable =="total_sv"],
# x = "adjusted_value", add = "median", rug = T, color = "label", fill = "label", palette = "npg" )
# ggsave("density_total_sv_sc100.png", p)

# levels(sc100$variable) # There are 7 levels.

save_boxplot_compare_sv <- function(p.data, plot.file) {
    # Save the box plots together with p value for the given plot data.
    p <- ggboxplot(p.data,
        x = "label", y = "adjusted_value",
        color = "label", palette = "npg",
        alpha = 0.7,
        ylab = "SV per million",
        add = c("jitter", "median_iqr")
    ) +
        stat_compare_means(aes(group = label),
            method = "wilcox.test",
            label = "..p.format..",
            # label = "..p.signif..",
            # paired = T,
            size = 3
        ) + facet_wrap(~variable, nrow = 3)
    ggsave(plot.file, p)
}

save_violin_compare_sv <- function(p.data, plot.file, plot.fun = ggviolin) {
    # With the given p.data, plot for rloop comparision for each tumor type
    # Loading ICGC tumor sample and type information
    icgc_tumor_type <- fread("input/nature2020_panCancerWGS_suppl_table1.txt")
    icgc <- icgc_tumor_type[, .(tumour_specimen_aliquot_id, histology_abbreviation)]
    names(icgc) <- c("sample", "tumor")

    # Select only total SV for the violin plot
    p.data.total_sv <- p.data[variable == "total_sv"]
    plot.data <- merge(p.data.total_sv, icgc, by = "sample")
    # plot.data <- plot.data[tumor %in% short_tumor_list] # only limiting to the plot tumor list
    # plot.data$tumor <- factor(plot.data$tumor, levels = short_tumor_list)
    plot.data$label <- sapply(strsplit(plot.data$label, "_", fixed = T), `[`, 1)
    plot.data$label <- factor(plot.data$label, levels = c("rloop", "norloop"))

    # Save the box plots together with p value for the given plot data.
    p <- plot.fun(plot.data,
        # x = "tumor", y = "adjusted_value",
        x = "label", y = "adjusted_value",
        color = "label", palette = "npg",
        alpha = 0.7,
        ylab = "SV per million BP",
        add = c("jitter", "median_iqr")
    ) +
        stat_compare_means(aes(group = label),
            method = "wilcox.test",
            label = "..p.format..",
            # label = "..p.signif..",
            # paired = T,
            size = 3
        ) + facet_wrap(~region, ncol = 2)

    mkdirp("output/violin_plot/")
    pfile <- paste0("output/violin_plot/", plot.file)
    ggsave(pfile, p, width = 10, height = 5)
}

get_icgc_tumor_types <- function() {
    # Return a list tumor types of ICGC table
    tumor_types <- list()

    # Loading cancer information
    icgc <- fread("input/nature2020_panCancerWGS_suppl_table1.txt")
    # table(icgc$project_code)
    # table(icgc$histology_abbreviation)
    tumor_types$icgc <- icgc

    tumor_types$breast <- icgc[histology_abbreviation == "Breast-AdenoCA"]
    tumor_types$ov <- icgc[histology_abbreviation == "Ovary-AdenoCA"]
    tumor_types$panc <- icgc[histology_abbreviation == "Panc-AdenoCA"]
    tumor_types$prostate <- icgc[histology_abbreviation == "Prost-AdenoCA"]
    tumor_types$melanoma <- icgc[histology_abbreviation == "Skin-Melanoma"]
    tumor_types$liver <- icgc[histology_abbreviation == "Liver-HCC"]

    # Return the tumor types
    return(tumor_types)
}
icgc_tumors <- get_icgc_tumor_types()

# Filter tumor type specific data
filter_tumor_type <- function(x, tumor) {
    # filter from  the x, for only tumor$tumour_specimen_aliquot_id
    sel <- x[sample %in% tumor$tumour_specimen_aliquot_id]
    return(sel)
}

save_tumor_boxplots <- function(tumor.data, tumor.name,
                                svdata = svdata, figpath_prefix = "plots/all_tumors/boxplot_") {
    # Using the specific tumor types data, and save the plot with tumor.name
    sc100.name <- paste0(figpath_prefix, tumor.name, "_sv_sc100.png")
    consensus.name <- paste0(figpath_prefix, tumor.name, "_sv_consensus.png")
    genebody.name <- paste0(figpath_prefix, tumor.name, "_sv_genebody.png")
    tss.name <- paste0(figpath_prefix, tumor.name, "_sv_tss_score100.png")
    tts.name <- paste0(figpath_prefix, tumor.name, "_sv_tts_score100.png")
    tss200.name <- paste0(figpath_prefix, tumor.name, "_sv_tss_score200.png")
    tts200.name <- paste0(figpath_prefix, tumor.name, "_sv_tts_score200.png")


    save_boxplot_compare_sv(
        filter_tumor_type(svdata$sc100_data, tumor.data),
        sc100.name
    )

    save_boxplot_compare_sv(
        filter_tumor_type(svdata$consensus_data, tumor.data),
        consensus.name
    )

    save_boxplot_compare_sv(
        filter_tumor_type(svdata$genebody_data, tumor.data),
        genebody.name
    )
    save_boxplot_compare_sv(
        filter_tumor_type(svdata$tss_data_score100, tumor.data),
        tss.name
    )
    save_boxplot_compare_sv(
        filter_tumor_type(svdata$tts_data_score100, tumor.data),
        tts.name
    )
    save_boxplot_compare_sv(
        filter_tumor_type(svdata$tss_data_score200, tumor.data),
        tss200.name
    )
    save_boxplot_compare_sv(
        filter_tumor_type(svdata$tts_data_score200, tumor.data),
        tts200.name
    )
}


test_muttype_region <- function(testdata, flabel = "tss", paired.test = FALSE) {
    # Compare rloop vs norloop relevant regions
    # Return a list/result object from the testdata MutationType by comparing
    # each type in the test data

    meano <- function(x) mean(na.omit(x)) + 5e-324
    mediano <- function(x) median(na.omit(x)) + 5e-324

    allrows <- data.table() # keep all rows as value to return

    for (i in unique(testdata$variable)) {
        # Obtain each type of SV data for testing/calculations
        subdata <- testdata[variable == i]
        rloop <- subdata[label == paste0("rloop_", flabel)]
        norloop <- subdata[label == paste0("norloop_", flabel)]
        rowd.pvalue <- wilcox.test(adjusted_value ~ label,
            data = subdata
        )$p.value
        # alternative = "less")$p.value
        mean_rloop <- meano(rloop$adjusted_value)
        mean_norloop <- meano(norloop$adjusted_value)
        rowd.fold <- mean_rloop / mean_norloop

        # Obtain the median count/ratio
        median_rloop <- mediano(rloop$adjusted_value)
        median_norloop <- mediano(norloop$adjusted_value)

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

save_tumor_sv_test_result <- function(tumor.data, tumor.name) {
    # Save the given tumor.data, tumor.name into the output/sv_test_results/
    # folder.
    mkdirp("output/sv_test_results")

    sc100.name <- paste0("output/sv_test_results/", tumor.name, "_sv_sc100.csv")
    sc200.name <- paste0("output/sv_test_results/", tumor.name, "_sv_sc200.csv")
    consensus.name <- paste0("output/sv_test_results/", tumor.name, "_sv_consensus.csv")

    # For the following use tss-sc100 so that when split strings, tss-sc100
    # will tbe taken as an unit together
    tss100.name <- paste0("output/sv_test_results/", tumor.name, "_sv_tss-sc100.csv")
    tts100.name <- paste0("output/sv_test_results/", tumor.name, "_sv_tts-sc100.csv")
    tss200.name <- paste0("output/sv_test_results/", tumor.name, "_sv_tss-sc200.csv")
    tts200.name <- paste0("output/sv_test_results/", tumor.name, "_sv_tts-sc200.csv")

    filter_tumor_type(svdata$sc100_data, tumor.data) %>%
        test_muttype_region("sc100") %>%
        fwrite(sc100.name)

    filter_tumor_type(svdata$sc200_data, tumor.data) %>%
        test_muttype_region("sc200") %>%
        fwrite(sc200.name)

    filter_tumor_type(svdata$consensus_data, tumor.data) %>%
        test_muttype_region("consensus") %>%
        fwrite(consensus.name)

    filter_tumor_type(svdata$tss_data_score100, tumor.data) %>%
        test_muttype_region("tss") %>%
        fwrite(tss100.name)

    filter_tumor_type(svdata$tts_data_score100, tumor.data) %>%
        test_muttype_region("tts") %>%
        fwrite(tts100.name)

    filter_tumor_type(svdata$tss_data_score200, tumor.data) %>%
        test_muttype_region("tss") %>%
        fwrite(tss200.name)

    filter_tumor_type(svdata$tts_data_score200, tumor.data) %>%
        test_muttype_region("tts") %>%
        fwrite(tts200.name)
}

###########################################################
main_boxplot_different_tumors <- function() {
    # Main execution part of the script
    mkdirp("plots/all_tumors/")
    save_boxplot_compare_sv(svdata$sc100_data, "plots/all_tumors/boxplot_sv_sc100.png")
    save_boxplot_compare_sv(svdata$sc200_data, "plots/all_tumors/boxplot_sv_sc200.png")
    save_boxplot_compare_sv(svdata$consensus_data, "plots/all_tumors/boxplot_sv_consensus.png")
    save_boxplot_compare_sv(svdata$genebody_data, "plots/all_tumors/boxplot_sv_genebody.png")
    save_boxplot_compare_sv(svdata$tss_data_score100, "plots/all_tumors/boxplot_sv_tss_score100.png")
    save_boxplot_compare_sv(svdata$tts_data_score100, "plots/all_tumors/boxplot_sv_tts_score100.png")
    save_boxplot_compare_sv(svdata$tss_data_score200, "plots/all_tumors/boxplot_sv_tss_score200.png")
    save_boxplot_compare_sv(svdata$tts_data_score200, "plots/all_tumors/boxplot_sv_tts_score200.png")

    save_tumor_boxplots(icgc_tumors$breast, "breast", svdata = svdata)
    save_tumor_boxplots(icgc_tumors$ov, "ov", svdata = svdata)
    save_tumor_boxplots(icgc_tumors$panc, "pancreatic", svdata = svdata)
    save_tumor_boxplots(icgc_tumors$prostate, "prostate", svdata = svdata)
    save_tumor_boxplots(icgc_tumors$melanoma, "melanoma", svdata = svdata)
}

main_generate_tumors_test_result <- function() {
    # Testing result using different cutoff
    # a1 <- test_muttype_region(svdata$sc100_data, "sc100")
    # a2 <- test_muttype_region(svdata$sc200_data, "sc200")
    # b1 <- test_muttype_region(svdata$tts_data_score100, "tts")
    # b2 <- test_muttype_region(svdata$tts_data_score200, "tts")

    save_tumor_sv_test_result(icgc_tumors$breast, "breast")
    save_tumor_sv_test_result(icgc_tumors$ov, "ov")
    save_tumor_sv_test_result(icgc_tumors$panc, "pancreatic")
    save_tumor_sv_test_result(icgc_tumors$prostate, "prostate")
    save_tumor_sv_test_result(icgc_tumors$melanoma, "melanoma")
    save_tumor_sv_test_result(icgc_tumors$liver, "liver")
    save_tumor_sv_test_result(icgc_tumors$icgc, "icgc")
}

# main_violin_plot_selected_tumors <- function() {
main <- function() {
    # save_violin_compare_sv(svdata$tss_data_score200, "violin_tss_score200_tumors.png")
    # save_violin_compare_sv(svdata$tts_data_score200, "violin_tts_score200_tumors.png")
    # save_violin_compare_sv(svdata$sc200_data,
    #     "violin_consensus_score200_tumors.png",
    #     plot.fun = ggboxplot
    # )
    tss <- svdata$tss_data_score200
    tts <- svdata$tts_data_score200
    # adding region type field
    tss$region <- "TSS"
    tts$region <- "TTS"

    # Combine tss and tts
    tss_and_tts <- rbindlist(list(tss, tts))
    save_violin_compare_sv(tss_and_tts,
        "boxplot_tss_and_tts_score200_tumors.png",
        plot.fun = ggboxplot
    )

    # Loading tumor info
    icgc_tumor_type <- fread("input/nature2020_panCancerWGS_suppl_table1.txt")
    icgc_tumor_type$tumor_type <- icgc_tumor_type$histology_tier2
    t.type <- icgc_tumor_type[, .(tumour_specimen_aliquot_id, tumor_type)]

    tss_and_tts <- merge(tss_and_tts, t.type, by.x = "sample", by.y = "tumour_specimen_aliquot_id", all.x = T)

    return(tss_and_tts)
}


save_plot_ov_biallelic <- function(input, figfile, tumor.type = "OV") {
    # Plot biallelic vs. control for ICGC ovarian tumors only
    ov <- fread("../../input/biallelic_sample_lists/icgc_ov_biallelic_cat.csv")
    ov.biallelic <- ov[cat == "Biallelic"]
    ov.control <- ov[cat == "control"]
    fd <- input
    # fd = icgc_tumors$ov
    # fd <- svdata$sc200_data
    fd$allelic_status <- ifelse(fd$sample %in% ov.biallelic$tumour_specimen_aliquot_id,
        "biallelic", ifelse(fd$sample %in% ov.control$tumour_specimen_aliquot_id,
            "control", "unknown"
        )
    )
    fd <- fd[allelic_status != "unknown"]
    fd$rloop <- sapply(strsplit(fd$label, "_", fixed = T), `[`, 1)
    fd$rloop <- factor(fd$rloop, levels = c("rloop", "norloop"))

    p1 <- ggboxplot(fd,
        x = "allelic_status", y = "adjusted_value",
        color = "rloop", palette = "npg",
        size = 0.3,
        alpha = 0.5,
        ylab = "Adjusted SNV per million",
        title = paste0(tumor.type, " Biallelic vs Control"),
        add = c("median_iqr", "jitter"),
    ) + font("x", size = 8) +
        # yscale("log10", .format = T) +
        # ylim(y.lim[1], y.lim[2]) +
        font("legend.text", size = 6) +
        # theme(axis.text.x = element_text(angle = 90)) +
        stat_compare_means(aes(group = rloop),
            method = "wilcox.test",
            label = "..p.format..",
            # label = "..p.signif..",
            # paired = T,
            size = 3
        ) + facet_wrap(~variable, ncol = 1, scales = "free")

    p2 <- ggboxplot(fd,
        x = "rloop", y = "adjusted_value",
        color = "allelic_status",
        palette = "jco",
        size = 0.3,
        alpha = 0.5,
        ylab = "Adjusted SNV per million",
        title = paste0(tumor.type, " Rloops vs Non-Rloops"),
        add = c("median_iqr", "jitter"),
    ) + font("x", size = 8) +
        # yscale("log10", .format = T) +
        # ylim(y.lim[1], y.lim[2]) +
        font("legend.text", size = 6) +
        # theme(axis.text.x = element_text(angle = 90)) +
        stat_compare_means(aes(group = allelic_status),
            method = "wilcox.test",
            label = "..p.format..",
            # label = "..p.signif..",
            # paired = T,
            size = 3
        ) + facet_wrap(~variable, ncol = 1, scales = "free")

    p <- ggarrange(p1, p2, ncol = 2)

    # Make the folder and save tests result
    mkdirp("output/ov_biallelic_tests/")
    ggsave(paste0("output/ov_biallelic_tests/", figfile), p, width = 10, height = 10)
}

main_biallelic <- function() {
    # main_generate_tumors_test_result()
    # main_violin_plot_selected_tumors()

    ## Save the biallelic comparison plot for OV tumors
    save_plot_ov_biallelic(svdata$tss_data_score200, "OV_tss_score200.png", "OV TSS")
    save_plot_ov_biallelic(svdata$tts_data_score200, "OV_tts_score200.png", "OV TTS")
    save_plot_ov_biallelic(svdata$sc200_data, "OV_sc200.png", "OV sc200")
}
# main_violin_plot_selected_tumors()
# a <- main()
# fwrite(a, "/data/projects/peix/rloop_project/output/fig_plot_data/SV_counts_all-tumors.tsv")

## Save the SV breakpoints data of consensus Rloop vs non-Rloop region for analysis with RAD52 comparison
# fwrite(svdata$consensus_data, "output/rloop-consensus_sv-data_for_RAD52_comparison.tsv")
# fwrite(svdata$sc100_data, "output/rloop-sc100_sv-data_for_RAD52_comparison.tsv")
# fwrite(svdata$sc200_data, "output/rloop-sc200_sv-data_for_RAD52_comparison.tsv")


save_plot_pcawg_hrd <- function(input, figfile, tumor.type = "OV") {
    ## Get nature 2020 HRD paper HRD status info
    h <- fread("../../input/2020_nature_com_hrd_chord_paper_data/CHORD_HRD_status_2020_Nature_Comm.txt")
    h <- h[group == "PCAWG"]
    h.hrd <- h[hr_status == "HR_deficient"]
    h.hrp <- h[hr_status == "HR_proficient"]


    # Plot biallelic vs. control for ICGC ovarian tumors only
    ov <- fread("../../input/biallelic_sample_lists/icgc_ov_biallelic_cat.csv")
    ov.biallelic <- ov[cat == "Biallelic"]
    ov.control <- ov[cat == "control"]
    fd <- input
    # fd = icgc_tumors$ov
    # fd <- svdata$sc200_data

    fd$HR_status <- ifelse(fd$sample %in% h.hrd$sample,
        "HR_deficient", ifelse(fd$sample %in% h.hrp$sample,
            "HR_proficient", "unknown"
        )
    )
    fd <- fd[HR_status != "unknown"]
    fd$rloop <- sapply(strsplit(fd$label, "_", fixed = T), `[`, 1)
    fd$rloop <- factor(fd$rloop, levels = c("rloop", "norloop"))

    p1 <- ggboxplot(fd,
        x = "HR_status", y = "adjusted_value",
        color = "rloop", palette = "npg",
        size = 0.3,
        alpha = 0.5,
        ylab = "SV breakpoint density",
        title = paste0(tumor.type, " HR proficient tumors vs deficient"),
        add = c("median_iqr", "jitter"),
    ) + font("x", size = 8) +
        # yscale("log10", .format = T) +
        # ylim(y.lim[1], y.lim[2]) +
        font("legend.text", size = 6) +
        # theme(axis.text.x = element_text(angle = 90)) +
        stat_compare_means(aes(group = rloop),
            method = "wilcox.test",
            label = "..p.format..",
            # label = "..p.signif..",
            # paired = T,
            size = 3
        ) + facet_wrap(~variable, ncol = 1, scales = "free")

    p2 <- ggboxplot(fd,
        x = "rloop", y = "adjusted_value",
        color = "HR_status",
        palette = "jco",
        size = 0.3,
        alpha = 0.5,
        ylab = "SV breakpoint density",
        title = paste0(tumor.type, " Rloops vs Non-Rloops"),
        add = c("median_iqr", "jitter"),
    ) + font("x", size = 8) +
        # yscale("log10", .format = T) +
        # ylim(y.lim[1], y.lim[2]) +
        font("legend.text", size = 6) +
        # theme(axis.text.x = element_text(angle = 90)) +
        stat_compare_means(aes(group = HR_status),
            method = "wilcox.test",
            label = "..p.format..",
            # label = "..p.signif..",
            # paired = T,
            size = 3
        ) + facet_wrap(~variable, ncol = 1, scales = "free")

    p <- ggarrange(p1, p2, ncol = 2)

    # Make the folder and save tests result
    mkdirp("output/PCAWG_HR_status_tests/")
    ggsave(paste0("output/PCAWG_HR_status_tests/", figfile), p, width = 10, height = 10)
    ggsave(paste0("output/PCAWG_HR_status_tests/rloop_hr_compare_", figfile), p2, width = 10, height = 10)
}

main_pcwag_hrd <- function() {
    ## Save the HRD comparision within PCAWG tumors
    save_plot_pcawg_hrd(svdata$tss_data_score200, "PCAWG_tss_score200.png", "PCAWG TSS")
    save_plot_pcawg_hrd(svdata$tts_data_score200, "PCAWG_tts_score200.png", "PCAWG TTS")
    save_plot_pcawg_hrd(svdata$sc200_data, "PCAWG_sc200.png", "PCAWG sc200")
}

# main_pcwag_hrd()

## Save boxplots of SV comparision for specific tumors using non-HRD tumors
main_boxplot_nohrd_tumors <- function() {
    # Using the hrd/hrp icgc paper to limit the SVData to only HRP tumors

    ## Get nature 2020 HRD paper HRD status info
    h <- fread("../../input/2020_nature_com_hrd_chord_paper_data/CHORD_HRD_status_2020_Nature_Comm.txt")
    h <- h[group == "PCAWG"]
    h.hrd <- h[hr_status == "HR_deficient"]
    h.hrp <- h[hr_status == "HR_proficient"]

    ## Subset svdata to svdata.hrp (HRP tumors)
    svdata.hrp <- svdata
    for (i in names(svdata)) {
        svdata.hrp[[i]] <- svdata[[i]][sample %in% h.hrp$sample, ]
    }
    # set svdata to svdata.hrp
    svdata <- svdata.hrp

    mkdirp("plots")
    save_boxplot_compare_sv(svdata$sc100_data, "plots/hrp_boxplot_sv_sc100.png")
    save_boxplot_compare_sv(svdata$sc200_data, "plots/hrp_boxplot_sv_sc200.png")
    save_boxplot_compare_sv(svdata$consensus_data, "plots/hrp_boxplot_sv_consensus.png")
    save_boxplot_compare_sv(svdata$genebody_data, "plots/hrp_boxplot_sv_genebody.png")
    save_boxplot_compare_sv(svdata$tss_data_score100, "plots/hrp_boxplot_sv_tss_score100.png")
    save_boxplot_compare_sv(svdata$tts_data_score100, "plots/hrp_boxplot_sv_tts_score100.png")
    save_boxplot_compare_sv(svdata$tss_data_score200, "plots/hrp_boxplot_sv_tss_score200.png")
    save_boxplot_compare_sv(svdata$tts_data_score200, "plots/hrp_boxplot_sv_tts_score200.png")

    # Save plots using updated function
    save_tumor_boxplots(icgc_tumors$breast, "breast", svdata = svdata.hrp, figpath_prefix = "plots/boxplot_hrp_")
    save_tumor_boxplots(icgc_tumors$ov, "ov", svdata = svdata.hrp, figpath_prefix = "plots/boxplot_hrp_")
    save_tumor_boxplots(icgc_tumors$panc, "pancreatic", svdata = svdata.hrp, figpath_prefix = "plots/boxplot_hrp_")
    save_tumor_boxplots(icgc_tumors$prostate, "prostate", svdata = svdata.hrp, figpath_prefix = "plots/boxplot_hrp_")
    save_tumor_boxplots(icgc_tumors$melanoma, "melanoma", svdata = svdata.hrp, figpath_prefix = "plots/boxplot_hrp_")

    return(svdata.hrp)
}

# a <- main_boxplot_nohrd_tumors()
main_boxplot_different_tumors()
