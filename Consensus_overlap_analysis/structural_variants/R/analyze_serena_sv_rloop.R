# Analyze the structural variants for rloop regions
# -- Adapted for Serena ERpos and TrpNeg tumors
rm(list = ls())
library(data.table)
library(magrittr)
library(ipfun)
library(ggpubr)

get_serena_samples <- function() {
    # Loading Serena breast tumor type info
    sample <- fread("../../input/bed-version-2021-05-08/BRCA-EU/serena_nature_2016_560_breast_tumor_table1_clindata.csv")
    sample$sampleid <- gsub("(a|b).*$", "", sample$sample_name)

    # This is the cohort we are using for final testing 320 ER+/HER2-
    er_pos_sample <- sample[final.ER == "positive" & final.HER2 == "negative"]
    erPosSamples <- unique(er_pos_sample$sampleid)
    trp_neg_sample <- sample[final.ER == "negative" &
        final.HER2 == "negative" &
        final.PR == "negative"]
    trpNegSamples <- unique(trp_neg_sample$sampleid)

    rl <- list()
    rl$erPosSamples <- erPosSamples
    rl$trpNegSamples <- trpNegSamples
    return(rl)
}
serena <- get_serena_samples()

load_assemble_sv_count <- function(samplelist,
                                   sv.count.file,
                                   label = "no_consensus",
                                   region_length) {
    # Assemble the total sv counts per sample and different category of sv
    # count
    # -- Adapted for Serena's SV input, input is the list of sampleids
    # samplelist = serena$erPosSamples
    a <- data.table(sample = samplelist)

    # Cast the different category of SV counts
    d <- fread(sv.count.file)
    # d = fread("output/intersect_serena_erpos_rloops/no_consensus.bedpe")
    d2 <- d[, .N, by = .(V13, V11)] # V13 is the sample, V11 is the SV type
    d3 <- dcast(d2, V13 ~ V11, value.var = "N", fill = 0)

    # Combine total count and count of each category
    # -- These change slightly from the icgc consensus data as only 1 cat. inv
    b <- merge(a, d3, by.x = "sample", by.y = "V13", all.x = T)
    b$DEL <- ifelse(is.na(b$DEL), 0, b$DEL)
    b$DUP <- ifelse(is.na(b$DUP), 0, b$DUP)
    b$TRA <- ifelse(is.na(b$TRA), 0, b$TRA)
    b$INV <- ifelse(is.na(b$INV), 0, b$INV)
    # Further transformation and melt
    b$total_sv <- as.numeric(b$DEL + b$DUP + b$TRA + b$INV)

    m <- melt(b, id.vars = c(1), measure.vars = c(2, 3, 4, 5, 6))
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


get_all_test_data <- function(tumor_sample_list = serena$erPosSamples,
                              tumor_intersected_folder = "output/intersect_serena_erpos_rloops/") {
    # Return a list of differeent groups of testdata
    testdata <- list()

    local_load_assemble_sv_count <- function(bedpe_name, label, adjust_region_length) {
        # wrapper function to return the sv data with region adjustment and
        # tag
        fd <- load_assemble_sv_count(
            tumor_sample_list,
            paste0(tumor_intersected_folder, bedpe_name),
            label,
            adjust_region_length
        )
        return(fd)
    }

    # Get the sv intersected data from different  regions
    no_consensus <- local_load_assemble_sv_count(
        "no_consensus.bedpe", "norloop_consensus",
        all_region_length$no_consensus_region
    )
    consensus <- local_load_assemble_sv_count(
        "consensus_no-cutoff.bedpe", "rloop_consensus",
        all_region_length$consensus_region
    )
    sc100 <- local_load_assemble_sv_count(
        "consensus_sc100.bedpe", "rloop_sc100",
        all_region_length$sc100_region
    )
    sc200 <- local_load_assemble_sv_count(
        "consensus_sc200.bedpe", "rloop_sc200",
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
        local_load_assemble_sv_count(
            "rloop_genebody-sc100.bedpe", "rloop_genebody",
            all_region_length$genebody
        ),
        local_load_assemble_sv_count(
            "norloop_genebody-sc100.bedpe", "norloop_genebody",
            all_region_length$no_genebody
        )
    )

    testdata$tss_data_score100 <- rbind(
        local_load_assemble_sv_count(
            "rloop_tss-sc100.bedpe", "rloop_tss",
            all_region_length$tss
        ),
        local_load_assemble_sv_count(
            "norloop_tss-sc100.bedpe", "norloop_tss",
            all_region_length$no_tss
        )
    )

    testdata$tts_data_score100 <- rbind(
        local_load_assemble_sv_count(
            "rloop_tts-sc100.bedpe", "rloop_tts",
            all_region_length$tts
        ),
        local_load_assemble_sv_count(
            "norloop_tts-sc100.bedpe", "norloop_tts",
            all_region_length$no_tts
        )
    )

    testdata$tss_data_score200 <- rbind(
        local_load_assemble_sv_count(
            "rloop_tss-sc200.bedpe", "rloop_tss",
            all_region_length$tss_sc200
        ),
        local_load_assemble_sv_count(
            "norloop_tss-sc200.bedpe", "norloop_tss",
            all_region_length$no_tss_sc200
        )
    )

    testdata$tts_data_score200 <- rbind(
        local_load_assemble_sv_count(
            "rloop_tts-sc200.bedpe", "rloop_tts",
            all_region_length$tts_sc200
        ),
        local_load_assemble_sv_count(
            "norloop_tts-sc200.bedpe", "norloop_tts",
            all_region_length$no_tts_sc200
        )
    )

    # Return all the assembled data
    return(testdata)
}

# Obtain ER and triple negative SV intersected data for different groups
er.svdata <- get_all_test_data()
tri.svdata <- get_all_test_data(
    serena$trpNegSamples,
    "output/intersect_serena_trpneg_rloops/"
)


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

save_tumor_boxplots <- function(tumor.data, tumor.name) {
    # Using the specific tumor types data, and save the plot with tumor.name
    sc100.name <- paste0("output/serena_tumor_boxplots/boxplot_", tumor.name, "_sv_sc100.png")
    consensus.name <- paste0("output/serena_tumor_boxplots/boxplot_", tumor.name, "_sv_consensus.png")
    genebody.name <- paste0("output/serena_tumor_boxplots/boxplot_", tumor.name, "_sv_genebody.png")
    tss100.name <- paste0("output/serena_tumor_boxplots/boxplot_", tumor.name, "_sv_tss_sc100.png")
    tts100.name <- paste0("output/serena_tumor_boxplots/boxplot_", tumor.name, "_sv_tts_sc100.png")
    tss200.name <- paste0("output/serena_tumor_boxplots/boxplot_", tumor.name, "_sv_tss_sc200.png")
    tts200.name <- paste0("output/serena_tumor_boxplots/boxplot_", tumor.name, "_sv_tts_sc200.png")

    save_boxplot_compare_sv(
        tumor.data$sc100_data,
        sc100.name
    )

    save_boxplot_compare_sv(
        tumor.data$consensus_data,
        consensus.name
    )

    save_boxplot_compare_sv(
        tumor.data$genebody_data,
        genebody.name
    )
    save_boxplot_compare_sv(
        tumor.data$tss_data_score100,
        tss100.name
    )
    save_boxplot_compare_sv(
        tumor.data$tts_data_score100,
        tts100.name
    )
    save_boxplot_compare_sv(
        tumor.data$tss_data_score200,
        tss200.name
    )
    save_boxplot_compare_sv(
        tumor.data$tts_data_score200,
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
    allrows$test_type <- "rloop-vs-norloop"
    return(allrows)
}

test_muttype_allelic <- function(testdata, flabel = "tss", paired.test = FALSE) {
    # Compare the allelic status for the input testdata
    # -- the input data is assumed to be the same type: such as only rloop region breakpoint counts
    # Return a list/result object from the testdata MutationType by comparing
    # each type in the test data

    # Assign biallelic vs control status
    fd <- testdata
    fd$allelic_status <- ifelse(fd$sample %in% davies$biallelic,
        "biallelic", ifelse(fd$sample %in% davies$control, "control", "unknown")
    )
    fd <- fd[allelic_status != "unknown"]

    meano <- function(x) mean(na.omit(x)) + 5e-324
    mediano <- function(x) median(na.omit(x)) + 5e-324

    allrows <- data.table() # keep all rows as value to return

    for (i in unique(fd$variable)) {
        # Obtain each type of SV data for testing/calculations
        subdata <- fd[variable == i]
        biallelic <- subdata[allelic_status == "biallelic"]
        control <- subdata[allelic_status == "control"]
        rowd.pvalue <- wilcox.test(adjusted_value ~ allelic_status,
            data = subdata
        )$p.value
        # alternative = "less")$p.value
        mean_biallelic <- meano(biallelic$adjusted_value)
        mean_control <- meano(control$adjusted_value)
        rowd.fold <- mean_biallelic / mean_control

        # Obtain the median count/ratio
        median_biallelic <- mediano(biallelic$adjusted_value)
        median_control <- mediano(control$adjusted_value)

        rowd.median.fold <- median_biallelic / median_control
        rowdt <- data.table(
            MutationType = i, pval = rowd.pvalue, foldChange = rowd.fold,
            biallelic_mean = mean_biallelic, control_mean = mean_control,
            biallelic_median = median_biallelic, control_median = median_control,
            foldChangeMedian = rowd.median.fold
        )
        allrows <- rbind(allrows, rowdt)
    }
    allrows$test_type <- "biallelic-vs-control"
    return(allrows)
}

test_muttype_allelic_compare_rloop_region <- function(testdata, label = "tss", group = "tss_data") {
    # with the given testdata, separate it into two groups of biallelic vs. control
    # For each of the groups: generate the test results, and then label with either biallelic or control
    # Return the combined two groups results

    # Assign biallelic vs control status
    fd <- testdata
    fd$allelic_status <- ifelse(fd$sample %in% davies$biallelic,
        "biallelic", ifelse(fd$sample %in% davies$control, "control", "unknown")
    )
    fd <- fd[allelic_status != "unknown"]

    fd.biallelic <- fd[allelic_status == "biallelic"]
    fd.control <- fd[allelic_status == "control"]

    fr.biallelic <- test_muttype_region(fd.biallelic, label)
    ## test_group is the field indicating which group the comparion is done with
    fr.biallelic$test_group <- "biallelic"
    fr.control <- test_muttype_region(fd.control, label)
    fr.control$test_group <- "control"
    fr <- rbind(fr.biallelic, fr.control)
    fr$group <- group
    return(fr)
}


test_muttype_rloop_region_compare_allelic <- function(testdata, flabel = "tss", group = "tss_data_score200") {
    # With the given testdata, separate testdata into rloop and nonrloop,
    # -- and then in each region type, compare biallelic vs. control
    # -- then combined the final results with right label and return

    # Separate data into rloop and nonrloop
    fd <- testdata
    rloop <- fd[label == paste0("rloop_", flabel)]
    norloop <- fd[label == paste0("norloop_", flabel)]

    fr.rloop <- test_muttype_allelic(rloop, label)
    ## test_group indicates within which group the test is done
    fr.rloop$test_group <- "rloop"
    fr.norloop <- test_muttype_allelic(norloop, label)
    fr.norloop$test_group <- "norloop"

    fr <- rbind(fr.rloop, fr.norloop)
    fr$group <- group
    return(fr)
}


save_tumor_sv_test_result <- function(tumor.data, tumor.name) {
    # With the serena's ER and triple negative breast tumor data as separte
    # data object, here the function utilization is different:
    # - tumor.data is the data directly provided
    # - tumor.name is just for labeling output data.
    # Save the given tumor.data, tumor.name into the output/serena_sv_test_results/
    # folder.
    mkdirp("output/serena_sv_test_results")

    sc100.name <- paste0("output/serena_sv_test_results/", tumor.name, "_sv_sc100.csv")
    sc200.name <- paste0("output/serena_sv_test_results/", tumor.name, "_sv_sc200.csv")
    consensus.name <- paste0("output/serena_sv_test_results/", tumor.name, "_sv_consensus.csv")

    # For the following use tss-sc100 so that when split strings, tss-sc100
    # will tbe taken as an unit together
    tss100.name <- paste0("output/serena_sv_test_results/", tumor.name, "_sv_tss-sc100.csv")
    tts100.name <- paste0("output/serena_sv_test_results/", tumor.name, "_sv_tts-sc100.csv")
    tss200.name <- paste0("output/serena_sv_test_results/", tumor.name, "_sv_tss-sc200.csv")
    tts200.name <- paste0("output/serena_sv_test_results/", tumor.name, "_sv_tts-sc200.csv")

    tumor.data$sc100_data %>%
        test_muttype_region("sc100") %>%
        fwrite(sc100.name)

    tumor.data$sc200_data %>%
        test_muttype_region("sc200") %>%
        fwrite(sc200.name)

    tumor.data$consensus_data %>%
        test_muttype_region("consensus") %>%
        fwrite(consensus.name)

    tumor.data$tss_data_score100 %>%
        test_muttype_region("tss") %>%
        fwrite(tss100.name)

    tumor.data$tts_data_score100 %>%
        test_muttype_region("tts") %>%
        fwrite(tts100.name)

    tumor.data$tss_data_score200 %>%
        test_muttype_region("tss") %>%
        fwrite(tss200.name)

    tumor.data$tts_data_score200 %>%
        test_muttype_region("tts") %>%
        fwrite(tts200.name)
}

###########################################################
main_boxplots <- function() {
    # This block is the old boxplot for icgc tumors- all tumors; NA for
    # Serena's breast tumor data
    # save_boxplot_compare_sv(sc100_data, "boxplot_sv_sc100.png")
    # save_boxplot_compare_sv(consensus_data, "boxplot_sv_consensus.png")
    # save_boxplot_compare_sv(svdata$genebody_data, "boxplot_sv_genebody.png")
    # save_boxplot_compare_sv(svdata$tss_data, "boxplot_sv_tss.png")
    # save_boxplot_compare_sv(svdata$tts_data, "boxplot_sv_tts.png")


    save_tumor_boxplots(er.svdata, "ER-Pos")
    save_tumor_boxplots(tri.svdata, "Trp-Neg")
}

main_save_test_results <- function() {
    # Testing result using different cutoff

    save_tumor_sv_test_result(er.svdata, "ER-Pos")
    save_tumor_sv_test_result(tri.svdata, "Trp-Neg")
}

## Using the biallelic information to plot
davies <- readRDS("../../input/biallelic_sample_lists/Davies_Serena_biallelic_samples.rds")
# davies <- readRDS("input/biallelic_sample_lists/Davies_Serena_biallelic_samples.rds")
# fd <- er.svdata$tss_data_score200

save_plot_compare_biallelic <- function(input, figfile, tumor.type = "ER+") {
    fd <- input
    fd$allelic_status <- ifelse(fd$sample %in% davies$biallelic,
        "biallelic", ifelse(fd$sample %in% davies$control, "control", "unknown")
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
        theme(axis.text.x = element_text(angle = 90)) +
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
        theme(axis.text.x = element_text(angle = 90)) +
        stat_compare_means(aes(group = allelic_status),
            method = "wilcox.test",
            label = "..p.format..",
            # label = "..p.signif..",
            # paired = T,
            size = 3
        ) + facet_wrap(~variable, ncol = 1, scales = "free")

    p <- ggarrange(p1, p2, ncol = 2)
    mkdirp("output/serena_biallelic_tests/")
    ggsave(paste0("output/serena_biallelic_tests/", figfile), p, width = 10, height = 10)
}

save_all_biallelic_comparison_plots <- function() {
    # Save different version of the biallelic comparison plots
    save_plot_compare_biallelic(er.svdata$tss_data_score200, "er_tss_score200.png", "ER+ TSS")
    save_plot_compare_biallelic(er.svdata$tts_data_score200, "er_tts_score200.png", "ER+ TTS")
    save_plot_compare_biallelic(er.svdata$sc200_data, "er_sc200.png", "ER+ sc200")
    save_plot_compare_biallelic(tri.svdata$tss_data_score200, "tri_tss_score200.png", "Trp- TSS")
    save_plot_compare_biallelic(tri.svdata$tts_data_score200, "tri_tts_score200.png", "Trp- TTS")
    save_plot_compare_biallelic(tri.svdata$sc200_data, "tri_sc200.png", "Trp- sc200")
}

get_tumor_allelic_compare_rloop <- function(tumor.data) {
    a <- list()
    a$tts_data_score200 <- test_muttype_allelic_compare_rloop_region(
        tumor.data$tts_data_score200,
        "tts",
        group = "tts_data_score200"
    )
    a$tss_data_score200 <- test_muttype_allelic_compare_rloop_region(
        tumor.data$tss_data_score200,
        "tss",
        group = "tss_data_score200"
    )
    a$sc200 <- test_muttype_allelic_compare_rloop_region(
        tumor.data$sc200_data,
        "sc200",
        group = "sc200_data"
    )
    return(rbindlist(a))
}

get_tumor_rloop_compare_allelic <- function(tumor.data) {
    a <- list()
    a$tts_data_score200 <- test_muttype_rloop_region_compare_allelic(
        tumor.data$tts_data_score200,
        "tts",
        group = "tts_data_score200"
    )
    a$tss_data_score200 <- test_muttype_rloop_region_compare_allelic(
        tumor.data$tss_data_score200,
        "tss",
        group = "tss_data_score200"
    )
    a$sc200 <- test_muttype_rloop_region_compare_allelic(
        tumor.data$sc200_data,
        "sc200",
        group = "sc200_data"
    )
    return(rbindlist(a))
}

save_all_serena_biallelic_test_result <- function() {
    # Save all Serena's biallelic different groups:
    # 1. among rloop/nonrloop regon, compare biallelic vs. control
    # 2. Among biallelic/control samples, compare rloop vs nonrloop regions

    a <- list()
    a$er.1 <- get_tumor_allelic_compare_rloop(er.svdata)
    a$er.1$tumor <- "ER+"
    a$tri.1 <- get_tumor_allelic_compare_rloop(tri.svdata)
    a$tri.1$tumor <- "Trp-"

    a$er.2 <- get_tumor_rloop_compare_allelic(er.svdata)
    a$er.2$tumor <- "ER+"
    a$tri.2 <- get_tumor_rloop_compare_allelic(tri.svdata)
    a$tri.2$tumor <- "Trp-"
    rbindlist(a) %>%
        fwrite("output/serena_biallelic_tests/all_serena_biallelic_test_results.csv")
}

######################################################################
main <- function() {
    # main_boxplots()
    # main_save_test_results()
    save_all_biallelic_comparison_plots()

    # save_all_serena_biallelic_test_result()
}
# main()
# a <- fread("output/serena_biallelic_tests/all_serena_biallelic_test_results.csv")
# table(a$group)
