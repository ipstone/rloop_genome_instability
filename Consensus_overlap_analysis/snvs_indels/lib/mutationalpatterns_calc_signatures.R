# To use mutationalPattern to fit each sample's mutation to cosmic signatures
# - use the fit_to_sigantures function in mutationalPatterns
rm(list = ls())
library(data.table)
library(magrittr)
library(ipfun)
library(MutationalPatterns) # conda activate MutationalPatterns
library(ggpubr)
source("lib/lib_load_dataprep.R")

# load cosmic signature
cosmic_sigs <- get_known_signatures()
rloop_types <- get_rloop_type()
rloop_types$path <- paste0("ICGC_OV/", rloop_types$all)

# ----------------------------------------------
# Below is the reference code previously using mutationalPatterns:
# -- the input mut_mat.tsv is the format needed for mutationalPatterns:
# ----------------------------------------------
# mutatonalPatterns signatures sections
arrange_mutMat_to_mutationalPatternFormat <- function(x) {
    # Order the input x with rownames (96 pattern) to mutationalPatterns convensions
    mut_mat <- as.matrix(read.table("input/mutationalPatterns_mut_mat.tsv", sep = "\t", header = T))
    new_order <- match(rownames(mut_mat), rownames(x))
    y <- x[new_order, ]
    return(y)
}

## Coding plans
# 1. generate fit, fitstrict, heatmap plot for all 4x2+1 conditions for the different rloops
# 2. Select the sigantures to compare (such as siganture 6, 40) within each  rloop genome type
#       - using best math fit
#       - using strict fit

calc_rloop_region_sigfitres <- function(rloop_path) {
    # Loading the SBS96 profile from the rloop_path
    # Return an object with fr$fit.best = all sig fit; fr$fit.strict = strict fit to cosmic signatures
    dfile <- paste0("analysis/", rloop_path, "/output/SBS/icgc_ov.SBS96.all")
    d <- fread(dfile)
    d[["="]] <- NULL
    dm <- as.matrix(d, rownames = "MutationType")
    mutmat <- arrange_mutMat_to_mutationalPatternFormat(dm)
    # Calculate math best fit
    fit.res <- fit_to_signatures(mutmat, cosmic_sigs)
    # Calculate strict fitting critera  -- this takes sometime to calculate
    fit.strict <- fit_to_signatures_strict(mutmat, cosmic_sigs,
        max_delta = 0.02
    )
    fit.strict.res <- fit.strict$fit_res # Fit result with stricter removing step

    # function return / result
    fr <- list()
    fr$fit.best <- fit.res
    fr$fit.strict <- fit.strict.res
    return(fr)
}

calc_all_rloop_regions_fitsigres <- function() {
    # Calculate all the signature fit for different regions/types, save the combined object in a rds file
    fr <- list()
    fr$all <- calc_rloop_region_sigfitres("icgc_ov_analysis")
    for (i in rloop_types$all) {
        writeLines(paste("Calculate", i))
        path <- paste0("ICGC_OV/", i)
        fr[[i]] <- calc_rloop_region_sigfitres(path)
    }
    return(fr)
}

# save all regons signature fits
save_all_sigfitres <- function() {
    fr <- calc_all_rloop_regions_fitsigres()
    saveRDS(
        fr,
        "output/mutationalpatterns_signatures/all_fit_signatures.rds"
    )
    return(fr)
}
# r <- save_all_sigfitres()
r <- readRDS("output/mutationalpatterns_signatures/all_fit_signatures.rds")

# Save signature contribution heatmaps for all cases
saveplot_heatmap_sigcontrib <- function(allsigfit) {
    # Save all the heatmap for signatures decomposed
    prf <- "output/mutationalpatterns_signatures/siganture_contribution_heatmaps/"

    for (i in names(allsigfit)) {
        # obtain fit object
        fitobj <- allsigfit[[i]]

        # setup path
        allf <- paste0(prf, "all/")
        mkdirp(allf)
        allf.png <- paste0(allf, i, ".png")
        # setup strict sig heatmap path
        strictf <- paste0(prf, "strict/")
        mkdirp(strictf)
        strictf.png <- paste0(strictf, i, ".png")

        # setup strict sig heatmap path
        combined <- paste0(prf, "combined/")
        mkdirp(combined)
        combined.png <- paste0(combined, i, ".png")

        # Save actual all plot
        writeLines(paste0("save png plot for ", allf.png))
        # print(head(fitobj$fit.best$contribution)) # Debug
        # -- for mutationalpatterns: the plots are ggplot object, so need to print to evaluate or use ggsave etc.
        p1 <- plot_contribution_heatmap(fitobj$fit.best$contribution,
            cluster_samples = T,
            cluster_sigs = T
        )
        ggsave(allf.png, p1, width = 12, height = 12)

        writeLines(paste0("save png plot for ", strictf.png))
        p2 <- plot_contribution_heatmap(fitobj$fit.strict$contribution,
            cluster_samples = T,
            cluster_sigs = T
        )
        ggsave(strictf.png, p2, width = 12, height = 12)

        writeLines(paste0("save png plot for ", combined.png))
        p <- ggarrange(p1, p2, nrow = 1)
        ggsave(combined.png, p, width = 24, height = 12)
    }
}
# saveplot_heatmap_sigcontrib(r)

saveplot_heatmap_rloop_side_by_side <- function(r) {
    # r is the signature fit object for all cases
    # Save all the combined plots for different rloop vs norloop for 4 type of genome regions

    prf <- "output/mutationalpatterns_signatures/siganture_contribution_heatmaps/rloop_types_combined/"
    mkdirp(prf)

    # plot the 4 regions types
    region_types <- c("exon", "genebody", "tss", "tts")

    for (i in region_types) {
        rloop_name <- paste0("rloop_", i)
        norloop_name <- paste0("norloop_", i)
        p.best.png <- paste0(prf, "best_fit_", i, ".png")
        p.strict.png <- paste0(prf, "strict_fit_", i, ".png")
        print(p.best.png)
        print(p.strict.png)

        # rloop fit contribution
        rloop.fit.best <- r[[rloop_name]]$fit.best$contribution
        rloop.fit.strict <- r[[rloop_name]]$fit.strict$contribution
        # norloop fit contribution
        norloop.fit.best <- r[[norloop_name]]$fit.best$contribution
        norloop.fit.strict <- r[[norloop_name]]$fit.strict$contribution

        # mathematically best fit
        p.best.rloop <- plot_contribution_heatmap(rloop.fit.best,
            cluster_samples = T,
            # cluster_sigs = T
        ) + ggtitle(rloop_name)
        p.best.norloop <- plot_contribution_heatmap(norloop.fit.best,
            cluster_samples = T,
            # cluster_sigs = T
        ) + ggtitle(norloop_name)
        p.best <- ggarrange(p.best.rloop, p.best.norloop, nrow = 1)
        ggsave(p.best.png, p.best, width = 24, height = 12)

        # stricter fit
        p.strict.rloop <- plot_contribution_heatmap(rloop.fit.strict,
            cluster_samples = T,
            # cluster_sigs = T
        ) + ggtitle(rloop_name)
        p.strict.norloop <- plot_contribution_heatmap(norloop.fit.strict,
            cluster_samples = T,
            # cluster_sigs = T
        ) + ggtitle(norloop_name)
        p.strict <- ggarrange(p.strict.rloop, p.strict.norloop, nrow = 1)
        ggsave(p.strict.png, p.strict, width = 24, height = 12)
    }
}
# saveplot_heatmap_rloop_side_by_side(r)


## Work on plot tests for each siganture category for comparion between rloop vs non-rloop for contributions of each individual samples decomposed

get_loop_region_fit_data <- function(data = r, loop_type = "rloop", region = "exon", fit = "fit.best") {
    # With the given input, return a data.table with the relevant contributions
    rloop_region <- paste0(loop_type, "_", region)
    # Select from the result fit contributions
    a <- data[[rloop_region]][[fit]]$contribution
    # Transform into right formats
    a <- t(a) # transpose
    b <- a / rowSums(a)
    b <- data.table(b, keep.rownames = T)
    bm <- melt(b, id.vars = "rn")
    bm$region <- region
    bm$type <- loop_type
    return(bm)
}

# Take the fit results calculated, and return  combined fit percentages for testing
get_all_fitdata_for_test <- function(data = r, fit.type = "fit.best") {
    # Combine all the different rloop types and regions, of the fit kind
    fr <- list()
    for (lloop in c("rloop", "norloop")) {
        for (lregion in c("exon", "genebody", "tss", "tts")) {
            loop_region <- paste0(lloop, "_", lregion)
            fr[[loop_region]] <- get_loop_region_fit_data(data, lloop, lregion, fit.type)
        }
    }
    fr2 <- rbindlist(fr)
    return(fr2)
}

saveplot_all_sigfit_compare_rloop_region <- function() {
    # Save the plots for comparing every signatures from the decomposed

    # Start with best fit
    pdf("output/mutationalpatterns_signatures/rloop_signature_comparisons/bestfit_rloop_sig_comparison.pdf")
    b <- get_all_fitdata_for_test(data = r, fit.type = "fit.best")
    for (i in unique(b$variable)) {
        print(paste0("Working on siganture ", i))
        # Selecting for the particular signature
        d <- b[variable == i]
        p <- ggboxplot(d,
            x = "type", y = "value",
            color = "type", palette = "npg",
            # size = 0.1
            alpha = 0.7,
            ylab = "Signature contribution",
            title = i,
            add = c("jitter", "median_iqr"),
        ) + stat_compare_means(aes(group = type),
            method = "wilcox.test",
            label = "..p.format.."
        ) + facet_wrap(~region, nrow = 1)
        print(p)
    }
    dev.off()

    # Start with best fit
    pdf("output/mutationalpatterns_signatures/rloop_signature_comparisons/strictfit_rloop_sig_comparison.pdf")
    b <- get_all_fitdata_for_test(data = r, fit.type = "fit.strict")
    for (i in unique(b$variable)) {
        print(paste0("Working on siganture ", i))
        # Selecting for the particular signature
        d <- b[variable == i]
        p <- ggboxplot(d,
            x = "type", y = "value",
            color = "type", palette = "npg",
            # size = 0.1
            alpha = 0.7,
            ylab = "Signature contribution",
            title = i,
            add = c("jitter", "median_iqr"),
        ) + stat_compare_means(aes(group = type),
            method = "wilcox.test",
            label = "..p.format.."
        ) + facet_wrap(~region, nrow = 1)
        print(p)
    }
    dev.off()
}
# saveplot_all_sigfit_compare_rloop_region()

## from rloop_path, load, return mut matrix 96 profiles
load_rloop_region_mutmat <- function(rloop_path) {
    # Loading the SBS96 profile from the rloop_path
    # Return mut profile matrix (96 profile for mutationalpatterns format)
    dfile <- paste0("analysis/ICGC_OV/", rloop_path, "/output/SBS/icgc_ov.SBS96.all")
    d <- fread(dfile)
    d[["="]] <- NULL
    dm <- as.matrix(d, rownames = "MutationType")
    mutmat <- arrange_mutMat_to_mutationalPatternFormat(dm)
    return(mutmat)
}

plot_region_mutmat <- function(region) {
    # Return a ggplot on mutprofile of mut matrix
    rloop <- load_rloop_region_mutmat(paste0("rloop_", region))
    norloop <- load_rloop_region_mutmat(paste0("norloop_", region))
    both <- cbind(rloop, norloop)
    types <- c(rep("rloop", ncol(rloop)), rep("norloop", ncol(norloop)))
    p <- plot_profile_heatmap(both, by = types)
    return(p)
}
## Save all rloop, region types mut matrix heatmap plot
saveplot_all_mutprofile_heatmap <- function() {
    plot_region <- function(x) {
        pfile <- paste("output/mutationalpatterns_signatures/mutprofile_heatmap_by_rloop_", x, "_regions.png")

        png(pfile, width = 960, height = 960)
        print(plot_region_mutmat(x))
        dev.off()
    }
    # previously there is some problem when try to save the plots into 1 pdf file
    lapply(c("exon", "genebody", "tss", "tts"), plot_region)
}
# saveplot_all_mutprofile_heatmap()


## Plot 6 mutation profiles - trying to use plot_spectrum function
get_spectrum_data <- function(rpath) {
    d.fp <- paste0("analysis/", rpath, "/output/SBS/icgc_ov.SBS6.all")
    # d <- fread("analysis/icgc_ov_analysis/output/SBS/icgc_ov.SBS6.all")
    d <- fread(d.fp)
    d[["="]] <- NULL
    dm <- as.matrix(d, rownames = "MutationType")
    dmt <- t(dm)
    dmt.df <- as.data.frame(dmt)
    return(dmt.df)
    # p <- plot_spectrum(dmt.df)
    # ggsave("output/mutationalpatterns_signatures/SBS6_spectrum/all_regions.png", p)
}

saveplot_spectrum_all_regions <- function() {
    # Save all different regons rloop vs norloop 6 profile spectrum
    plot_region <- function(x) {
        pfile <- paste("output/mutationalpatterns_signatures/SBS6_spectrum/spectrum_", x, "_regions.png")
        png(pfile, width = 960, height = 960)
        r.fp <- paste0("ICGC_OV/rloop_", x)
        nor.fp <- paste0("ICGC_OV/norloop_", x)
        p.r <- plot_spectrum(get_spectrum_data(r.fp)) + ylim(0, 0.45) + ggtitle("rloop")
        p.nor <- plot_spectrum(get_spectrum_data(nor.fp)) + ylim(0, 0.45) + ggtitle("norloop")
        p <- ggarrange(p.r, p.nor, nrow = 1)
        print(p)
        dev.off()
    }
    # previously there is some problem when try to save the plots into 1 pdf file
    lapply(c("exon", "genebody", "tss", "tts"), plot_region)
}
saveplot_spectrum_all_regions()








#------------------------------------------------------------
experiment_mutationalpatterns <- function() {
    # Loading matrix to see if fit with the format
    d <- fread("analysis/icgc_ov_analysis/output/SBS/icgc_ov.SBS96.all")
    d[["="]] <- NULL
    dm <- as.matrix(d, rownames = "MutationType")
    mutmat <- arrange_mutMat_to_mutationalPatternFormat(dm)
    fit.res <- fit_to_signatures(mutmat, cosmic_sigs)

    # Calculate strict fitting critera  -- this takes sometime to calculate
    fit.strict <- fit_to_signatures_strict(mutmat, cosmic_sigs,
        max_delta = 0.02
    )
    fit.strict.res <- fit.strict$fit_res # Fit result with stricter removing step

    # Save the signature plots etc.
    png("output/mutationalpatterns_signatures/icgc_ov_sigfitres_contribution_absolute.png")
    plot_contribution(fit.res$contribution, mode = "absolute")
    dev.off()
    png("output/mutationalpatterns_signatures/icgc_ov_sigfitres_contribution_relative.png")
    plot_contribution(fit.res$contribution, mode = "relative")
    dev.off()

    # Saving the heatmap plots
    png("output/mutationalpatterns_signatures/icgc_ov_sigfitres_heatmap.png",
        width = 960, height = 960
    )
    plot_contribution_heatmap(fit.res$contribution,
        cluster_samples = TRUE,
        cluster_sigs = TRUE
    )
    dev.off()

    # Plot profile heatmap
    png("output/mutationalpatterns_signatures/icgc_ov_sigfitres_profile-heatmap.png",
        width = 960, height = 960
    )
    plot_profile_heatmap(mutmat)
    dev.off()

    tissues <- c(rep("tissue_1", 40), rep("tissue_2", 31))
    png("output/mutationalpatterns_signatures/icgc_ov_sigfitres_profile-heatmap_testing_groups.png",
        width = 960, height = 960
    )
    plot_profile_heatmap(mutmat, by = tissues)
    dev.off()
    ## ------------------------------------------------------
    ## plot for stricter fit
    # Save the signature plots etc.
    png("output/mutationalpatterns_signatures/icgc_ov_sigfitres-strict_absolute.png")
    plot_contribution(fit.strict.res$contribution, mode = "absolute")
    dev.off()
    png("output/mutationalpatterns_signatures/icgc_ov_sigfitres-strict_relative.png")
    plot_contribution(fit.strict.res$contribution, mode = "relative")
    dev.off()

    png("output/mutationalpatterns_signatures/icgc_ov_sigfitres-strict_heatmap.png",
        width = 960, height = 960
    )
    plot_contribution_heatmap(fit.strict.res$contribution,
        cluster_samples = TRUE,
        cluster_sigs = TRUE
    )
    dev.off()
}

## Test plot_profile_heatmap on SBS6 etc.
not_working <- function() {
    # This is not working: the plot_profile_heatmap seems only working with 96 profiles
    d <- fread("analysis/icgc_ov_analysis/output/SBS/icgc_ov.SBS6.all")
    d[["="]] <- NULL
    dm <- as.matrix(d, rownames = "MutationType")
    png("output/mutationalpatterns_signatures/_explore_mutpattern_plots/SBS6_profile_all.png")
    print(plot_profile_heatmap(dm))
    dev.off()
}

## Try plot_spectrum to see how it's done
experiment_plot_spectrum <- function() {
    library(BSgenome)
    ref_genome <- "BSgenome.Hsapiens.UCSC.hg19"
    library(ref_genome, character.only = TRUE)

    vcf_files <- list.files(system.file("extdata", package = "MutationalPatterns"),
        pattern = "sample.vcf", full.names = TRUE
    )
    sample_names <- c(
        "colon1", "colon2", "colon3",
        "intestine1", "intestine2", "intestine3",
        "liver1", "liver2", "liver3"
    )
    grl <- read_vcfs_as_granges(vcf_files, sample_names, ref_genome)
    tissue <- c(rep("colon", 3), rep("intestine", 3), rep("liver", 3))

    muts <- mutations_from_vcf(grl[[1]])
    types <- mut_type(grl[[1]])
    # type_context <- type_context(grl[[1]], ref_genome)
    type_occurrences <- mut_type_occurrences(grl, ref_genome)

    ##########################################
    # > type_occurrences
    #            C>A C>G C>T T>A T>C T>G C>T at CpG C>T other
    # colon1      28   5 109  12  30  12         59        50
    # colon2      77  29 345  36  90  21        209       136
    # colon3      79  19 243  25  61  23        165        78
    # intestine1  19   8  74  19  26   4         33        41
    # intestine2 118  49 423  57 126  27        258       165
    # intestine3  54  27 298  32  67  22        192       106
    # liver1      43  22  94  30  77  34         18        76
    # liver2     144  93 274 103 209  73         48       226
    # liver3      39  28  61  15  32  23          7        54
}