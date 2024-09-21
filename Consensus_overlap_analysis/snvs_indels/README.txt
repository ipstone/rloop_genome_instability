These are the codes to analyze SNVs and Indels mutations in the context of R-loop regions in multiple tumor types.

At the top level, the mutation data (SNVs or indels) from TCGA or ICGC were
intersected with either R-loop regions or the control non-R-loop regions (also
in the context of genes, such TSS/TES regions) as the comparison datasets. Then
these mutation counts/frequencies were further compared and visualized in R.

Below is the note for each sub-folder for these analysis:

- The convention to run these codes are at this folder level, such that for example:
    snakemake -s smk/run_gen_mutmatrix.smk
    Rscript lib/j



* smk: Snakemake and Make files:
    These are the snakemake or make pipelines to generate intersected mutation
    data with the R-loop related bed files and making relevant plots.


    run_gen_mutmatrix.smk
    : Generate all mutation matrix through sigProfiler Matrix Gen for all type
        of tumors, for both SNVs and Indels

    run_plot_indel_all-tumors.smk
    : Run tests and make plots for indels.

    run_plot_SBS6_SBS18_all_tumors.smk
    : Analysze the SBS6 and SBS18 signatures (just an example, can be used to analyze other sigantures).

    intersect_snv_beds.mk
    : This is an earlier version Make setup to intersect mutations by R-loop relatd bed files.

    generate_all_sigProfExtract_sigs.mk
    : Extract all signatures through SigProfilerExtractor

    run_clean_intersected_bed_snv.smk
    : Cleaning up script - to clean up previous calculation and archive them.

* lib: R and python code for analysis and visualization

  Additionally, two packages' source code were used for our analysis, which are cloned into this folder:
        SigProfilerMatrixGenerator
        MutationalPatterns

* output: caluclation output and numbers for the paper.
