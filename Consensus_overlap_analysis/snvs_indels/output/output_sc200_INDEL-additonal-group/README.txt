2021-10-01 Generating additonal INDEL category data for paper

* Generate the long_indel (long ins, del, and MH), and long_simple_indel
  category for testing/results for paper
  - This calculation is through the ./plot_figure_heatmap_INDEL_paper_data.R script
    in the folder 
      This script also generate the heatmaps for all indel category and tumor
      types in the folder.

* The indel test results folder contained here are generated with the
  analyze-indel scriopt in the lib folder:
      The code is updated to include additonal categories:
            adj_long_indel (long ins, del, + MH)
            adj_long_simple_indel  (long ins, del)
