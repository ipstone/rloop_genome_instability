from SigProfilerExtractor import sigpro as sig

# Using provided data folder to specify SBS96 matrix and output folrer
from sys import argv

input_data_folder = argv[1]

# Specify the rloop_region type input file
rtype = input_data_folder.strip().split("/")[-1]
snv_suffix = ".region"
if rtype == "all":
    snv_suffix = ".all"


def main_function():
    # to get input from table format (mutation catalog matrix)
    path_to_example_table = sig.importdata("matrix")
    data = input_data_folder + "/output/SBS/icgc_ov.SBS96" + snv_suffix
    # data = input_data_folder + "/output/SBS/icgc_ov.SBS96.region"
    output_folder = input_data_folder + "/extracted_signatures"
    # data = "analysis/icgc_ov_analysis/output/SBS/icgc_ov.SBS96.all"  # you can put the path to your tab delimited file containing the mutational catalog matrix/table

    sig.sigProfilerExtractor(
        "matrix",
        output_folder,
        # "analysis/icgc_ov_analysis/extracted_sigantures",
        data,
        # opportunity_genome="GRCh37",
        minimum_signatures=1,
        maximum_signatures=7,
        # max_nmf_iterations=500000,
    )


if __name__ == "__main__":
    main_function()

