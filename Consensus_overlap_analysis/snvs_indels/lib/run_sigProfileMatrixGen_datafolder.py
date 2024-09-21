# Originally planned to run on the data folder, together with bed file
# -- but found that using the bed file taking about 30 mins to run
# -- decided to use bedtools for the intersection and then run in each folder without bedfiles

from SigProfilerMatrixGenerator.scripts import SigProfilerMatrixGeneratorFunc as matGen

# There are two input file path: snv file path, bed file path
from sys import argv

input_data_folder = argv[1]
argv_len = len(argv)

# Config whether there's bedfile given as 2nd argument
no_bed_file = True
bed_file_fullpath = None

if argv_len > 2:
    no_bed_file = False
    bed_file_fullpath = argv[2]


def main():
    project = "icgc_ov"
    genome = "GRCh37"
    # vcfFiles = "analysis/icgc_ov_analysis_bed_gene_body"
    vcfFiles = input_data_folder
    # -- bed file need to use the full path
    # bedFile = (
    #     "/data/projects/peix/rloop_project/input/rloop_bed_files/rloop_genebody.bed"
    # )
    # -- bed file need to use the full path

    matrices = matGen.SigProfilerMatrixGeneratorFunc(
        project,
        genome,
        vcfFiles,
        exome=False,
        bed_file=bed_file_fullpath,
        # bed_file=None,
        # bed_file=bedFile,
        chrom_based=False,
        plot=False,  # Selecting no plot to make it faster
        # plot=True,
        tsb_stat=False,
        seqInfo=False,
    )


if __name__ == "__main__":
    main()

