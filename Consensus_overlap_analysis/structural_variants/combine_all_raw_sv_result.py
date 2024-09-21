#!/usr/bin/env python3
## use this simple script to combine all the intersected sv file together

from sys import argv
from glob import glob
import os

input_folder = argv[1]
all_bed_files = glob(input_folder + "/*.bedpe")
first_ever_line = True

for input_file in all_bed_files:
    input_file_name = os.path.basename(input_file)
    sampleid = input_file_name.split(".")[0]

    with open(input_file) as f:
        header = f.readline()
        if first_ever_line:
            print("sampleid" + "\t" + header.strip().replace("#", "") + "\tinputfile")
            # print("this is working")
            # print(header)
            first_ever_line = False
        for l in f:
            print(sampleid + "\t" + l.strip() + "\t" + input_file)

