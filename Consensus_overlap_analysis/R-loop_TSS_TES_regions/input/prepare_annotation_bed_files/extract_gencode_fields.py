#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# @Author: peix
# @Date:  2021-06-08  
# Extract the annotation fields from gencodeV19 for gene/transcript
# information


# Extract the useful info fields for SV features for analysis
from sys import argv

def split_for_dict(x):
    "Split x by =, if not splittable, the use empty string"
    xs = x.strip().split()
    # print(xs)
    if len(xs)==2:
        return(xs)
    else:
        xs[1] ="no-split-2-item:" + x
        return(xs[0:1])


with open(argv[1]) as del_file:
    header = del_file.readline()
    output_header = "chrom\tstart\tend\tstrand\tsource\tfield_type\tgene_type\tgene_status\tgene_name\ttranscript_type\ttranscript_status\ttranscript_name\tgene_id\ttranscript_id"
    print(output_header)

    for l in del_file:
        # Go through each line to generate dictionary
        fields = l.strip().replace("\"", "").split("\t")
        infos = fields[8].strip().split(";") # The information field
        # print(infos)
        fd = dict( split_for_dict(item) for item in infos if item )

        out_fields = ["gene_type", "gene_status", "gene_name",
                     "transcript_type", "transcript_status",
                     "transcript_name","gene_id", "transcript_id"]

        for okey in out_fields:
            if okey not in fd:
                fd[okey] = ""

        out_values = [fd[i] for i in out_fields]

        chrom = fields[0]
        source = fields[1]
        field_type = fields[2]
        start = fields[3]
        end = fields[4]
        strand = fields[6]
        # -- there are " quotes around text which needs to be removed

        output = '\t'.join([chrom, start, end, strand, source, field_type] + out_values).replace("\"", "")
        print(output)




