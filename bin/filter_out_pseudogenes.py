#!/usr/bin/env python3 
# 4-space indented, v0.0.1
# File name: filter_out_pseudogenes.py
# Description: Filters out reads for which any alignment maps to a pseudogene.
# Author: Robert Linder
# Date: 2024-03-20

import argparse
import pysam 
import re
import pandas as pd
import os

def parse_args():
    """this enables users to provide input as command line arguments to minimize the need for directly modifying the script; please enter '-h' to see the full list of options"""
    parser = argparse.ArgumentParser(description="Load file for filtring duplicate entries")
    parser.add_argument("gencDNA_bam", type=str, help="gencDNA bam file to filter")
    parser.add_argument("ps_file", type=str, help="text file of reads that overlap known pseudogenes")
    args = parser.parse_args()
    return args

def filter_ps_entries(bam, ps, output_name):
    bamfile = pysam.AlignmentFile(bam, "rb")
    output_file = pysam.AlignmentFile(output_name + ".pseudogeneless.bam", "wb", template=bamfile)
    # Need to check to ensure the file is not empty before attempting to read it in 
    if os.stat(ps).st_size > 0:
        ps_reads = pd.read_csv(ps, header=None)
        read_names = ps_reads[0].values
        just_read_names = [r.split('-')[0] for r in read_names]
        for read in bamfile.fetch():
            just_read = read.query_name.split('-')[0]
            if not just_read in just_read_names:
                output_file.write(read)
    else:
        for read in bamfile.fetch():
            output_file.write(read)
    output_file.close()
    bamfile.close()

def main():
    inputs = parse_args()
    bam =  inputs.gencDNA_bam
    pysam.index(bam)
    ps = inputs.ps_file
    output = bam.split('/')[-1].split('.')[:-1]
    output_name = '.'.join(output)
    filter_ps_entries(bam, ps, output_name)

if __name__ =="__main__":
    main()