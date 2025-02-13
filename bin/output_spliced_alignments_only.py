#!/usr/bin/env python3 
# 4-space indented, v0.0.1
# File name: output_spliced_alignments_only.py
# Description: Filters out reads that do not have a 30+ bp deletion or any length of a spliced out intron (N).
# Author: Robert Linder
# Date: 2024-03-20

import argparse
import pysam 
import re

def parse_args():
	"""this enables users to provide input as command line arguments to minimize the need for directly modifying the script; please enter '-h' to see the full list of options"""
	parser = argparse.ArgumentParser(description="Load file for filtring duplicate entries")
	parser.add_argument("gencDNA_bam", type=str, help="gencDNA bam file to filter")
	args = parser.parse_args()
	return args

def filter_unspliced_entries(bam, output_name):
    bamfile = pysam.AlignmentFile(bam, "rb")
    output_file = pysam.AlignmentFile(output_name + ".deletions.bam", "wb", template=bamfile)
    for read in bamfile.fetch():
        find_dels = re.compile("[3][0-9]D")
        if 'N' in read.cigarstring or find_dels.match(read.cigarstring):
            output_file.write(read)
    output_file.close()
    bamfile.close()

def main():
    inputs = parse_args()
    bam =  inputs.gencDNA_bam
    output = bam.split('/')[-1].split('.')[:-1]
    output_name = '.'.join(output)
    filter_unspliced_entries(bam, output_name)
    pysam.index(output_name + ".deletions.bam")

if __name__ =="__main__":
    main()