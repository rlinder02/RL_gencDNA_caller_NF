#!/usr/bin/env python3 
# 4-space indented, v0.0.1
# File name: filter_duplicate_entries.py
# Description: Filters out duplicate entries from the list of filtered, potential gencDNAs.
# Author: Robert Linder
# Date: 2023-11-29

import argparse
import pysam 

def parse_args():
	"""this enables users to provide input as command line arguments to minimize the need for directly modifying the script; please enter '-h' to see the full list of options"""
	parser = argparse.ArgumentParser(description="Load file for filtring duplicate entries")
	parser.add_argument("gencDNA_bam", type=str, help="gencDNA bam file to filter")
	args = parser.parse_args()
	return args

def filter_duplicate_entries(bam, output_name):
    bamfile = pysam.AlignmentFile(bam, "rb")
    output_file = pysam.AlignmentFile(output_name + ".nodups.bam", "wb", template=bamfile)
    reads_processed = []
    for read in bamfile.fetch():
        gene = bamfile.get_reference_name(read.reference_id)
        position = read.reference_end
        cigar = read.cigarstring
        identifier = f"{gene}_{position}_{cigar}"
        if identifier not in reads_processed:
            output_file.write(read)
            reads_processed.append(identifier)
    output_file.close()
    bamfile.close()

def main():
    inputs = parse_args()
    bam =  inputs.gencDNA_bam
    output = bam.split('/')[-1].split('.')[:-1]
    output_name = '.'.join(output)
    filter_duplicate_entries(bam, output_name)
    pysam.index(output_name + ".nodups.bam")

if __name__ =="__main__":
    main()