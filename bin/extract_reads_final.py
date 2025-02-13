#!/usr/bin/env python3 
# 4-space indented, v0.0.1
# File name: extract_reads_final.py
# Description: Extracts reads harboring potential gencDNA insertions from the original Minimap2 alignment.
# Author: Robert Linder
# Date: 2024-01-09

import pandas as pd
import re
import os
import argparse
from subprocess import Popen, PIPE, check_output
from itertools import combinations
import pysam

def parse_args():
    """this enables users to provide input as command line arguments to minimize the need for directly modifying the script; please enter '-h' to see the full list of options"""
    parser = argparse.ArgumentParser(description="Load files for final filtration step of the gencDNA caller")
    parser.add_argument("gencDNA_file", type=str, help="potential gencDNA bed file to filter")
    parser.add_argument("bam_file", type=str, help="Sorted and indexed bam file generated from the original minimap2 alignment")
    args = parser.parse_args()
    return args

def subset_bam(genc, bam, file_name):
    """This extracts the read names that will get kept in the bam file"""
    # need to check to ensure there are any hits at all
    if os.stat(genc).st_size > 0:
        genc_df = pd.read_csv(genc, sep='\t', header=None)
        if genc_df.iloc[1].str.contains("NativeHit").any():
            read_id = genc_df[3].str.split('_').values
            read_name = ['_'.join(r[0:4]) for r in read_id]
            read_type = [r[4] for r in read_id]
            alignment = genc_df[29].values
            alignment_type = genc_df[30].values
            as_score = genc_df[31].values
            all_id_info = zip(read_name, alignment, alignment_type, as_score, read_type)
            to_id_reads = [f"{r[0]}_{r[1]}_{r[2]}_{r[3]}_{r[4]}" for r in all_id_info]
        elif genc_df.iloc[1].str.contains("Sniffles2").any():
            read_names = genc_df[3].values
            to_id_reads = read_names
        with pysam.AlignmentFile(bam, "rb") as bamfile:
            with pysam.AlignmentFile(file_name + ".filtered.bam", "wb", template=bamfile) as output_file:
                counter = 0
                for read in bamfile.fetch():
                    if "Sniffles2" in read.query_name:
                        read_id = read.query_name.split('-')
                        read_id = read_id[1]
                    else:
                        read_id = read.query_name
                    if read_id in to_id_reads:
                        counter += 1
                        output_file.write(read)
    else:
        # if there are no hits left, create a bam file as output
        header = { 'HD': {'VN': '1.0'},
            'SQ': [{'LN': 1575, 'SN': 'chr1'},
                   {'LN': 1584, 'SN': 'chr2'}] }
        with pysam.AlignmentFile(file_name + ".filtered.bam", "wb", header=header) as outf:
            a = pysam.AlignedSegment()
            a.query_name = "fake_read"
            a.query_sequence="AGCTTAGCTAGCTACCTATATCTTGGTCTTGGCCG"
            a.flag = 99
            a.reference_id = 0
            a.reference_start = 32
            a.mapping_quality = 20
            a.cigar = ((0,10), (2,1), (0,25))
            a.next_reference_id = 0
            a.next_reference_start=199
            a.template_length=167
            a.query_qualities = pysam.qualitystring_to_array("<<<<<<<<<<<<<<<<<<<<<:<9/,&,22;;<<<")
            a.tags = (("NM", 1),
                    ("RG", "L1"))
            outf.write(a)
def main():
    inputs = parse_args()
    genc =  inputs.gencDNA_file
    bam = inputs.bam_file
    pysam.index(bam)
    file_name = bam.split('/')[-1].split('.')[:-1]
    file_name = '.'.join(file_name)
    subset_bam(genc, bam, file_name)

if __name__ =="__main__":
    main()