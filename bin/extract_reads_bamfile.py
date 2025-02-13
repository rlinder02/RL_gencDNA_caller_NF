#!/usr/bin/env python3 
# 4-space indented, v0.0.1
# File name: extract_reads_bamfile.py
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
    parser.add_argument("vcf_file", type=str, help="Sniffles2 generated vcf file to pull read names from")
    parser.add_argument("bam_file", type=str, help="Sorted and indexed bam file generated from the original minimap2 alignment")
    args = parser.parse_args()
    return args

def find_read_names(genc, vcf):
    """This extracts the read names from the Sniffles2 vcf as a dictionary with the insertion id as the key and read names supporting it as the values"""
    genc_df = pd.read_csv(genc, sep='\t', header=None)
    vcf_df = pd.read_csv(vcf, sep='\t', header=None)
    ins_names = set(genc_df[3])
    ins_df = vcf_df.loc[vcf_df[2].isin(ins_names), [2,5]]
    finding_reads = re.compile("RNAMES=.*")
    info = ins_df[5].str.split(';').tolist()
    all_reads = []
    for item in info:
        find_reads = [r.split('=')[1].split(',') for r in item if finding_reads.match(r)]
        all_reads.append(find_reads)
    all_reads2 = [x for xs in all_reads for x in xs]
    reads = ins_df[5].str.split(';').str[5].str.split('=').str[1].str.split(',')
    reads.index = reads.index.astype(int)
    insertions = vcf_df.iloc[reads.index,2]
    ins_chroms= vcf_df.iloc[reads.index,0]
    # Need some way of differentiating when two reads map to different insertions- need the cigar string- or one other parameter to use, as just the insertion ID and read name aren't by themselves enough to differentiate
    ins_reads = list(zip(insertions, all_reads2, ins_chroms))
    return ins_reads 

def subset_bam(bam, read_names, file_name):
    """Extract reads matching read names in the dictionary above from the bam file, renaming reads to include the insertion id detected within them (eg. read_id-insertion-id)"""
    ins_identifiers = [x[0] for x in read_names]
    read_identifiers = [x[1] for x in read_names]
    all_ins_reads = [x for xs in read_identifiers for x in xs]
    with pysam.AlignmentFile(bam, "rb") as bamfile:
        with pysam.AlignmentFile(file_name + ".filtered.bam", "wb", template=bamfile) as output_file:
            for read in bamfile.fetch():
                read_id = read.query_name
                if read_id in all_ins_reads:
                    counter = 0
                    for record in read_names:
                        counter += 1
                        if read_id in record[1]:
                            #print("It is here")
                            #print(record[0])
                            read.query_name = read_id + "-" + record[0] + '-' + record[2]
                            #print(read.query_name)
                            output_file.write(read)
                    #print(counter)

def main():
    inputs = parse_args()
    genc =  inputs.gencDNA_file
    vcf = inputs.vcf_file
    bam = inputs.bam_file
    pysam.index(bam)
    read_names = find_read_names(genc, vcf)
    file_name = bam.split('/')[-1].split('.')[:-1]
    file_name = '.'.join(file_name)
    subset_bam(bam, read_names, file_name)

if __name__ =="__main__":
    main()