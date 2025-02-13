#!/usr/bin/env python3 
# 4-space indented, v0.0.1
# File name: delineate_aligned_insertions.py
# Description: Compile list of where the insertions align.
# Author: Robert Linder
# Date: 2024-03-12

import pandas as pd
import re
import os
import argparse
from subprocess import Popen, PIPE, check_output
from itertools import combinations
import pysam
from Bio import SeqIO

def parse_args():
    """this enables users to provide input as command line arguments to minimize the need for directly modifying the script; please enter '-h' to see the full list of options"""
    parser = argparse.ArgumentParser(description="Load files for final filtration step of the gencDNA caller")
    parser.add_argument("bam_file", type=str, help="Sorted and indexed bam file generated from the original minimap2 alignment that has been filtered")
    parser.add_argument("bed_file", type=str, help="Final pseudogeneless bed file")
    args = parser.parse_args()
    return args

def extract_mapped_read(bam, bed, file_name):
    all_insertions = []
    insertions_ann = []
    with open(bed, 'r') as bedfile:
        for line in bedfile:
            fields = line.split('\t')
            align_type = fields[30]
            as_score = int(float(fields[31]))
            if not "ccs-" in line:
                info_field = fields[20].split(';')
                genes = fields[25]
                read_id = fields[3]
                chrom = fields[0]
                starting = int(fields[1])
            exon_finder = re.compile("exon_id=.*|ID=exon-")
            exons = [e.split('=')[1] for e in info_field if exon_finder.match(e)]
            exons = ';'.join(exons)
            all_insertions.append([read_id, chrom, int(starting), align_type, as_score])
            insertions_ann.append([read_id, chrom, starting, align_type, as_score, exons, genes])
    ins_info = []
    bamfile = pysam.AlignmentFile(bam, "rb")
    for read in bamfile.fetch():
        if "NativeHit" in read.query_name:
            name_split = read.query_name.split('_')
            name1 = name_split[0:4] + [name_split[7]]
            name = '_'.join(name1)
            alignment_type = name_split[5]
            as_score = name_split[6]
            c = bamfile.get_reference_name(read.reference_id)
            start_pos = read.reference_start
        else:
            name = read.query_name.split('_')[0]
            alignment_type = read.query_name.split('_')[2]
            as_score = int(float(read.query_name.split('_')[3]))
            c = bamfile.get_reference_name(read.reference_id)
            start_pos = read.reference_start
        read_lst = [name, c, int(start_pos), alignment_type, as_score]
        if read_lst in all_insertions:
            all_exons = [ex[5] for ex in insertions_ann if name in ex and c in ex and start_pos in ex and alignment_type in ex and as_score in ex]
            all_exons = ';'.join(all_exons)
            all_genes = [g[6] for g in insertions_ann if name in g and c in g and start_pos in g and alignment_type in g and as_score in g and len(g[6]) > 0]
            all_genes = ';'.join(all_genes)
            c = f"{c};"
            covered_pos = read.get_blocks()
            start_positions = [str(pos[0]) for pos in covered_pos]
            end_positions = [str(pos[1]) for pos in covered_pos]
            c_length = len(start_positions)
            start_positions = ';'.join(start_positions)
            end_positions = ';'.join(end_positions)
            chroms = c*c_length
            chroms = chroms[:-1]
            read_info = [name, chroms, start_positions, end_positions, all_exons, all_genes]
            ins_info.append(read_info)
    ins_df = pd.DataFrame(ins_info, columns=['InsertionID', 'Chromosome', 'Start', 'End', 'Exons', 'Genes'])
    # check if anything is in this dataframe 
    if len(ins_df) > 0:
        if "NativeHit" in ins_df['InsertionID'].values[0]:
            ins_df.insert(loc=1, column ='Read', value = ins_df['InsertionID'])
            ins_df['InsertionID'] = ['' for i in range(ins_df.shape[0])]
            grouped_df_lst = [df for _, df in ins_df.groupby('Exons')]
            counter = 0
            for read_df in grouped_df_lst:
                counter += 1
                find_rows = read_df.index
                ins_df.loc[ins_df.index[find_rows], 'InsertionID'] = counter
            read_df_lst = [df for _, df in ins_df.groupby('Read')]
            for read_df in read_df_lst:
                insertion_id = read_df.iloc[0, 0:1]
                insertion_id = insertion_id.values[0]
                find_rows = read_df.index
                ins_df.loc[ins_df.index[find_rows], 'InsertionID'] = insertion_id
            ins_df = ins_df.sort_values('InsertionID')
            ins_df.to_csv(f"{file_name}.native_locus.report.txt", index=False, sep = '\t', header=True)
        else:
            ins_df.to_csv(f"{file_name}.ins.report.txt", index=False, sep = '\t', header=True)
    else:
        ins_df.to_csv(f"{file_name}.ins.report.txt", index=False, sep = '\t', header=True)

def main():
    inputs = parse_args()
    bam = inputs.bam_file
    pysam.index(bam)
    bed = inputs.bed_file
    file_name = bam.split('/')[-1].split('.')[0]
    extract_mapped_read(bam, bed, file_name)

if __name__ =="__main__":
    main()