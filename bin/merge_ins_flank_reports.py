#!/usr/bin/env python3 
# 4-space indented, v0.0.1
# File name: merge_ins_flank_reports.py
# Description: Merge the insertion and flanking reports 
# Author: Robert Linder
# Date: 2024-03-29

import pandas as pd
import re
import os
import argparse
from subprocess import Popen, PIPE, check_output
from itertools import combinations
import pysam
from Bio import SeqIO
import itertools
import functools
from operator import itemgetter

#pd.set_option("display.width",2000)
#pd.set_option("display.max_columns",None)

def parse_args():
    """this enables users to provide input as command line arguments to minimize the need for directly modifying the script; please enter '-h' to see the full list of options"""
    parser = argparse.ArgumentParser(description="Load files for final filtration step of the gencDNA caller")
    parser.add_argument("ins_file", type=str, help="Text file of insertions called")
    parser.add_argument("flank_file", type=str, help="Text file of flanks supporting insertion sites")
    parser.add_argument('-o', '--output', type=str, help="The file-name prefix")
    args = parser.parse_args()
    return args

def merge_files(ins_file, flank_file):
    """This merges the insertion and flank files into a final report"""
    ins_df = pd.read_csv(ins_file, header = 0, sep = "\t")
    flank_df = pd.read_csv(flank_file, header = 0, sep = "\t")
    if len(ins_df) > 0 and len(flank_df) > 0:
        flank_df = flank_df.astype(str)
        # Only insertions with complete flanking information are kept
        merged_df = ins_df.merge(flank_df, how='inner', on='InsertionID')
        read_id_count = merged_df['ReadID'].values
        read_counts = []
        for r in read_id_count:
            split_r = r.split(';')
            read_counts.append(len(split_r))
        merged_df['Read_count'] = read_counts
    else:
        data = []
        merged_df = pd.DataFrame(data)
    return(merged_df)
    
def filter_merged_file(merge_df, file_name):
    """This removes rows for duplicate insertions that overlap the same set of genes and at least one flanking chromosome dependent on which row has the most supporting reads"""
    if len(merge_df) > 0:
        genes = merge_df['Genes'].str.split(';')
        # for gene in genes:
        unique_genes = list(map(set, genes.values))
        unique_genes = [';'.join(list(g)) for g in unique_genes]
        merge_df['Unique_genes'] = unique_genes
        grouped_df_lst = [df for _, df in merge_df.groupby('Unique_genes')]
        for read_df in grouped_df_lst:
            read_df['Read_type'] = read_df['ReadID'].str.split(';').str[0].str.split('_').str[1]
            unique_lf_chroms = read_df['LF_chr'].str.split(';').values
            unique_rf_chroms = read_df['RF_chr'].str.split(';').values
            lf_chroms = set([item for row in unique_lf_chroms for item in row])
            rf_chroms = set([item for row in unique_rf_chroms for item in row])
            lf_chrom_lst = []
            for chrom in lf_chroms:
                lf_row_lst = []
                for i,row in enumerate(unique_lf_chroms):
                    if chrom in row:
                        lf_row_lst.append(i)
                lf_chrom_lst.append(lf_row_lst)
            lf_order = sorted(lf_chrom_lst, key=itemgetter(0))
            lf_common_elements = functools.reduce(lambda a,b: a+b if len(list(set(a).intersection(b))) > 0 else a, lf_order)
            lf_row_group = list(set(sorted(lf_common_elements)))
            rf_chrom_lst = []
            for chrom in rf_chroms:
                rf_row_lst = []
                for i,row in enumerate(unique_rf_chroms):
                    if chrom in row:
                        rf_row_lst.append(i)
                rf_chrom_lst.append(rf_row_lst)
            rf_order = sorted(rf_chrom_lst, key=itemgetter(0))
            rf_common_elements = functools.reduce(lambda a,b: a+b if len(list(set(a).intersection(b))) > 0 else a, rf_order)
            rf_row_group = list(set(sorted(rf_common_elements)))
            combined_row_group = list(set(lf_row_group + rf_row_group))
            combined_grp_df = read_df.iloc[combined_row_group]
            most_reads_df = combined_grp_df[combined_grp_df['Read_count'] == max(combined_grp_df['Read_count'])]
            not_most_reads_df = combined_grp_df[combined_grp_df['Read_count'] != max(combined_grp_df['Read_count'])]
            not_most_reads_df_idx = not_most_reads_df.index.tolist()
            merge_df.drop(not_most_reads_df_idx, axis=0, inplace=True)
            merge_df.to_csv(file_name, index=False, sep = '\t', header=True)
    else:
        merge_df.to_csv(file_name, index=False, header=False)
def main():
    inputs = parse_args()
    ins_file =  inputs.ins_file
    flank_file = inputs.flank_file
    prefix = inputs.output
    file_name = f"{prefix}.ins_flank.report.txt"
    merged_file = merge_files(ins_file, flank_file)
    filter_merged_file(merged_file, file_name)

if __name__ =="__main__":
    main()