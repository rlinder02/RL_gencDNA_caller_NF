#!/usr/bin/env python3 
# 4-space indented, v0.0.1
# File name: extract_ins_reads.py
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
    parser.add_argument("ins_pos_file", type=str, help="Text file of insertion locations")
    parser.add_argument("bam_file", type=str, help="Sorted and indexed bam file generated from the original minimap2 alignment")
    args = parser.parse_args()
    return args

def find_ins_reads(bam, ins_locs, file_name):
    """Extract reads from which the insertions calls originated from"""
    ins_loc_df = pd.read_csv(ins_locs, header = None, sep = "\t")
    rf_start_finder = re.compile("END=.*")
    strand_finder = re.compile("STRAND=.*")
    info_field = ins_loc_df[3].str.split(';').values
    rf_locs_lst = []
    strand_info = []
    for info in info_field:
        rf_locs = [loc.split('=')[1] for loc in info if rf_start_finder.match(loc)]
        rf_locs_lst.append(rf_locs[0])
        strand = [s.split('=')[1] for s in info if strand_finder.match(s)]
        strand_info.append(strand)
    ins_loc_df['end'] = rf_locs_lst
    ins_loc_df['strand'] = strand_info
    bamfile = pysam.AlignmentFile(bam, "rb")
    reads_to_keep = {}
    for read in bamfile.fetch():
        key = read.query_name
        # check if incoming bam file is a placeholder
        if key != "fake_read":
            name = read.query_name.split('-')[1]
            read_chrom = bamfile.get_reference_name(read.reference_id)
            tags = read.get_tags()
            cigar = read.cigarstring
            SA_tag = [t for t in tags if 'SA' in t[0]]
            tp_tag = [t[1] for t in tags if 'tp' in t[0]]
            ins = re.findall(r"(\d+)I", read.cigarstring)
            introns = re.findall(r"N", read.cigarstring)
            proper_ins = [int(i) for i in ins if int(i) > 30]
            if name in ins_loc_df[2].values:
                ins_df = ins_loc_df[ins_loc_df[2] == name]
                ins_strand = ins_df['strand'].values[0][0]
                same_strand = False
                if ins_strand == '-' and read.is_reverse:
                    same_strand = True
                elif ins_strand == '+' and not read.is_reverse:
                    same_strand = True
                elif ins_strand == '+-':
                    same_strand = True
                if read_chrom == ins_df[0].values[0] and same_strand:
                    aligned_pairs = read.get_aligned_pairs()         
                    ins_start_ref = int(ins_df[1].values[0])
                    ins_end_ref = int(ins_df['end'].values[0]) + 1
                    ins_start_query = [s[0] for s in aligned_pairs if s[1] == ins_start_ref]
                    ins_end_query = [s[0] for s in aligned_pairs if s[1] == ins_end_ref]
                    # Keep the read alignment if can find the insertion start and end positions in the alignment; as there may be multiple alignments that fit these parameters, prioritize the higher quality alignments, write the others out labelled as -SA; need to change so put into dictionary then loop through again
                    if ins_start_query and ins_end_query:
                        if len(proper_ins) > 0 and len(introns) > 0 and tp_tag[0] == 'P':
                            reads_to_keep[key] = [cigar, 5]
                        elif len(proper_ins) > 0 and tp_tag[0] == 'P' and reads_to_keep.get(key, [0,0])[1] < 5:
                            reads_to_keep[key] =  [cigar, 4]
                        elif len(proper_ins) > 0 and len(introns) > 0 and reads_to_keep.get(key, [0,0])[1] < 4:
                            reads_to_keep[key] = [cigar, 3]
                        elif len(proper_ins) > 0 and reads_to_keep.get(key, [0,0])[1] < 3:
                            reads_to_keep[key] =  [cigar, 2]
                        elif len(SA_tag) > 0 and reads_to_keep.get(key, [0,0])[1] < 2:
                            reads_to_keep[key] =  [cigar, 1]
                        elif reads_to_keep.get(key, [0,0])[1] < 1:
                            reads_to_keep[key] =  [cigar, 0]
    bamfile.close()
    if len(reads_to_keep) > 0:
        bamfile = pysam.AlignmentFile(bam, "rb")
        output_file = pysam.AlignmentFile(file_name + ".ins.reads.bam", "wb", template=bamfile)
        reads_to_keep_cigar = {k:v[0] for (k, v) in reads_to_keep.items()}
        for read in bamfile.fetch():
            matching_cigar = reads_to_keep_cigar.get(read.query_name)
            tags = read.get_tags()
            SA_tag = [t for t in tags if 'SA' in t[0]]
            if read.cigarstring == matching_cigar:
                output_file.write(read)
            elif len(SA_tag) > 0:
                read.query_name = f"{read.query_name}-SA"
                output_file.write(read)
        bamfile.close()
        output_file.close()
    # if there are no insertions left after filtering, continue to propagate a place-holder bam file
    else:
        header = { 'HD': {'VN': '1.0'},
            'SQ': [{'LN': 1575, 'SN': 'chr1'},
                   {'LN': 1584, 'SN': 'chr2'}] }
        with pysam.AlignmentFile(file_name + ".ins.reads.bam", "wb", header=header) as outf:
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
    bam = inputs.bam_file
    pysam.index(bam)
    ins_locs = inputs.ins_pos_file
    file_name = bam.split('/')[-1].split('.')[:-1]
    file_name = '.'.join(file_name)
    find_ins_reads(bam, ins_locs, file_name)

if __name__ =="__main__":
    main()