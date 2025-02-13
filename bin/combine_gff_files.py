#!/usr/bin/env python3 
# 4-space indented, v0.0.1
# File name: combine_gff_files.py
# Description: Reformats the bed files.
# Author: Robert Linder
# Date: 2024-03-21

import pandas as pd
import re
import os
import argparse
from subprocess import Popen, PIPE, check_output
from itertools import combinations

def parse_args():
    """this enables users to provide input as command line arguments to minimize the need for directly modifying the script; please enter '-h' to see the full list of options"""
    parser = argparse.ArgumentParser(description="Load file for final filtration step of the gencDNA caller")
    parser.add_argument("ps_gff", type=str, help="The additional pseudogene annotation files in gff format")
    parser.add_argument("species", type=str, help="The species being interogatted")
    args = parser.parse_args()
    return args

def replace_str(s, r):
    return ''.join([r.get(x, x) for x in s])

def reformat_gff(ps_gff, file_name, species):
    with open(f"{file_name}.additional.gff", 'w') as outfile:
        with open(ps_gff, 'r') as ps:
            # Edit the info field for downstream processing of the appropriate primary annotation file for the reference genome used 
            if species == "human":
                for line in ps:
                    if not line.startswith("#"):
                        if 'pseudo' in line or 'exon' in line:
                            fields = line.split('\t')
                            mapping_table = {'%3B': ';', '%3D': '=', 'gene_name': 'gene', 'gene_type': 'product', 'thickStart': 'transcript_id'}
                            for key, value in mapping_table.items():
                                fields[8] = fields[8].replace(key, value)
                            # change the third column to exons so make it through the overlap filter
                            fields[2] = 'exon'
                            new_line = '\t'.join(fields)
                            outfile.write(new_line)
                        # for the repeat gffs generated from the same genome as used for alignment
                        elif 'RepeatMasker' in line:
                            fields = line.split('\t')
                            fields[2] = 'exon'
                            new_line = '\t'.join(fields)
                            outfile.write(new_line)
            elif species == "mouse":
                for line in ps:
                    if not line.startswith("#"):
                        if 'pseudo' in line or 'IG_' in line or 'exon' in line:
                            fields = line.split('\t')
                            mapping_table = {'%3B': ';', '%3D': '=', '"': '', 'blockCount=gene_id': 'gene_id=', 'gene_biotype': 'gene_type=', 'transcript_id': 'Parent=transcript:', 'gene_name': 'common_name='}
                            for key, value in mapping_table.items():
                                fields[8] = fields[8].replace(key, value)
                            if 'common_name' not in fields[8]:
                                gene_id = re.compile("gene_id=.*")
                                info_split = fields[8].split(';')
                                gene_match = [gene.split('=')[1] for gene in info_split if gene_id.match(gene)]
                                fields[8] = f"{fields[8].rstrip()};common_name={gene_match[0]}\n"
                            else:
                                fields[8] = fields[8]
                            # change the third column to exons so make it through the overlap filter
                            fields[2] = 'exon'
                            new_line = '\t'.join(fields)
                            # only keep lines which contain actual exons for the immunoglobulin genes and exons of known genes; not for pseudogenes
                            if 'IG_' in line or 'exon' in line:
                                if "exon_id" in new_line:
                                    outfile.write(new_line)
                            else:
                                outfile.write(new_line)
                        # for the repeat gffs generated from the same genome as used for alignment
                        elif 'RepeatMasker' in line:
                            fields = line.split('\t')
                            fields[2] = 'exon'
                            new_line = '\t'.join(fields)
                            outfile.write(new_line)

def main():
    inputs = parse_args()
    ps_gff = inputs.ps_gff
    species = inputs.species
    file_name = ps_gff.split('/')[-1].split('.')[:-1]
    file_name = '.'.join(file_name)
    reformat_gff(ps_gff, file_name, species)

if __name__ =="__main__":
    main()
