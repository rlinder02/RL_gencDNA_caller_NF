#!/usr/bin/env python3 
# 4-space indented, v0.0.1
# File name: filter_insertion_fasta.py
# Description: Filters the insertion fasta file for the final list of insertions keeping.
# Author: Robert Linder
# Date: 2024-01-10

import pandas as pd
import re
import os
import argparse
from subprocess import Popen, PIPE, check_output
from itertools import combinations
from Bio import SeqIO


def parse_args():
    """this enables users to provide input as command line arguments to minimize the need for directly modifying the script; please enter '-h' to see the full list of options"""
    parser = argparse.ArgumentParser(description="Load files for final filtration step of the gencDNA caller")
    parser.add_argument("ins_names", type=str, help="file with list of insertions keeping")
    parser.add_argument("fasta_file", type=str, help="fasta file of insertion sequences")
    args = parser.parse_args()
    return args

def extract_fasta_entries(genc, fasta, file_name):
    """This extracts a filtered list of fasta entries"""
    ins_df = pd.read_csv(genc, header = None)
    insertions = ins_df[0].to_string(index=False)
    with open(f"{file_name}.filtered.fasta", 'w') as output:
        pass
    for record in SeqIO.parse(fasta, "fasta"):
        if record.description in insertions:
            with open(f"{file_name}.filtered.fasta", 'a') as output:
                output.write(record.format("fasta"))
                
def main():
    inputs = parse_args()
    genc =  inputs.ins_names
    fasta = inputs.fasta_file
    file_name = fasta.split('/')[-1].split('.')[:-1]
    file_name = '.'.join(file_name)
    extract_fasta_entries(genc, fasta, file_name)

if __name__ =="__main__":
    main()