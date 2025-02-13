#!/usr/bin/env python3 
# 4-space indented, v0.0.1
# File name: mapping_insertions.py
# Description: Formats insertions from Sniffles2 output into Fasta files, which can then be aligned to the ref genome using minimap2 in splice-aware mode.
# Author: Robert Linder
# Date: 2023-11-18

import argparse

def parse_args():
    """this enables users to provide input as command line arguments to minimize the need for directly modifying the script; please enter '-h' to see the full list of options"""
    parser = argparse.ArgumentParser(description="Load text file for converting insertion file into fasta format")
    parser.add_argument("ins_file", type=str, help="insertion file to process, which must be a tab-delimited txt file")
    args = parser.parse_args()
    return args 

def convert_to_fasta(file, output_name):
    with open(f"{output_name}.fasta", 'w') as output:
        with open(file, 'r') as f:
            for line in f:
                name = '>' + line.split('\t')[0]
                seq = line.split('\t')[1]
                if not "INS" in seq: # very large insertions are not written out but have a place-holder INS - ask how could retrieve this sequence 
                    output.write(name + "\n" + seq)


def main():
    inputs = parse_args()
    file = inputs.ins_file
    output = file.split('/')[-1].split('.')[:-1]
    output_name = '.'.join(output)
    convert_to_fasta(file, output_name)

if __name__ == "__main__":
    main()