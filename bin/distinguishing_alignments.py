#!/usr/bin/env python3 
# 4-space indented, v0.0.1
# File name: distinguishing_alignments.py
# Description: Add a distinguishing mark to reads that are mapped to multiple locations to enable downstream filtering.
# Author: Robert Linder
# Date: 2024-01-04

import argparse
import pysam
import re 

def parse_args():
	"""this enables users to provide input as command line arguments to minimize the need for directly modifying the script; please enter '-h' to see the full list of options"""
	parser = argparse.ArgumentParser(description="Load file for marking primary vs secondary alignments")
	parser.add_argument("bam_file", type=str, help="bam file to process")
	args = parser.parse_args()
	return args

def mark_alignments(bam, output_name):
    """Iterate through the sorted/indexed bam file to mark multiple alignments so they can be distinguished in downstream processes"""
    bamfile = pysam.AlignmentFile(bam, "rb")
    output_file = pysam.AlignmentFile(output_name + ".marked.bam", "wb", template=bamfile)
    reads_processed = []
    for read in bamfile.fetch():
            reads_processed.append(read.query_name)
            #print(read.query_name)
            number = reads_processed.count(read.query_name)
            tags = read.get_tags()
            tp = [t[1] for t in tags if t[0] == "tp"]
            as_score = [t[1] for t in tags if t[0] == "AS"]
            if re.match("Sniffles2", read.query_name): 
                read.query_name = f"{read.query_name}_{str(number)}_{tp[0]}_{as_score[0]}"
            else:
                read.query_name = f"{read.query_name}_{str(number)}_{tp[0]}_{as_score[0]}_NativeHit"
            output_file.write(read)
    output_file.close()
    bamfile.close()

def main():
    inputs = parse_args()
    bam =  inputs.bam_file
    output = bam.split('/')[-1].split('.')[:-1]
    output_name = '.'.join(output)
    mark_alignments(bam, output_name)
    pysam.sort("-m", "4G", "-@", "4", "-o", output_name + ".marked.sorted.bam",  output_name + ".marked.bam")

if __name__ =="__main__":
    main()
