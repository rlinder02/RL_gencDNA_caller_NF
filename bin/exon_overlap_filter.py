#!/usr/bin/env python3 
# 4-space indented, v0.0.1
# File name: exon_overlap_filter.py
# Description: Filter by the amount of exon overlapped by reads.
# Author: Robert Linder
# Date: 2024-03-07

import argparse

def parse_args():
	"""this enables users to provide input as command line arguments to minimize the need for directly modifying the script; please enter '-h' to see the full list of options"""
	parser = argparse.ArgumentParser(description="Load bed file for filtering by amount of exon overlap")
	parser.add_argument("bed_file", type=str, help="bed file to process")
	parser.add_argument('-o', "--overlap", type=float, help="The insertion must overlap at least this fraction of an exon.")
	args = parser.parse_args()
	return args

def exon_overlap_filter(bed, min_overlap, output_name):
	"""Iterate through the sorted/indexed bam file to mark multiple alignments so they can be distinguished in downstream processes"""
	with open(f"{output_name}.overlap.bed", 'w') as output_file:
		with open(bed ,'r') as bedfile:
			for line in bedfile:
				fields = line.split('\t')
				overlap_length = int(fields[-1])
				exon_length = int(fields[-6]) - int(fields[-7]) + 1
				overlap_frac = overlap_length/exon_length
				# Keep hits that overlap exons as well as hits that overlap repeats (without specifying the length of the overlap- this is done later in the repeat filter module)
				if overlap_frac >= min_overlap or fields[13] == "RepeatMasker":
					output_file.write(line) 

def main():
	inputs = parse_args()
	bed = inputs.bed_file
	min_overlap = inputs.overlap
	output = bed.split('/')[-1].split('.')[:-1]
	output_name = '.'.join(output)
	exon_overlap_filter(bed, min_overlap, output_name)

if __name__=="__main__":
	main()