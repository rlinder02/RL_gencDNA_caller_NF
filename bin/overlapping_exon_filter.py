#!/usr/bin/env python3 
# 4-space indented, v0.0.1
# File name: overlapping_exon_filter.py
# Description: Filters out exon-exon joins that span divergent transcripts that overlap in genome. This will let pass exons overlapping different genes from different chromosomes on the same insertion, so long as the alignment is the same - multiply aligned insertions would be filtered out if solo.
# Author: Robert Linder
# Date: 2023-11-29

import pandas as pd
import re
import os
import argparse
from subprocess import Popen, PIPE, check_output
from itertools import combinations
import networkx as nx
import polars as pl
import pysam
from operator import itemgetter
import pyarrow

def parse_args():
    """this enables users to provide input as command line arguments to minimize the need for directly modifying the script; please enter '-h' to see the full list of options"""
    parser = argparse.ArgumentParser(description="Load file for final filtration step of the gencDNA caller")
    parser.add_argument("gencDNA_file", type=str, help="potential gencDNA bed file to filter")
    parser.add_argument("bam_file", type=str, help="Insertion bam file")
    parser.add_argument("repeat_bed", type=str, help="Bed file of overlaps with repeats")
    parser.add_argument('-g', '--gap', type=float, default=0.2, help="The fraction of the distance between adjacent exons on the read vs in the reference allowed")
    parser.add_argument('-r', '--repeat_overlap', type=float, default=0.2, help="The fraction of the overlap between a repeat and the length of exons within an alignement to filter out")
    args = parser.parse_args()
    return args

def filter_same_exon_overlaps(file):
    """This filters out exon-exon junctions that span annotated overlapping exons from different transcripts of the same gene as well as removing the additional pseudogene annotations to avoid downstream conflicts as hits at known pseudogenes should have been filtered out by this step"""
    genc_df = pd.read_csv(file, sep='\t', header=None)
    genc_df = genc_df[genc_df[13] != "data"]
    genc_df = genc_df.reset_index(drop=True)
    genc_df.rename(columns={22: 'transcript', 23: 'gene', 24: 'gene_type', 25: 'gene_name', 26: 'chr', 27: 'overlap_rows_1', 28: 'nonoverlap_rows_1', 29: 'alignment', 30: 'tp', 31: 'AS'}, inplace=True)
    # create a list of dataframes grouped by column 3
    grouped_df_lst = [df for _, df in genc_df.groupby(3)]
    # Make all pairwise comparisons between the individual overlaps within each group
    for read_df in grouped_df_lst:
        align_dfs = [df for _, df in read_df.groupby('alignment')]
        for align_df in align_dfs:
            gene_dfs = [df for _, df in align_df.groupby('gene')]
            for gene_df in gene_dfs:
                chr_dfs = [df for _, df in gene_df.groupby('chr')]
                for chr_df in chr_dfs: # add this to distinguish b/w alt haplotypes
                # check if different transcripts from the same gene are present
                    if chr_df.transcript.nunique() != 1:
                        # By default, rows that don't meet the following criteria will be marked as 0; rows marked as 1 will be deleted from the final .bed file output 
                        exon_intervals = chr_df.apply(lambda row: pd.Interval(row[15], row[16], closed='both'), axis=1)
                        exon_combos = [e for e in combinations(exon_intervals, 2)]
                        row_combos = [r for r in combinations(chr_df.index, 2)]
                        counter = -1
                        overlap_list = []
                        for combo in exon_combos:
                            counter += 1
                            if combo[0].overlaps(combo[1]):
                                overlapping_rows = row_combos[counter]
                                rows_overlapping = [int(overlapping_rows[0]), int(overlapping_rows[1])]
                                overlap_list.append(rows_overlapping)
                        # if "Sniffles2.INS.380S12.1" in chr_df[3].values:
                        G=nx.Graph()
                        G.add_edges_from(overlap_list)
                        connected_rows = [n for n in nx.connected_components(G)]
                        grp_counter = 0
                        for grp in connected_rows:
                            grp_counter += 1
                            genc_df.loc[genc_df.index[list(grp)], 'nonoverlap_rows_1'] = grp_counter
                            genc_df.loc[genc_df.index[list(grp)], 'overlap_rows_1'] = 1
    return(genc_df)
    
def keep_nonoverlapping_intervals(genc_df):
    """This keeps the first row for each nonoverlapping group and reverts the overlapping rows value to 0 so it's not deleted"""
    grouped_df_lst = [df for _, df in genc_df.groupby(3)]
    filtered_data = []
    for read_df in grouped_df_lst:
        align_dfs = [df for _, df in read_df.groupby('alignment')]
        for align_df in align_dfs:
            gene_dfs = [df for _, df in align_df.groupby('gene')]
            for gene_df in gene_dfs:
                chr_dfs = [df for _, df in gene_df.groupby('chr')]
                for chr_df in chr_dfs:
                    if any(chr_df.nonoverlap_rows_1 > 0):
                        subset_dfs = [df for _, df in chr_df.groupby('nonoverlap_rows_1')]
                        for sub_df in subset_dfs:
                            first_row = sub_df.iloc[:1]
                            filtered_data.append(first_row)
                    else:
                        filtered_data.append(chr_df)
    filtered_df = pd.concat(filtered_data)
    filtered_df = filtered_df.reset_index(drop=True)
    filtered_df.loc[filtered_df['nonoverlap_rows_1'] > 0, 'overlap_rows_1'] = 0
    return(filtered_df)

def correct_overlapping_rows(filtered_df):
    """Revert the overlap_rows_1 value to 0 if the insertion overlaps multiple different genes"""
    filtered_df['chr_char'] = filtered_df['chr'].astype(str).str[0]
    grouped_df_lst = [df for _, df in filtered_df.groupby(3)]
    for read_df in grouped_df_lst:
        align_dfs = [df for _, df in read_df.groupby('alignment')]
        for align_df in align_dfs:
            chr_dfs = [df for _, df in align_df.groupby('chr_char')]
            for chr_df in chr_dfs:
                if chr_df.gene.nunique() != 1:
                        for gene in chr_df.gene.unique():
                            first_row = chr_df.loc[chr_df['gene'] == gene].index[0]
                            filtered_df.loc[filtered_df.index[first_row], 'overlap_rows_1'] = 0 
    return(filtered_df)

def filter_divergent_transcripts(same_exon_df):
     """This filters out exon-exon junctions that span known adjacent exons from divergent transcripts coded on different strands"""
     grouped_df_lst = [df for _, df in same_exon_df.groupby(3)]
     for read_df in grouped_df_lst:
        align_dfs = [df for _, df in read_df.groupby('alignment')]
        for align_df in align_dfs:
            chr_dfs = [df for _, df in align_df.groupby('chr')]
            for chr_df in chr_dfs:
                if chr_df.gene.nunique() > 1:
                    exon_intervals = chr_df.apply(lambda row: pd.Interval(row[15], row[16], closed='both'), axis=1) 
                    exon_combos = [e for e in combinations(exon_intervals, 2)]
                    row_combos = [r for r in combinations(chr_df.index, 2)]
                    gene_combos = [g for g in combinations(chr_df.gene, 2)]
                    counter = -1
                    overlap_list = []
                    for combo in exon_combos:
                        counter += 1
                        if combo[0].overlaps(combo[1]):
                            if gene_combos[counter][0] != gene_combos[counter][1]:
                                overlapping_rows = row_combos[counter]
                                rows_overlapping = [int(overlapping_rows[0]), int(overlapping_rows[1])]
                                overlap_list.append(rows_overlapping)
                    G=nx.Graph()
                    G.add_edges_from(overlap_list)
                    connected_rows = [n for n in nx.connected_components(G)]
                    for grp in connected_rows:
                        same_exon_df.loc[same_exon_df.index[list(grp)], 'overlap_rows_1'] = 1
     return(same_exon_df)

def correct_overlapping_divergent_rows(filtered_df):
    """Revert the overlap_rows_1 value to 0 if the there are multiple exons for divergent transcripts"""
    grouped_df_lst = [df for _, df in filtered_df.groupby(3)]
    for read_df in grouped_df_lst:
        align_dfs = [df for _, df in read_df.groupby('alignment')]
        for align_df in align_dfs:
            gene_dfs = [df for _, df in align_df.groupby('gene')]
            for gene_df in gene_dfs:
                chr_dfs = [df for _, df in gene_df.groupby('chr')]
                for chr_df in chr_dfs:
                    if len(chr_df) > 1:
                        filtered_df.loc[filtered_df.index[chr_df.index], 'overlap_rows_1'] = 0 
    return(filtered_df)

def single_exon_filter(filtered_df):
    # Mark for deletion rows in which there is only a single exon overlap detected (eg, a deletion within an exon)
    grouped_df_lst = [df for _, df in filtered_df.groupby(3)]
    for read_df in grouped_df_lst:
        align_dfs = [df for _, df in read_df.groupby('alignment')]
        for align_df in align_dfs:
            chr_dfs = [df for _, df in align_df.groupby('chr_char')]
            for chr_df in chr_dfs:
                if len(chr_df[chr_df['overlap_rows_1'] == 0]) == 1:
                    find_rows = chr_df.index
                    filtered_df.loc[filtered_df.index[find_rows], 'overlap_rows_1'] = 1
    return(filtered_df)

def filter_different_alignments(filtered_df):
    """This removes rows in which the reads supporting the hit are different alignments"""
    grouped_df_lst = [df for _, df in filtered_df.groupby(3)]
    for read_df in grouped_df_lst:
        primary_aligns = read_df[read_df['tp'] == 'P']
        if len(primary_aligns) > 0:
            secondary_aligns = read_df[read_df['tp'] == 'S']
            secondary_idxs = secondary_aligns.index
            filtered_df.loc[filtered_df.index[secondary_idxs], 'overlap_rows_1'] = 1
        else:
            best_as_score = max(read_df['AS'])
            best_as_rows = read_df[read_df['AS'] == best_as_score]
            not_best_rows = read_df[read_df['AS'] != best_as_score]
            not_best_idxs = not_best_rows.index
            filtered_df.loc[filtered_df.index[not_best_idxs], 'overlap_rows_1'] = 1
            if best_as_rows.alignment.nunique() > 1:
                default_choice = min(best_as_rows['alignment'])
                not_default_choice = best_as_rows[best_as_rows['alignment'] != default_choice]
                not_default_idxs = not_default_choice.index
                filtered_df.loc[filtered_df.index[not_default_idxs], 'overlap_rows_1'] = 1
        align_dfs = [df for _, df in read_df.groupby('alignment')]
        for align_df in align_dfs:
            if len(align_df[align_df['overlap_rows_1'] == 0]) == 1:
                find_row = align_df.index
                filtered_df.loc[filtered_df.index[find_row], 'overlap_rows_1'] = 1
    return(filtered_df)

def single_exon_filter2(filtered_df):
    # Mark for deletion rows in which there is only a single exon overlap detected after accounting for multiple alignments
    grouped_df_lst = [df for _, df in filtered_df.groupby(3)]
    for read_df in grouped_df_lst:
        align_dfs = [df for _, df in read_df.groupby('alignment')]
        for align_df in align_dfs:
            chr_dfs = [df for _, df in align_df.groupby('chr_char')]
            for chr_df in chr_dfs:
                if len(chr_df[chr_df['overlap_rows_1'] == 0]) == 1:
                    find_rows = chr_df.index
                    filtered_df.loc[filtered_df.index[find_rows], 'overlap_rows_1'] = 1
    return(filtered_df)

def remove_flagged_overlaps(different_alignments_df):
    final_transcripts_df = different_alignments_df[different_alignments_df.overlap_rows_1 == 0]
    #final_transcripts_df.to_csv("test.bed", sep = "\t", header = True)
    return(final_transcripts_df)
    
def exon_exon_gaps(genc, bam, exon_gap):
    with pysam.AlignmentFile(bam, "rb") as bamfile:
        genc_df = pl.from_pandas(genc)
        genc_df = genc_df.with_columns(
                pl.concat_str(
                    [pl.col("3"), pl.col("alignment"), pl.col("tp"), pl.col("AS")], separator="_"
                ).alias("Insertion_ID")
            )
        for read in bamfile.fetch():
            if 'NativeHit' in read.query_name:
                name1, name2, name3, name4, alignment, tp, score, read_type = read.query_name.split("_")
                new_query_name = '_'.join([name1, name2, name3, name4, read_type, alignment, tp, score])
                read.query_name = new_query_name
            if read.query_name in genc_df.get_column("Insertion_ID").to_list():
                ins_df = genc_df.filter(pl.col("Insertion_ID") == read.query_name)
                aligned_pairs = read.get_aligned_pairs()
                ins_df = ins_df.with_columns(pl.col("15").map_elements(lambda x: min([tup[0] for tup in aligned_pairs if tup[0] and tup[1] and tup[1] >= x]), return_dtype=pl.Int64).alias('Overlap_start'))
                ins_df = ins_df.with_columns(pl.col("16").map_elements(lambda x: max([tup[0] for tup in aligned_pairs if tup[0] and tup[1] and tup[1] <= x]), return_dtype=pl.Int64).alias('Overlap_end')).sort("Overlap_start")
                ins_df = ins_df.with_columns((pl.col("Overlap_start").shift(-1) - pl.col("Overlap_end")).alias("Distance"))                
                ins_df = ins_df.with_columns((pl.col("15").shift(-1) - pl.col("16")).alias("Exon_Distance"))
                ins_df = ins_df.with_row_index("Index")
                for i, row in enumerate(ins_df.iter_rows(named=True)):
                    if row['Distance']:
                        if row['Distance'] > 10:
                            old_distance = row['Distance']
                            next_start = ins_df.row(i+1, named=True)
                            next_start = next_start['Overlap_start']
                            current_end = row['Overlap_end']
                            between_exons = []
                            query_length = 0
                            for o,l in read.cigartuples:
                                if o in [1,4,7,8]:
                                    query_length += l
                                # Have to add 1 to current end to ensure not getting the operations before or after the split
                                if query_length > (current_end+1) and query_length < (next_start):
                                    between_exons.append([l,o])
                            between_exons = between_exons[:-1]
                            find_n = [n for n in between_exons if n[1] == 3]
                            # If there are multiple introns before the next exon, there was likely a mapping issue; if there is more than a single intron between two adjacent exons, then proceed; a single spliced intron between two exons is normal 
                            if len(find_n) > 1:
                                insertion_lengths = [n[0] for n in between_exons if n[1] == 1]
                                match_lengths = [n[0] for n in between_exons if n[1] == 7]
                                mismatch_lengths = [n[0] for n in between_exons if n[1] == 8]
                                total_length = sum(insertion_lengths + match_lengths + mismatch_lengths)
                                new_distance = old_distance - total_length
                                ins_df = ins_df.with_columns(
                                    pl.when(pl.col("Index") == i).then(new_distance).otherwise(pl.col("Distance")).alias("Distance")
                                )
                ins_df = ins_df.head(-1)
                ins_df = ins_df.with_columns((pl.col("Distance") / pl.col("Exon_Distance")).alias("Distance_frac"))
                max_distance = ins_df.select(pl.max("Distance")).to_numpy()[0][0]
                #max_distance_frac = ins_df.select(pl.max("Distance_frac")).to_numpy()[0][0]
                min_distance_frac = ins_df.select(pl.min("Distance_frac")).to_numpy()[0][0]
                min_exon_distance = ins_df.select(pl.min("Exon_Distance")).to_numpy()[0][0]
                # if any annotated exons are less than 20 bp apart, there is likely an error in the annotation file for this gene
                if min_exon_distance < 20:
                    insertion = list(set(ins_df.get_column("3").to_list()))[0]
                    genc_df = genc_df.filter(pl.col("3") != insertion)
                elif max_distance > 10:
                    # If any of the distances between exons on the read goes above 10 bp, then only filter if the distance on the read is more than 50% of the distance in the genome, as this likely represents a real splice event still; even if there are multiple exons, but there is only a splice event between just two exons, keep this as a potential hit
                    if min_distance_frac > exon_gap:
                        insertion = list(set(ins_df.get_column("3").to_list()))[0]
                        genc_df = genc_df.filter(pl.col("3") != insertion)
        #genc_df.write_csv("test_after_gap_filter.bed", separator = '\t', include_header = True)
        return(genc_df)
        

def repeat_filter(genc_df, repeats, repeat_overlap, file_name):
    """Filters out insertions in which there is an alignment in which the fraction of the total length of repeats overlapping exons over the total length of exons is greater than or equal to the repeat_overlap parameter"""
    repeats_df = (
         pl.read_csv(repeats, has_header=False, separator='\t')

    )
    exon_lens_df = genc_df.with_columns((pl.col("16")-pl.col("15")).alias("Exon_Length"))
    joined_df = exon_lens_df.join(repeats_df, left_on = "Insertion_ID", right_on = "column_4", how = "inner").filter(pl.col("0") == pl.col("column_1")).with_columns(overlap_min=pl.min_horizontal("16", "column_17")).with_columns(overlap_max=pl.max_horizontal("15", "column_16")).with_columns((pl.col("overlap_min") - pl.col("overlap_max") + 1).alias("Repeat_overlap")).with_columns(pl.when(pl.col("Repeat_overlap") < 0).then(0).otherwise(pl.col("Repeat_overlap")).alias("Repeat_overlap_corrected")).group_by(["Insertion_ID", "column_21"]).agg(pl.col("Repeat_overlap_corrected").sum().alias("Total_repeat_overlap"), pl.col("Exon_Length").sum().alias("Total_exons_length")).group_by("Insertion_ID").agg(pl.col("Total_repeat_overlap").sum(), pl.col("Total_exons_length").first()).with_columns((pl.col("Total_repeat_overlap")/pl.col("Total_exons_length")).alias("Fraction_overlap")).filter(pl.col("Fraction_overlap") >= repeat_overlap).with_columns(pl.col("Insertion_ID").str.extract(r"(Sniffles2.INS.\w+.\d+)").alias("Insertion"))
    #print(joined_df.filter(pl.col("Insertion_ID") == "Sniffles2.INS.86S6.2_1_P_431"))
    #joined_df.write_csv("Test.filtered2.bed", separator = '\t', include_header = True)
    ins_to_filter = list(set(joined_df.get_column("Insertion").to_list()))
    genc_df = genc_df.filter(~pl.col("3").is_in(ins_to_filter))
    genc_df = genc_df.drop(["Insertion_ID"])
    genc_df.write_csv(f"{file_name}.filtered.bed", separator = '\t', include_header=False)
    return(repeats_df)
    
def line1_ins(repeat_df, file_name):
    line1_df = repeat_df.filter(pl.col("column_21").str.contains("LINE/L1"))
    line1_df.write_csv(f"{file_name}.LINE1.bed", separator = '\t', include_header=False)

def main():
    inputs = parse_args()
    file =  inputs.gencDNA_file
    bam = inputs.bam_file
    repeats = inputs.repeat_bed
    exon_gap = inputs.gap
    repeat_overlap = inputs.repeat_overlap
    pysam.index(bam)
    file_name = file.split('/')[-1].split('.')[:-1]
    file_name = '.'.join(file_name)
    same_exon_df = filter_same_exon_overlaps(file)
    nonoverlap_df = keep_nonoverlapping_intervals(same_exon_df)
    corrected_rows_df = correct_overlapping_rows(nonoverlap_df)
    diverged_transcripts_df = filter_divergent_transcripts(corrected_rows_df)
    corrected_diverged_rows_df = correct_overlapping_divergent_rows(diverged_transcripts_df)
    multiple_exons_df = single_exon_filter(corrected_diverged_rows_df)
    different_alignments_df = filter_different_alignments(multiple_exons_df)
    single_exon_filtered_df = single_exon_filter2(different_alignments_df)
    removed_overlaps_df = remove_flagged_overlaps(single_exon_filtered_df)
    gapped_exons = exon_exon_gaps(removed_overlaps_df, bam, exon_gap)
    filtered_repeats = repeat_filter(gapped_exons, repeats, repeat_overlap, file_name)
    line1_ins(filtered_repeats, file_name)
if __name__ =="__main__":
    main()