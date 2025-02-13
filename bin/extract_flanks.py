#!/usr/bin/env python3 
# 4-space indented, v0.0.1
# File name: extract_flanks.py
# Description: Extract flanking regions of insertions.
# Author: Robert Linder
# Date: 2024-01-11

import pandas as pd
import re
import os
import argparse
from subprocess import Popen, PIPE, check_output
from itertools import combinations
import pysam
from Bio import SeqIO
from operator import itemgetter
import numpy as np
import copy

def parse_args():
    """this enables users to provide input as command line arguments to minimize the need for directly modifying the script; please enter '-h' to see the full list of options"""
    parser = argparse.ArgumentParser(description="Load files for final filtration step of the gencDNA caller")
    parser.add_argument("gencDNA_fasta", type=str, help="potential gencDNA fasta file to filter")
    parser.add_argument("ins_pos_file", type=str, help="Text file of insertion locations")
    parser.add_argument("bam_file", type=str, help="Sorted and indexed bam file generated from the original minimap2 alignment that has been filtered")
    #parser.add_argument('-r', '--rlen', type=int, default=1000, help="The amount of the flanking region to extract")
    args = parser.parse_args()
    return args

def left_flank_calc(lf_end, ref_cigar_lengths, ref_cigar_operations, sa_start_ref):
    interval = lf_end
    lf_ref_cigar_lengths = copy.deepcopy(ref_cigar_lengths)
    lf_deletion_idces = [idx for idx in range(len(ref_cigar_operations)) if ref_cigar_operations[idx] == "D"]
    orig_del_lengths = [lf_ref_cigar_lengths[idx] for idx in range(len(lf_ref_cigar_lengths)) if idx in lf_deletion_idces]
    orig_del_zip = zip(lf_deletion_idces, orig_del_lengths)
    for i in lf_deletion_idces:
        lf_ref_cigar_lengths[i] = 0
    lf_ref_cigar_lengths_cumsum = np.cumsum(lf_ref_cigar_lengths)
    find_idx = min([i for i in range(len(lf_ref_cigar_lengths_cumsum)) if lf_ref_cigar_lengths_cumsum[i] >= interval])
    find_diff_array = lf_ref_cigar_lengths_cumsum[:(find_idx+1)]
    find_diffs = interval - find_diff_array
    if any(find_diffs > 0):
        leftover_idx = max([i for i in range(len(find_diffs)) if find_diffs[i] >= 0])
        leftover = find_diffs[leftover_idx]
    else:
        leftover_idx = -1
        leftover = interval
    if leftover > 0:
        lf_ref_cigar_lengths[leftover_idx+1] = leftover
    for i,v in orig_del_zip:
        lf_ref_cigar_lengths[i] = v
    ref_add = 0
    lf_ref_end_find = zip(lf_ref_cigar_lengths[:(leftover_idx+2)], ref_cigar_operations[:(leftover_idx+2)])
    for l,o in lf_ref_end_find:
        if o in ['M', 'D', 'X', '=']:
            ref_add += l
    sa_end_ref = int(sa_start_ref) + int(ref_add)  
    return sa_end_ref

def right_flank_calc(rf_ref_start, ref_cigar_lengths, ref_cigar_operations, sa_end_ref):
    interval = rf_ref_start
    rf_ref_cigar_lengths = copy.deepcopy(ref_cigar_lengths)
    rf_deletion_idces = [idx for idx in range(len(ref_cigar_operations)) if ref_cigar_operations[idx] == "D"]
    orig_del_lengths = [rf_ref_cigar_lengths[idx] for idx in range(len(rf_ref_cigar_lengths)) if idx in rf_deletion_idces]
    orig_del_zip = zip(rf_deletion_idces, orig_del_lengths)
    for i in rf_deletion_idces:
        rf_ref_cigar_lengths[i] = 0
    rf_ref_cigar_lengths_cumsum = np.cumsum(rf_ref_cigar_lengths[::-1]) 
    find_idx = min([i for i in range(len(rf_ref_cigar_lengths_cumsum)) if rf_ref_cigar_lengths_cumsum[i] >= interval])
    find_diff_array = rf_ref_cigar_lengths_cumsum[:(find_idx+1)]
    find_diffs = interval - find_diff_array
    if any(find_diffs > 0):
        leftover_idx = max([i for i in range(len(find_diffs)) if find_diffs[i] >= 0])
        leftover = find_diffs[leftover_idx]
    else:
        leftover_idx = 0
        leftover = interval
    if leftover > 0:
        rf_ref_cigar_lengths[leftover_idx] = leftover
    for i,v in orig_del_zip:
        rf_ref_cigar_lengths[i] = v
    ref_add = 0
    rf_ref_end_find = zip(rf_ref_cigar_lengths[(leftover_idx):], ref_cigar_operations[(leftover_idx):])
    for l,o in rf_ref_end_find:
        if o in ['M', 'D', 'X', '=']:
            ref_add += l
    sa_start_ref = int(sa_end_ref) - int(ref_add)
    return sa_start_ref 

def find_intervals(flank_df, flank):
    interval_df_lst = []
    grouped_df_lst = [df for _, df in flank_df.groupby('InsertionID')]
    for ins_df in grouped_df_lst:
        chr_dfs = [df for _, df in ins_df.groupby(f"{flank}_chr")]
        for chr_df in chr_dfs:
            chrom = chr_df[f"{flank}_chr"].values[0]
            min_start = min(list(map(int, chr_df[f"{flank}_start"].values)))
            max_end = max(list(map(int, chr_df[f"{flank}_end"].values)))
            align_dict = {'InsertionID': chr_df['InsertionID'].values[0], 'ReadID': ','.join(chr_df['ReadID'].values), f"{flank}_chr": chrom, f"{flank}_start": str(min_start), f"{flank}_end": str(max_end), f"{flank}_aln_order": ','.join(chr_df[f"{flank}_aln_order"].values)}
            align_dict_df = pd.DataFrame(data=align_dict, index = [0])
            interval_df_lst.append(align_dict_df)
    interval_df_merged = pd.concat(interval_df_lst)
    interval_df_merged = interval_df_merged.reset_index(drop=True)
    collapsed_flank_df = interval_df_merged.groupby(interval_df_merged['InsertionID']).agg({'ReadID': ';'.join, f"{flank}_chr": ';'.join, f"{flank}_start": ';'.join, f"{flank}_end": ';'.join, f"{flank}_aln_order": ';'.join})
    unique_reads = collapsed_flank_df['ReadID'].values
    unique_read_lst = []
    for s in unique_reads:
        read_split = re.split("[;,]+", s)
        read_split_dict = list(dict.fromkeys(read_split))
        unique_read_lst.append(';'.join(read_split_dict))
    collapsed_flank_df['ReadID'] = unique_read_lst
    return collapsed_flank_df

def extract_flanks(fasta, bam, ins_locs, file_name):
    """This extracts the read names from the Sniffles2 vcf as a dictionary with the insertion id as the key and read names supporting it as the values"""
    ins_loc_df = pd.read_csv(ins_locs, header = None, sep = "\t")
    rf_start_finder = re.compile("END=.*")
    info_field = ins_loc_df[3].str.split(';').values
    rf_locs_lst = []
    for info in info_field:
        rf_locs = [loc.split('=')[1] for loc in info if rf_start_finder.match(loc)]
        rf_locs_lst.append(rf_locs[0])
    ins_loc_df['end'] = rf_locs_lst
    bamfile = pysam.AlignmentFile(bam, "rb")
    fasta_dict = SeqIO.to_dict(SeqIO.parse(fasta, "fasta"))
    flanks_lst = []
    lf_dict = {}
    rf_dict = {}
    # Initialize a list that will contain the insertion ID, read name, tsd, and whether or not there is a polyA tail present
    tsd_polya_dict = {}
    with open(f"{file_name}.tsd_find.fasta", 'w') as tsd_output:
        for read in bamfile.fetch():
            if read.query_name != "fake_read":
                read_name = read.query_name.split('-')[0]
                name = read.query_name.split('-')[1]
                read_name_id = '-'.join([name, read_name])
                if name in fasta_dict and not "-SA" in read.query_name:
                    aligned_pairs = read.get_aligned_pairs()
                    sequence = read.query_sequence
                    read_length = len(sequence)            
                    ins_df = ins_loc_df[ins_loc_df[2] == name]
                    ins_start_ref = int(ins_df[1].values[0])
                    ins_chrom = ins_df[0].values[0]
                    ins_end_ref = int(ins_df['end'].values[0]) + 1
                    # need to correct the ref start as a couple of bases back from what Sniffles2 says
                    ins_start_ref_corrected = ins_start_ref - 2
                    ins_start_query = [s[0] for s in aligned_pairs if s[0] and s[1] and s[1] == ins_start_ref_corrected]
                    ins_end_query = [s[0] for s in aligned_pairs if s[0] and s[1] and s[1] == ins_end_ref]
                    # if the insertion start position is in an insertion, deletion, or intron, need to shift leftwards
                    if len(ins_start_query) == 0 and ins_start_ref_corrected >= read.reference_start:
                        ins_start_query = [max([tup[0] for tup in aligned_pairs if tup[0] and tup[1] and tup[1] < ins_start_ref_corrected])]
                    elif len(ins_start_query) == 0:
                        ins_start_query = [0]
                    # if the insertion end position is in an insertion, deletion, or intron, need to shift rightwards
                    if len(ins_end_query) == 0 and ins_end_ref < read.reference_end:
                        ins_end_query = [min([tup[0] for tup in aligned_pairs if tup[0] and tup[1] and tup[1] > ins_end_ref])]
                    elif len(ins_end_query) == 0:
                        ins_end_query = [max((tup[0] for tup in aligned_pairs if tup[0] and tup[1]))]
                    chrom = bamfile.get_reference_name(read.reference_id)

                    ##### Can add in here a module to get a fasta which contains for each read, the 100 bases flanking the insertion and the insertion to find the TSD in a downstream process 
                    if ins_start_query[0] > 100 and ins_end_query[0] < (len(read.query_sequence)- 100):
                        seq_name = f">{name}-{read_name}"
                        full_seq = read.query_sequence[(ins_start_query[0]-100):(ins_end_query[0]+101)]
                        # sanity check to ensure only include results that have the insertion present as well 
                        if len(full_seq) > 250:
                            ### Find a TSD if present and polyA tail; find TSD by searching for the longest common substring searching 50bp into the insertions and 100bp into the flank on each side 
                            lf = full_seq[:151]
                            rf = full_seq[-150:]
                            lst_common_substrings = [lf[i:j] for i in range(len(lf)) for j in range(len(lf)) if lf[i:j] in rf]
                            lcs = max(lst_common_substrings, key = len)
                            polya = False
                            if 'A'*10 in rf or 'T'*10 in lf:
                                polya = True
                            if len(lcs) >5:
                                # If the lcs is greater than 5bp in length, write it out along with the source sequence
                                if name in tsd_polya_dict:
                                    tsd_polya_dict[name].append([name, lcs, polya])
                                else:
                                    tsd_polya_dict[name] = [[name, lcs, polya]]
                                tsd_output.write(seq_name + "\n" + full_seq + "\n" + seq_name + "_tsd" + "\n" + lcs + "\n")
                            else:
                                if name in tsd_polya_dict:
                                    tsd_polya_dict[name].append([name, '', polya])
                                else:
                                    tsd_polya_dict[name] = [[name, '', polya]]
                    first_soft_read = list(map(int, re.findall(r"^(\d+)S", read.cigarstring)))
                    last_soft_read = list(map(int, re.findall(r"(\d+)S$", read.cigarstring)))
                    if first_soft_read:
                        if first_soft_read[0] > ins_start_query[0]:
                            lf_end = ins_start_query[0]
                            lf_soft_clip = ins_start_query[0]
                        else:
                            lf_end = first_soft_read[0]
                            lf_soft_clip = first_soft_read[0]
                    if not first_soft_read:
                        lf_soft_clip = 0
                        lf_end = ins_start_query[0]
                    lf_aln_start = read.query_alignment_start
                    lf_aln_end = ins_start_query[0]
                    if lf_aln_start <= lf_aln_end:
                        lf_aln_sequence = sequence[lf_aln_start:lf_aln_end]
                    else:
                        lf_aln_sequence = ''
                    aln_start_ref = int(read.reference_start)
                    if lf_aln_end > 0:
                        aln_end_ref = [s[1] for s in aligned_pairs if s[0] == lf_aln_end]
                    else: 
                        aln_end_ref = [aln_start_ref + 1] 
                    aln_entry = f"{chrom}:{aln_start_ref}-{aln_end_ref[0]}"
                    read_name = f"{read.query_name}-LF"
                    if read_name_id in lf_dict:
                        lf_dict[read_name_id].append([aln_entry, lf_soft_clip, lf_aln_sequence, read.cigarstring])
                    else:
                        lf_dict[read_name_id] = [[aln_entry, lf_soft_clip, lf_aln_sequence, read.cigarstring]]
                    if last_soft_read:
                        if (read_length - last_soft_read[0]) < ins_end_query[0]:
                            rf_start = ins_end_query[0]
                            rf_from_end = read_length - rf_start
                        else:
                            rf_start = read_length - last_soft_read[0]
                            rf_from_end = last_soft_read[0]
                            rf_soft_clip = last_soft_read[0]
                    if not last_soft_read:
                        rf_soft_clip = 0
                    rf_aln_end = read.query_alignment_end
                    rf_aln_start = ins_end_query[0]
                    rf_aln_sequence = sequence[rf_aln_start:rf_aln_end]
                    aln_start_ref = [s[1] for s in aligned_pairs if s[0] == rf_aln_start]
                    aln_end_ref = int(read.reference_end)
                    aln_entry = f"{chrom}:{aln_start_ref[0]}-{aln_end_ref}"
                    read_name = f"{read.query_name}-RF"
                    if read_name_id in rf_dict:
                        rf_dict[read_name_id].append([aln_entry, rf_soft_clip, rf_aln_sequence, read.cigarstring])
                    else:
                        rf_dict[read_name_id] = [[aln_entry, rf_soft_clip, rf_aln_sequence, read.cigarstring]]
                    if read.has_tag("SA"):
                        sa = [tag for tag in read.get_tags() if "SA" in tag]
                        sa_lst = sa[0][1].split(';')[:-1]
                        for s in sa_lst:
                            sa_split = s.split(',')
                            s_chr = sa_split[0]
                            s_start = int(sa_split[1])
                            strand = sa_split[2]
                            cigar = sa_split[3]
                            find_m = list(map(int, re.findall(r"(\d+)M", cigar)))
                            find_d = list(map(int, re.findall(r"(\d+)D", cigar)))
                            find_i = list(map(int, re.findall(r"(\d+)I", cigar)))
                            first_soft_sa = list(map(int, re.findall(r"^(\d+)S", cigar)))
                            last_soft_sa = list(map(int, re.findall(r"(\d+)S$", cigar)))
                            ref_cigar_lengths = list(map(int, re.findall(r"(\d+)[IMSD]", cigar)))
                            ref_cigar_operations = list(re.findall(r"\d+([IMSD])", cigar))
                            query_length = sum(find_m) + sum(find_i)
                            reference_length = sum(find_m) + sum(find_d)
                            s_end = s_start + reference_length
                            # If left flank info could not be found in the primary alignment, see if it can be located in the SAs
                            if lf_end == 0 and ins_chrom == s_chr and ins_start_ref > s_start and ins_start_ref < s_end:
                                lf_reference_end = ins_start_ref_corrected
                                left_flank_ref_length = lf_reference_end - s_start
                                if left_flank_ref_length <= 0:
                                    left_flank_ref_length = ins_start_ref - s_start
                                lf_reference_cigar_lengths = copy.deepcopy(ref_cigar_lengths)
                                lf_insertion_idces = [idx for idx in range(len(ref_cigar_operations)) if ref_cigar_operations[idx] in ["I", "S"]]
                                orig_ins_lengths = [lf_reference_cigar_lengths[idx] for idx in range(len(lf_reference_cigar_lengths)) if idx in lf_insertion_idces]
                                orig_ins_zip = zip(lf_insertion_idces, orig_ins_lengths)
                                for i in lf_insertion_idces:
                                    lf_reference_cigar_lengths[i] = 0
                                lf_reference_cigar_lengths_cumsum = np.cumsum(lf_reference_cigar_lengths)
                                find_idx = min([i for i in range(len(lf_reference_cigar_lengths_cumsum)) if lf_reference_cigar_lengths_cumsum[i] >= left_flank_ref_length])
                                find_diff_array = lf_reference_cigar_lengths_cumsum[:(find_idx+1)]
                                find_diffs = left_flank_ref_length - find_diff_array
                                if any(find_diffs > 0):
                                    leftover_idx = max([i for i in range(len(find_diffs)) if find_diffs[i] >= 0])
                                    leftover = find_diffs[leftover_idx]
                                else:
                                    leftover_idx = -1
                                    leftover = left_flank_ref_length
                                if leftover > 0:
                                    lf_reference_cigar_lengths[leftover_idx+1] = leftover
                                for i,v in orig_ins_zip:
                                    lf_reference_cigar_lengths[i] = v
                                query_add = 0
                                lf_query_end_find = zip(lf_reference_cigar_lengths[:(leftover_idx+2)], ref_cigar_operations[:(leftover_idx+2)])
                                for l,o in lf_query_end_find:
                                    if o in ['M', 'I', 'X', '=', 'S']:
                                        query_add += l
                                lf_end = query_add
                            same_strand = False
                            if strand == '-' and read.is_reverse:
                                same_strand = True
                            elif strand == '+' and not read.is_reverse:
                                same_strand = True
                            if same_strand:
                                # Ensure that the SA is on the same strand as the primary alignment, otherwise this can cause issues getting the flanking info correct downstream
                                # Need to overwrite lf_end and rf_start if the SA is on a different chromosome from the insertion site
                                if s_chr != chrom or s_start > ins_start_ref or s_end < ins_start_ref:
                                    if first_soft_read:
                                        lf_end = first_soft_read[0]
                                    if last_soft_read:
                                        rf_start = read_length - last_soft_read[0]
                                # if there is soft-clipping at the start of both the primary and supplementary alignments that is at least 100 bases long
                                if first_soft_sa and first_soft_read:
                                    if first_soft_sa[0] < first_soft_read[0] and first_soft_read[0] >= 100:
                                        # only interested in left flank sequence that doesn't overlap the primary alignment
                                        lf_sequence = sequence[(first_soft_sa[0]+1):lf_end]
                                        sa_start_ref = s_start
                                        sa_end_ref = left_flank_calc(lf_end, ref_cigar_lengths, ref_cigar_operations, sa_start_ref)
                                        sa_entry = f"{s_chr}:{sa_start_ref}-{sa_end_ref}"
                                        if s_chr == chrom and sa_start_ref < ins_start_ref and sa_end_ref >= ins_start_ref:
                                            continue
                                        else:
                                            read_name = f"{read.query_name}-LF"
                                            # append to dictionary if key exists, which is read + lf, add new key/value pair if not
                                            if read_name_id in lf_dict:
                                                lf_dict[read_name_id].append([sa_entry, first_soft_sa[0], lf_sequence, cigar])
                                            else:
                                                lf_dict[read_name_id] = [[sa_entry, first_soft_sa[0], lf_sequence, cigar]]
                                # if there is soft-clipping at the start of just the primary alignment
                                elif first_soft_read:
                                    if first_soft_read[0] >= 100:
                                        lf_sequence = sequence[0:lf_end]
                                        sa_start_ref = s_start
                                        sa_end_ref = left_flank_calc(lf_end, ref_cigar_lengths, ref_cigar_operations, sa_start_ref)
                                        sa_entry = f"{s_chr}:{sa_start_ref}-{sa_end_ref}"
                                        if s_chr == chrom and sa_start_ref < ins_start_ref and sa_end_ref >= ins_start_ref:
                                            continue
                                        else:
                                            read_name = f"{read.query_name}-LF"
                                            if read_name_id in lf_dict:
                                                lf_dict[read_name_id].append([sa_entry, 0, lf_sequence, cigar])
                                            else:
                                                lf_dict[read_name_id] = [[sa_entry, 0, lf_sequence, cigar]]
                                # if there is soft-clipping at the end of both the primary and supplementary alignments
                                if last_soft_sa and last_soft_read:
                                    if last_soft_sa[0] < last_soft_read[0] and last_soft_read[0] >= 100:
                                        rf_sequence = sequence[rf_start:-last_soft_sa[0]]
                                        sa_end_ref = s_start + reference_length
                                        rf_ref_start = rf_from_end + 1
                                        sa_start_ref = right_flank_calc(rf_ref_start, ref_cigar_lengths, ref_cigar_operations, sa_end_ref)
                                        sa_entry = f"{s_chr}:{sa_start_ref}-{sa_end_ref}"
                                        read_name = f"{read.query_name}-RF"
                                        if read_name_id in rf_dict:
                                            rf_dict[read_name_id].append([sa_entry, last_soft_sa[0], rf_sequence, cigar])
                                        else:
                                            rf_dict[read_name_id] = [[sa_entry, last_soft_sa[0], rf_sequence, cigar]]
                                # if there is soft-clipping at the end of just the primary alignment
                                elif last_soft_read:
                                    if last_soft_read[0] >= 100:
                                        rf_sequence = sequence[rf_start:]
                                        sa_end_ref = s_start + reference_length
                                        rf_ref_start = rf_from_end + 1
                                        sa_start_ref = right_flank_calc(rf_ref_start, ref_cigar_lengths, ref_cigar_operations, sa_end_ref)
                                        sa_entry = f"{s_chr}:{sa_start_ref}-{sa_end_ref}"
                                        read_name = f"{read.query_name}-RF"
                                        if read_name_id in rf_dict:
                                            rf_dict[read_name_id].append([sa_entry, 0, rf_sequence, cigar])
                                        else:
                                            rf_dict[read_name_id] = [[sa_entry, 0, rf_sequence, cigar]]
    return [lf_dict, rf_dict, tsd_polya_dict]

def order_flanks(flank_lst, file_name):
    lf,rf = flank_lst[0:2] 
    lf_list = []
    rf_list = []
    if len(lf) > 0:
        with open(f"{file_name}.lf.fasta", 'w') as output:
            read_counter = 0
            for key,value in lf.items():
                read_counter += 1 
                insertion_id,read_id = key.split('-')[0:2]
                # if there are multiple LF alignments, need to order them and trim overlapping sequences from sequential SAs (this has already been done for the primary alignment for both flanks); for now, leave as is, can add this in later if needed 
                if len(value) > 1:
                    lf_order = sorted(value, key=itemgetter(1))
                    # only need to correct the left flank coordinates before the final left flank coordinate - that one stays as is
                    lf_ends = [sc[1]+1 for (i,sc) in enumerate(lf_order) if i > 0]
                    ref_cigar_lengths = [list(map(int, re.findall(r"(\d+)[IMSDX=]", cigar[3]))) for (i,cigar) in enumerate(lf_order) if i < (len(lf_order)-1)]
                    ref_cigar_operations = [list(re.findall(r"\d+([IMSDX=])", cigar[3])) for (i,cigar) in enumerate(lf_order) if i < (len(lf_order)-1)]
                    ref_starts = [c[0].split(':')[1].split('-')[0] for (i,c) in enumerate(lf_order) if i < (len(lf_order)-1)]
                    ref_ends = list(map(left_flank_calc, lf_ends, ref_cigar_lengths, ref_cigar_operations, ref_starts))
                    final_lf_end = [len(sc[2]) for (i,sc) in enumerate(lf_order) if i == (len(lf_order)-1)]
                    lf_ends.append(final_lf_end[0])
                    lf_seqs = [s[2] for s in lf_order]
                    final_ref_start = [c[0].split(':')[1].split('-')[0] for (i,c) in enumerate(lf_order) if i == (len(lf_order)-1)]
                    ref_starts.append(final_ref_start[0])
                    final_ref_end = [c[0].split(':')[1].split('-')[1] for (i,c) in enumerate(lf_order) if i == (len(lf_order)-1)]
                    ref_ends.append(final_ref_end[0])
                    chroms = [c[0].split(':')[0] for c in lf_order]
                    info_zip = zip(lf_ends, lf_seqs, ref_starts, ref_ends)
                    counter = 0
                    for i,flank in enumerate(info_zip):
                        sequence = flank[1][:flank[0]]
                        if len(sequence) >= 100:
                            counter += 1
                            name = f">{key}-LF{read_counter}_{counter}"
                            output.write(name + "\n" + sequence + "\n")
                            chrom = chroms[i]
                            ref_start = int(flank[2])
                            ref_end = int(flank[3])
                            lf_list.append([insertion_id, read_id, chrom, ref_start, ref_end, f"{read_counter}_{counter}"])
                else:
                    counter = 1
                    name = f">{key}-LF{read_counter}_{counter}"
                    sequence = value[0][2]
                    output.write(name + "\n" + sequence + "\n")
                    chrom = value[0][0].split(':')[0]
                    ref_start = value[0][0].split(':')[1].split('-')[0]
                    ref_end = value[0][0].split(':')[1].split('-')[1]
                    lf_list.append([insertion_id, read_id, chrom, ref_start, ref_end, f"{read_counter}_{counter}"])
    else:
        with open(f"{file_name}.lf.fasta", 'w') as output:
            output.write("Nothing to report")
    lf_df = pd.DataFrame(lf_list, columns=['InsertionID', 'ReadID', 'LF_chr', 'LF_start', 'LF_end', 'LF_aln_order'])
    if len(rf) > 0:
        with open(f"{file_name}.rf.fasta", 'w') as output:
            read_counter = 0
            for key,value in rf.items():
                read_counter += 1
                insertion_id,read_id = key.split('-')[0:2]
                # if there are multiple RF alignments, need to order them and trim overlapping sequences from sequential SAs (this has already been done for the primary alignment for both flanks); for now, leave as is, can add this in later if needed 
                if len(value) > 1: 
                    rf_order = sorted(value, key=itemgetter(1))
                    rf_starts = [(sc[1] - rf_order[i-1][1]) for (i,sc) in enumerate(rf_order) if i > 0]
                    ref_cigar_lengths = [list(map(int, re.findall(r"(\d+)[IMSDX=]", cigar[3]))) for (i,cigar) in enumerate(rf_order) if i < (len(rf_order)-1)]
                    ref_cigar_operations = [list(re.findall(r"\d+([IMSDX=])", cigar[3])) for (i,cigar) in enumerate(rf_order) if i < (len(rf_order)-1)]
                    ref_ends = [c[0].split(':')[1].split('-')[1] for (i,c) in enumerate(rf_order) if i < (len(rf_order)-1)]
                    ref_starts = list(map(right_flank_calc, rf_starts, ref_cigar_lengths, ref_cigar_operations, ref_ends))
                    final_rf_start = [len(sc[2]) for (i,sc) in enumerate(rf_order) if i == (len(rf_order)-1)]
                    rf_starts.append(final_rf_start[0])
                    rf_seqs = [s[2] for s in rf_order]
                    final_ref_start = [c[0].split(':')[1].split('-')[0] for (i,c) in enumerate(rf_order) if i == (len(rf_order)-1)]
                    ref_starts.append(final_ref_start[0])
                    final_ref_end = [c[0].split(':')[1].split('-')[1] for (i,c) in enumerate(rf_order) if i == (len(rf_order)-1)]
                    ref_ends.append(final_ref_end[0])
                    chroms = [c[0].split(':')[0] for c in rf_order]
                    info_zip = zip(rf_starts, rf_seqs, ref_starts, ref_ends)
                    counter = 0
                    for i,flank in enumerate(info_zip):
                        sequence = flank[1][-flank[0]:]
                        if len(sequence) >= 100:
                            counter += 1
                            name = f">{key}-RF{read_counter}_{counter}"
                            output.write(name + "\n" + sequence + "\n")
                            chrom = chroms[i]
                            ref_start = int(flank[2])
                            ref_end = int(flank[3])
                            rf_list.append([insertion_id, read_id, chrom, ref_start, ref_end, f"{read_counter}_{counter}"])
                else:
                    counter = 1
                    name = f">{key}-RF{read_counter}_{counter}"
                    sequence = value[0][2]
                    output.write(name + "\n" + sequence + "\n")
                    chrom = value[0][0].split(':')[0]
                    ref_start = value[0][0].split(':')[1].split('-')[0]
                    ref_end = value[0][0].split(':')[1].split('-')[1]
                    rf_list.append([insertion_id, read_id, chrom, ref_start, ref_end, f"{read_counter}_{counter}"])
    else:
        with open(f"{file_name}.rf.fasta", 'w') as output:
            output.write("Nothing to report")
    rf_df = pd.DataFrame(rf_list, columns=['InsertionID', 'ReadID', 'RF_chr', 'RF_start', 'RF_end', 'RF_aln_order'])
    return [lf_df, rf_df]

def flank_intervals(flank_dfs, flank_lst, file_name):
    tsd_polya = flank_lst[2]
    for k,v in tsd_polya.items():
        if len(v) > 1:
            ins_name = [tsd[0] for tsd in v][0]
            tsds = [tsd[1] for tsd in v]
            polya = [tsd[2] for tsd in v]
            tsd_counts = dict()
            for i in tsds:
                tsd_counts[i] = tsd_counts.get(i, 0) + 1
            most_common = max(tsd_counts.values())
            all_most_common = [k for k, v in tsd_counts.items() if v == most_common]
            longest_most_common = max(all_most_common, key=len)
            for i in tsds:
                if len(i) > len(longest_most_common):
                    longest_most_common = i
            if True in polya:
                tsd_polya[k] = [[ins_name, longest_most_common, True]]
            else:
                tsd_polya[k] = [[ins_name, longest_most_common, False]]
    tsd_polya_lst = []
    for entry in tsd_polya.values():
        tsd_polya_lst.append(entry[0])
    tsd_polya_df = pd.DataFrame(tsd_polya_lst, columns=['InsertionID', 'TSD', 'PolyA_tail'])
    lf_df,rf_df = flank_dfs[0:2]
    if len(lf_df) > 0 and len(rf_df) > 0:
        lf_max_intervals = find_intervals(lf_df, "LF")
        rf_max_intervals = find_intervals(rf_df, "RF")
        merged_df = lf_max_intervals.merge(rf_max_intervals, how='outer', on='InsertionID')
        merged_df = merged_df.dropna()
        merged_df['ReadID'] = merged_df[['ReadID_x', 'ReadID_y']].agg(';'.join, axis = 1)
        unique_merged_read = merged_df['ReadID'].values
        unique_merged_read_lst = []
        for s in unique_merged_read:
            read_split = re.split("[;,]+", s)
            read_split_dict = list(dict.fromkeys(read_split))
            unique_merged_read_lst.append(';'.join(read_split_dict))
        merged_df['ReadID'] = unique_merged_read_lst
        merged_df.drop(['ReadID_x', 'ReadID_y'], axis=1, inplace=True)
        merged_df = merged_df.reset_index()
        merged_df = merged_df.merge(tsd_polya_df, how="inner", on="InsertionID")
        merged_df = merged_df.dropna()
        merged_df.to_csv(f"{file_name}.flank.report.txt", index=False, sep = '\t', header=True)
    else:
        with open(f"{file_name}.flank.report.txt", 'w') as outf:
            outf.write("Nothing to report")
def main():
    inputs = parse_args()
    fasta =  inputs.gencDNA_fasta
    ins_locs = inputs.ins_pos_file
    bam = inputs.bam_file
    pysam.index(bam)
    file_name = bam.split('/')[-1].split('.')[0]
    get_flanks = extract_flanks(fasta, bam, ins_locs, file_name)
    flank_dfs = order_flanks(get_flanks, file_name)
    flank_intervals(flank_dfs, get_flanks, file_name)

if __name__ =="__main__":
    main()