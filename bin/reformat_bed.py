#!/usr/bin/env python3 
# 4-space indented, v0.0.1
# File name: reformat_bed.py
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
    parser.add_argument("gff_file", type=str, help="The annotation file")
    parser.add_argument("gencDNA_file", type=str, help="potential gencDNA bed file to filter")
    args = parser.parse_args()
    return args

def reformat_bed(file, gff, file_name):
    balb_transcript_lst = ["C_gene_segment","D_gene_segment","J_gene_segment","V_gene_segment","gene_segment","lnc_RNA","mRNA","miRNA","ncRNA","pseudogenic_transcript","rRNA","scRNA","snRNA","snoRNA","tRNA","unconfirmed_transcript"]
    genc_df = pd.read_csv(file, sep='\t', header=None)
    if 'gene_biotype' in ''.join(genc_df[20].values) and 'gene_name' in ''.join(genc_df[20].values):
        # This is for the mm39 gff
        gene_name = genc_df[20].str.split(';')
        r = re.compile(".*gene_id.*")
        gt = re.compile(".*gene_biotype.*")
        gcn = re.compile(".*gene_name.*")
        transc = re.compile(".*transcript_id.*")
        gene_name_list = gene_name.values.tolist()
        gene_list = []
        common_name_list = []
        gene_type_list = []
        transc_list = []
        for entry in gene_name_list:
            entry = [gene.replace("'", "").replace('"', '').strip() for gene in entry]
            matches = [gene.split(' ')[1] for gene in entry if r.match(gene)]
            gt_matches = [gene.split(' ')[1] for gene in entry if gt.match(gene)]
            transc_matches = [gene.split(' ')[1] for gene in entry if transc.match(gene)]
            common_matches = [gene.split(' ')[1] for gene in entry if gcn.match(gene)]
            if len(common_matches) == 0:
                common_name_list.append(matches)
            else:
                common_name_list.append(common_matches)
            gene_list.append(matches)
            gene_type_list.append(gt_matches)
            transc_list.append(transc_matches)
        gene_list_flattened = [item for row in gene_list for item in row]
        gene_type_list_flattened = [item for row in gene_type_list for item in row]
        transc_list_flattened = [item for row in transc_list for item in row]
        common_list_flattened = [item for row in common_name_list for item in row]
        genc_df['transcript'] = transc_list_flattened
        genc_df['gene'] = gene_list_flattened
        genc_df['gene_type'] = gene_type_list_flattened
        # The common name is listed in the same gene= part of the info field as the LOC gene name if there is no common name
        genc_df['gene_name'] = common_list_flattened
    elif 'gene_name' in ''.join(genc_df[20].values):
        # This is for the hg38 ensembl gff
        gene_name = genc_df[20].str.split(';')
        genc_df['transcript'] = [gene[7] for gene in gene_name]
        genc_df['gene'] = [gene[5] for gene in gene_name]
        genc_df['gene_type'] = [gene[4] for gene in gene_name]
        # Need to add gene_name column for common name, so have same number of columns for all output files  
    elif 'gene=' in ''.join(genc_df[20].values):
         # This is for the T2T gff
        gene_name = genc_df[20].str.split(';')
        r = re.compile("gene=.*")
        gt = re.compile("product=.*")
        transc = re.compile("transcript_id=.*")
        gene_name_list = gene_name.values.tolist()
        gene_list = []
        gene_type_list = []
        transc_list = []
        for entry in gene_name_list:
            matches = [gene.split('=')[1] for gene in entry if r.match(gene)]
            gt_matches = [gene.split('=')[1] for gene in entry if gt.match(gene)]
            transc_matches = [gene.split('=')[1] for gene in entry if transc.match(gene)]
            gene_list.append(matches)
            gene_type_list.append(gt_matches)
            transc_list.append(transc_matches)
        gene_list_flattened = [item for row in gene_list for item in row]
        gene_type_list_flattened = [item for row in gene_type_list for item in row]
        transc_list_flattened = [item for row in transc_list for item in row]
        genc_df['transcript'] = transc_list_flattened
        genc_df['gene'] = gene_list_flattened
        genc_df['gene_type'] = gene_type_list_flattened
        # The common name is listed in the same gene= part of the info field as the LOC gene name if there is no common name
        genc_df['gene_name'] = gene_list_flattened
    elif 'Parent' in ''.join(genc_df[20].values):
        # This is for the balb v3 gff
        gff_file = pd.read_csv(gff, sep='\t', header=None, comment='#')
        gff_file = gff_file[gff_file[1] != 'data']
        # Now need to split the dataframe such that any gffs that got added to the primary gff can be accounted for; split by the second column (contains data for additional gffs added on), then concat them back together at the end; can leave out the repetitive element annotations
        primary_gff_genc_df = genc_df[~genc_df[13].isin(['data', 'RepeatMasker']) & ~genc_df[20].str.contains("ID=gene:")]
        secondary_gff_genc_df = genc_df[genc_df[13] == 'data']
        tertiary_gff_genc_df = genc_df[genc_df[20].str.contains("ID=gene:")]
        # The tertiary is for pseudogenes in the balb annotation that are genes (not exons)
        if len(tertiary_gff_genc_df) > 0:
            ter_gene_name = tertiary_gff_genc_df[20].str.split(';')
            ter_gene_name_list = ter_gene_name.values.tolist()
            r = re.compile("gene_id=.*")
            gt = re.compile("biotype=.*")
            gcn = re.compile(".*parent_gene_display_xref.*")
            transc = re.compile("ID=gene.*")
            add_gene_list = []
            add_gene_type_list = []
            add_common_name_list = []
            transc_list = []
            for entry in ter_gene_name_list:
                matches = [gene.split('=')[1] for gene in entry if r.match(gene)]
                gt_matches = [gene.split('=')[1] for gene in entry if gt.match(gene)]
                # Giving common names default of gene ids as strings aren't consistent for parsing
                common_matches = [gene.split('=')[1] for gene in entry if r.match(gene)]
                transc_matches = [gene.split('=')[1] for gene in entry if transc.match(gene)]
                add_gene_list.append(matches)
                add_gene_type_list.append(gt_matches)
                add_common_name_list.append(common_matches)
                transc_list.append(transc_matches)
            add_gene_list_flattened = [item for row in add_gene_list for item in row]
            add_gene_type_list_flattened = [item for row in add_gene_type_list for item in row]
            add_common_name_list_flattened = [item for row in add_common_name_list for item in row]
            ter_transc_list_flattened = [item for row in transc_list for item in row]
            tertiary_gff_genc_df['transcript'] = ter_transc_list_flattened
            tertiary_gff_genc_df['gene'] = add_gene_list_flattened
            tertiary_gff_genc_df['gene_type'] = add_gene_type_list_flattened
            tertiary_gff_genc_df['gene_name'] = add_common_name_list_flattened
        
        if len(secondary_gff_genc_df) > 0:
            sec_gene_name = secondary_gff_genc_df[20].str.split(';')
            sec_gene_name_list = sec_gene_name.values.tolist()
            r = re.compile("gene_id=.*")
            gt = re.compile("gene_type=.*")
            gcn = re.compile("common_name=.*")
            transc = re.compile(r'Parent=transcript.*|thickEnd=gene')
            add_gene_list = []
            add_gene_type_list = []
            add_common_name_list = []
            transc_list = []
            for entry in sec_gene_name_list:
                matches = [gene.split('=')[1] for gene in entry if r.match(gene)]
                gt_matches = [gene.split('=')[1] for gene in entry if gt.match(gene)]
                common_matches = [gene.split('=')[1] for gene in entry if gcn.match(gene)]
                transc_matches = [gene.split('=')[1] for gene in entry if transc.match(gene)]
                add_gene_list.append(matches)
                add_gene_type_list.append(gt_matches)
                add_common_name_list.append(common_matches)
                transc_list.append(transc_matches)
            add_gene_list_flattened = [item for row in add_gene_list for item in row]
            add_gene_type_list_flattened = [item for row in add_gene_type_list for item in row]
            add_common_name_list_flattened = [item for row in add_common_name_list for item in row]
            sec_transc_list_flattened = [item for row in transc_list for item in row]
            secondary_gff_genc_df['transcript'] = sec_transc_list_flattened
            secondary_gff_genc_df['gene'] = add_gene_list_flattened
            secondary_gff_genc_df['gene_type'] = add_gene_type_list_flattened
            secondary_gff_genc_df['gene_name'] = add_common_name_list_flattened
        
        gene_name = primary_gff_genc_df[20].str.split(';')
        primary_gff_genc_df['transcript'] = [gene[0] for gene in gene_name]
        transcripts = primary_gff_genc_df['transcript'].str.split(':')
        transcripts = [tr[1] for tr in transcripts]
        parent_transcript_list = gff_file.loc[gff_file[2].isin(balb_transcript_lst), 8].values.tolist()
        parent_genes_df = gff_file[gff_file[2] == "gene"]
        gene_list = []
        common_name_list = []
        gene_type_list = []
        for mRNA in parent_transcript_list:
            transcript = mRNA.split(';')[0].split(':')[1]
            if transcript in transcripts:
                #print(transcript)
                #print(mRNA)
                info_list = mRNA.split(';')
                gene_id_grep = re.compile(".*=gene:.*")
                gene_id = [g.split('=')[1].split(':')[1] for g in info_list if gene_id_grep.match(g)]
                gene_id = ';'.join(gene_id)
                gene_type = mRNA.split(';')[2].split('=')[1]
                find_common_name = parent_genes_df.loc[parent_genes_df[8].str.contains(gene_id), 8].values.tolist()
                if find_common_name:
                    split_common_names = find_common_name[0].split(';')
                    common_name = re.compile("Name=.*")
                    matches = [g.split('=')[1] for g in split_common_names if common_name.match(g)]
                    matches = ';'.join(matches)
                    if matches:
                        common_name_list.append([transcript, matches])
                    else:
                        # editing so if there is no common name, get at least the ensembl gene id 
                        common_name_list.append([transcript, gene_id])
                else:
                    common_name_list.append([transcript, gene_id])
                gene_list.append([transcript, gene_id])
                gene_type_list.append([transcript, gene_type])
        reordered_common_gene_list = []
        reordered_gene_list = []
        reordered_gene_type_list = []
        for trans in transcripts:
            new_gene_list = [tr[1] for tr in gene_list if tr[0] == trans]
            new_gene_type_list = [tr[1] for tr in gene_type_list if tr[0] == trans]
            new_gene_name_list = [tr[1] for tr in common_name_list if tr[0] == trans]
            reordered_gene_list.append(new_gene_list)
            reordered_gene_type_list.append(new_gene_type_list)
            reordered_common_gene_list.append(new_gene_name_list)
        gene_list_flattened = [item for row in reordered_gene_list for item in row]
        gene_type_list_flattened = [item for row in reordered_gene_type_list for item in row]
        common_gene_name_list_flattened = [item for row in reordered_common_gene_list for item in row]
        primary_gff_genc_df['gene'] = gene_list_flattened
        primary_gff_genc_df['gene_type'] = gene_type_list_flattened
        primary_gff_genc_df['gene_name'] = common_gene_name_list_flattened
        # now concatenate the primary and secondary dataframes 
        if len(secondary_gff_genc_df) > 0:
            genc_df = pd.concat([primary_gff_genc_df, secondary_gff_genc_df], ignore_index=True)
        else:
            genc_df = primary_gff_genc_df
        if len(tertiary_gff_genc_df) > 0:
            genc_df = pd.concat([genc_df, tertiary_gff_genc_df], ignore_index=True)
        else:
            genc_df = genc_df
    # now fill in the gene, gene_type, and gene_name columns for the additional gff entries 
    genc_df['chr'] = genc_df[12]
    genc_df['overlap_rows_1'] = 0
    genc_df['nonoverlap_rows_1'] = 0
    alignment = genc_df[3].str.split('_')
    # This is to process insertions separately from potential hits at the native locus
    if re.match("Sniffles2", genc_df[3].values[0]):
        genc_df['alignment'] = [aligned[1] if len(aligned) > 1 else 0 for aligned in alignment]
        genc_df['tp'] = [aligned[2] if len(aligned) > 1 else 0 for aligned in alignment]
        genc_df['AS'] = [aligned[3] if len(aligned) > 1 else 0 for aligned in alignment]
        genc_df[3] = [aligned[0] for aligned in alignment]
    else:
        genc_df['alignment'] = [aligned[4] if len(aligned) > 1 else 0 for aligned in alignment]
        genc_df['tp'] = [aligned[5] if len(aligned) > 1 else 0 for aligned in alignment]
        genc_df['AS'] = [aligned[6] if len(aligned) > 1 else 0 for aligned in alignment]
        genc_df[3] = ['_'.join(aligned[0:4] + [aligned[7]]) if len(aligned) > 1 else 0 for aligned in alignment]
    genc_df = genc_df.sort_values(by=[3, 'chr', 'gene'])
    genc_df = genc_df.reset_index(drop=True)
    genc_df.to_csv(f"{file_name}.reformatted.bed", index=False, sep = '\t', header=False)

def main():
    inputs = parse_args()
    file =  inputs.gencDNA_file
    gff = inputs.gff_file
    file_name = file.split('/')[-1].split('.')[:-1]
    file_name = '.'.join(file_name)
    reformat_bed(file, gff, file_name)

if __name__ =="__main__":
    main()
