#!/usr/bin/env python3 
# 4-space indented, v0.0.1
# File name: pseudogene_readnames.py
# Description: Filters out known pseudogenes.
# Author: Robert Linder
# Date: 2024-01-10

import pandas as pd
import polars as pl
import re
import os
import argparse
from subprocess import Popen, PIPE, check_output
from itertools import combinations
from Bio import SeqIO

def parse_args():
    """this enables users to provide input as command line arguments to minimize the need for directly modifying the script; please enter '-h' to see the full list of options"""
    parser = argparse.ArgumentParser(description="Load file for filtering known pseudogenes")
    parser.add_argument("gencDNA_file", type=str, help="potential gencDNA bed file to filter")
    args = parser.parse_args()
    return args

def extract_read_names(genc, file_name):
    """This extracts the read names from reads at which both the original alignment and the mapped insertion overlap the same annotated pseudogene and outputs a bed file with same pseudogene-containing entries removed"""
    bed = (
         pl.read_csv(genc, has_header=False, separator='\t', dtypes={"column_11": pl.Utf8, "column_12": pl.Utf8})
        .with_columns(pl.when(pl.col("column_4").str.contains("Sniffles2"))
                    .then(pl.col("column_4").str.extract(r"(Sniffles2.INS.\w+.\d+)"))
                    .otherwise(pl.col("column_4"))
                    .alias("Insertion_ID")
                    )
        .with_columns(pl.when((pl.col("column_4").str.contains("ccs-Sniffles2")) & (pl.col("column_14") != "RepeatMasker"))
                    .then(0)
                    .when(pl.col("column_4").str.contains("NativeHit"))
                    .then(2)
                    .when(pl.col("column_14") == "RepeatMasker")
                    .then(3)
                    .otherwise(1)
                    .alias("Read_type")
                    )
        .with_columns(pl.when(pl.col("column_21").str.contains("pseudo") )
                    .then(1)
                    .when(pl.col("column_25").str.contains("pseudo") )
                    .then(1)
                    .otherwise(0)
                    .alias("Pseudogene_status")
                    )
        .rename({"column_24": "Gene_name"})
        .with_columns(pl.when(pl.col("Read_type") == 2)
                    .then(pl.col("column_4").str.extract(r"(.*)_NativeHit"))
                    .when(pl.col("Read_type") == 3)
                    .then(pl.lit("Repeat"))
                    .otherwise(pl.lit("Insertion"))
                    .alias("Read_origin")
                    )
        .with_columns(pl.when(pl.col("Read_type") == 2)
                    .then(pl.concat_str([pl.col("Read_origin"), pl.col("column_30"), pl.col("column_31"), pl.col("column_32"), pl.lit("NativeHit")], separator="_"))
                    .otherwise(pl.col("column_4"))
                    .alias("Full_read_ID")
                    )
        .with_columns(pl.when(pl.col("column_4").str.contains("ccs-Sniffles2"))
                    .then(pl.col("column_4").str.extract(r".*Sniffles2.INS.\w+.\d+-(\w+\d+)"))
                    .otherwise(pl.col("column_1"))
                    .alias("True_Flank_Chr")
                    )
    )
    bed2 = bed.filter(pl.col("Read_origin") != "Repeat")
    # Need to exclude flanking rows in which the read was mapped to the insertion (likely a supplementary alignment); for now, this is implemented by excluding flanking rows in which the chromosome mapped to is the same as for the insertion and a gene mapped to the insertion site is mapped to the flanking region as well (so same chromosome and same gene)
    ins_chr_df = bed2.filter(pl.col("Read_type") == 1).group_by("Insertion_ID").agg(pl.col("column_1"), pl.col("column_26"), pl.col("Gene_name")).with_columns(pl.col("column_1").list.drop_nulls().list.unique()).with_columns(pl.col("column_26").list.drop_nulls().list.unique()).with_columns(pl.col("Gene_name").list.drop_nulls().list.unique()).with_columns(all_genes1=pl.col("Gene_name").list.concat("column_26")).drop(['Gene_name', 'column_26'])
    bed2 = bed2.join(ins_chr_df, on = "Insertion_ID", how="inner")
    bed2 = bed2.with_columns(pl.col('column_21').str.to_uppercase().alias("Info")).with_row_index()
    rows_to_keep1 = []
    for i, row in enumerate(bed2.iter_rows(named=True)):
        if any(val.upper() in row['Info'] for val in row['all_genes1']) and row['Read_type'] == 0:
            rows_to_keep1.append(i)
    bed2 = bed2.filter(~((pl.col("Read_type") == 0) & (pl.col("column_1").is_in(pl.col("column_1_right"))) & (pl.col("index").is_in(rows_to_keep1))))
    bed2 = bed2.drop(['index', 'all_genes1', 'Info', "column_1_right"])
    #bed.write_csv("J.ps_text.bed", separator='\t', include_header=True)
    native_bed = bed.filter(pl.col("Read_type") == 2)
    find_ins_pseudo = bed2.filter(pl.col("Pseudogene_status") == 1).get_column("Insertion_ID").to_list()
    unique_ins_pseudo = list(set(find_ins_pseudo))
    # the arg_true function returns the index of the values in a list column that meet the specified condition
    bed2 = bed2.filter(pl.col("Insertion_ID").is_in(unique_ins_pseudo)).group_by("Insertion_ID").agg(pl.col("Gene_name"), pl.col("Read_type"), pl.col("column_21"), pl.col("Full_read_ID"), pl.col("column_26"), pl.col("True_Flank_Chr"), pl.col("column_1"), pl.col("Pseudogene_status")).filter(pl.col("Read_type").list.contains(1)).filter(pl.col("Read_type").list.contains(0)).with_columns(pl.col("Read_type")
                    .list.eval((pl.element() == 0).arg_true())
                    .alias("Flanking_index")
                ).with_columns(pl.col("Read_type")
                    .list.eval((pl.element() == 1).arg_true())
                    .alias("Insertion_index")
    )  
    if not bed2.is_empty():
        bed3 = bed2.with_columns(
            pl.col("Read_type")
            .list.sum()
            .alias("Read_types_sum")
        ).with_columns(
            pl.col("Read_type")
            .list.len()
            .alias("Read_type_count")
        ).with_columns(
            (pl.col("Read_types_sum")/pl.col("Read_type_count"))
            .alias("Read_type_fraction")
        ).filter((pl.col("Read_type_fraction") > 0) & (pl.col("Read_type_fraction") < 1)
        ).with_columns(pl.col("Gene_name")
                    .list.drop_nulls()
                    .list.unique()
                    .alias("Unique_genes")
        ).with_columns(pl.col("column_26")
                    .list.drop_nulls()
                    .list.unique()
                    .alias("Unique_common_genes")
        )
        bed3 = bed3.with_columns(all_genes=pl.col("Unique_genes").list.concat("Unique_common_genes"))
        # Index column 21 values by the values in the read type index column
        # exclusion_df = bed3.with_columns(pl.col('column_21').list.gather(pl.col("Flanking_index"))).explode("all_genes").with_columns(pl.col('column_21').list.join(",")).filter(pl.col("column_21").str.contains(r"(?i)"+pl.col("all_genes")))  
        exclusion_df = bed3.with_columns(pl.col('column_21').list.gather(pl.col("Flanking_index"))).with_columns(pl.col('column_21').list.join(",")).with_columns(pl.col('column_21').str.to_uppercase()).with_row_index()
        rows_to_keep = []
        for i, row in enumerate(exclusion_df.iter_rows(named=True)):
            if any(val.upper() in row['column_21'] for val in row['all_genes']):
                rows_to_keep.append(row['index'])
        exclusion_df = exclusion_df.filter(pl.col("index").is_in(rows_to_keep))
        # print(exclusion_df.filter(pl.col("Insertion_ID") == "Sniffles2.INS.ADCS5.1").get_column("column_21").to_list())
        # Don't exclude insertions in which a pseudogene was detected in both the insertion and flanking region if they both map to the same chromosome, but the actual flank was detected on another chromsome from the insertion as per Sniffles2, as this implies secondary alignment of the flanking region by minimap2, which can lead to erroneous false negatives; also added in a bit to ensure instances in which the pseudogene is detected in the true flanking region as well as the insertion are discarded
        #print(exclusion_df.with_columns(pl.col('column_1').list.gather(pl.col("Insertion_index")).alias("Insertion_chr")))
        if len(exclusion_df) > 0:
            inclusion_lst = exclusion_df.with_columns(pl.col('column_1').list.gather(pl.col("Insertion_index")).alias("Insertion_chr")).with_columns(pl.col('column_1').list.gather(pl.col("Flanking_index")).alias("Flanking_chr")).with_columns(pl.col('True_Flank_Chr').list.gather(pl.col("Flanking_index")).alias("True_Flanking_chr")).with_columns(pl.col('Pseudogene_status').list.gather(pl.col("Flanking_index")).alias("Flanking_ps_status")).with_columns(pl.col("Insertion_chr").list.drop_nulls().list.unique()).with_columns(pl.col("Flanking_chr").list.drop_nulls().list.unique()).with_columns(pl.col("True_Flanking_chr").list.drop_nulls().list.unique()).with_columns(pl.col("Flanking_ps_status").list.drop_nulls().list.unique()).explode("Insertion_chr").with_columns(pl.col('Flanking_chr').list.join(",")).filter(pl.col("Flanking_chr").str.contains(pl.col("Insertion_chr"))).explode("True_Flanking_chr").filter(~pl.col("Insertion_chr").str.contains(pl.col("True_Flanking_chr"))).filter(~pl.col("Flanking_ps_status").list.contains(1)).get_column("Insertion_ID").to_list()
            unique_inclusion_lst = list(set(inclusion_lst))
            exclusion_df = exclusion_df.filter(~pl.col("Insertion_ID").is_in(unique_inclusion_lst))
            all_exclusion_lst = exclusion_df.get_column("Insertion_ID").to_list()
            unique_exclusion_lst = list(set(all_exclusion_lst)) # list of insertion IDs to exclude
            # rejoin the repeat dataframe with the 
            insertion_bed_df = bed.filter(~pl.col("Insertion_ID").is_in(unique_exclusion_lst), ~pl.col("Full_read_ID").str.contains("-")).drop(["Insertion_ID", "Read_type", "Pseudogene_status", "Read_origin", "Full_read_ID"])
            insertion_bed_df.write_csv(f"{file_name}.pseudogeneless.bed", separator = '\t', include_header=False)
            excluded_read_lst = exclusion_df.get_column("Full_read_ID").to_list()
            excluded_read_lst = [x for xs in excluded_read_lst for x in xs]
            excluded_read_lst = [r for r in excluded_read_lst if '-' in r]
            unique_excluded_read_lst = list(set(excluded_read_lst))
            unique_excluded_read_df = pl.DataFrame({"reads": unique_excluded_read_lst})
            unique_excluded_read_df.write_csv(f"{file_name}.pseudogene.reads.txt", separator = '\t', include_header=False)
            included_read_lst = bed.filter(~pl.col("Insertion_ID").is_in(unique_exclusion_lst)).filter(~pl.col("Full_read_ID").str.contains('-')).get_column("Full_read_ID").to_list()
            unique_included_read_lst = list(set(included_read_lst))
            unique_included_read_df = pl.DataFrame({"reads": unique_included_read_lst})
            unique_included_read_df.write_csv(f"{file_name}.pseudogeneless.reads.txt", separator = '\t', include_header=False)
        else:
            bed_df = bed.drop(["Insertion_ID", "Read_type", "Pseudogene_status", "Read_origin", "Full_read_ID"])
            bed_df.write_csv(f"{file_name}.pseudogeneless.bed", separator = '\t', include_header=False)
            included_read_lst = bed.filter(~pl.col("Full_read_ID").str.contains('-')).get_column("Full_read_ID").to_list()
            unique_included_read_lst = list(set(included_read_lst))
            unique_included_read_df = pl.DataFrame({"reads": unique_included_read_lst})
            unique_included_read_df.write_csv(f"{file_name}.pseudogeneless.reads.txt", separator = '\t', include_header=False)
            unique_excluded_read_df = pl.DataFrame()
            unique_excluded_read_df.write_csv(f"{file_name}.pseudogene.reads.txt", separator = '\t', include_header=False)
    # For when the native locus hits are being filtered
    elif not native_bed.is_empty():
        exclusion_df = native_bed.filter(pl.col("Pseudogene_status") == 1)
        exclusion_lst = exclusion_df.get_column("Full_read_ID").to_list()
        unique_exclusion_lst = list(set(exclusion_lst))
        native_bed_df = native_bed.filter(~pl.col("Full_read_ID").is_in(unique_exclusion_lst)).drop(["Insertion_ID", "Read_type", "Pseudogene_status", "Read_origin", "Full_read_ID"])
        native_bed_df.write_csv(f"{file_name}.pseudogeneless.bed", separator = '\t', include_header=False)
        excluded_read_lst = exclusion_df.get_column("Full_read_ID").to_list()
        unique_excluded_read_lst = list(set(excluded_read_lst))
        unique_excluded_read_df = pl.DataFrame({"reads": unique_excluded_read_lst})
        unique_excluded_read_df.write_csv(f"{file_name}.pseudogene.reads.txt", separator = '\t', include_header=False)
        included_read_lst = native_bed.filter(~pl.col("Full_read_ID").is_in(unique_exclusion_lst)).get_column("Full_read_ID").to_list()
        unique_included_read_lst = list(set(included_read_lst))
        unique_included_read_df = pl.DataFrame({"reads": unique_included_read_lst})
        unique_included_read_df.write_csv(f"{file_name}.pseudogeneless.reads.txt", separator = '\t', include_header=False)
    # if there are no psuedogenes to filter out
    else:
        bed_df = bed.drop(["Insertion_ID", "Read_type", "Pseudogene_status", "Read_origin", "Full_read_ID"])
        bed_df.write_csv(f"{file_name}.pseudogeneless.bed", separator = '\t', include_header=False)
        included_read_lst = bed.filter(~pl.col("Full_read_ID").str.contains('-')).get_column("Full_read_ID").to_list()
        unique_included_read_lst = list(set(included_read_lst))
        unique_included_read_df = pl.DataFrame({"reads": unique_included_read_lst})
        unique_included_read_df.write_csv(f"{file_name}.pseudogeneless.reads.txt", separator = '\t', include_header=False)
        unique_excluded_read_df = pl.DataFrame()
        unique_excluded_read_df.write_csv(f"{file_name}.pseudogene.reads.txt", separator = '\t', include_header=False)
        
def main():
    inputs = parse_args()
    genc =  inputs.gencDNA_file
    file_name = genc.split('/')[-1].split('.')[:-1]
    file_name = '.'.join(file_name)
    extract_read_names(genc, file_name)

if __name__ =="__main__":
    main()