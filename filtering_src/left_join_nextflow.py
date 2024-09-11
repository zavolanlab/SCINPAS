# -*- coding: utf-8 -*-
"""

Created on Sun Nov 26 22:55:32 2023

@author: Youngbin Moon (y.moon@unibas.ch)
"""











import argparse
import itertools
from pathlib import Path
from typing import Tuple
import numpy as np
import pandas as pd
import os

def write_result(df, out):
    df.to_csv(out, sep = '\t', header = True, index = False) 

def filter_rows(df, filter_dictionary):
    
    filtered_final_df = df.query(create_filter_query_for_df(filter_dictionary)).copy()
    # print('filtered_df: ' + str(filtered_final_df))
    return filtered_final_df

def create_filter_query_for_df(filter_dict):
    # What if my column names have whitespace, or other weird characters?
    # From pandas 0.25, you can wrap your column name in backticks so this works:
    list_pairs = [f'`{k}`>{v}' for k, v in filter_dict.items()]
    query = '|'.join(list_pairs)
    '(col20 > 0) | (col25 > 0)......'
    # print(query)
    return query
    
def left_join(modified_cs_dir, beds_list):
    
    sample_columns = {}
    modified_cs = pd.read_csv(modified_cs_dir, delimiter = '\t', header = 0)
    # need to drop supp column, because this is support across all samples. You want support for each sample, and then aggregate for same organ.
    # all_samples_modified_unique_cs.bed has 'seqid', 'start', 'end', 'id', 'score', 'strand', 'supp'
    modified_cs = modified_cs[['seqid', 'start', 'end', 'id', 'score', 'strand']]
    original_modified_cs_length = len(modified_cs)
    for file in beds_list:
        df = pd.read_csv(file, delimiter = '\t', header = 0)
        sample = df.columns[-1]
        
        # right df has Chrom, start, end, id, score(rpm), strand, supp, sample(=# reads)
        # need to subset right df to avoid duplicate columns of 'id', 'score'. also you dont need sample (# of reads) (8th column)
        # you dont need supp (7th column) either
        # What you need is 'score' (5th column) of each sample because that indicates sample level expression of a particular cleavage site
        # You will use this to aggregate over organs, and then get organ-specific score (rpm)
        subset_df = df.iloc[:, [0, 1, 2, 5, 4]]
        print('subset_df:' + str(subset_df))
        # you need to change sample score column of individual sample into 'sample' column because you want to merge scores over same organs.
        subset_df.rename(columns = {'score': sample}, inplace = True)
        # here you won't have duplicate rows because you did exact match of seqid, start, end and strand. 
        # not like when you did merge between PAS and unique_cs in terms of 'id' 
        modified_cs = pd.merge(left = modified_cs, right = subset_df, how = 'left', on = ['seqid', 'start', 'end', 'strand'])
       
        sample_columns[sample] = 0
        
    final_df = modified_cs.fillna(0, inplace = False)
    print('length of original df: ' + str(original_modified_cs_length))
    print('length of added df (before filtering): ' + str(len(final_df)))
    assert(original_modified_cs_length == len(final_df))
    
    return final_df, sample_columns

def get_args():        
    parser = argparse.ArgumentParser(description="left join modified_unique_cs and unique_cs beds of a particular organ to generate 1. modified_unique_cs.bed of the organ and 2. PAS of the organ ")

    parser.add_argument('in_bed',
      nargs='*', help='input original unique cleavage site bed files of a particular organ.')

    parser.add_argument('--modified_unique_cs', dest = 'modified_unique_cs',
                        required = True,
                        help = 'a modified unique cleavage site bed file')
    
    parser.add_argument('--out_cs', dest = 'out_cs',
                        required = True,
                        help = 'modified_unique_cs.bed of a particular organ')

    parser.add_argument('--out_pas', dest = 'out_pas',
                        required = True,
                        help = 'PAS of a particular organ')
    
    parser.add_argument('--chrom', dest = 'chrom',
                        required = True,
                        help = 'chromosome')

    parser.add_argument('--direction', dest = 'direction',
                        required = True,
                        help = 'direction of cs')

    parser.add_argument('--organ', dest = 'organ',
                        required = True,
                        help = 'target organ')

    parser.add_argument('--n_cores', dest = 'n_cores',
                        required = True,
                        help = 'the number of cores')
    
    args = parser.parse_args()
    
    bed_in = args.in_bed
    modified_unique_cs = args.modified_unique_cs
    out_cs_template = args.out_cs
    out_pas_template = args.out_pas
    
    chrom = args.chrom
    direction = args.direction
    organ = args.organ
    n_cores = int(args.n_cores)
    
    return bed_in, modified_unique_cs, out_cs_template, out_pas_template, chrom, direction, organ, n_cores

def run_process():

    bed_in, modified_unique_cs, out_cs_template, out_pas_template, chrom, direction, organ, n_cores = get_args()
    
    print('successfully got arguments')
    
    first_elem = bed_in.pop(0)
    assert(first_elem == 'in_bed')
    
    print('bed_in: ' + str(bed_in))
    
    final_df, sample_columns = left_join(modified_unique_cs, bed_in)    
    print('successfully done left join')
        
    filtered_final_df = filter_rows(final_df, sample_columns)
    print('successfully got filtered modified unique cs with sample columns')
    
    out_cs_filtered = out_cs_template + '_filtered_' + chrom + '_' + direction + '.bed'
    write_result(filtered_final_df, out_cs_filtered)
    print('successfully got 2nd output')
        
if __name__ == "__main__":
    run_process()
    print("success")