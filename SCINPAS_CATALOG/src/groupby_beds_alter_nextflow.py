# -*- coding: utf-8 -*-
"""

Created on Sun Dec 19 22:55:32 2021

@author: Youngbin Moon (y.moon@unibas.ch)
"""











import argparse
from collections import Counter
import pandas as pd
import numpy as np
import os
import multiprocessing as mp

def write_result(df, out):
    df.to_csv(out, sep = '\t', header = True, index = False) 

def merge_all_samples_bed(beds_list, chromosome):

    # os.chdir(r"C:/Users/Geniu/Downloads/new_folder")
        
    # rule for grouping
    # exclude id here because you will group by id
    aggregation_functions = {'seqid': 'first', 'start': 'first',\
                             'end': 'first', 'score': 'sum', 'strand': 'first',\
                             'supp': 'sum'}   
    
    print('beds list before pop: ' + str(beds_list))
    first_bed_dir = beds_list.pop(0)
    
    first_df = pd.read_csv(first_bed_dir, delimiter = '\t', header = 0)
    total_df = first_df.iloc[:, 0:7].copy()
    
    n = 1
    for file in beds_list:
        df = pd.read_csv(file, delimiter = '\t', header = 0)
        
        df_subset = df.iloc[:, 0:7]
        # concat bed files and convert NA values to 0
        total_df = pd.concat([total_df, df_subset]).copy()
        # Group by cleavage site id and then aggregate by the rules we defined.
        # When you aggregate within the same group, score, supp and samples are summed
        # whilst, other columns, we use first one. (because columns like chr, start, end, direction are identical within the group)
        # Since cleavage site id = chr20:6112552:- for example. same id means same chr, same start, same end and same direction.
        
        total_df = total_df.groupby(by = 'id', as_index = False).agg(aggregation_functions).reindex(columns = total_df.columns).copy()
        n += 1
        print('successfully merged ' + str(n) + 'th bed file of chromosome' + str(chromosome))
    
    print(total_df)
    return total_df
  
def get_args():        
    parser = argparse.ArgumentParser(description="concat all unique cleavage sites bed files from all samples")

    parser.add_argument('in_bed',
      nargs='*', help='input bed files.')
    
    parser.add_argument('--bed_out', dest = 'bed_out',
                        required = True,
                        help = 'merged bed output')

    parser.add_argument('--chrom', dest = 'chrom',
                        required = True,
                        help = 'chromosome')

    parser.add_argument('--n_cores', dest = 'n_cores',
                        required = True,
                        help = 'the number of cores')

    parser.add_argument('--num_samples', dest = 'num_samples',
                        required = True,
                        help = 'the number of samples')

    parser.add_argument('--direction', dest = 'direction',
                        required = True,
                        help = 'direction of cs')
    
    args = parser.parse_args()
    
    bed_in = args.in_bed
    bed_out = args.bed_out
    chrom = args.chrom
    
    n_cores = int(args.n_cores)
    n_samples = int(args.num_samples)
    
    direction = args.direction
    
    return bed_in, bed_out, chrom, n_cores, n_samples, direction

def run_process():

    bed_in, bed_out, chrom, n_cores, n_samples, direction = get_args()
    print('successfully got arguments')
    
    first_elem = bed_in.pop(0)
    assert(first_elem == 'in_bed')
    print('bed_in: ' + str(bed_in))
    assert(len(bed_in) == n_samples)
    
    total_df = merge_all_samples_bed(bed_in, chrom)
    print('successfully merged bed files for chromosome' + str(chrom))
    
    modified_bed_out = bed_out + '_' + chrom + '_' + direction + '.bed'

    write_result(total_df, modified_bed_out)
    print('successfully saved the result for chromosome' + str(chrom))

    print('successfully done')
    
if __name__ == "__main__":
    run_process()
    print("success")