# -*- coding: utf-8 -*-
"""
Created on Mon Jul 18 19:25:43 2024

@author: Youngbin Moon (y.moon@unibas.ch)
"""











import argparse
import pandas as pd

def write_to_bed(final_df, out_file):
    
    final_df.to_csv(out_file, sep = '\t', header = False, index = False)
    
def concat_all(beds_list):
    all_dfs_list = []
    for file in beds_list:
        df = pd.read_csv(file, delimiter = '\t', header = 0)
        all_dfs_list.append(df)
    
    final_df = pd.concat(all_dfs_list)
    return final_df

def get_args():        
    parser = argparse.ArgumentParser(description="concat all bed files over all chromosomes and directions")

    parser.add_argument('in_bed',
      nargs='*', help='pas bed file for each chromosome and direction')             

    parser.add_argument('--out', dest = 'out',
                        required = True,
                        help = 'output name')
    
    args = parser.parse_args()
    
    bed_in = args.in_bed
    first_elem = bed_in.pop(0)
    
    out = args.out
    
    return bed_in, out

def run_process_alter():
    bed_in, out = get_args()
    print('successfully got inputs')
    
    final_df = concat_all(bed_in)
    print('successfully got final df')
    
    write_to_bed(final_df, out)
    print('successfully saved the final df')
    
if __name__ == "__main__":
    run_process_alter()
    print('success')


