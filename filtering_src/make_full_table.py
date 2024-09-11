# -*- coding: utf-8 -*-
"""
Created on Tue Aug  6 17:01:34 2024

@author: Youngbin Moon (y.moon@unibas.ch)
"""











import argparse
import pandas as pd

def write_to_bed(final_df, out_file):
    final_df.to_csv(out_file, sep = '\t', header = True, index = False)
    
def add_presence_column(full_df, subset_df, column_name):
    subset_df[column_name] = 1
    merged_df = pd.merge(full_df, subset_df[['seqid', 'start', 'end', 'id', 'score', 'strand', column_name]],
                         on = ['seqid', 'start', 'end', 'id', 'score', 'strand'], how = 'left')
    merged_df[column_name] = merged_df[column_name].fillna(0).astype(int)    
    return merged_df

def left_join(original_pas, beds_list):
    
    for file in beds_list:
        threshold = int(file.split('/')[-1].split('_')[0])
        print(file.split('/')[-1].split('_')[0])
        column_name = f"{threshold}_MP"
        subset_df = pd.read_csv(file, delimiter = '\t', header = 0)
        original_pas = add_presence_column(original_pas, subset_df, column_name)       
    
    return original_pas

def get_args():        
    parser = argparse.ArgumentParser(description="generate total polyAsite table with column info about filtered pas at different motif percentage threshold")

    parser.add_argument('in_bed',
      nargs='*', help='a list of filtered pas at different motif percentage threshold')
           
    parser.add_argument('--original_pas', dest = 'original_pas',
                        required = True,
                        help = 'total pas bed')

    parser.add_argument('--out', dest = 'out',
                        required = True,
                        help = 'output name')
    
    args = parser.parse_args()
    
    bed_in = args.in_bed
    original_pas_dir = args.original_pas
    out = args.out  
    
    # bed_in = ['in_bed', r'C:/Users/Geniu/Downloads/80_human_filtered_pas_v3_random_test.bed', r'C:/Users/Geniu/Downloads/20_human_filtered_pas_v3_random_test.bed']
    # original_pas_dir = r'C:/Users/Geniu/Downloads/kendall_reassigned_gene_pas_8_+.bed'
    # out = 'C:/Users/Geniu/Downloads/test_merged_pas_240807_2.bed'
    first_elem = bed_in.pop(0)
    assert(first_elem == 'in_bed')
    
    pas = pd.read_csv(original_pas_dir, delimiter = '\t', header = 0)
    
    return bed_in, pas, out

def run_process_alter():
    in_bed, pas, out = get_args()
    print('successfully got inputs')
    
    final_df = left_join(pas, in_bed)
    print('successfully make full dataframe')
    
    write_to_bed(final_df, out)
    print('successfully saved the output')
    
    print(final_df[['id', '20_MP', '80_MP']])
    
if __name__ == "__main__":
    run_process_alter()
    print('success')