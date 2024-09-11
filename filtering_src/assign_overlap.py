# -*- coding: utf-8 -*-
"""
Created on Mon Jul 18 19:25:43 2024

@author: Youngbin Moon (y.moon@unibas.ch)
"""











import argparse
import pandas as pd

def write_to_bed(final_df, out_file):
    
    final_df.to_csv(out_file, sep = '\t', header = True, index = False)

def add_col(df, col_val):
    df['overlap'] = col_val
    return df
    
def get_args():        
    parser = argparse.ArgumentParser(description="add a column to pas bed file indicating whether pas is mapping to overlapping, nonoverlapping or intergenic")

    parser.add_argument('--non_overlapping', dest = 'non_overlapping',
                        required = True,
                        help = 'pas mapping to non_overlapping regions')             

    parser.add_argument('--overlapping', dest = 'overlapping',
                        required = True,
                        help = 'pas mapping to overlapping region')

    parser.add_argument('--intergenic', dest = 'intergenic',
                        required = True,
                        help = 'pas mapping to intergenic region')
    
    parser.add_argument('--out', dest = 'out',
                        required = True,
                        help = 'output name')
    
    args = parser.parse_args()
    
    non_overlapping = args.non_overlapping
    overlapping = args.overlapping
    intergenic = args.intergenic
    out = args.out
    
    non_overlapping_pas = pd.read_csv(non_overlapping, delimiter = '\t', names = ['seqid', 'start', 'end', 'id', 'score', 'strand'])
    overlapping_pas = pd.read_csv(overlapping, delimiter = '\t', names = ['seqid', 'start', 'end', 'id', 'score', 'strand'])
    intergenic_pas = pd.read_csv(intergenic, delimiter = '\t', names = ['seqid', 'start', 'end', 'id', 'score', 'strand'])
    
    return non_overlapping_pas, overlapping_pas, intergenic_pas, out

def run_process_alter():
    non_overlapping_pas, overlapping_pas, intergenic_pas, out = get_args()
    print('successfully got inputs')
    
    non_overlapping_pas_modified = add_col(non_overlapping_pas, 'non_overlapping')
    print('successfully add new column to nonoverlapping pas')
    
    overlapping_pas_modified = add_col(overlapping_pas, 'overlapping')
    print('successfully add new column to overlapping pas')

    intergenic_pas_modified = add_col(intergenic_pas, 'intergenic')
    print('successfully add new column to intergenic pas')
    
    final_df = pd.concat([non_overlapping_pas_modified, overlapping_pas_modified, intergenic_pas_modified], ignore_index = True)
    print('successfully got full dataframe')
    
    write_to_bed(final_df, out)
    print('successfully saved the final df')
    
if __name__ == "__main__":
    run_process_alter()
    print('success')


