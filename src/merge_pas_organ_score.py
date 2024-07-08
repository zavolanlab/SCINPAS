# -*- coding: utf-8 -*-
"""
Created on Fri May 31 16:58:34 2024

@author: Youngbin Moon (y.moon@unibas.ch)
"""

import numpy as np
import pandas as pd
import argparse

def write_to_bed(df, out_file):
    df.to_csv(out_file, sep = '\t', header = True, index = False)

# pas has: seqid, start, end, id, score, strand, class, gene_id
# organ_bed has: id and avg score
# From 9th column, it is organ specific column
def concat_organ_score(df, bed_list):
    
    for bed_dir in bed_list:
        organ_bed = pd.read_csv(bed_dir, delimiter = '\t', header = 0)
        df = pd.merge(df, organ_bed, on = ['id'], how = 'left')
    
    copy_df = df.copy()
    # Fill NaN values with 0 in organ columns only (from the 9th column onwards)
    if len(df.columns) > 8:
        # from 9th column, columns are organ columns
        organ_columns = df.columns[8:]
        copy_df[organ_columns] = df[organ_columns].fillna(0)
    
    return copy_df
    
def get_args():        
    
    parser = argparse.ArgumentParser(description="compute average score of a PAS per organ")

    parser.add_argument('in_bed',
      nargs='*', help='organ specific df of id and average organ score.')
    
    parser.add_argument('--pas', dest = 'pas',
                        required = True,
                        help = 'pas bed file')

    parser.add_argument('--out_name', dest = 'out_name',
                        required = True,
                        help = 'bed file out name')

    args = parser.parse_args()
    
    in_bed_dir = args.in_bed
    pas_dir = args.pas
    out_name = args.out_name
    
    print('in_bed_dir: ' + str(in_bed_dir))
    in_bed_dir.pop(0)
    print('in_bed_dir: ' + str(in_bed_dir))
    print('length of in_bed_dir: ' + str(len(in_bed_dir)))
    
    pas = pd.read_csv(pas_dir, delimiter = '\t', header = 0)
    return  pas, in_bed_dir, out_name

def run_process():
    
    pas, in_bed_dir, out_name = get_args()
    print('successfully got inputs')
    
    final_df = concat_organ_score(pas, in_bed_dir)
    print('successfully got the final df')
    
    write_to_bed(final_df, out_name)
    print('successfully wrote the result')
    
if __name__ == "__main__":
    run_process()
    print("success")



