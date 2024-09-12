# -*- coding: utf-8 -*-
"""
Created on Mon Jul 18 19:25:43 2024

@author: Youngbin Moon (y.moon@unibas.ch)
"""











import argparse
import pandas as pd

def write_to_bed(final_df, out_file):
    
    final_df.to_csv(out_file, sep = '\t', header = False, index = False)

def get_args():        
    parser = argparse.ArgumentParser(description="change into real BED format so that it is compatible with bedtools")

    parser.add_argument('--bed_in', dest = 'bed_in',
                        required = True,
                        help = 'bed_in')               
    
    args = parser.parse_args()
    
    bed_in = args.bed_in
        
    return bed_in

def run_process_alter():
    bed_in = get_args()
      
    input_df = pd.read_csv(bed_in, delimiter = '\t', header = 0)
    
    subset_df = input_df.iloc[:, [0,1,2,3,4,5]]
    print(subset_df)
    
    bed_template = 'modified_v2_full_' + bed_in.split('.')[0]
    out_file = bed_template + '.bed'
    write_to_bed(subset_df, out_file)
    
    print('successfully saved the merged bed file')
        
if __name__ == "__main__":
    run_process_alter()
    print('success')