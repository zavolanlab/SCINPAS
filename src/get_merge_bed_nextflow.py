# -*- coding: utf-8 -*-
"""

Created on Sun Dec 19 22:55:32 2021

@author: Youngbin Moon (y.moon@unibas.ch)
"""











import pysam
import argparse
import pandas as pd
import os

def write_to_bed(final_df, out_file):
    
    final_df.to_csv(out_file, sep = '\t', header = False, index = False,\
    columns = ['seqid', 'start', 'end', 'id', 'score', 'strand'])

def get_args():        
    parser = argparse.ArgumentParser(description="concatenate unique cleavage site bed files for all chromosomes")

    parser.add_argument('in_bed',
      nargs='*', help='Input BED files.')

    parser.add_argument('--bed_template', dest = 'bed_template',
                        required = True,
                        help = 'output bed file name')
               
    args = parser.parse_args()
    
    bed_in = args.in_bed   
    bed_template = args.bed_template
        
    return bed_in, bed_template

def run_process_alter():
    bed_in, bed_template = get_args()
    
    first_elem = bed_in.pop(0)
    assert(first_elem == 'in_bed')
  
    full_dataframe = pd.DataFrame(columns = ['seqid', 'start', 'end', 'id', 'score', 'strand']) 
        
    for file in bed_in:
        if file.endswith(".bed"):
            
            input_bed = pd.read_csv(file, delimiter = '\t', header = 0)

            input_bed.columns = ['seqid', 'start', 'end', 'id', 'score', 'strand']    
                
            print('successfully got 1 input bed file')           
            print('input_bed: ' + str(input_bed))
            
            full_dataframe = pd.concat([full_dataframe, input_bed])
            print('successfully merged 1 bed file')
        
    full_dataframe.sort_values(by=['seqid', 'start', 'end', 'id'], inplace = True)
    
    out_file = bed_template + '.bed'
    write_to_bed(full_dataframe, out_file)   
    print('successfully saved the merged bed file')
        
if __name__ == "__main__":
    run_process_alter()
    print('success')