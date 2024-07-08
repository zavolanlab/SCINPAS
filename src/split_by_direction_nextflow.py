# -*- coding: utf-8 -*-
"""

Created on Sun Dec 19 22:55:32 2021

@author: Youngbin Moon (y.moon@unibas.ch)
"""











import argparse
from collections import Counter
from pathlib import Path
import pandas as pd
import numpy as np
import os

def write_BED(out_path: Path, df: pd.DataFrame) -> None:
    df.to_csv(out_path, sep="\t", header=True, index=False)
    
def split_by_direction(df, dir_mrna):
    print('df strand: ' + str(df['strand']))
    print('dir_mRNA: ' + str(dir_mrna))
    subset = df[df['strand'] == dir_mrna].copy()
    print('subset: ' + str(subset))
    
    subset.sort_values(by=['seqid', 'start', 'end', 'id'], inplace = True)
    print('subset after: ' + str(subset))
    return subset
   
def read_BED(file: Path) -> pd.DataFrame:
    df = pd.read_csv(file, sep='\t', header = 0)
    return(df)
    
def get_args():        
    parser = argparse.ArgumentParser(description="divide unique cleavage sites by strandness")

    parser.add_argument('--bed_in', dest = 'bed_in',
                        required = True,
                        help = 'input bed file')
    
    parser.add_argument('--direction', dest = 'direction',
                        required = True,
                        help = 'direction of cs')

    parser.add_argument('--chromo', dest = 'chromo',
                        required = True,
                        help = 'chromosome')
    
    args = parser.parse_args()
    
    bed_in = args.bed_in
    direction = str(args.direction)
    
    chromo = str(args.chromo)
    return bed_in, direction, chromo

def run_process():

    bed_in, direction, chromo = get_args()
    print('successfully got arguments')
    print('bed_in: ' + str(bed_in))
    
    input_df = read_BED(bed_in)
    print('successfully read the df')
    
    subset_df = split_by_direction(input_df, direction)
    print('successfully got subset')
    
    print('subset_df: ' + str(subset_df))
    
    file_name = os.path.basename(bed_in)
    
    number = file_name.split('.')[0].split('_')[-1]
    template = '_'.join(file_name.split('.')[0].split('_')[0:-1])
    
    bed_out = template + 'ForClustering_' + number + '_' + direction + '.bed'
    
    assert(number == chromo)
    
    write_BED(bed_out, subset_df)
    print('successfully saved output')
    
    bed_out2 = template + 'ForOrganGrouping_' + number + '_' + direction + '.bed'

    write_BED(bed_out2, subset_df)
    print('successfully saved output2')
    
    print('successfully done')
    
if __name__ == "__main__":
    run_process()
    print("success")