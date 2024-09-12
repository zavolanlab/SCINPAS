# -*- coding: utf-8 -*-
"""
Created on Mon Jul 22 13:44:47 2024

@author: Youngbin Moon (y.moon@unibas.ch)
"""











import argparse
import pandas as pd

def write_to_bed(final_df, out_file):
    final_df.to_csv(out_file, sep = '\t', header = True, index = False)
    
def merge_all_df(beds_list, intergenic):
    all_dfs = []
    for file in beds_list:
        df = pd.read_csv(file, delimiter = '\t', header = 0)
        all_dfs.append(df)
    
    all_dfs.append(intergenic)
    merged_df = pd.concat(all_dfs)
    
    final_df = merged_df.sort_values(by=['seqid', 'start', 'end', 'id'])
    return final_df

def get_args():        
    parser = argparse.ArgumentParser(description="merge all genic PAS.bed across chromosomes and directions and then concat with intergenic PAS")

    parser.add_argument('in_bed',
      nargs='*', help='genic PAS with reassigned gene_id and class for each chromosome and direction')
           
    parser.add_argument('--intergenic', dest = 'intergenic',
                        required = True,
                        help = 'intergenic pas that is assigned to a gene')
    
    parser.add_argument('--out', dest = 'out',
                        required = True,
                        help = 'output name')
    
    args = parser.parse_args()
    
    bed_in = args.in_bed
    intergenic_dir = args.intergenic
    out = args.out  
    
    first_elem = bed_in.pop(0)
    assert(first_elem == 'in_bed')
    
    intergenic_pas = pd.read_csv(intergenic_dir, delimiter = '\t', header = 0)

    return bed_in, intergenic_pas, out

def run_process():
    bed_in, intergenic_pas, out = get_args()
    print('successfully got inputs')
    
    final_df = merge_all_df(bed_in, intergenic_pas)
    print('successfully computed merged df')
    
    write_to_bed(final_df, out)
    print('successfully saved the output')
    
if __name__ == "__main__":
    run_process()
    print('success')