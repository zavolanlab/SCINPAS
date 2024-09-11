# -*- coding: utf-8 -*-
"""
Created on Fri Aug  2 19:21:43 2024

@author: Youngbin Moon (y.moon@unibas.ch)
"""











import pandas as pd
import argparse

def write_to_bed(final_df, out_file):
    final_df.to_csv(out_file, sep = '\t', header = False, index = False)
    
def select_row(group):
    if len(group) == 1:
        return group.iloc[0]
    # Check for the highest priority group
    highest_priority = group['priority'].min()
    highest_priority_genes = group[group['priority'] == highest_priority].reset_index(drop = True)
    # Prioritize ENSG genes if available in the highest priority group
    ensg_rows = highest_priority_genes[highest_priority_genes['id'].str.startswith('ENSG')]
    if not ensg_rows.empty:
        return ensg_rows.sample(n=1).iloc[0]
    else:
        return highest_priority_genes.sample(n=1).iloc[0]
    
def filter_all_duplicates(df):
    df['priority'] = df['gene_type'].apply(lambda x: 1 if x=='protein_coding' else 2 if x=='lncRNA' else 3)
    # Group by the specified columns
    grouped = df.groupby(['seqid', 'start', 'end', 'strand'])
    
    # Apply the function to each group and create a new DataFrame
    result = grouped.apply(select_row).reset_index(drop = True)
    
    return result
def get_args():        
    parser = argparse.ArgumentParser(description="remove duplicates in genes with same chr, start, end, strands")
           
    parser.add_argument('--genes', dest = 'genes',
                        required = True,
                        help = 'filtered_genes.bed')
    
    parser.add_argument('--out', dest = 'out',
                        required = True,
                        help = 'output name')
    
    args = parser.parse_args()
    
    genes_dir = args.genes
    out = args.out  
    
    genes = pd.read_csv(genes_dir, delimiter = '\t', names = ['seqid', 'start', 'end', 'id', 'score', 'strand', 'gene_type'])
    
    return genes, out

def run_process_alter():
    genes, out = get_args()
    print('successfully got inputs')
    
    intermediate_df = filter_all_duplicates(genes)
    print('successfully removed duplicates')

    final_df = intermediate_df[['seqid', 'start', 'end', 'id', 'score', 'strand', 'gene_type']]
    
    write_to_bed(final_df, out)
    print('successfully saved the result')
    
if __name__ == "__main__":
    run_process_alter()
    print('success')