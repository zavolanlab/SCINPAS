# -*- coding: utf-8 -*-
"""

Created on Sun Dec 19 22:55:32 2021

@author: Youngbin Moon (y.moon@unibas.ch)
"""
import pandas as pd
import argparse
import csv






"""
Aim : get total number of genes in the species.
"""
def write_out(row_data, o_file):
    """
    Parameters
    ----------
    row_data : list
        1st element = sample name
        2nd element = number of genes covered by sample
        3rd element = % of genes overlapped between sample and gtf
            
    o_file : string
        output file name    
       
    Returns
    -------        
    returns nothing but writes row_data at the csv file with output name: o_file.
    """        
    w = csv.writer(open(o_file, "w"))
    w.writerow(row_data)
    
def get_total_number_genes(bed, Species):
    """
    Parameters
    ----------
    bed : bed file
        a bed file containing all terminal exons of the species of interest.
        (generated from gtf file)
        
    Species : str
        species. e.g. mouse, human......
       
    Returns
    -------        
    total_num_genes : int
        total number of genes in that species.
    """ 
    
    if Species == 'mouse':
        numbers = list(range(1, 20))
        numbers += ['X', 'Y']
        
    elif Species == 'human':
        numbers = list(range(1, 23))
        numbers += ['X', 'Y']
    
    chromosomes = ['chr'+ str(number) for number in numbers]
    unique_genes = set()
    total_num_genes = 0
    
    bed.columns = ['seqid', 'start', 'end', 'id', 'score', 'strand']
    
    for chromosome in chromosomes:
        rows = bed[bed['seqid'] == chromosome]
        filtered_rows = rows[rows['score'] <= 3]
        unique_genes.update(filtered_rows['id'])
    
    total_num_genes = len(list(unique_genes))
    
    return total_num_genes

def get_args():
        
    parser = argparse.ArgumentParser(description="get number of genes covered by more than 1 read")

    parser.add_argument('--terminal_exons', dest = 'terminal_exons',
                        required = True,
                        help = 'terminal exons.bed input file')
            
    parser.add_argument('--out_name', dest = 'out_name',
                        required = True,
                        help = 'output csv file name')

    parser.add_argument('--species', dest = 'species',
                        required = True,
                        help = 'species of the data')
    
    args = parser.parse_args()

    terminal_exons_dir = args.terminal_exons
    terminal_exons = pd.read_csv(terminal_exons_dir, delimiter = '\t', header = None) 
        
    out_name = args.out_name
    species = args.species
        
    return terminal_exons, out_name, species

def run_process():

    terminal_exons, out_name, species = get_args()
    print('successfully got inputs')
    
    tot_num_genes = get_total_number_genes(terminal_exons, species)
    print('total number of genes expressed')
    
    sample_name = 'total # of genes'
    final_data = [sample_name, tot_num_genes]
    
    write_out(final_data, out_name)
    print('successfully saved the result')
    
if __name__ == "__main__":
    run_process()
    print("success")