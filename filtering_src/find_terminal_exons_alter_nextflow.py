# -*- coding: utf-8 -*-

"""
Created on Wed Jan  5 09:37:59 2022

@author: Youngbin Moon (y.moon@unibas.ch)
"""
import pandas as pd
from gtfparse import read_gtf
import argparse
import numpy as np
import multiprocessing as mp
from collections import Counter
import time
import itertools
from multiprocessing import Pool
"""
Aim: From genes.gtf retrieve only terminal exons
"""

# Note: since you did merge overlapping exons, you wont have duplicate terminal exons. No need to check them.

def write_bed_file (final_df, out_file):
    """
    Parameters
    ----------
    final_df : dataframe
        dataframe containing terminal exons but removed terminal exons that are exactly identical.
    
    out_file : string
        output file name. This output will be in the bed format.

    Returns
    -------
    returns nothing but writes the output in the bed format.
    """
    
    final_df.to_csv(out_file, sep = '\t', header = False, index = False)

def check_and_remove_duplicate(grouped_df):
    
    assert(len(grouped_df) >= 1)
    if len(grouped_df) == 1:
        return grouped_df
    
    elif len(grouped_df) > 1:
        print('duplicate terminal exons exist')
        index = np.argmin(np.asarray(grouped_df['transcript_support_level']))
        print(grouped_df.iloc[index, :].to_frame().T)
        return grouped_df.iloc[index, :].to_frame().T
    
def remove_duplicated_exons (df, n_cores):
    
    groupedBy_df = df.groupby(['seqname', 'start', 'end', 'gene_id', 'strand'])
    pool = mp.Pool(n_cores)
    contents = [contents_of_group for name_of_group, contents_of_group in groupedBy_df]
    print(df)
    with Pool(n_cores) as pool:
        result = pool.starmap(check_and_remove_duplicate, zip(contents))
        
    final_df = pd.concat(result)
    
    final_df.sort_values(by=['seqname', 'start', 'end', 'gene_id'], inplace = True)      
    final_df['transcript_support_level'].fillna(2300, inplace =  True)
    final_df['transcript_support_level'].replace(['', 'NA', ' '], 2300, inplace = True)
       
    return final_df
      
def get_last_exon(group_df):
    strands = list(set(group_df['strand']))
    assert(len(strands) == 1)
    direction = strands[0]
    # print('group_df: ' + str(group_df))
    last_exons_df = pd.DataFrame(columns = group_df.columns)
    if direction == '+':
        # index is integer array. have to use iloc. array is length of 1.
        index = np.argmax(np.asarray(group_df['end']))
        # print(group_df.iloc[index, :].to_frame().T)
        last_exons_df = pd.concat([last_exons_df, group_df.iloc[index, :].to_frame().T])

    
    elif direction == '-':
        # index is integer array. have to use iloc. array is length of 1.
        index = np.argmin(np.asarray(group_df['start']))
        # print(group_df.iloc[index, :].to_frame().T)
        last_exons_df = pd.concat([last_exons_df, group_df.iloc[index, :].to_frame().T])
    
    # print(last_exons_df)
    return last_exons_df

def get_last_exons(exons_df, n_cores):
    """
    Parameters
    ----------
    df : dataframe
        dataframe that contains genes.gtf infomation.

    Returns
    -------
    final_df : dataframe
        dataframe that contains only terminal exons. (a subset of genes.gtf)
    """    
    groupedBy_geneId_transcriptId = exons_df.groupby(['gene_id', 'transcript_id'])
    pool = mp.Pool(n_cores)
    contents = [contents_of_group for name_of_group, contents_of_group in groupedBy_geneId_transcriptId]
    names = [name_of_group[0] for name_of_group, contents_of_group in groupedBy_geneId_transcriptId]

    with Pool(n_cores) as pool:
        result = pool.starmap(get_last_exon, zip(contents))
    
    all_last_exons_df = pd.concat(result)
    
    all_last_exons_df.sort_values(by=['seqname', 'start', 'end', 'gene_id'], inplace = True)
    # print(all_last_exons_df)
    
    all_last_exons_df = all_last_exons_df[['seqname', 'start', 'end', 'gene_id', 'transcript_support_level', 'strand']]
    
    return all_last_exons_df

def read_file (input_dir, species):
    """
    Parameters
    ----------
    input_dir : str
        directory towards input genes.gtf file.
                
    Returns
    -------
    input_df : dataframe
        dataframe that contains genes.gtf infomation.
    """
            
    intermediate_df = read_gtf(input_dir)
    exon_df = intermediate_df[intermediate_df['feature'] == 'exon']
    
    if species == 'worm':
        input_df = exon_df[['seqname', 'source', 'feature', 'start', 'end', 'strand', 'frame', 'gene_id', 'transcript_id', 'gene_biotype', 'gene_name', 'exon_number', 'exon_id']]
        # worm gtf does not have transcript_support_level
        input_df['transcript_support_level'] = 2300
        
    else:
        input_df = exon_df[['seqname', 'source', 'feature', 'start', 'end', 'transcript_support_level', 'strand', 'frame', 'gene_id', 'transcript_id', 'gene_type', 'gene_name', 'exon_number', 'exon_id']]
    
    return input_df
    
def get_terminal_args():
    parser = argparse.ArgumentParser(description="get terminal exons bed file")
    parser.add_argument('--gtf_file', dest = 'gtf_file',
                        required = True,
                        help = 'input gtf file')
    
    parser.add_argument('--bed_out', dest = 'bed_out',
                        required = True,
                        help = 'output bed file')    
    
    parser.add_argument('--n', dest = 'n',
                        required = True,
                        help = 'number of cores')

    parser.add_argument('--species', dest = 'species',
                        required = True,
                        help = 'species')
    
    args = parser.parse_args()
    return args

def get_inputs():
    args = get_terminal_args()
    input_dir = args.gtf_file
    output_file = args.bed_out
    n = int(args.n)
    species = args.species
    # input_dir = r"C:/Users/Geniu/Desktop/success/NEXT_PROJECT_CATALOG/SCINPAS_ALL_SAMPLES/temp_results/240126/extended_merged_noNA_recovered_gtf.gtf"  
    # output_file = r"C:/Users/Geniu/Downloads/sample_TE_removed_duplicate_final.bed"
    # do_custom = True
    # n = 4
    return input_dir, output_file, n, species

def run_process():
    i_dir, o_file, n, species = get_inputs()
    print('successfully got inputs')
    
    dataframe = read_file(i_dir, species)
    print('successfully read the file')
    
    result_df = get_last_exons(dataframe, n)
    print('successfully got terminal exons and starting deduplication of exons......')
    
    final_df = remove_duplicated_exons(result_df, n)
    print('successfully remvoed duplicate terminal exon')
    
    write_bed_file (final_df, o_file)
    print("successfully wrote the bed file")
        
if __name__ == "__main__":
    run_process()
    print('success')