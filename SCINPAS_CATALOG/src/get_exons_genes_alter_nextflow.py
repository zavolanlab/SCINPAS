# -*- coding: utf-8 -*-
"""

Created on Fri May  6 19:05:00 2022

@author: Youngbin Moon (y.moon@unibas.ch)
"""
from gtfparse import read_gtf
import pandas as pd
import argparse
import numpy as np
"""
Aim 1: get bed file that contains genes only.
Aim 2: get bed file that contains exons (including terminal exons) only. 
"""
def write_as_bed_file (final_df, out_file):
    """
    Parameters
    ----------
    final_df : dataframe
        dataframe containing either genes only (or exons only)
    
    out_file : string
        Output file name. This output will be in bed format.

    Returns
    -------
    returns nothing but writes the output in the bed format.
    """
    # chr, start, end are compulsory for bed file. others are accessory. adjust according to your needs
    final_df.to_csv(out_file, sep = '\t', header = False, index = False)
    
def filter_genes(df, species):
    
    g_df = df[df['feature'] == 'gene']
    # genes do not have transcript support level by default.
    g_df['transcript_support_level'] = [0]*len(g_df)
    
    if species == 'worm':
        unsorted_genes_df = g_df[['seqname', 'start', 'end', 'gene_id', 'transcript_support_level', 'strand', 'gene_biotype']]
        
    else:
        unsorted_genes_df = g_df[['seqname', 'start', 'end', 'gene_id', 'transcript_support_level', 'strand', 'gene_type']]
        
    genes_df = unsorted_genes_df.sort_values(by=['seqname', 'start', 'end'], inplace = False)  
    
    return genes_df    

def filter_exons(df, species):
    
    pre_e_df = df[df['feature'] == 'exon']
    
    # worm does not have transcript_support_level. Need to add them.
    if species == 'worm':
        pre_e_df['transcript_support_level'] = 2300
    
    # Replace NaNs with 2300
    pre_e_df['transcript_support_level'].fillna(2300, inplace=True)
    # Replace empty strings with 2300
    pre_e_df['transcript_support_level'].replace(['', 'NA', ' '], 2300, inplace=True)
    
    # subset
    unsorted_exons_df = pre_e_df[['seqname', 'start', 'end', 'gene_id', 'transcript_support_level', 'strand', 'transcript_id']]
    exons_df = unsorted_exons_df.sort_values(by=['seqname', 'start', 'end'], inplace = False)  
    
    return exons_df

def read_and_filter(input_dir, species):
    """
    Parameters
    ----------
    input_dir : str
        directory towards input genes.gtf file.
        
    Returns
    -------
    genes_df : dataframe
        gtf dataframe that contains only genes

    exons_df : dataframe
        gtf dataframe that contains only exons       
    """
        
    input_df = read_gtf(input_dir)
    
    genes_df = filter_genes(input_df, species)
    exons_df = filter_exons(input_df, species)
        
    return genes_df, exons_df
    
def get_input_noGX():
    
    parser = argparse.ArgumentParser(description="get bed files that either contain genes or exons exclusively.")
    
    parser.add_argument('--gtf_dir', dest = 'gtf_dir',
                        required = True,
                        help = 'genes.gtf file directory')
      
    parser.add_argument('--genes_bed_out', dest = 'genes_bed_out',
                        required = True,
                        help = 'bed file that only contains genes')

    parser.add_argument('--exons_bed_out', dest = 'exons_bed_out',
                        required = True,
                        help = 'bed file that only contains exons')  

    parser.add_argument('--species', dest = 'species',
                        required = True,
                        help = 'species')
        
    args = parser.parse_args()
    
    gtf_dir = args.gtf_dir
    
    genes_bed_out = args.genes_bed_out
    exons_bed_out = args.exons_bed_out
    species = args.species
    # gtf_dir = r"C:/Users/Geniu/Desktop/success/NEXT_PROJECT_CATALOG/SCINPAS_ALL_SAMPLES/temp_results/240126/extended_merged_noNA_recovered_gtf.gtf"
    # genes_bed_out = r"C:/Users/Geniu/Downloads/gene_out_custom.bed"
    # exons_bed_out = r"C:/Users/Geniu/Downloads/exon_out_custom.bed"
    # do_custom = True
    
    return gtf_dir, genes_bed_out, exons_bed_out, species

def run_process():
    
    gtf_dir, genes_bed_out, exons_bed_out, species = get_input_noGX()
    print('successfully got inputs')
    
    filtered_genes_df, filtered_exons_df = read_and_filter(gtf_dir, species)
    print('successfully got filtered exons and filtered genes')
    
    write_as_bed_file(filtered_genes_df, genes_bed_out)
    print('successfully saved genes_bed file')
        
    write_as_bed_file(filtered_exons_df, exons_bed_out)
    print('successfully saved exons_bed file')
         
if __name__ == "__main__":
    run_process()
    print('success')