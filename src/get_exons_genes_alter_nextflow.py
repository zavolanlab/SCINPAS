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
    
def filter_genes(df, custom):
    
    if custom:
        g_df = modify_format(df, 'gene')
    
    else:
        g_df = df[df['feature'] == 'gene']
    
    g_df['transcript_support_level'] = [0]*len(g_df)
    unsorted_genes_df = g_df[['seqname', 'start', 'end', 'gene_id', 'transcript_support_level', 'strand', 'gene_type']]
    genes_df = unsorted_genes_df.sort_values(by=['seqname', 'start', 'end'], inplace = False)  
    
    return genes_df    

def filter_exons(df, custom):
    if custom:
        pre_e_df = modify_format(df, 'exon')
    
    else:
        pre_e_df = df[df['feature'] == 'exon']
    
    # Replace NaNs with 2300
    pre_e_df['transcript_support_level'].fillna(2300, inplace=True)
    # Replace empty strings with 2300
    pre_e_df['transcript_support_level'].replace(['', 'NA', ' '], 2300, inplace=True)
    
    # subset
    unsorted_exons_df = pre_e_df[['seqname', 'start', 'end', 'gene_id', 'transcript_support_level', 'strand', 'transcript_id']]
    exons_df = unsorted_exons_df.sort_values(by=['seqname', 'start', 'end'], inplace = False)  
    
    return exons_df

def modify_format(df, feature_type):
    """
    Parameters
    ----------
    df : dataframe
        dataframe that contains our custom genes.gtf infomation.
    
    feature_type : str
        either exon or gene
        
    Returns
    -------
    new_df : dataframe
        dataframe that contains either only exons or only genes. This is to make it easier to use our own custom gtf df.
    """    
    gene_exon_df = df[df['feature'] == feature_type]
    
    attributes = np.asarray(gene_exon_df['attribute'])
    gene_ids = []
    gene_types = []
    gene_names = []    
    transcript_ids = []
    exon_numbers = []
    exon_ids = []
    
    new_df = gene_exon_df.copy()
    for elem in attributes:
        # 1 element of attributes array
        each_attributes = elem.split(';')[0:-1]
        
        if feature_type == 'gene':
            # within each element, you have several attributes (e.g. gene_id)
            for each_attribute in each_attributes:
                if 'gene_id' in each_attribute:
                    gene_id = each_attribute.split(' ')[1]
                    gene_ids.append(gene_id)
                    
                elif 'gene_type' in each_attribute:
                    gene_type = each_attribute.split(' ')[2]
                    gene_types.append(gene_type)
                
                elif 'gene_name' in each_attribute:
                    gene_name = each_attribute.split(' ')[2]
                    gene_names.append(gene_name)

        elif feature_type == 'exon':
            # within each element, you have several attributes (e.g. gene_id)
            for each_attribute in each_attributes:
                if 'gene_id' in each_attribute:
                    gene_id = each_attribute.split(' ')[1]
                    gene_ids.append(gene_id)

                elif 'transcript_id' in each_attribute:
                    transcript_id = each_attribute.split(' ')[2]
                    transcript_ids.append(transcript_id)
                    
                elif 'gene_type' in each_attribute:
                    gene_type = each_attribute.split(' ')[2]
                    gene_types.append(gene_type)
                
                elif 'gene_name' in each_attribute:
                    gene_name = each_attribute.split(' ')[2]
                    gene_names.append(gene_name)

                elif 'exon_number' in each_attribute:
                    exon_number = each_attribute.split(' ')[2]
                    exon_numbers.append(exon_number)
                
                elif 'exon_id' in each_attribute:
                    exon_id = each_attribute.split(' ')[2]
                    exon_ids.append(exon_id)            
                              
    if feature_type == 'gene':            
        new_df['gene_id'] = gene_ids
        new_df['gene_type'] = gene_types
        new_df['gene_name'] = gene_names

    elif feature_type == 'exon':            
        new_df['gene_id'] = gene_ids
        new_df['transcript_id'] = transcript_ids
        new_df['gene_type'] = gene_types
        new_df['gene_name'] = gene_names
        new_df['exon_number'] = exon_numbers
        new_df['exon_id'] = exon_ids
    
    new_df.drop('attribute', inplace=True, axis=1)
    return new_df

def read_and_filter(input_dir, custom):
    """
    Parameters
    ----------
    input_dir : str
        directory towards input genes.gtf file.
        
    custom : bool
        whether to use our custom gtf or not
        
    Returns
    -------
    genes_df : dataframe
        gtf dataframe that contains only genes

    exons_df : dataframe
        gtf dataframe that contains only exons       
    """
    
    print('custom: ' + str(custom))
    if custom:
        input_df = pd.read_csv(input_dir, sep='\t',\
                              names = ['seqname', 'source', 'feature', 'start', 'end', 'transcript_support_level', 'strand', 'frame', 'attribute'], low_memory = False)
        
    else:
        input_df = read_gtf(input_dir)
    
    genes_df = filter_genes(input_df, custom)
    exons_df = filter_exons(input_df, custom)
        
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

    parser.add_argument('--custom', dest = 'custom',
                        required = True,
                        help = 'whether to use custom gtf or not')   
        
    args = parser.parse_args()
    
    gtf_dir = args.gtf_dir
    
    genes_bed_out = args.genes_bed_out
    exons_bed_out = args.exons_bed_out
    do_custom = bool(int(args.custom))
    
    # gtf_dir = r"C:/Users/Geniu/Desktop/success/NEXT_PROJECT_CATALOG/SCINPAS_ALL_SAMPLES/temp_results/240126/extended_merged_noNA_recovered_gtf.gtf"
    # genes_bed_out = r"C:/Users/Geniu/Downloads/gene_out_custom.bed"
    # exons_bed_out = r"C:/Users/Geniu/Downloads/exon_out_custom.bed"
    # do_custom = True
    
    return gtf_dir, do_custom, genes_bed_out, exons_bed_out

def run_process():
    
    gtf_dir, do_custom, genes_bed_out, exons_bed_out = get_input_noGX()
    print('successfully got inputs')
    
    filtered_genes_df, filtered_exons_df = read_and_filter(gtf_dir, do_custom)
    print('successfully got filtered exons and filtered genes')
    
    write_as_bed_file(filtered_genes_df, genes_bed_out)
    print('successfully saved genes_bed file')
        
    write_as_bed_file(filtered_exons_df, exons_bed_out)
    print('successfully saved exons_bed file')
         
if __name__ == "__main__":
    run_process()
    print('success')