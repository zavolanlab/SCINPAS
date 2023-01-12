# -*- coding: utf-8 -*-
"""

Created on Fri May  6 19:05:00 2022

@author: Youngbin Moon (y.moon@unibas.ch)
"""
from gtfparse import read_gtf
import pandas as pd
import argparse

"""
Aim 1: get bed file that contains genes only.
Aim 2: get bed file that contains exons (including terminal exons) only. 
"""

def read_csvfile (input_dir):
    """
    Parameters
    ----------
    input_dir : str
        directory towards input genes.gtf file.

    Returns
    -------
    df : dataframe
        dataframe that contains genes.gtf infomation.
    """
    
    df = read_gtf (input_dir)
    return df

def filter_genes(df):
    """
    Parameters
    ----------
    df : dataframe
        dataframe that contains genes.gtf infomation.

    Returns
    -------
    result_df : dataframe
        dataframe that only contains genes.     
    """
    
    filtered_rows = []
    for index, row in df.iterrows():
        # if feature == 'gene', there is no transcript support level
        if row['feature'] == 'gene' and row['gene_id'] != '':
            if row['gene_type'] == 'protein_coding' or row['gene_type'] == 'lncRNA':
                filtered_rows.append((row['seqname'], row['start'], row['end'], row['gene_id'], 0, row['strand']))
            
    result_df = pd.DataFrame(filtered_rows, columns = ['seqid', 'start', 'end', 'id', 'score', 'strand'])
    result_df.sort_values(by=['seqid', 'start', 'end'], inplace = True)  
    return result_df

def filter_exons(df):
    """
    Parameters
    ----------
    df : dataframe
        dataframe that contains genes.gtf infomation.

    Returns
    -------
    result_df : dataframe
        dataframe that only contains exons (including terminal exons).
    """
    filtered_rows = []
    for index, row in df.iterrows():
        if row['feature'] == 'exon' and row['gene_id'] != '' and row['transcript_id'] != '' and row['exon_id'] != '':
            # newly added line
            if row['transcript_support_level'] != 'NA' and row['transcript_support_level'] != '' and int(float(row['transcript_support_level'])) <= 3:
                # newly added line
                if row['gene_type'] == 'protein_coding' or row['gene_type'] == 'lncRNA':
                    filtered_rows.append((row['seqname'], row['start'], row['end'], \
                                          row['gene_name'], row['transcript_support_level'], row['strand']))
    
    result_df = pd.DataFrame(filtered_rows, columns = ['seqid', 'start', 'end', 'id', 'score', 'strand'])
    result_df.sort_values(by=['seqid', 'start', 'end'], inplace = True)       
    return result_df

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
    final_df.to_csv(out_file, sep = '\t', header = False, index = False,
    columns = ['seqid', 'start', 'end', 'id', 'score', 'strand'])

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
        
    args = parser.parse_args()
    
    gtf_dir = args.gtf_dir
    gtf = read_csvfile(gtf_dir)
    
    genes_bed_out = args.genes_bed_out
    exons_bed_out = args.exons_bed_out

    return gtf, genes_bed_out, exons_bed_out

def run_process():
    
    gtf, genes_bed_out, exons_bed_out = get_input_noGX()
    print('successfully got inputs')
    
    genes_df = filter_genes(gtf)
    print('successfully got genes bed dataframe')
    
    write_as_bed_file(genes_df, genes_bed_out)
    print('successfully saved genes_bed file')
    
    exons_df = filter_exons(gtf)
    print('successfully got exons_bed dataframe')
    
    write_as_bed_file(exons_df, exons_bed_out)
    print('successfully saved exons_bed file')
            
if __name__ == "__main__":
    run_process()
    print('success')