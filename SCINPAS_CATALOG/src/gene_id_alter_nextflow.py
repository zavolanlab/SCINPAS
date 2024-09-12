# -*- coding: utf-8 -*-
"""
Created on Fri May 31 16:58:34 2024

@author: Youngbin Moon (y.moon@unibas.ch)
"""













import numpy as np
import pandas as pd
import pybedtools
import argparse
import matplotlib.pyplot as plt
import seaborn as sns

def write_to_bed(df, out_file):
    # print(df)
    df.to_csv(out_file, sep = '\t', header = True, index = False)

def change_header(df):
    df_final = df.rename(columns = {
        'pas_seqid': 'seqid',
        'pas_start': 'start',
        'pas_end': 'end',
        'pas_id': 'id',
        'pas_score': 'score',
        'pas_strand': 'strand',
        'pas_class': 'class',
        'gene_id': 'gene_id'
    })
    
    return df_final

def compute_priority(row):
    
    if row['gene_type'] == 'protein_coding':
        priority = 1
    
    elif row['gene_type'] == 'lncRNA':
        priority = 2
    
    else:
        priority = 3
    
    return priority

def intersect(pas_dir, genes_dir):
    
    pas_bed = pybedtools.BedTool(pas_dir)
    genes_bed = pybedtools.BedTool(genes_dir)
    
    # Intersect PAS with genes on the same strand, retaining information from both
    # it will have both columns from pas_bed and columns from genes_bed for the pas that overlaps
    intersection = pas_bed.intersect(genes_bed, s = True, wa = True, wb = True)
    
    # Create a DataFrame from the intersections
    # Adjust column names based on the exact format of your files
    columns = ['pas_seqid', 'pas_start', 'pas_end', 'pas_id', 'pas_score', 'pas_strand', 'pas_class',\
               'gene_chr', 'gene_start', 'gene_end', 'gene_id', 'gene_score', 'gene_strand', 'gene_type']
    
    # intersection is a pybedTools.BedTool object. 
    # .fn property gets the filename of the termporary file where it contains saved output of intersection
    df = pd.read_table(intersection.fn, header = None, names = columns)
    
    # compute priority along the rows using 'apply'. axis =1 means do it for all rows
    # You can have a pas that overlaps with multiple genes and intersection reports all overlaps
    # Hence you will have duplicates of the same pas
    # In such case you remove duplicates by prioritizing protein_coding, lncRNA first
    df['priority'] = df.apply(compute_priority, axis = 1)
    
    # sort by priority so that you can keep the first (highest priority)
    # subset = [......] defines what is the duplicate
    df = df.sort_values(by = ['pas_seqid', 'pas_start', 'pas_end', 'pas_id', 'pas_strand', 'priority'], ascending=[True, True, True, True, True, True])\
        .drop_duplicates(subset = ['pas_seqid', 'pas_start', 'pas_end', 'pas_id', 'pas_strand'], keep = 'first')
    
    intersect_df = df[['pas_seqid', 'pas_start', 'pas_end', 'pas_id', 'pas_score', 'pas_strand', 'pas_class', 'gene_id']]
    
    final_intersect_df = change_header(intersect_df)
    print('successfully changed the header')
    
    no_intersection = pas_bed.intersect(genes_bed, s= True, wa = True, v = True)
    no_intersect_df = pd.read_table(no_intersection.fn, header = None, names = ['seqid', 'start', 'end', 'id', 'score', 'strand', 'class'])
    no_intersect_df['gene_id'] = 'not_available'
    
    df_final = pd.concat([final_intersect_df, no_intersect_df], ignore_index = True)
    df_final.sort_values(by = ['seqid', 'start', 'end', 'id', 'strand'], inplace = True)
    
    return df_final

def get_args():        
    parser = argparse.ArgumentParser(description="assign gene name that PAS belongs to ")

    parser.add_argument('--pas', dest = 'pas',
                        required = True,
                        help = 'pas bed file')

    parser.add_argument('--genes', dest = 'genes',
                        required = True,
                        help = 'bed file format of gtf containing genes')
    
    parser.add_argument('--out', dest = 'out',
                        required = True,
                        help = 'output bed file name')
       
    args = parser.parse_args()
    
    pas = args.pas
    genes = args.genes
    out_template = args.out
    
    # pas = r"/scicore/home/zavolan/moon0000/AdditionalColAllsamples_polyA_cluster_out_1_+.bed"
    # genes = r"/scicore/home/zavolan/moon0000/filtered_genes.bed"
    # out_template = r"/scicore/home/zavolan/moon0000/intersection_pas_1_+"   
    
    return pas, genes, out_template

def run_process():

    pas, genes, out_template = get_args()
    print('successfully got arguments')
    
    final_df = intersect(pas, genes)
    print('successfully got final df')
        
    out_name = out_template + '_intersect_out.bed'
    write_to_bed(final_df, out_name)
    print('successfully saved the result')
    
if __name__ == "__main__":
    run_process()
    print("success")

