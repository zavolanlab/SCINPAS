# -*- coding: utf-8 -*-
"""
Created on Mon Mar 18 11:36:39 2024

@author: Youngbin Moon (y.moon@unibas.ch)
"""













import argparse
from collections import Counter
import pandas as pd
import numpy as np
import multiprocessing as mp
from multiprocessing import Pool

def write_to_bed(df, out_file):
    
    df.to_csv(out_file, sep = '\t', header = False, index = False)
    
def get_intronic_region_per_gene_transcript(df):
    
    print('df: ' + str(df))
    chroms = df['seqid'].values[:-1] # Adjusted to have one fewer element
    starts = df['start'].values
    ends = df['end'].values
    # It does not matter that much because we use this to do bedtools intersect 
    # but it is better to use transcript_id to identify transcripts within the same gene
    ids = df['transcript_id'].values[:-1] # Adjusted to have one fewer element
    scores = df['score'].values[:-1] # Adjusted to have one fewer element
    strands = df['strand'].values[:-1] # Adjusted to have one fewer element
        
    # from first element to 2nd last element
    # bed file is 0 indexed (start, end]
    intronic_starts = ends[:-1]
    intronic_ends = starts[1:] -1
    
    assert(len(chroms) == len(intronic_starts) == len(intronic_ends) == len(ids) == len(scores) == len(strands))    
    
    # if it is same transcript from same gene: chromosome, ids (gene_id or transcript_id), tsl (score), strands should be the same 
    data = np.column_stack((chroms, intronic_starts, intronic_ends, ids, scores, strands))
    final_df = pd.DataFrame(data, columns=['seqid', 'start', 'end', 'id', 'score', 'strand'])
    print('final_df: ' + str(final_df))
    
    return final_df

# chunk is a list of tuple (list of df). 2nd element of a tuple is df (grouped by gene id and transcript id)
# e.g. [(name21, df21), (name22, df22), (name23, df23)]
def process_chunk(chunk):
    
    # group = df21 or df22 or df23
    sub_results = [get_intronic_region_per_gene_transcript(group) for _, group in chunk]    
    
    # sub_total = pd.concat(sub_results)
    # sub_total is partial merged output for df21, df22, df23.
    print('successfully processed 1 chunk')
    return pd.concat(sub_results)

def chunkify(data, n):

    chunk_size = max(1, len(data) // n)
    # data[i:i+chunk_size] returns subset list
    chunks = [data[i:i+chunk_size] for i in range(0, len(data), chunk_size)]
    return chunks

def get_all_intronic_regions(df, n_cores):

    groupedBy_gene_transcript = df.groupby(['id', 'transcript_id'])
    
    # each element is df
    grouped_list = list(groupedBy_gene_transcript)
    # list of list. Each element is a list of df
    # e.g. [[(name18, df18), (name19, df19), (name20, df20)], [(name21, df21), (name22, df22), (name23, df23)], .....]
    chunks = chunkify(grouped_list, n_cores)
    print('successfully chunked')
    
    with Pool(n_cores) as pool:
        results = pool.map(process_chunk, chunks)
    
    # get total merged output for df18, df19, df20, df21, df22 and df23
    total_intronic_df = pd.concat(results)
    
    total_intronic_df.sort_values(by=['seqid', 'start', 'end'], inplace=True, ignore_index=True)
    print('total_intronic_df: ' + str(total_intronic_df))
    
    # check the number of gene id and transcript id pairs
    # num_transcripts_original >> num_transcripts_final becauase many transcripts have 1 exon only.
    # If you have 1 exon, there is no intron
    num_transcripts_original = len(grouped_list)
    num_transcripts_final = len(list(total_intronic_df.groupby(['id'])))
    
    print('num_transcripts original: ' + str(num_transcripts_original))
    print('num_transcripts final: ' + str(num_transcripts_final))
    
    return total_intronic_df

def get_args():        
    
    parser = argparse.ArgumentParser(description="get intronic positions.bed")

    parser.add_argument('--exon_bed', dest = 'exon_bed',
                        required = True,
                        help = 'bed file containing exons only')
       
    parser.add_argument('--out_name', dest = 'out_name',
                        required = True,
                        help = 'bed file name containing intronic_region')

    parser.add_argument('--n', dest = 'n',
                        required = True,
                        help = 'number of cores')
    
    args = parser.parse_args()
    
    exon_bed_dir = args.exon_bed
    out_name = args.out_name
    n = int(args.n)
    
    # exon_bed_dir = r"C:/Users/Geniu/Downloads/filtered_exons.bed"
    # out_name = r"C:/Users/Geniu/Downloads/test_intronic.bed"
    # n = 1
    # Specify data types for each column
    dtype_dict = {
        'seqid': 'category',
        'start': 'int32',
        'end': 'int32',
        'id': 'category',
        'score': 'float32',
        'strand': 'category',
        'transcript_id': 'category'
    }
    exons_df = pd.read_csv(exon_bed_dir, sep='\t',\
                          names = ['seqid', 'start', 'end', 'id', 'score', 'strand', 'transcript_id'], dtype = dtype_dict, low_memory = False)
    
    exons_df['score'].fillna(2300, inplace = True)
    exons_df['score'].replace(['', 'NA', ' '], 2300, inplace = True)
    
    return exons_df, out_name, n

def run_process():

    exons_df, out_name, n = get_args()
    print('successfully got arguments')
    
    total_intronic_df = get_all_intronic_regions(exons_df, n)
    print('successfully got the total intronic df')
    
    write_to_bed(total_intronic_df, out_name)
    print('successfully saved the output')
    
if __name__ == "__main__":
    run_process()
    print("success")