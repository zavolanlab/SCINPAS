# -*- coding: utf-8 -*-
"""
Created on Mon Jul 22 13:44:47 2024

@author: Youngbin Moon (y.moon@unibas.ch)
"""











import argparse
import pandas as pd
from multiprocessing import Pool
import os
import numpy as np

def write_to_bed(final_df, out_file):
    final_df.to_csv(out_file, sep = '\t', header = True, index = False)

def get_pas_info(pas):
    chrom = pas['seqid']
    start = pas['start']
    end = pas['end']
    gene_id = pas['reassigned_g']
    strand = pas['strand']
    return chrom, start, end, gene_id, strand

def find_class(pas, te, exons, introns):
    chrom, start, end, gene_id, strand = get_pas_info(pas)
    print(f"chrom: {chrom}, start: {start}, end: {end}, gene_id: {gene_id}, strand: {strand}")
    # Extend start and end positions by 1 bp
    start -= 1
    end += 1

    # Check for intersection with introns
    # First, find the transcripts related to the gene id
    transcripts = exons[exons['id'] == gene_id]['transcript_id'].unique()
    print('transcripts: ' + str(transcripts))
    # introns['id'].isin(transcripts) checks each id in the introns DataFrame to see if it is present in the transcripts array. 
    # The isin method returns True for each id in introns that is found in transcripts and False otherwise.
    in_intersect = introns[(introns['id'].isin(transcripts)) & (introns['seqid'] == chrom) & 
                           (introns['start'] <= end) & (introns['end'] >= start) & (introns['strand'] == strand)]
    
    print('in_intersect: ' + str(in_intersect))
    if not in_intersect.empty:
        return 'intronic'

    # Check for intersection with exons
    ex_intersect = exons[(exons['id'] == gene_id) & (exons['seqid'] == chrom) & (exons['start'] <= end) & 
                         (exons['end'] >= start) & (exons['strand'] == strand)]
    
    print('ex_intersect: ' + str(ex_intersect))
    if not ex_intersect.empty:
        return 'exonic'
    
    # check for intersection with terminal exons
    te_intersect = te[(te['id'] == gene_id) & (te['seqid'] == chrom) & (te['start'] <= end) & 
                      (te['end'] >= start) & (te['strand'] == strand)]
    
    print('te_intersect: ' + str(te_intersect))
    if not te_intersect.empty:
        return 'TE'
    
    print('you reached here......')
    print('it is because you reassigned gene')
    return 'ambiguous'

def process_chunk(chunk, te, exons, introns):
    print("Processing the following chunk:")
    print(chunk)
    copy_df = chunk.copy()
    copy_df['re_class'] = chunk.apply(lambda row: find_class(row, te, exons, introns), axis = 1)
    return copy_df

def get_all_classes(df, te, exons, introns):
    num_cores = os.cpu_count()    
    print('num_cores: ' + str(num_cores))
    # Split the DataFrame into chunks based on the number of cores
    # It also works for dataframe
    # return value: a list of dataframe
    num_chunks = min(num_cores, len(df))
    chunks = np.array_split(df, num_chunks)
    
    with Pool(num_cores) as pool:
        results = pool.starmap(process_chunk, [(chunk, te, exons, introns) for chunk in chunks])    
    return pd.concat(results)

def get_args():        
    parser = argparse.ArgumentParser(description="assign PAS gene ID according to 2 criterion and change class of PAS")
           
    parser.add_argument('--pas', dest = 'pas',
                        required = True,
                        help = 'total pas bed')
    
    parser.add_argument('--exons', dest = 'exons',
                        required = True,
                        help = 'exons.bed')

    parser.add_argument('--te', dest = 'te',
                        required = True,
                        help = 'terminal_exons.bed')
    
    parser.add_argument('--introns', dest = 'introns',
                        required = True,
                        help = 'introns.bed')

    parser.add_argument('--out', dest = 'out',
                        required = True,
                        help = 'output name')
    
    parser.add_argument('--chrom', dest = 'chrom',
                        required = True,
                        help = 'chromosome')

    parser.add_argument('--strand', dest = 'strand',
                        required = True,
                        help = 'direction of PAS')   

    parser.add_argument('--species', dest = 'species',
                        required = True,
                        help = 'species for example human, mouse etc')
    
    args = parser.parse_args()
    
    pas_dir = args.pas
    exons_dir = args.exons
    te_dir = args.te
    introns_dir = args.introns
    out = args.out
    chrom = args.chrom
    strand = args.strand
    species = args.species
    
    if species == 'worm':
        chromosome = chrom
    else:    
        chromosome = 'chr' + chrom
        
    pas = pd.read_csv(pas_dir, delimiter = '\t', header = 0)
    exons = pd.read_csv(exons_dir, delimiter = '\t', names = ['seqid', 'start', 'end', 'id', 'score', 'strand', 'transcript_id'])
    te = pd.read_csv(te_dir, delimiter = '\t', names = ['seqid', 'start', 'end', 'id', 'score', 'strand'])
    introns = pd.read_csv(introns_dir, delimiter = '\t', names = ['seqid', 'start', 'end', 'id', 'score', 'strand'])
    
    exons_subset = exons[(exons['seqid'] == chromosome) & (exons['strand'] == strand)]
    te_subset = te[(te['seqid'] == chromosome) & (te['strand'] == strand)]
    introns_subset = introns[(introns['seqid'] == chromosome) & (introns['strand'] == strand)]
    
    out_name = out + '_' + chrom + '_' + strand + '.bed'

    return pas, exons_subset, te_subset, introns_subset, out_name

def run_process_alter():
    pas, exons_subset, te_subset, introns_subset, out_name = get_args()
    print('successfully got inputs')
    
    final_df = get_all_classes(pas, te_subset, exons_subset, introns_subset)
    print('successfully got classes for all pas')
    
    write_to_bed(final_df, out_name)
    print('successfully saved the output')
    
if __name__ == "__main__":
    run_process_alter()
    print('success')