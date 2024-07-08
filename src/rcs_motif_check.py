# -*- coding: utf-8 -*-
"""
Created on Sun Jan 26 22:55:32 2024

@author: Youngbin Moon (y.moon@unibas.ch)
"""













import argparse
import pandas as pd
import pysam
import multiprocessing
import itertools

def write_to_bed(df, out_file):
    df.to_csv(out_file, sep = '\t', header = True, index = False)  
    
def all_motives_presence(df, motives):
    motives_list = list(motives['motif'])
    
    copy_df =  df.copy()
    # sum over all motif columns per each row
    copy_df['all_motif'] = df.loc[:, motives_list].sum(axis = 1).astype(int)
    # Change value of all_motif column where 'all_motif' > 0
    # This is to use mean as an estimator in pointplot 
    # mean would be the motif_percentage and we get CI from pointplot    
    copy_df.loc[copy_df['all_motif'] > 0, 'all_motif'] = 1
    
    # get 'supp' column which means the number of samples supporting that PAS
    # Apply the lambda function to the 'id' column (for each row = x)
    copy_df['supp'] = copy_df['id'].apply(lambda x: x.split(':')[-1])
    return copy_df

def change_dna_to_rna(sequence, direction):
    """
    Parameters
    ----------    
    sequence : string
        a current genome DNA sub-sequence with same length as motif 
        (always + strand because reference genome is always + strand).
        
        we need to convert this DNA into RNA so that we can decide
        whether this sub-sequence is identical to the motif or not.
        
    direction : character
        direction of DNA in which a read maps to.
        
    Returns
    -------
    corrected_string : string
        RNA version (5' -> 3') of the current genome DNA sub-sequence.
        Now you can directly compare it with the motif.
        
        i.e. change a subesequence of DNA into RNA so that it becomes compatible with motif (5' -> 3')
    """     
    corrected_sequence = []
    # if a read is mapping to - strand,
    # revert the DNA sequence and then make a complementary
    if direction == '-':
        reverted_subsequence = list(reversed(sequence))
        for elem in reverted_subsequence:
            if elem == 'A':
                corrected_sequence.append('U')
            elif elem == 'T':
                corrected_sequence.append('A')
            elif elem == 'G':
                corrected_sequence.append('C')
            elif elem == 'C':
                corrected_sequence.append('G')
    
    # if a read is mapping to + strand
    # only change T in the DNA -> U in RNA. other nucleotides stay the same
    elif direction == '+':
        for elem in sequence:
            if elem == 'A':
                corrected_sequence.append('A')
            elif elem == 'T':
                corrected_sequence.append('U')
            elif elem == 'G':
                corrected_sequence.append('G')
            elif elem == 'C':
                corrected_sequence.append('C')
                
    # convert a list of characters into a single string
    corrected_string = ''.join(corrected_sequence)
    return corrected_string

def get_extra_col(group, motives_df, fasta_path):
    fasta_f = pysam.FastaFile(fasta_path)
    name = group[0]
    pas_df = group[1]
    rcs_list = list(pas_df['id'])

    for index, row in motives_df.iterrows():
        motif = row['motif']
        # upper and lower boundaries are negative values
        # e.g. upper = -35, lower = -10
        upper = int(row['upper'])
        lower = int(row['lower'])
        print(f'Motif: {motif}, Upper: {upper}, Lower: {lower}')
        values = []
        for rcs_id in rcs_list:
            chrom = rcs_id.split(':')[0]
            rcs = int(rcs_id.split(':')[1])
            strand = rcs_id.split(':')[2]
            if strand == '-':
                # rcs + 10
                sequence_start = rcs - lower
                # rcs + 35
                sequence_end = rcs - upper
            
            elif strand == '+':
                # rcs - 35
                sequence_start = rcs + upper
                # rcs - 10
                sequence_end = rcs + lower

            # Add debug statements to verify the computed positions
            print(f'Processing {rcs_id}:')
            print(f'Chromosome: {chrom}')
            print(f'Strand: {strand}')
            print(f'Sequence start: {sequence_start}')
            print(f'Sequence end: {sequence_end}')
            
            dna_subsequence = fasta_f.fetch(reference=chrom, start=sequence_start, end=sequence_end + 1)
            print(f'fetching sequence for {rcs_id} with chrom {chrom}, start {sequence_start}, end {sequence_end}')
            # correct DNA subsequence into RNA so that it becomes compatible with motif
            corrected_subsequence = change_dna_to_rna(dna_subsequence, strand)
            if motif in corrected_subsequence:
                print('corrected_subsequence match : ' + str(corrected_subsequence))
                print('motif match: ' + str(motif))                
                value = 1
            
            else:
                print('corrected_subsequence no match : ' + str(corrected_subsequence))
                print('motif no match: ' + str(motif))
                value = 0
            
            values.append(value)
        
        assert(len(values) == len(rcs_list))
        # append a new column
        pas_df[motif] = values
    
    return pas_df

def add_columns_serial(pas_df, motives_df, fasta_path):
    final_dfs = []
    groups = pas_df.groupby(['seqid', 'strand'])
    
    for group in groups:
        final_df = get_extra_col(group, motives_df, fasta_path)
        final_dfs.append(final_df)
    
    final_df = pd.concat(final_dfs)
    return final_df

def get_args():
    parser = argparse.ArgumentParser(description="add extra column to the representative cleavage site BED file which designates the presence of a motif in these cleavage sites")

    parser.add_argument('--rcs_dir', dest = 'rcs_dir',
                        required = True,
                        help = 'input pas bed file')
    
    parser.add_argument('--motif_dir', dest = 'motif_dir',
                        required = True,
                        help = 'directory that contains motif info')

    parser.add_argument('--fasta_dir', dest = 'fasta_dir',
                        required = True,
                        help = 'directory towards fasta file')
    
    parser.add_argument('--rcs_out', dest = 'rcs_out',
                        required = True,
                        help = 'representative cleavage sites in BED format having extra column')
    
    args = parser.parse_args()
    
    pas_dir = args.rcs_dir
    motif_dir = args.motif_dir
    fasta_dir = args.fasta_dir
    rcs_out = args.rcs_out
    
    # pas_dir = r"/scicore/home/zavolan/moon0000/pas_to_rcs_21_+.bed"
    # motif_dir = r"/scicore/home/zavolan/moon0000/intergenic_analysis/motif_info_2.csv"
    # fasta_dir = r"/scicore/home/zavolan/moon0000/CATALOG/data/human/genome.fa"
    # rcs_out = r"/scicore/home/zavolan/moon0000/rsc_out_ex_col_21_+.bed"
    pas = pd.read_csv(pas_dir, delimiter = '\t', header = 0)
    motives = pd.read_csv(motif_dir, delimiter = ',', header = 0)
    
    return pas, motives, fasta_dir, rcs_out

def run_process():
    pas, motives, fasta_dir, rcs_out = get_args()
    print('successfully got inputs')
    
    intermediate_df = add_columns_serial(pas, motives, fasta_dir)
    print('successfully added motif columns')
    
    final_df = all_motives_presence(intermediate_df, motives)
    print('final_df: ' + final_df.to_string())
    
    write_to_bed(final_df, rcs_out)
    print('successfully saved the result')
    
if __name__ == '__main__':
    run_process()
    print('success')  