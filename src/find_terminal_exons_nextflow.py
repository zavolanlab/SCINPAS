# -*- coding: utf-8 -*-

"""
Created on Wed Jan  5 09:37:59 2022

@author: Youngbin Moon (y.moon@unibas.ch)
"""
import pandas as pd
from gtfparse import read_gtf
import argparse

"""
Aim: From genes.gtf retrieve only terminal exons
"""

def read_file (input_dir):
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


def get_last_exons(df):
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
    final_dict = {}
    last_exons = []
    for index, row in df.iterrows():
        if row['feature'] == 'exon' and row['gene_id'] != '' and row['transcript_id'] != '' and row['exon_id'] != '':
            # in cellranger version gtf, key is gene_type rather than gene_biotype
            if row['gene_type'] == 'protein_coding' or row['gene_type'] == 'lncRNA':
                if row['gene_id'] not in final_dict.keys():
                    final_dict[row['gene_id']] = {}
                # check if 2nd key exists
                # 1st key = gene_id
                # 2nd key = transcript_id
                if row['transcript_id'] not in final_dict[row['gene_id']].keys():
                    final_dict[row['gene_id']][row['transcript_id']] = []
                # save as gene_name for visualization in IGV use gene_id for practical use.
                chromosome_id = row['seqname']
                
                if row['transcript_support_level'] == 'NA':
                    score = 2100
                
                elif row['transcript_support_level'] == '':
                    score = 2000
                   
                else:
                    score = row['transcript_support_level']
                    
                final_dict[row['gene_id']][row['transcript_id']].append((chromosome_id, row['start'], row['end'], 
                      row['gene_id'], score, row['strand']))
    
    # get the most distal exon.
    for first_key in final_dict.keys():
        for second_key in final_dict[first_key].keys():
            starts = [elem[1] for elem in final_dict[first_key][second_key]]
            ends = [elem[2] for elem in final_dict[first_key][second_key]]
            strands = [elem[5] for elem in final_dict[first_key][second_key]]
            if strands[0] == '+':
                index = ends.index(max(ends)) 
                last_exon = final_dict[first_key][second_key][index]
            if strands[0] == '-':
                index = starts.index(min(starts))
                last_exon = final_dict[first_key][second_key][index]
            last_exons.append(last_exon)
    final_df = pd.DataFrame(last_exons, columns = ['seqid', 'start', 'end', 'id', 'score', 'strand'])
    final_df.sort_values(by=['seqid', 'id', 'start', 'end'], inplace = True)       
    return final_df


def check_same_exon (current_row, next_row):
    """
    Parameters
    ----------
    current_row : 1 row of the dataframe
        current row of the dataframe containing terminal exon info.
    
    next_row : 1 row of the dataframe
        next row of the dataframe containing terminal exon info.

    Returns
    -------
    bool
        Return true if current row (terminal exon info) is identical to the next row (terminal exon info).
        Since dataframe is sorted, only need to check nearby row to see if terminal exon is "Exactly" identical or not.
    """    
    current_chr = current_row['seqid']
    next_chr = next_row['seqid']
    
    current_start = current_row['start']
    next_start = next_row['start']
    
    current_end = current_row['end']
    next_end = next_row['end']
    
    current_gene = current_row['id']
    next_gene = next_row['id']
    
    current_dir = current_row['strand']
    next_dir = next_row['strand']        
    
    if current_chr == next_chr and current_start == next_start and current_end == next_end and current_gene == next_gene and current_dir == next_dir:
        return True
    else:
        return False

# input: df containing last exons 
# output: remove terminal exons that are identical.
# final_df sorted by chromosome, start, end
def remove_duplicated_exons (final_df):
    """
    Parameters
    ----------
    final_df : dataframe
        dataframe containing terminal exons

    Returns
    -------
    final_deduplicated_df : dataframe
        dataframe containing terminal exons but removed terminal exons that are exactly identical.
    """
    final_deduplicated_df = pd.DataFrame(columns = ['seqid', 'start', 'end', 'id', 'score', 'strand'])
    for index, row in final_df.iterrows():
        if index <= len(final_df)-2: 
            current_row = final_df.iloc[index].copy()
            next_row = final_df.iloc[index+1].copy()
            same = check_same_exon (current_row, next_row)
            
            if same:
                continue
            
            else:
                final_deduplicated_df = final_deduplicated_df.append(current_row, ignore_index = True)
                
    final_deduplicated_df.sort_values(by=['seqid', 'id', 'start', 'end'], inplace = True)          
    
    return final_deduplicated_df     

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
    
    final_df.to_csv(out_file, sep = '\t', header = False, index = False,
    columns = ['seqid', 'start', 'end', 'id', 'score', 'strand'])
    
def get_terminal_args():
    parser = argparse.ArgumentParser(description="get terminal exons bed file")
    parser.add_argument('--gtf_file', dest = 'gtf_file',
                        required = True,
                        help = 'input gtf file')
    
    parser.add_argument('--bed_out', dest = 'bed_out',
                        required = True,
                        help = 'output bed file')    
    
    args = parser.parse_args()
    return args

def get_inputs():
    args = get_terminal_args()
    input_dir = args.gtf_file
    output_file = args.bed_out
    
    return input_dir, output_file

def run_process():
    i_dir, o_file = get_inputs()
    dataframe = read_file(i_dir)
    result_df = get_last_exons(dataframe)
    print('successfully got terminal exons and starting deduplication of exons......')
    
    final_deduplicated_df = remove_duplicated_exons (result_df)
    print('successfully got deduplicated terminal exons and writing to BED file......')
    
    write_bed_file (final_deduplicated_df, o_file)
    print("successfully wrote the bed file")
    
if __name__ == "__main__":
    run_process()
    print('success')