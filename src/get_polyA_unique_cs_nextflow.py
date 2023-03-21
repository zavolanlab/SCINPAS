# -*- coding: utf-8 -*-
"""

Created on Sun Dec 19 22:55:32 2021

@author: Youngbin Moon (y.moon@unibas.ch)
"""

import pysam
import argparse
from collections import Counter
import pandas as pd
import numpy as np




"""
Aim : write a bed file of unique cleavage sites.
This is to make a cluster. Either by single linkage clustering or 
"""
def write_to_bed(final_df, out_file):
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
    
def convert_dict_to_df(cleavage_dictionary):
    """
    Parameters
    ----------
    cleavage_dictionary : dictionary
        key = cleavage site id = chromosome:fixed_cleavage_site:direction    
        value = how many times that cleavage site occured
        
    Returns
    -------
    final_df : dataframe
        a dataframe where each row = chromosome, start position (fixed cleavage site - 1), end position (fixed_cleavage site), cleavage_site_id, frequency, direction
    """    
    cleavage_list = []
    for cleavage_id in cleavage_dictionary.keys():
        chrom = cleavage_id.split(':')[0]
        end = cleavage_id.split(':')[1]
        start = str(int(end) - 1)
        score = cleavage_dictionary[cleavage_id]
        direction = cleavage_id.split(':')[2]
        
        cleavage_list.append((chrom, start, end, cleavage_id, score, direction))
        
    final_df = pd.DataFrame(cleavage_list, columns = ['seqid', 'start', 'end', 'id', 'score', 'strand'])
    final_df.sort_values(by=['seqid', 'id', 'start', 'end'], inplace = True)           
    
    return final_df

def get_num_unique_cleavage(sam, use_FC):
    """
    Parameters
    ----------
    bam : bam file
        a bam file containing polyA reads.
    
    use_FC : bool
        use fixed cleavage site or not.
        
    Returns
    -------
    unique_cleavage_site : dictionary
        key = cleavage site id = chromosome:fixed_cleavage_site:direction    
        value = how many times that cleavage site occured    
    """      
    unique_cleavage_site = Counter()
    for read in sam.fetch():
        chrom = read.reference_name
        rev = read.is_reverse
        alignedRefPositions = read.get_reference_positions()
        refStart = alignedRefPositions[0]
        refEnd = alignedRefPositions[-1]
        if rev == True:
            rev = '-'
            if not use_FC:
                # assert(refStart == int(read.get_tag('OC')))
                assert(refStart + 1 == int(read.get_tag('XO')))
                end = int(read.get_tag('XO'))
            
            elif use_FC:
                end = int(read.get_tag('XF'))      
        else:
            rev = '+'
            if not use_FC:
                # assert(refEnd == int(read.get_tag('OC')))
                assert(refEnd + 1 == int(read.get_tag('XO')))
                end = int(read.get_tag('XO'))
            
            elif use_FC:
                end = int(read.get_tag('XF'))            
                
        cleavage_site_id = str(chrom) + ":" + str(end) + ":" + str(rev)
        unique_cleavage_site[cleavage_site_id] += 1    
    
    return unique_cleavage_site
        
def get_args():        
    parser = argparse.ArgumentParser(description="get number of polyA sites")

    parser.add_argument('--bam', dest = 'bam',
                        required = True,
                        help = 'bam file containing polyA reads')
       
    parser.add_argument('--bed_out', dest = 'bed_out',
                        required = True,
                        help = 'output number of polyA sites csv file name')

    parser.add_argument('--use_fc', dest = 'use_fc',
                        required = True,
                        help = 'use fixed cleavage site or not')
    
    args = parser.parse_args()
    
    bam_dir = args.bam
    bam = pysam.AlignmentFile(bam_dir, "rb")  
    
    bed_out = args.bed_out
    use_fc = args.use_fc
    return bam, bed_out, use_fc

def run_process():

    bam, bed_out, use_fc = get_args()
    print('successfully got arguments')
    
    unique_cleavage_site = get_num_unique_cleavage(bam, use_fc)
    print('successfully got unique cleavage sites dictionary')
    
    final_df = convert_dict_to_df(unique_cleavage_site)
    print('successfully converted the dictionary into dataframe format')
    
    write_to_bed(final_df, bed_out)
    print('successfully saved the result')
    
if __name__ == "__main__":
    run_process()
    print("success")