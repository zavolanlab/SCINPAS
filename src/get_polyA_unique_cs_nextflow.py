# -*- coding: utf-8 -*-
"""

Created on Sun Oct 22 22:55:32 2023

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
def write_to_bed(final_df, out_file, all_mode, sample):
    """
    Parameters
    ----------
    final_df : dataframe
        dataframe containing terminal exons but removed terminal exons that are exactly identical.
    
    out_file : string
        output file name. This output will be in the bed format.

    all_mode : bool
        whether use all samples mode or not (i.e. use RPM and # of experiments that support this cs)

    sample : str
        sample name
                                             
    Returns
    -------
    returns nothing but writes the output in the bed format.
    """
    if all_mode:
        
        final_df.to_csv(out_file, sep = '\t', header = True, index = False,\
                        columns = ['seqid', 'start', 'end', 'id', 'score', 'strand', 'supp', sample])

    else:
        final_df.to_csv(out_file, sep = '\t', header = True, index = False,\
                        columns = ['seqid', 'start', 'end', 'id', 'score', 'strand'])
            
def convert_dict_to_df(cleavage_dictionary, all_mode, num_reads, sample):
    """
    Parameters
    ----------
    cleavage_dictionary : dictionary
        key = cleavage site id = chromosome:fixed_cleavage_site:direction    
        value = how many times that cleavage site occured
    
    all_mode : bool
        whether use all samples mode or not (i.e. use RPM and # of experiments that support this cs)
    
    num_reads : int
        the number of pA reads in the sample
    
    sample : str
        sample name
        
    Returns
    -------
    final_df : dataframe
        a dataframe where each row = chromosome, start position (fixed cleavage site - 1), end position (fixed_cleavage site), cleavage_site_id, frequency, direction
        if all_mode == True,
        score = RPM and you have additional column, supp_exp, which is the number of samples supporting this cleavage site.
    """    
    cleavage_list = []
    for cleavage_id in cleavage_dictionary.keys():
        chrom = cleavage_id.split(':')[0]
        cleavage_site = cleavage_id.split(':')[1]
        start = str(int(cleavage_site) - 1)
        end = str(int(cleavage_site))
        # start = str(int(cleavage_site))
        # end = str(int(cleavage_site) + 1) 
        
        direction = cleavage_id.split(':')[2]
        
        print('all_mode: ' + str(all_mode))
        if all_mode:
            score = cleavage_dictionary[cleavage_id]
            # RPM
            rpm = cleavage_dictionary[cleavage_id] * (10**6)/num_reads
            if rpm > 0:
                supp = 1
            else:
                supp = 0
                
            cleavage_list.append((chrom, start, end, cleavage_id, rpm, direction, supp, score))
            
        else:
            score = cleavage_dictionary[cleavage_id]
            cleavage_list.append((chrom, start, end, cleavage_id, score, direction))
    
    if all_mode:
        final_df = pd.DataFrame(cleavage_list, columns = ['seqid', 'start', 'end', 'id', 'score', 'strand', 'supp', sample])
   
    else:
        final_df = pd.DataFrame(cleavage_list, columns = ['seqid', 'start', 'end', 'id', 'score', 'strand'])
    
    final_df.sort_values(by=['seqid', 'id', 'start', 'end'], inplace = True)
    final_df['start'] = final_df['start'].astype(int)
    final_df['end'] = final_df['end'].astype(int)
    final_df['supp'] = final_df['supp'].astype(int)
    final_df[sample] = final_df[sample].astype(int)
    
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

    parser.add_argument('--split', dest = 'split',
                        required = True,
                        help = 'whether use chr-splitted bam or full bam')

    parser.add_argument('--multiple_samples', dest = 'multiple_samples',
                        default = 0,
                        help = 'whether use all samples mode or not (i.e. use RPM and # of experiments that support a particular cs')

    parser.add_argument('--count_dir', dest = 'count_dir', default = 0,
                        help = 'directory towards csv file containing the number of pA reads in 1 sample')
    
    args = parser.parse_args()
    
    bam_dir = args.bam
    
    bed_out = args.bed_out
    use_fc = args.use_fc
    
    split = args.split
    
    multiple_samples = args.multiple_samples
    count_dir = args.count_dir
    
    if multiple_samples:
        count_df = pd.read_csv(count_dir, delimiter = ',', low_memory=False, names = ['sample', 'type', 'count'])
        pA_count = count_df.loc[0, 'count']
        print('pA_count: ' + str(pA_count))
        print('multiple samples true')
        
    else:
        pA_count = 0
        print('multiple samples false')
        
    return bam_dir, bed_out, use_fc, split, multiple_samples, pA_count

def run_process():

    bam_dir, bed_out, use_fc, split, multiple_samples, pA_count = get_args()
    print('successfully got arguments')
    
    bam = pysam.AlignmentFile(bam_dir, "rb")
    
    unique_cleavage_site = get_num_unique_cleavage(bam, use_fc)
    print('successfully got unique cleavage sites dictionary')
    
    sample_name = '_'.join(bed_out.split('.')[0].split('_')[0:3])
    
    final_df = convert_dict_to_df(unique_cleavage_site, multiple_samples, pA_count, sample_name)
    print('successfully converted the dictionary into dataframe format')
    
    chrom_number = split
    out_template = bed_out.split('.')[0]
    bed_out = out_template + '_' + chrom_number + '.bed'
        
    write_to_bed(final_df, bed_out, multiple_samples, sample_name)
    print('successfully saved the result')
    print('successfully done')
    
if __name__ == "__main__":
    run_process()
    print("success")