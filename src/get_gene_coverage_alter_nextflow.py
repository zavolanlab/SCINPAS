# -*- coding: utf-8 -*-
"""

Created on Sun Dec 19 22:55:32 2021

@author: Youngbin Moon (y.moon@unibas.ch)
"""
import pysam
import argparse
from collections import Counter
import csv
import pandas as pd




"""
Aim : Use all polyA reads mapping to terminal exons in order to get the number of genes covered by more than k read(s).
where k is a parameter.
i.e. number of genes covered by sample.

Aim2: Using deduplicated bam file, estimate the number of genes expressed in the sample.

Note: gene should exist in both read (GX_tag) and in gtf (terminal_exons.bed)
"""

def write_out(row_data, o_file):
    """
    Parameters
    ----------
    row_data : list
        1st element = sample name
        2nd element = number of genes covered by sample
        3rd element = % of genes overlapped between sample and gtf
            
    o_file : string
        output file name    
     
    Returns
    -------        
    returns nothing but writes row_data at the csv file with output name: o_file.
    """        
    w = csv.writer(open(o_file, "w"))
    w.writerow(row_data)
        
def check_GX_overlap(tag, bed):
    """
    Parameters
    ----------
    tag : string
        current GX tag (e.g. ENSMUSG00000022995).
    
    bed : bed file
        a bed file containing all terminal exons from gtf.
        
    Returns
    -------
    True if GX tag from the read also exist in the bed file.
    False if GX tag from the read does not exist in the bed file.
    """         
    bed.columns = ['seqid', 'start', 'end', 'id', 'score', 'strand']
    filtered_exons = bed[bed['id'] == tag]
    
    if len(filtered_exons) >= 1:
        return True
    
    else:
        return False
    
def get_num_genes_covered(GX_counter, k, bedFile):
    """
    Parameters
    ----------
    GX_counter : counter object
        key = GX tag. (e.g. ENSMUSG00000022995)
        value = count how many times that GX tag (e.g.ENSMUSG00000022995) appeared.
    
    k : int
        a threshold of how many reads should at least support a gene 
        in order to consider this gene as covered.

    bedFile : bed file
        a bed file containing all terminal exons from gtf.  
        
    Returns
    -------
    num_genes_covered : int
        the number of genes that are covered by more than k reads.
    """
    num_genes_covered = 0
    for GX_tag in GX_counter.keys():
        count = GX_counter[GX_tag]

        if count >= k:
            is_overlap = check_GX_overlap(GX_tag, bedFile)
            
            if is_overlap:    
                num_genes_covered += 1    
    
    return num_genes_covered
 
def convert_dict(list_of_tuples):
    """
    Parameters
    ----------
    list_of_tuples : a list of tuples
        an object that contains infomation about 1 read from BAM file.
        1st element of a tuple: tag name. For example, 'GX' tag
        2nd element of a tuple: value of the tag. if it is 'GX' tag,
        then value is gene name(s) that a read maps to. For example, ENSMUSG00000022995.
        
    Returns
    -------
    tag_dict : dictionary
        key = tag name. e.g. 'GX' tag.
        value = value that the tag has. e.g. if key is 'GX', value is gene name.
    """
    tag_dict = dict()
    for key, tag in list_of_tuples:
        tag_dict.setdefault(key, []).append(tag)
    return tag_dict

def get_GX_dicts(sam):
    """
    Parameters
    ----------
    sam : bam file
        an input bam file which contains either:
        1) all polyA reads mapping to terminal exon within a sample.
        2) or dedup bam file
    
    Returns
    -------
    GX_dict : counter object
        key = GX tag. (e.g. ENSMUSG00000022995)
        value = count how many times that GX tag (e.g.ENSMUSG00000022995) appeared.
    """    
    GX_dicts = Counter()   
    
    for read in sam.fetch():
        list_tuples = read.tags
        tag_dict = convert_dict(list_tuples)  
        
        if 'GX' in tag_dict.keys():
            # GX_tags is a list of tags
            GX_tags = tag_dict['GX']
            # GX_tags = [ENSMUSG00000022995;ENSMUSG00000025779]            
            if GX_tags[0] != GX_tags[-1] or len(GX_tags) > 1:
                print('you have more than 1 element')
                
            # GX_tags = [ENSMUSG00000022995, ENSMUSG00000025779]    
            GX_tags = GX_tags[0].split(';')      
            # length of GX_tags in general is 1 (sometimes you can have more than 2 gene tags for a given read)
            # for each gene
            for GX_tag in GX_tags:
                # For Tabula muris, GX tag is written ENSMUSG00000023143.10
                # but algorithm uses ENSMUSG00000023143
                # if a tag doesnt have ".", it will just output the same thing.
                GX_tag = GX_tag.split(".")[0]
                GX_dicts[GX_tag] += 1
    
    return GX_dicts
    
def get_args():
        
    parser = argparse.ArgumentParser(description="get number of genes covered by more than 1 read")

    parser.add_argument('--bam', dest = 'bam',
                        required = True,
                        help = 'bam file containing all polyA reads + terminal exon or deduplicated bam file')
            
    parser.add_argument('--out_name', dest = 'out_name',
                        required = True,
                        help = 'output csv file name')

    parser.add_argument('--read_count_threshold', dest = 'read_count_threshold',
                        required = True,
                        help = 'a threshold of how many reads should at least support a gene in order to consider this gene as covered.')

    parser.add_argument('--terminal_exons', dest = 'terminal_exons',
                        required = True,
                        help = 'terminal exons.bed input file')

    parser.add_argument('--type_read', dest = 'type_read',
                        required = True,
                        help = 'type of reads. either polyA + terminal exon or dedup')
    
    args = parser.parse_args()
    
    bam_dir = args.bam
    bam = pysam.AlignmentFile(bam_dir, "rb")  
    
    out_name = args.out_name
    read_count_threshold = int(args.read_count_threshold)
    
    terminal_exons_dir = args.terminal_exons
    terminal_exons = pd.read_csv(terminal_exons_dir, delimiter = '\t', header = None)
    
    type_read = args.type_read
    return bam, out_name, read_count_threshold, terminal_exons, type_read

def run_process():

    bam, out_name, read_count_threshold, terminal_exons, type_read = get_args()
    print('successfully got inputs')

    # e.g. 10X_P4_7
    sample_name = '_'.join(out_name.split('_')[0:3])
    
    if "UmiDedup" in sample_name:
        sample_name = 'negative_control'
        
    GX_dicts = get_GX_dicts(bam)
    print('successfully got GX_dicts')
        
    if type_read == 'polyA':
        
        num_genes_covered = get_num_genes_covered(GX_dicts, read_count_threshold, terminal_exons)
        print('successfully got number of genes covered from polyA reads + terminal exons')
                
        final_data = [sample_name, num_genes_covered]
    
    elif type_read == 'dedup':
        
        num_genes_covered = get_num_genes_covered(GX_dicts, 1, terminal_exons)
        print('successfully got number of genes expressed in the dedup bam file')        
    
        final_data = [sample_name, num_genes_covered]
    
    write_out(final_data, out_name)
    print('successfully saved the result')
    
if __name__ == "__main__":
    run_process()
    print("success")