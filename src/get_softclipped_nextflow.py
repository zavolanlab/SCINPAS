# -*- coding: utf-8 -*-
"""

Created on Thu Dec  8 15:19:40 2022

@author: Youngbin Moon (y.moon@unibas.ch)
"""

import pysam
import argparse
import numpy as np
import matplotlib.pyplot as plt
import csv

"""
Aim1: plot a histograms of length of softclipped region.
(before fixing alignment of soft clipped region)


Aim2: save as a csv, % of reads with length of softclipped region between 1 and 3.
"""

def write_out(row_data, o_file):
    """
    Parameters
    ----------
    row_data : list
        A list that contains sample, total # of softclipped reads, # of reads with slength of softclipped region between 1 and 3., and its percentage
    
    o_file : string
        output csv file name    

    Returns
    -------
    returns nothing but saves the output in csv format.
    """
    
    w = csv.writer(open(o_file, "w"))
    w.writerow(row_data)
    
def get_read_properties(each_read):
    """
    Parameters
    ----------
    each_read : pysam object
        current read of our interest.
            
    Returns
    -------
    rev : bool
        whether a read is mapped to + or - strand of a genome.
    
    chromosome : string
        chromosome that a read maps to
    
    refStart : int
        start position of the mapped part of a read
    
    refEnd : int
        end position of the mapped part of a read
    
    left_end : tuple
        contains whether a read has a soft clipped region and if yes, how long?
        This is for a read mapping to (-) strand of the genome.
        
        if left_end[0] == 4 -> there is a soft clipped region in the left side of a read.
        left_end[1] -> gives you length of the soft clipped region in the left side of a read.
        
    right_end : tuple
        contains whether a read has a soft clipped region and if yes, how long?
        This is for a read mapping to (+) strand of the genome.   
        
        if right_end[0] == 4 -> there is a soft clipped region in the right side of a read.
        right_end[1] -> gives you length of the soft clipped region in the right side of a read.
    """      
    rev = each_read.is_reverse
    chromosome = each_read.reference_name
    alignedRefPositions = each_read.get_reference_positions()
    refStart = alignedRefPositions[0]
    refEnd = alignedRefPositions[-1]
    tuples = each_read.cigartuples
    left_end = tuples[0]
    right_end = tuples[-1]
    
    return rev, chromosome, refStart, refEnd, left_end, right_end

def get_softclipped(read, use_fc):
    """
    Parameters
    ----------
    read : pysam object
        current read of our interest.
        
    use_Fc : bool
        whether use fixed cleavage site or not.
    
    Returns
    -------
    size_softclip : int
        get length of a softclipped region of a given read. 
    """      
    rev, chromosome, refStart, refEnd, left_end, right_end = get_read_properties(read)    

    # if mapping to reverse direction and have a soft-clipped region    
    if rev == True and left_end[0] == 4:
        
        if not use_fc:
            size_softclip = left_end[1]
        
        elif use_fc:
            OCS = int(read.get_tag('XO'))
            FCS = int(read.get_tag('XF'))
            difference = OCS - FCS
              
            size_softclip = left_end[1] - difference
            
        if size_softclip == 0:
            return "no_softclipped"
        
        elif size_softclip > 0:
            return size_softclip
        
    # if mapping to foward direction and have a soft-clipped region
    elif rev == False and right_end[0] == 4:
        
        if not use_fc:
            size_softclip = right_end[1]

        elif use_fc:
            OCS = int(read.get_tag('XO'))
            FCS = int(read.get_tag('XF'))
            difference = FCS - OCS
                    
            size_softclip = right_end[1] - difference
            
        if size_softclip == 0:
            return "no_softclipped"
        
        elif size_softclip > 0:
            return size_softclip
        
    # in this case, it does not have a soft-clipped region. should not consider these reads
    else:
        return "no_softclipped"

def get_softclipped_distribution(sam, use_Fc):
    """
    Parameters
    ----------
    sam : bam file
        deduplicated bam file in a specific sample.
        
    use_Fc : bool
        whether use fixed cleavage site or not.
    
    Returns
    -------
    soft_clipped_list : list
        a list of length of soft clipped region for all reads.
        (It can be either before or after fixing alignment of soft clipped region.)
    """        
    soft_clipped_list = []
    for Read in sam.fetch():
        size_softclip = get_softclipped(Read, use_Fc)
        
        if size_softclip != "no_softclipped":
            soft_clipped_list.append(size_softclip)
        
        else:
            continue
        
    return soft_clipped_list

def get_mean(scores_list):
    """
    Parameters
    ----------
    scores_list : list
        a list of length of soft clipped region for all reads. 
        (It can be either before or after fixing alignment of soft clipped region.)
    
    Returns
    -------
    mean_val : float
        mean length of the soft-clipped regions.
    """       
    mean_val = np.mean(scores_list) 
    return mean_val

def plot_histogram_alone(file_template, distribution_unfixed):
    """
    Parameters
    ----------
    file_template : string
        output file template
        
    distribution_unfixed : list
        a list of length of soft clipped regions for all reads (before fixing alignment of soft clipped region). 
    
    Returns
    -------
    Returns nothing but plots a histograms of length of softclipped region.
    (before fixing alignment of soft clipped region)
    """       
    plt.figure()
    plt.xticks(fontsize='x-large')
    plt.yticks(fontsize='x-large')
    plt.xlabel('size of softclipped region', fontsize = 'x-large')
    plt.ylabel('frequency', fontsize = 'x-large')
            
    bins = np.arange(min(distribution_unfixed), max(distribution_unfixed), 1)
                     
    plt.hist(distribution_unfixed, bins = bins, histtype='bar', color ='b', ec = 'blue', fc = 'None') 
                              
    plt.yscale("log") 
    
    plt.rcParams['font.family'] = "DejaVu Sans"
    
    file1 = file_template + ".png"
    file2 = file_template + ".svg"
    
    plt.savefig(file1, bbox_inches='tight')
    plt.savefig(file2, bbox_inches='tight')
    
def plot_histogram(file, distribution_fixed, distribution_unfixed):
    """
    Parameters
    ----------
    file : string 
        output histogram file name
    
    distribution_fixed : list
        a list of length of soft clipped regions for all reads (after fixing alignment of soft clipped region). 
        
    distribution_unfixed : list
        a list of length of soft clipped regions for all reads (before fixing alignment of soft clipped region). 
    
    Returns
    -------
    returns nothing but plots 2 histograms of length of softclipped region.
    1 for before fixing alignment of soft clipped region.
    1 for after fixing alignment of soft clipped region.
    """       
    plt.figure()
    plt.xticks(fontsize='x-large')
    plt.yticks(fontsize='x-large')
    plt.xlabel('size of softclipped region', fontsize = 'x-large')
    plt.ylabel('frequency', fontsize = 'x-large')
    
    plt.rcParams['font.family'] = "Arial"
        
    bins = np.arange(min(min(distribution_fixed), min(distribution_unfixed)), max(max(distribution_fixed), max(distribution_unfixed)), 1)
                     
    plt.hist(distribution_unfixed, bins = bins, histtype='bar', color ='b', ec = 'blue', fc = 'None') 
    plt.hist(distribution_fixed, bins = bins, histtype='bar', color = 'gold', ec= 'gold', fc = 'None')
                              
    plt.legend(loc="upper right")
    plt.yscale("log") 
    
    plt.savefig(file, bbox_inches='tight')
    
def get_inputs():
    parser = argparse.ArgumentParser(description="get a distribution of softclipped region")
    parser.add_argument('--dedup_bam', dest = 'dedup_bam',
                        required = True,
                        help = 'dedup bam input')
    
    parser.add_argument('--file_name', dest = 'file_name',
                        required = True,
                        help = 'output file_name')    

    args = parser.parse_args()
    
    bamFile = args.dedup_bam
    file_name = args.file_name
    
    samFile = pysam.AlignmentFile(bamFile, "rb")
    
    return samFile, file_name

def run_process():
    samFile, file_name =  get_inputs()
    print('successfully got inputs')
    
    # soft_clipped_fixed_list = get_softclipped_distribution(samFile, True)
    # print('successfully got a distribution of soft clipped region with a fixed cleavage site')
    # print("soft_clipped_fixed_list: " + str(soft_clipped_fixed_list))
    
    soft_clipped_unfixed_list = get_softclipped_distribution(samFile, False)
    print('successfully got a distribution of soft clipped region with an unfixed cleavage site')
    print("soft_clipped_unfixed_list: " + str(soft_clipped_unfixed_list))
    
    # assert(len(soft_clipped_fixed_list) == len(soft_clipped_unfixed_list))
    
    # plot_histogram(file_name, soft_clipped_fixed_list, soft_clipped_unfixed_list)
    # print('successfully got plot histogram')
    
    plot_histogram_alone(file_name, soft_clipped_unfixed_list)
    
    total = len(soft_clipped_unfixed_list)
    # collect the number of reads with length of softclipped region between 1 and 3.
    subset = [1 if elem >= 1 and elem <= 3 else 0 for elem in soft_clipped_unfixed_list]
    len_subset = sum(subset)
    
    percentage = (len_subset * 100)/total
    
    sample_name = '_'.join(file_name.split('_')[0:3])
    
    row = [sample_name, total, len_subset, percentage]
    
    file_name2 = sample_name + "_subset_percentage.csv"
    write_out(row, file_name2)
    print('successfully got percentage of reads which have softclipped length between 1 and 3')
    
if __name__ == "__main__":
    run_process()
    print('success')
    