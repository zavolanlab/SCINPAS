# -*- coding: utf-8 -*-
"""

Created on Sun Jan  9 07:02:13 2022

@author: Youngbin Moon (y.moon@unibas.ch)
"""

import pysam
import argparse
import statistics
import csv




"""
Aim: get an average A nucleotide frequency occuring in non-softclipped reads.
"""

def write_out_csv(row_data, o_file):
    """
    Parameters
    ----------
    row_data : list
        A list that contains sample, bam type and count. This will be 1 row in the csv.
        Later, this will be merged across all samples and all bam types.
    
    o_file : string
        output csv file name    

    Returns
    -------
    returns nothing but saves the output (sample, bam type and count) in csv format.
    """
    
    w = csv.writer(open(o_file, "w"))
    w.writerow(row_data)
    
def get_mean_value(a_freq_list):
    """
    Parameters
    ----------
    a_frequencies : list
        a list that contains average frequency of "A" appearing in a read.
    
    Returns
    -------        
    average_A_freq : float
        an average A nucleotide frequency occuring in non-softclipped reads.
    """              
    average_A_freq = statistics.mean(a_freq_list)
    
    return average_A_freq

def count_As(sequence):
    """
    Parameters
    ----------
    sequence : string
        a whole mapped sequence
    
    Returns
    -------        
    number_A : int
        The number of "A"s in the sequence.
    """       
    list_count = [1 if elem == 'A' else 0 for elem in sequence]
    number_A = sum(list_count)
    return number_A

def get_A_frequency(bam):
    """
    Parameters
    ----------
    bam : bam file
        a bam file containing reads that do not contain softclipped region.
    
    Returns
    -------        
    a_frequencies : list
        a list that contains average frequency of "A" appearing in a read.
    """          
    a_frequencies = []
    for read in bam.fetch():
      full_sequence = read.get_forward_sequence()
      num_A = count_As(full_sequence)
      a_freq = num_A/len(full_sequence)
      a_frequencies.append(a_freq)
    
    return a_frequencies 

def get_dedup_input():
    
    parser = argparse.ArgumentParser(description="get average frequency of nucleotide A appearing in non-softclipped reads")
    parser.add_argument('--bam_input', dest = 'bam_input',
                        required = True,
                        help = 'bam file containing non-softlcipped reads')
    
    parser.add_argument('--o_file', dest = 'o_file',
                        required = True,
                        help = 'output csv file')

    args = parser.parse_args()
    
    bamFile = args.bam_input
    sam = pysam.AlignmentFile(bamFile, "rb")
        
    o_file = args.o_file
    
    return sam, o_file

def run_process():
    
    sam, o_file = get_dedup_input()
    print('successfully got inputs')
    
    a_frequencies = get_A_frequency(sam)
    print('successfully got average frequency of A appearing for each read')
    
    average_A_freq = get_mean_value(a_frequencies)
    print('successfully got an average A nucleotide frequency occuring in non-softclipped reads.')
    
    sample_name = o_file.split('_')[0:3]
    if "UmiDedup" in sample_name:
        sample_name = "negative_control"    
    
    data = [sample_name, average_A_freq]
    
    write_out_csv(data, o_file)
    print('successfully saved the result')
    
if __name__ == "__main__":    
    run_process()
    print('success')