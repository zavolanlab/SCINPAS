# -*- coding: utf-8 -*-
"""

Created on Thu Jan 20 13:15:22 2022

@author: Youngbin Moon (y.moon@unibas.ch)
"""

import pysam
import numpy as np
import matplotlib.pyplot as plt
import argparse





"""
Aim : draw A/T/G/C frequency plot for a specific sample and a specific class of reads. (considers all chromosomes.)
class = annotated, unannotated, intronic, intergenic or exonic
"""
def draw_lineGraphs(As, Ts, Gs, Cs, file_template, window):
    """
    Parameters
    ----------
    As : list
        A list which contains frequency of "A" nucleotide from -99 to 100bp of a cleavage site.
    
    Ts : list
        A list which contains frequency of "T" nucleotide from -99 to 100bp of a cleavage site.

    Gs : list
        A list which contains frequency of "G" nucleotide from -99 to 100bp of a cleavage site.

    Cs : list
        A list which contains frequency of "C" nucleotide from -99 to 100bp of a cleavage site.
    
    file_template : string
        output file template
        
    Returns
    -------
    returns nothing but plots of A/T/G/C frequency plot of a specific class in a specific sample.
    """         
    plt.figure()
    x = np.arange(-(window/2) + 1, (window/2) + 1, 1)
    plt.plot(x, As, label = 'A frequency')
    plt.plot(x, Ts, label = 'T frequency')
    plt.plot(x, Gs, label = 'G frequency')
    plt.plot(x, Cs, label = 'C frequency')
        
    plt.legend(loc = 'upper right')
     
    plt.xticks(fontsize='x-large')
    plt.yticks(fontsize='x-large')
    plt.xlabel('distance', fontsize = 'x-large')
    plt.ylabel('relative frequency', fontsize = 'x-large')
        
    plt.rcParams['font.family'] = "DejaVu Sans"
    
    file1 = file_template + ".png"
    file2 = file_template + ".svg"
    
    plt.savefig(file1, bbox_inches='tight')
    plt.savefig(file2, bbox_inches='tight')
    
def convert_to_list(matrix):
    """
    Parameters
    ----------
    matrix : matrix
        A matrix that contains all frequency of a specific nucelotide from -99 to 100bp of a cleavage site.
        This matrix is for all chromosomes in a given sample and a given class.
        1st row : frequency of "A" from -99 to 100bp of a cleavage site.
        2nd row : frequency of "T" from -99 to 100bp of a cleavage site.
        3rd row : frequency of "G" from -99 to 100bp of a cleavage site.
        4th row : frequency of "C" from -99 to 100bp of a cleavage site.  
        
    Returns
    -------
    A_list : list
        A list which contains frequency of "A" nucleotide from -99 to 100bp of a cleavage site.
    
    T_list : list
        A list which contains frequency of "T" nucleotide from -99 to 100bp of a cleavage site.

    G_list : list
        A list which contains frequency of "G" nucleotide from -99 to 100bp of a cleavage site.

    C_list : list
        A list which contains frequency of "C" nucleotide from -99 to 100bp of a cleavage site.
    """        
    A_list = []
    T_list = []
    G_list = []
    C_list = []
    # shape[1] is column
    for i in range(matrix.shape[1]):
        total_freq = sum(matrix[:,i])
        A_list.append(matrix[0,i]/total_freq)
        T_list.append(matrix[1,i]/total_freq)
        G_list.append(matrix[2,i]/total_freq)
        C_list.append(matrix[3,i]/total_freq)
    return A_list, T_list, G_list, C_list

def update_frequency_matrix(sequence, winDow):
    """
    Parameters
    ----------
    sequence : string
        a genomic sequence from -99 to 100bp of a single true cleavage site.
        
    Returns
    -------
    matrix : matrix
        A matrix that contains all frequency of a specific nucelotide from -99 to 100bp of a cleavage site.
        This matrix is for a single cleavage site. Hence value is either 1 or 0.
        1st row : frequency of "A" from -99 to 100bp of a cleavage site.
        2nd row : frequency of "T" from -99 to 100bp of a cleavage site.
        3rd row : frequency of "G" from -99 to 100bp of a cleavage site.
        4th row : frequency of "C" from -99 to 100bp of a cleavage site.        
        
    """     
    matrix = np.zeros((4, winDow))
    for i in range(len(sequence)):
        nucleotide = sequence[i]
        if nucleotide == 'A':
            matrix[0][i] += 1
        elif nucleotide == 'T':
            matrix[1][i] += 1
        elif nucleotide == 'G':
            matrix[2][i] += 1
        elif nucleotide =='C':
            matrix[3][i] += 1
            
    return matrix

def make_frequency_matrix(bed_directory, fasta, window):
    """
    Parameters
    ----------    
    bed_directory: bed file
        a bed file of a specific class of a specific sample.
        each row is a polyA site cluster which contains start, end of a cluster
        and also representative cleavage site of a cluster.
    
    fasta : fasta file
        contains the reference genome sequence.
    
             
    Returns
    -------
    frequency_matrix : matrix
        A matrix that contains all frequency of a specific nucelotide from -99 to 100bp of a cleavage site.
        This matrix is for a single chromosome.
        1st row : frequency of "A" from -99 to 100bp of a cleavage site.
        2nd row : frequency of "T" from -99 to 100bp of a cleavage site.
        3rd row : frequency of "G" from -99 to 100bp of a cleavage site.
        4th row : frequency of "C" from -99 to 100bp of a cleavage site.
    """        
    frequency_matrix = np.zeros((4, window))
    # for all clusters in a bed file
    with open(bed_directory, "r") as t:
        for line in t:
            cluster_id = line.split()[3]
            chromosome = cluster_id.split(':')[0]
            # read_end = the most frequent cleavage sites = "true" cleavage site
            read_end = cluster_id.split(':')[1]
            
            start_ind = int(read_end) - (window/2) + 1
            end_ind = int(read_end) + window/2   
            # get genomic sequence
            # end should be + 1 to include the 100th nucleotide
            seq = fasta.fetch(reference = chromosome, start = start_ind, end = end_ind + 1)
            temp_matrix = update_frequency_matrix(seq, window)
            
            copy_matrix = frequency_matrix.copy()
            frequency_matrix = np.add(copy_matrix, temp_matrix)
        
    return frequency_matrix

def get_args():        
    parser = argparse.ArgumentParser(description="ATGC frequency plot")

    parser.add_argument('--bed', dest = 'bed',
                        required = True,
                        help = 'bed file in specific class')
        
    parser.add_argument('--fasta', dest = 'fasta',
                        required = True,
                        help = 'fasta file dir')
    
    parser.add_argument('--out_name', dest = 'out_name',
                        required = True,
                        help = 'ATGC plot')
        
    parser.add_argument('--window_size', dest = 'window_size',
                        required = True,
                        help = 'window size. It has to be even number')
    
    args = parser.parse_args()
    
    bed_dir = args.bed
    fasta_dir = args.fasta
    fasta_file = pysam.FastaFile(fasta_dir)
    out_name = args.out_name
    
    window_size = int(args.window_size)
    
    return bed_dir, fasta_file, out_name, window_size

def run_process():

    bed_dir, fasta_file, out_name, window_size = get_args()
             
    total_frequency_matrix = np.zeros((4, window_size))

    frequency_matrix = make_frequency_matrix(bed_dir, fasta_file, window_size)
    print("successfully got full frquency matrix")
    
    copy_total_matrix = total_frequency_matrix.copy()
    
    total_frequency_matrix = np.add(copy_total_matrix, frequency_matrix)
    print('successfully added matrix')
    
    A_list, T_list, G_list, C_list = convert_to_list(total_frequency_matrix)
    print("successfully got frquency list of A, T, G and C")
        
    draw_lineGraphs(A_list, T_list, G_list, C_list, out_name, window_size)
    print("successfully got frquency line graphs")
    
    print("success")
    
if __name__ == "__main__":
    run_process()
    print("success")

