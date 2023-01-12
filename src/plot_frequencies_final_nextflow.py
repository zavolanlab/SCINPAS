# -*- coding: utf-8 -*-
"""

Created on Thu Jan 20 13:15:22 2022

@author: Youngbin Moon (y.moon@unibas.ch)
"""

import pysam
import numpy as np
from collections import Counter
import matplotlib.pyplot as plt
import os
import argparse

"""
Aim : draw A/T/G/C frequency plot for a specific sample and a specific class of reads. (considers all chromosomes.)
class = annotated, unannotated, intronic, intergenic or exonic
"""
def draw_lineGraphs(As, Ts, Gs, Cs, file, total_num_reads):
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
    
    file : string
        file name of A/T/G/C frequency plot of a specific class in a specific sample.
        
    Returns
    -------
    returns nothing but plots of A/T/G/C frequency plot of a specific class in a specific sample.
    """         
    plt.figure()
    x = np.arange(-99, 101, 1)
    plt.plot(x, As, label = 'A frequency')
    plt.plot(x, Ts, label = 'T frequency')
    plt.plot(x, Gs, label = 'G frequency')
    plt.plot(x, Cs, label = 'C frequency')
        
    plt.legend(loc = 'upper right')
     
    plt.xticks(fontsize='x-large')
    plt.yticks(fontsize='x-large')
    plt.xlabel('distance', fontsize = 'x-large')
    plt.ylabel('relative frequency', fontsize = 'x-large')
    
    plot_title = file.split('.')[0] + ': (' + str(total_num_reads) + ')'
    plt.title(plot_title, fontsize = 'x-large')
    
    plt.rcParams['font.family'] = "Arial"
    plt.savefig(file, bbox_inches='tight')
    
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

def update_frequency_matrix(sequence):
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
    matrix = np.zeros((4, 200))
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

def make_frequency_matrix(dictionary, fasta, chromosome):
    """
    Parameters
    ----------    
    dictionary : dictionary of list.
            key : a tuple of chromID and cluster_ID.
            value : a list which contains all "potential" cleavage sites with same chromID and cluster_ID.
            This dictionary has all potential cleavage sites in a specific sample, specific class and specific chromosome.
            
            The most frequent read end (most frequent potential cleavage sites) 
            amongst read ends with same chromID and cluster_ID will be selected as a "true" cleavage site.
    
    fasta : fasta file
        contains the reference genome sequence.
    
    chromosome : string
        current chromosome.            
    
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
    frequency_matrix = np.zeros((4, 200))
    # elem = chromID and cluster_ID
    # for all clusters in a specific chromosome
    for elem in dictionary.keys():
        read_ends = dictionary[elem]
        # read_ends_keys:a list of potential cleavage sites
        # read_ends_counts: a list of counts how many times this potential cleavage site appears
        read_ends_keys = list(Counter(read_ends).keys())
        read_ends_counts = list(Counter(read_ends).values())
        # choose the most frequent cleavage sites as a "true" cleavage sites
        read_end = read_ends_keys[read_ends_counts.index(max(read_ends_counts))]
        start_ind = int(read_end) -99
        end_ind = int(read_end) + 100
        # get genomic sequence
        # end should be + 1 to include the 100th nucleotide
        seq = fasta.fetch(reference = chromosome, start = start_ind, end = end_ind + 1)
        temp_matrix = update_frequency_matrix(seq)
        frequency_matrix = np.add(frequency_matrix, temp_matrix)
    return frequency_matrix

def compute_n_reads(bam):
    """
    Parameters
    ----------    
    bam : bam file
        A bam file of a specific sample, a specific class and a specific chromosome
                
    Returns
    -------
    partial_num_reads : int
        number of a reads in a specific sample, a specific class and a specific chromosome
    """      
    partial_num_reads = 0
    for read in bam.fetch():
        partial_num_reads += 1
    return partial_num_reads

def make_input_slc(in_template, sam, chromosome, distance_param, use_FC):
    """
    Parameters
    ----------    
    in_template : string
        input directory template of the single linkage clustering script. 
        e.g. ../inputSLC_above_
        In this function, you make a full input directory of single linkage clustering.
        e.g. ../inputSLC_above_chr6.txt

    sam : bam file
        A bam file of a specific sample, a specific class and a specific chromosome
        
    chromosome : string
        chromosome. e.g. chr6
    
    distance_param : character. e.g. '4'
        distance parameter of a single linkage clustering.
        
    use_fc : bool
        whether you use fixed cleavage site or not.
        
    Returns
    -------
    output : dataframe
        A dataframe that contains the output of the single linkage clustering.
        This dataframe has a set of clusters of reads generated according to distance parameter.
    """             
    input_dir = in_template + chromosome + ".txt"
    with open(input_dir, "a+") as t:
        for read in sam.fetch():
            alignedRefPositions = read.get_reference_positions()
            rev = read.is_reverse
            chrom = read.reference_name
            if rev == True:
                rev = '-'
            else:
                rev = '+'
            Id = read.query_name
            count = read.mapping_quality
            if not use_FC:
                refStart = alignedRefPositions[0]
                refEnd = alignedRefPositions[-1]
                t.writelines(Id + " " + chrom + " " + rev + " " + str(refStart) + " " + str(refEnd) + " " + str(count) + "\n")
            elif use_FC:
                start = int(read.get_tag('FC')) -1
                end = int(read.get_tag('FC'))
                print('used fixed cleavage site')
                t.writelines(Id + " " + chrom + " " + rev + " " + str(start) + " " + str(end) + " " + str(count) + "\n")
                
        t.close()                    
    print("###### Creating output file of " + str(chromosome) + " ######")
    input_dir = str(input_dir)
    distance_param = str(distance_param)
    
    stream = os.popen(f"./single_linkage {input_dir} {distance_param}")
    output = stream.read()
    
    print("###### Output of " + chromosome + "######")
    print("###### successfully done C++ script " + "######")
    return output

def write_temp (tempfile, out_put):
    """
    Parameters
    ----------    
    tempfile : string
        output directory of the single linkage clustering script.
        e.g. ../slc_output_above_6.txt

    out_put : dataframe
        A dataframe that contains the output of the single linkage clustering.
        This dataframe has a set of clusters of reads generated according to distance parameter.
        
    Returns
    -------
    returns nothing but writes the out_put in tempfile.
    """         
    with open(tempfile, "a+") as t:
        t.writelines(str(out_put))           

def make_cluster(tempfile, use_FC):
    """
    Parameters
    ----------    
    tempfile : string
        output directory of the single linkage clustering script.
        This file contains a dataframe which has a set of clusters of reads 
        generated according to distance parameter.
        e.g. ../slc_output_above_6.txt.

    use_FC : bool
        whether you use fixed cleavage site or not.
        
    Returns
    -------
    clusters_reads_ends : dictionary of list.
            key : a tuple of chromID and cluster_ID.
            value : a list which contains all "potential" cleavage sites with same chromID and cluster_ID.
            In other words, we make clusters of potential cleavage sites. 
            This dictionary has all potential cleavage sites in a specific sample, specific class and specific chromosome.
            
            In the "make_frequency_matrix" function,
            the most frequent read end (most frequent potential cleavage sites) 
            will be selected as a "true" cleavage site.
    """     
    clusters_read_ends = {}
    with open(tempfile, "r") as t:
        for line in t:
            # clusterID is always encountered first within each cluster
            if len(line.split())==4:
                cluster_ID = line.split()[1]
            if len(line.split())==6:
                chromID = line.split()[1]
                direction = line.split()[2]
                read_start = line.split()[3]
                read_end = line.split()[4]
                
                if not use_FC:
                    if direction == '-':
                        cleavage_site = read_start
                    elif direction == '+':
                        cleavage_site = read_end
                # fixed cleavage site already considers the direction
                elif use_FC:
                    cleavage_site = read_end
                    
                # add read end to the appropriate cluster
                if (chromID, cluster_ID) not in clusters_read_ends.keys():
                    clusters_read_ends[(chromID, cluster_ID)] = []
                clusters_read_ends[(chromID,cluster_ID)].append(cleavage_site)
    
    return clusters_read_ends

def generate_input(number, type_of_reads, bam_template, in_slc_dir):
    """
    Parameters
    ----------
    number : character
        chromosome number (just number). e.g. '6'
    
    type_of_reads : string
        class of reads. it can be: annotated, unannotated, intronic, intergenic, exonic    
    
    bam_template : string
        input bam file template. e.g. 10X_P7_15_intergenic_reads_sorted.
        In this function, it will be converted into a full directory towards bam file of 
        a specific sample, specific class and specific chromosome.
        e.g. 10X_P7_15_intergenic_reads_sorted_6.bam
            
    in_slc_dir : string
        directory toward a folder that contains input and output of single linkage clustering script (C++ script).
        
    Returns
    -------
    sam : bam file
        A bam file of a specific sample, a specific class and a specific chromosome
    
    in_templ : string
        input directory template of the single linkage clustering script. 
        e.g. ../inputSLC_above_
    
    chrom : string
        chromosome. e.g. chr6
    
    temp_file : string
        output directory of the single linkage clustering script.
        e.g. ../slc_output_above_6.txt
    """     
    chrom = 'chr' + str(number)
    # annotated
    if type_of_reads == 4:
        bamFile = bam_template + str(number) + "_sorted.bam"   
        in_templ = in_slc_dir + "/inputSLC_below_"
        temp_dir = in_slc_dir + "/slc_output_below_"
    # unannotated
    elif type_of_reads == 5:
        bamFile = bam_template + str(number) + "_sorted.bam"
        in_templ = in_slc_dir + "/inputSLC_above_"
        temp_dir = in_slc_dir + "/slc_output_above_"
    # nonpolyA
    elif type_of_reads == 6:
        bamFile = bam_template + "_sorted_" + str(number) + ".bam"
        in_templ = in_slc_dir + "/inputSLC_nonpolyA_"
        temp_dir = in_slc_dir + "/slc_output_nonpolyA_"     
    # intronic
    elif type_of_reads == 7:
        bamFile = bam_template + "_sorted_" + str(number) + ".bam"
        in_templ = in_slc_dir + "/inputSLC_intronic_"
        temp_dir = in_slc_dir + "/slc_output_intronic_" 
    # intergenic
    elif type_of_reads == 8:
        bamFile = bam_template + "_sorted_" + str(number) + ".bam"
        in_templ = in_slc_dir + "/inputSLC_intergenic_"
        temp_dir = in_slc_dir + "/slc_output_intergenic_" 
    # internal priming
    elif type_of_reads == 9:
        bamFile = bam_template + "_sorted_" + str(number) + ".bam"
        in_templ = in_slc_dir + "/inputSLC_internal_priming_"
        temp_dir = in_slc_dir + "/slc_output_internal_priming_"
    # exonic polyA
    elif type_of_reads == 10:
        bamFile = bam_template + "_sorted_" + str(number) + ".bam"
        in_templ = in_slc_dir + "/inputSLC_exonicA_"
        temp_dir = in_slc_dir + "/slc_output_exonicA_"

    sam = pysam.AlignmentFile(bamFile, "rb")    
    temp_file = temp_dir + str(number) + ".txt"
    return sam, in_templ, chrom, temp_file

def get_args():        
    parser = argparse.ArgumentParser(description="ATGC frequency plot")

    parser.add_argument('--bam_template', dest = 'bam_template',
                        required = True,
                        help = 'specific class of reads')
        
    parser.add_argument('--fasta', dest = 'fasta',
                        required = True,
                        help = 'fasta file dir')
    
    parser.add_argument('--slc_distance', type = int, dest = 'slc_distance',
                        required = True,
                        help = 'single linkage cluster dsitance parameter')

    parser.add_argument('--out_name', dest = 'out_name',
                        required = True,
                        help = 'ATGC plot')
    
    parser.add_argument('--annotated', type = int, dest = 'annotated',
                        required = True,
                        help = 'which reads did you recieve?')    
    
    parser.add_argument('--in_slc', dest = 'in_slc',
                        required = True,
                        help = 'directory towards in_slc folder')   
    
    parser.add_argument('--use_fc', dest = 'use_fc',
                        required = True,
                        help = 'whether use fixed cleavage site or not')

    parser.add_argument('--species', dest = 'species',
                        required = True,
                        help = 'species of the data')
    
    args = parser.parse_args()
    
    bamTemplate = args.bam_template
    fasta = args.fasta
    slc_distance = args.slc_distance
    out_name = args.out_name
    annotated = args.annotated 
    in_slc_dir = args.in_slc
    use_fc = bool(int(args.use_fc))
    species = args.species
    
    return bamTemplate, fasta, slc_distance, out_name, annotated, in_slc_dir, use_fc, species

def run_process():

    bamTemplate, fasta_dir, slc_distance, out_name, annotated, in_slc_dir, use_fc, species = get_args()
     
    if species == 'mouse':
        list_of_chromosomes = list(range(1, 20))
        list_of_chromosomes += ['X', 'Y']
        
    elif species == 'human':
        list_of_chromosomes = list(range(1, 23))
        list_of_chromosomes += ['X', 'Y']
        
    fasta_file = pysam.FastaFile(fasta_dir)
    total_frequency_matrix = np.zeros((4, 200))
    n_reads = 0
    
    for elem in list_of_chromosomes:
        # annotated = 4 means polyA reads that have distance(read, closest terminal exon) < 200
        # unannotated = 5 means polyA reads that have distance(read, closest terminal exon) >= 200
        
        sam, in_templ, chrom, temp_file = generate_input(elem, annotated, bamTemplate, in_slc_dir)
        print("successfully initialized input for chr" + str(elem))
            
        slc_output = make_input_slc(in_templ, sam, chrom, slc_distance, use_fc)
        print("successfully generated output of slc for chr" + str(elem))
        
        write_temp(temp_file, slc_output)
        print("successful saved slc_output for chr" + str(elem))
        
        clusters_read_ends = make_cluster(temp_file, use_fc)
        print("successful made dictionary for chr" + str(elem))

        frequency_matrix = make_frequency_matrix(clusters_read_ends, fasta_file, chrom)
        print("successfully got frquency matrix for chr" + str(elem))
        
        total_frequency_matrix = np.add(total_frequency_matrix, frequency_matrix)
        print('successfully added matrix')
        
        partial_n_reads = compute_n_reads(sam)
        n_reads += partial_n_reads
                              
    A_list, T_list, G_list, C_list = convert_to_list(total_frequency_matrix)
    print("successfully got frquency list of A, T, G and C")
        
    draw_lineGraphs(A_list, T_list, G_list, C_list, out_name, n_reads)
    print("successfully got frquency line graphs")
        
if __name__ == "__main__":
    run_process()
    print("success")

