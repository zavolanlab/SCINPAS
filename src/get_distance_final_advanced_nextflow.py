# -*- coding: utf-8 -*-
"""

Created on Mon Jan 10 10:26:42 2022

@author: Youngbin Moon (y.moon@unibas.ch)
"""
import pandas as pd
import pysam
import csv
import argparse
import math

"""
Aim: to get distance between a pA site (representative cleavage site) and an end of a terminal exon
followed by saving the output in csv file.

Note: if a pA site maps to multiple terminal exons,
choose the terminal exon that has the shortest distance with the pA site
because we use 3'end sequencing.
"""

def write_out_to_csv(data, file_name):
    """
    Parameters
    ----------
    data : a list of tuples
        a list of tuples in which each tuple is (min_distance, min_is_minus, min_log_distance) where:
        min_distance : a distance between a pA site(representative cleavage site) and an end of a terminal exon.
        min_is_minus : direction of a PA site
        min_log_distance : log10(min_distance)
        
        Note: if a pA site maps to multiple terminal exons,
        choose the terminal exon that has the shortest distance with the pA site
        because we use 3'end sequencing.
    
    file_name : string 
        output csv file that contains a distribution of 'absolute' distance between end of a read and end of terminal exon.
    
    Returns
    -------
    returns nothing but it writes output to the csv file.
    """    
    w = csv.writer(open(file_name, 'w'))
    for row in data:
        w.writerow(list(row))
        
def get_distance(cleavage_site, current_terminal_exon, rev):
    """
    Parameters
    ----------
    cleavage_site : int
        representative cleavage site of a polyA + T cluster.
    
    current_terminal_exon : dataframe
        terminal exon in which current read is thought to map to.
        It is 1 row (subset) of the terminal_exons.bed file.
    
    rev : string
        direction of a strand a read maps to.
        
    Returns
    -------
    distance : int
        An 'absolute' distance between a representative cleavage site and an end of a terminal exon.
    
    """    
    # need to consider the direction of terminal exon as well.
    if rev == '-':
        terminal_exon_end = int(current_terminal_exon['start'])
    
    elif rev == '+':
        terminal_exon_end = int(current_terminal_exon['end'])
                
    distance = (terminal_exon_end - cleavage_site)
    
    if distance < 0 :
        is_minus = '-'
    
    elif distance >= 0 :
        is_minus = '+'
    
    return (abs(distance), is_minus)

def get_d_for_polyA_reads(bed_dir, terminal_exons):
    """
    Parameters
    ----------
    bed_dir : string
        a directory towards a bed file containing all pA sites mapping to terminal exons.
    
    terminal_exons : dataframe
        a dataframe containing all terminal exons from gtf.
        
    Returns
    -------
    min_distances_distribution : list
        a list of tuples in which each tuple is (min_distance, min_is_minus, min_log_distance) where:
        min_distance : a distance between a pA site(representative cleavage site) and an end of a terminal exon.
        min_is_minus : direction of a PA site
        min_log_distance : log10(min_distance)
        
        Note: if a pA site maps to multiple terminal exons,
        choose the terminal exon that has the shortest distance with the pA site
        because we use 3'end sequencing.
    """       
    min_distances_distribution = []
    with open(bed_dir, "r") as t:
        for line in t:
            # use cluster span to find which terminal exons that a cluster intersects with.
            chrom = line.split()[0]
            cluster_start = int(line.split()[1])
            cluster_end = int(line.split()[2])
            cluster_id = line.split()[3]
            score = int(line.split()[4])
            direction = line.split()[5]
            
            cs = int(cluster_id.split(':')[1])
            
            # subset terminal exon bed file (gtf) with overlapping span.
#            filtered_terminal_exons = terminal_exons[(terminal_exons['seqid'] == chrom) & (terminal_exons['start'] <= cluster_end)\
#                                                     & (terminal_exons['end'] >= cluster_start) & (terminal_exons['strand'] == direction)].copy()

            filtered_terminal_exons = terminal_exons[(terminal_exons['seqid'] == chrom) & (terminal_exons['start'] <= cluster_end + 1)\
                                                     & (terminal_exons['end'] >= cluster_start - 1) & (terminal_exons['strand'] == direction)].copy()            
            if len(filtered_terminal_exons) >= 2:
                print(filtered_terminal_exons)    
            distances = []
            are_minus = []
            # compute distance between terminal exon and current representative cleavage site
            # for each terminal exon
            # this 'for' phrase is needed because a read can be mapped to more than 1 terminal exons of the same gene (or rarely different genes).
            # this can happen even if we only keep reads that are primarily aligned.
            # because of how 3'end sequencing is designed, terminal exon that gives minimal distance to a representative cleavage site
            # is the true terminal exon that the representative read actually maps to.
            for index, terminal_exon in filtered_terminal_exons.iterrows():
                distance, is_minus = get_distance(cs, terminal_exon, direction)
                distances.append(distance)
                are_minus.append(is_minus)
            
            # This condition is for negative control which is at the deduplicated level.
            if len(distances) > 0:
                # According to definition of 3'end sequencing, the terminal exon that gives the shortest read-terminal exon distance
                # is the true terminal that a read maps to.
                # Hence min(distances) gives the distance between the true terminal exon and the read.
                # This holds for even if reads map to multiple genes. The gene in which it has the shortest read-terminal exon distance
                # is the true gene that a read maps to.
                min_d = min(distances)
                idx = distances.index(min(distances))
                
                min_d_is_minus = are_minus[idx]
                min_distance = distances[idx]
                
                assert(min_d == min_distance)
    
                # min_log_distance always >= 0
                # add a pseudocount so that distance > 0
                min_log_distance = math.log10(1 + min_distance)
                # compute distance at the level of reads rather than at the level of clusters. 
                # Create a list of distances as many as the number of reads that belong to the cluster.
                partial_distances_list = [(min_distance, min_d_is_minus, min_log_distance)] * score
                min_distances_distribution += partial_distances_list
            
    return min_distances_distribution

def get_inputs():
    
    parser = argparse.ArgumentParser(description="get distance between polyA read and its closest terminal exon")
    parser.add_argument('--bed_input', dest = 'bed_input',
                        required = True,
                        help = 'bed file containing all pA clusters mapping to terminal exons from a sample')
    
    parser.add_argument('--gtf_input', dest = 'gtf_input',
                        required = True,
                        help = 'bed file containing all terminal exons from gtf file')   
    
    parser.add_argument('--distance_out', dest = 'distance_out',
                        required = True,
                        help = 'distance distribution output csv file name')

    args = parser.parse_args()
    
    bed_input = args.bed_input
    gtf_dir = args.gtf_input
    distance_out = args.distance_out
         
    gtfFile = pd.read_csv(gtf_dir, delimiter = '\t', header = None)
    gtfFile.columns = ['seqid', 'start', 'end', 'id', 'score', 'strand']
    
    return (bed_input, gtfFile, distance_out)

def run_process():

    (bed_input, gtfFile, distance_out) = get_inputs()
    print('successfully got inputs')
    
    min_distances_distribution = get_d_for_polyA_reads(bed_input, gtfFile)
    print("successfully got distance distribution from read end to end of terminal exon")
    
    write_out_to_csv(min_distances_distribution, distance_out)
    print('successfully wrote distribution to csv')
        
if __name__ == "__main__":
    
    run_process()
    print("success")
    
    