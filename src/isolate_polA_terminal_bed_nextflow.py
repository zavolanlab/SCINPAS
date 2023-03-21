
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  5 14:10:04 2022

@author: Youngbin Moon (y.moon@unibas.ch)
"""

import matplotlib.pyplot as plt
import argparse
import numpy as np
from collections import Counter
import pandas as pd





"""
Aim : separate a bed file containing all pA sites (from sample) mapping to terminal exons
into annotated and unannotated class according to the distance threshold (default: 100bp)
"""
def write_to_bed(final_df, out_file):
    """
    Parameters
    ----------
    final_df : dataframe
        dataframe containing pA sites of either annotated or unannotated class in a given sample.
    
    out_file : string
        output file name. This output will be in the bed format.

    Returns
    -------
    returns nothing but writes the output in the bed format.
    """    
    final_df.to_csv(out_file, sep = '\t', header = False, index = False,
    columns = ['seqid', 'start', 'end', 'id', 'score', 'strand'])
    
def get_distance(cleavage_site, current_terminal_exon, rev):
    """
    Parameters
    ----------
    cleavage_site : int
        representative cleavage site of a polyA + T cluster.
    
    current_terminal_exon : dataframe
        a current terminal exon in which a current pA site is thought to map to.
        It is 1 row (subset) of the terminal_exons.bed file.
    
    rev : string
        direction of a strand in which a current pA site maps to.
        
    Returns
    -------
    distance : int
        An 'absolute' distance between representative cleavage site of a pA site and end of a terminal exon.
    
    """    
    # need to consider the direction of terminal exon as well.
    if rev == '-':
        terminal_exon_end = int(current_terminal_exon['start'])
    
    elif rev == '+':
        terminal_exon_end = int(current_terminal_exon['end'])
                
    distance = abs(terminal_exon_end - cleavage_site)
    
    return distance

def separate_cluster(bed_dir, terminal_exons, threshold):
    """
    Parameters
    ----------
    bed_dir : string
        a directory towards bed file that contains all pA sites (from sample) mapping to terminal exons.
        
    threshold : int
        distance threshold for splitting pA clusters (pA sites) mapping to terminal exon into 2 categories
        e.g. 100bp by default
        
    Returns
    -------        
    unannotated_df : dataframe
        dataframe in bed format containing all pA sites from sample that map to terminal exon with distance >= threshold.

    annotated_df : dataframe
        dataframe in bed format containing all pA sites from sample that map to terminal exon with distance < threshold.        
    """           
    unannotated_bed = []
    annotated_bed = []
    with open(bed_dir, "r") as t:
        for line in t:
            # use cluster span to find which terminal exons that a cluster intersects with.
            chrom = line.split()[0]
            cluster_start = int(line.split()[1])
            cluster_end = int(line.split()[2])
            cluster_id = line.split()[3]
            score = line.split()[4]
            direction = line.split()[5]
            
            cs = int(cluster_id.split(':')[1])
            
            # subset terminal exon bed file (gtf) with overlapping span.
            filtered_terminal_exons = terminal_exons[(terminal_exons['seqid'] == chrom) & (terminal_exons['start'] <= cluster_end)\
                                                     & (terminal_exons['end'] >= cluster_start) & (terminal_exons['strand'] == direction)]
            print(filtered_terminal_exons)
            distances = []
            # compute distance between terminal exon and current representative cleavage site
            # for each terminal exon
            # this 'for' phrase is needed because a read can be mapped to more than 1 terminal exons of the same gene.
            # this can happen even if we only keep reads that are primarily aligned.
            # because of how 3'end sequencing is designed, terminal exon that gives minimal distance to a representative cleavage site
            # is the true terminal exon that the representative read actually maps to.
            for index, terminal_exon in filtered_terminal_exons.iterrows():
                distance = get_distance(cs, terminal_exon, direction)
                distances.append(distance)

            min_d = min(distances)
            
            # check minimum distance with threshold.
            if min_d >= threshold:
                # have to add a read once and only once
                unannotated_bed.append((chrom, cluster_start, cluster_end, cluster_id, score, direction))
            
            elif min_d < threshold:
                annotated_bed.append((chrom, cluster_start, cluster_end, cluster_id, score, direction))

    unannotated_df = pd.DataFrame(unannotated_bed, columns = ['seqid', 'start', 'end', 'id', 'score', 'strand'])
    unannotated_df.sort_values(by=['seqid', 'id', 'start', 'end'], inplace = True)
    
    annotated_df = pd.DataFrame(annotated_bed, columns = ['seqid', 'start', 'end', 'id', 'score', 'strand'])
    annotated_df.sort_values(by=['seqid', 'id', 'start', 'end'], inplace = True)
    
    return unannotated_df, annotated_df

def get_args():        
    parser = argparse.ArgumentParser(description="searching for motives upstream of cleavage sites")

    parser.add_argument('--bed_input', dest = 'bed_input',
                        required = True,
                        help = 'bed file of a specific sample containing clusters that map to polyA + Terminal exon region')
        
    parser.add_argument('--terminal_exons_gtf', dest = 'terminal_exons_gtf',
                        required = True,
                        help = 'directory towards terminal_exons bed file from gtf file')  

    parser.add_argument('--annotated_out', dest = 'annotated_out',
                        required = True,
                        help = 'directory towards annotated class output bed file')

    parser.add_argument('--unannotated_out', dest = 'unannotated_out',
                        required = True,
                        help = 'directory towards unannotated class output bed file')  

    parser.add_argument('--distance_threshold', type = int, dest = 'distance_threshold',
                        required = True,
                        help = 'distance threshold for splitting polyA reads mapping to terminal exon into 2 categories')     
      
    args = parser.parse_args()
    
    bed_input_dir = args.bed_input
    terminal_exons_gtf_dir = args.terminal_exons_gtf

    terminal_exons_gtf = pd.read_csv(terminal_exons_gtf_dir, delimiter = '\t', header = None)
    terminal_exons_gtf.columns = ['seqid', 'start', 'end', 'id', 'score', 'strand']
    
    annotated_out = args.annotated_out
    unannotated_out = args.unannotated_out
    d_threshold = args.distance_threshold

    
    return bed_input_dir, terminal_exons_gtf, annotated_out, unannotated_out, d_threshold

def run_process():

    bed_input_dir, terminal_exons_gtf, annotated_out, unannotated_out, d_threshold = get_args()
    print('successfully got inputs')
    
    unannotated_df, annotated_df = separate_cluster(bed_input_dir, terminal_exons_gtf, d_threshold)
    print('successfully separated unannoated df and annotated df')
    
    write_to_bed(unannotated_df, unannotated_out)
    print('successfully saved the unannotated bed file')    

    write_to_bed(annotated_df, annotated_out)
    print('successfully saved the annotated bed file')   
    
    print('length of unannotated_df: ' + str(len(unannotated_df)))
    print('length of annotated_df: ' + str(len(annotated_df)))
    
    bed = pd.read_csv(bed_input_dir, delimiter = '\t', header = None)
    bed.columns = ['seqid', 'start', 'end', 'id', 'score', 'strand']

    print('length of input bed: ' + str(len(bed)))
     
if __name__ == "__main__":
    run_process()
    print("success")