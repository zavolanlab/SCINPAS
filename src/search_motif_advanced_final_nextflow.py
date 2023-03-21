
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 21 14:50:21 2022

@author: Youngbin Moon (y.moon@unibas.ch)
"""

import pysam
import matplotlib.pyplot as plt
import argparse
import csv
import numpy as np
from collections import Counter
import pandas as pd
import statistics
"""
Aim 1 : Using all pA clusteres in a particular sample and class, draw all 18 motif frequency plots.
Aim 2 : Using all pA clusteres in a particular sample and class, draw an overlaid motif frqeuncy plot which contains all 18 motif frequency plots.
Aim 3 : Using all pA clusteres in a particular sample and class, compute a motif score.(1 score per 1 sample)
Aim 4 : Using all pA clusteres in a particular sample and class, save order of motives 
        (from highest frequency to lowest frequency) in the csv.

Aim 5 : Using all reads in a particular sample and class, save exact peak positions for each motif. 
        (peak positions and score are saved together)
"""
def modify_class_name(class_type):
    """
    Parameters
    ----------        
    class_type : string
        class of reads. (either annotated(4), unannotated(5), intronic(7), intergenic(8) or exonic(10))
           
    Returns
    -------
    modified_class_name : string
        modified class name of reads in order to use same terminology as paper. 
        (either ATE, UTE, I, IG, NTE)
    """  
    if class_type == 4:
        modified_class_name = "ATE"

    elif class_type == 5:
        modified_class_name = "UTE"
 
    elif class_type == 7:
        modified_class_name = "I"

    elif class_type == 8:
        modified_class_name = "IG"
        
    elif class_type == 10:
        modified_class_name = "NTE"        
        
    return modified_class_name

def write_out(row_data, o_file):
    """
    Parameters
    ----------
    row_data : list
        it can be either :
        1) a list of motif orders (from highest frequency to lowest)
        2) [sample name, score]
            
    o_file : string
        output file name    
       
    Returns
    -------        
    returns nothing but writes row_data at the csv file with output name: o_file.
    """        
    w = csv.writer(open(o_file, "w"))
    w.writerow(row_data)
        
def get_smoothen_list(freq_list):
    """
    Parameters
    ----------
    freq_list : list
       frequency of occurence of a motif of interest at each position.
       This frequency is a float and computed by considering all representative cleavage sites involved.
       
    Returns
    -------        
    smoothen_list : list
        it is a smoothen frequency_list where each element is averaged by its nearby elements.        
    """        
    smoothen_list = []
    for i in range(len(freq_list)):
        
        if i == 0:
            mean_value = statistics.mean(freq_list[0 : 6])
        
        elif i > 0 and i < 5:
            mean_value = statistics.mean(freq_list[0 : i + 6])
        
        elif i >= 5 and i <= len(freq_list) - 6:
            mean_value = statistics.mean(freq_list[i - 5 : i + 6])
        
        elif i > len(freq_list) - 6 and i < len(freq_list) - 1:
            mean_value = statistics.mean(freq_list[i - 5 : len(freq_list)])
        
        elif i == len(freq_list) - 1:
            mean_value = statistics.mean(freq_list[len(freq_list) - 5 : len(freq_list)])
        
        smoothen_list.append(mean_value)
    
    return smoothen_list
       
def get_score(frequency_dictionary, up_bp, length_motif_interest, down_bp, peak_position):
    """
    Parameters
    ----------
    frequency_dictionary : dictionary
        key = position in bp
        value = frequency. (considering all fixed representative cleavage sites)
        
        e.g. if read_1 has the motif at -40, -30, -20 and if read_2 has motif at -30 and -20bp,
        frequency_dict: -40: 1/3, -30: 5/6, -20: 5/6
                
    up_bp : int (positive integer)
        how many number of base pairs do we go upstream of fixed representative cleavage site?
    
    length_motif_interest : int
        length of the most frequently appearing motif at that moment.
            
    down_bp : int (negative integer)
        how many number of base pairs do we go downstream of fixed representative cleavage site?
    
    peak_position : dictionary of tuple
        key = motif
        value = (lower boundary, upper boundary)      
        
    Returns
    -------
    score_peak : tuple
        1st elemet = score        
        score = 1 if a peak is within the boundary generated by gtf.
        score = 0 if a peak is not within the boundary generated by gtf.
        
        2nd element =  peak
        position of a peak of the motif frequency (in bp).
    """     
    # index is currently -5 ..... -40 or +20 ...... -40
    # you need to reverse it so that it becomes -40 ...... -5 or -40 ...... +20
    existing_up_bp = list(frequency_dictionary.keys())
    
    if down_bp == 0:
        base_pairs_upstream = list(range(-up_bp, -length_motif_interest + 2))
    else:
        base_pairs_upstream = list(range(-up_bp, -down_bp + 1))
    
    frequency_list = [frequency_dictionary[elem] if elem in existing_up_bp else 0 for elem in base_pairs_upstream]   

    # smoothen the graph
    smoothen_list = get_smoothen_list(frequency_list)
    peak = base_pairs_upstream[smoothen_list.index(max(smoothen_list))]
    min_boundary = int(peak_position[0])
    max_boundary = int(peak_position[1])
    
    if peak >= min_boundary and peak <= max_boundary:
        score_peak = (1, peak)
    else:
        score_peak = (0, peak)
    
    return score_peak

def plot_linegraph(frequency_dictionary, motif_of_interest, file, up_bp, length_motif_interest, down_bp, tot_freq):
    """
    Parameters
    ----------
    frequency_dictionary : dictionary
        key = position in bp
        value = frequency. (considering all fixed representative cleavage sites that has the motif)
        
        e.g. if read_1 has the motif at -40, -30, -20 and if read_2 has motif at -30 and -20bp,
        frequency_dict: -40: 1/3, -30: 5/6, -20: 5/6
        
    motif_of_interest : string
        a motif of interest (the most frequently appearing motif at that moment)
        
    file : string
        ouput figure file name template.
        
    up_bp : int (positive integer)
        how many number of base pairs do we go upstream of fixed representative cleavage site?
    
    length_motif_interest : int
        length of the most frequently appearing motif at that moment.
            
    down_bp : int (negative integer)
        how many number of base pairs do we go downstream of fixed representative cleavage site?
    
    tot_freq : int
        the maximum column sum of full_dataframe at that moment.
        i.e. the number of occurrences of the most frequently appearing motif at that moment.       
        
    Returns
    -------        
    returns nothing but plot 1 motif frequency plot of interest.
    """     
    plt.figure()
    full_out_name = file + '_' + motif_of_interest + ".png"

    # index is currently -5 ..... -40 or +20 ...... -40
    # you need to reverse it so that it becomes -40 ...... -5 or -40 ...... +20
    existing_up_bp = list(frequency_dictionary.keys())
    
    if down_bp == 0:
        base_pairs_upstream = list(range(-up_bp, -length_motif_interest + 2))
    else:
        base_pairs_upstream = list(range(-up_bp, -down_bp + 1))
      
    frequency_list = [frequency_dictionary[elem] if elem in existing_up_bp else 0 for elem in base_pairs_upstream]   
    
    # smoothen the graph
    smoothen_list = get_smoothen_list(frequency_list)
    plt.plot(base_pairs_upstream, smoothen_list, label = 'frequency graph for motif: ' + str(motif_of_interest) + 'with total count: ' + str(tot_freq))
    
    plt.legend(bbox_to_anchor=(1.04, 1), loc="upper left")
    plt.xticks(fontsize='x-large')
    plt.yticks(fontsize='x-large')
    plt.xlabel('bp upstream', fontsize = 'x-large')
    plt.ylabel('frequency', fontsize = 'x-large')
    
    plt.rcParams['font.family'] = "Arial"
    plt.savefig(full_out_name, bbox_inches='tight')
    
def plot_linegraph_overlay(frequency_dictionary, motif_of_interest, up_bp, length_motif_interest, down_bp, tot_freq, motif_color, width):
    """
    Parameters
    ----------
    frequency_dictionary : dictionary
        key = position in bp
        value = frequency. (considering all fixed representative cleavage sites that has the motif)
        
        e.g. if read_1 has the motif at -40, -30, -20 and if read_2 has motif at -30 and -20bp,
        frequency_dict: -40: 1/3, -30: 5/6, -20: 5/6
        
    motif_of_interest : string
        a motif of interest (the most frequently appearing motif at that moment)
                    
    up_bp : int (positive integer)
        how many number of base pairs do we go upstream of fixed representative cleavage site?
    
    length_motif_interest : int
        length of the most frequently appearing motif at that moment.
            
    down_bp : int (negative integer)
        how many number of base pairs do we go downstream of fixed representative cleavage site?
    
    tot_freq : int
        the maximum column sum of full_dataframe at that moment.
        i.e. the number of occurrences of the most frequently appearing motif at that moment.       
     
    motif_color_dioct : string
        a color that is assigned to a particular motif.
    
    width : float
        width of the line. gets thicker if you have high frequency (# of PAS).
        
    Returns
    -------        
    returns nothing but you will overlay 1 motif frequency plot of an interest at a time.
    """     
    # index is currently -5 ..... -40 or +20 ...... -40
    # you need to reverse it so that it becomes -40 ...... -5 or -40 ...... +20
    existing_up_bp = list(frequency_dictionary.keys())
    
    if down_bp == 0:
        base_pairs_upstream = list(range(-up_bp, -length_motif_interest + 2))
    else:
        base_pairs_upstream = list(range(-up_bp, -down_bp + 1))
      
    frequency_list = [frequency_dictionary[elem] if elem in existing_up_bp else 0 for elem in base_pairs_upstream]   
    
    # smoothen the graph
    smoothen_list = get_smoothen_list(frequency_list)
    print('smoothen list: ' + str(smoothen_list))
    print('smoothen list: ' + str(sum(smoothen_list)))
    smoothen_total = sum(smoothen_list)
    
    if smoothen_total != 0:
        normalized_smoothen_list = [elem/smoothen_total for elem in smoothen_list]
        plt.plot(base_pairs_upstream, normalized_smoothen_list, label = str(motif_of_interest) + ': ' + str(tot_freq), color = motif_color, linewidth = width)
    
def get_frequency_by_position(selected_row_names, selected_motif_only, up_until, length_motif, fasta_f, down_bp):
    """
    Parameters
    ----------
    selected_row_names : list
        a list of fixed representative cleavage sites which has the motif of interest 
        (the most frequently appearing motif at that moment) within the range specified.
        
    selected_motif_only : string
        a motif of interest (the most frequently appearing motif at that moment)
        
    up_until : int (positive integer)
        how many number of base pairs do we go upstream of fixed representative cleavage site?
    
    length_motif : int
        length of the motif
        
    fasta_f : a fasta flie
        contains the reference genome sequence.
    
    down_bp : int (negative integer)
        how many number of base pairs do we go downstream of fixed representative cleavage site?
         
    Returns
    -------        
    frequency_dict : dictionary
        key = position in bp
        value = frequency. (considering all fixed representative cleavage sites that has the motif)
        
        e.g. if read_1 has the motif at -40, -30, -20 and if read_2 has motif at -30 and -20bp,
        frequency_dict: -40: 1/3, -30: 5/6, -20: 5/6
    """       
    if down_bp == 0:    
        upstream_basepairs = list(range(length_motif - 1, up_until + 1))   
    else:
        upstream_basepairs = list(range(down_bp, up_until + 1))
        
    frequency_dict = Counter()
    
    for row_concatenate in selected_row_names:
        row = row_concatenate.split('_')
        chrom = row[0]
        read_direction = row[1]
        cleavage_site = int(row[2])
        
        num_appeared = 0
        # search ranges from 5 to 40 (i.e. -40 ~ -5) if down_bp == 0
        # search ranges from -20 to 40 (i.e. -40 ~ +20) if down_bp != 0
        for upstream_bp in upstream_basepairs:
            if read_direction == '-':
                start_point = cleavage_site + upstream_bp 
                end_point = cleavage_site + upstream_bp + length_motif -1   
                
            elif read_direction == '+':
                start_point = cleavage_site - upstream_bp
                end_point = cleavage_site - upstream_bp + length_motif - 1              
        
            dna_subsequence = fasta_f.fetch(reference = chrom, start = start_point, end = end_point + 1)        
            # correct subsequence so that it becomes compatible with motif
            rna_corrected_subsequence = change_dna_to_rna(dna_subsequence, read_direction)      
        
            if rna_corrected_subsequence == selected_motif_only:
                # if a motif appeared -30, -25 and -20 then assign count of 1/3 to each of positions
                num_appeared += 1           
            else:
                continue
            
        for upstream_bp in upstream_basepairs:
            if read_direction == '-':
                start_point = cleavage_site + upstream_bp 
                end_point = cleavage_site + upstream_bp + length_motif -1   
                
            elif read_direction == '+':
                start_point = cleavage_site - upstream_bp
                end_point = cleavage_site - upstream_bp + length_motif - 1              
        
            dna_subsequence = fasta_f.fetch(reference = chrom, start = start_point, end = end_point + 1)        
            # correct subsequence so that it becomes compatible with motif
            rna_corrected_subsequence = change_dna_to_rna(dna_subsequence, read_direction)      
        
            if rna_corrected_subsequence == selected_motif_only:
                # upstream_bp ranges from -20 ...... 40
                # but it is actually saved as 20 ....... -40
                
                # upstream_bp ranges from 5 ..... 40
                # or it is actually saved as  -5 ...... -40
                frequency_dict[-upstream_bp] += 1/num_appeared            
            else:
                continue
            
    return frequency_dict
    
def get_colsum_list(data, motives):
    """
    Parameters
    ----------
    data: dataframe
        a full dataframe which considers all chromosomes in a given sample and class.
        each row: fixed representative cleavage site where it is tuple of (chromId, direction, representative_fc)  
        each col: motif
        
        each cell value = 1 if a motif is present in the range defined by that particular fixed representative cleavage site. 
        (does not matter where you have the motif as long as it is within the range)
        
        each cell value = 0 if a motif is not present in the range defined by that particular fixed representative cleavage site. 
    
    motives : list
        a list of motif names.
        
    Returns
    -------
    col_sums_list : list
        a list which contains all column sums of the full dataframe (i.e. data)
    """     
    col_sums_list = []
    for motif in motives:
        column = list(data[motif])
        col_sum = sum(column)
        col_sums_list.append(col_sum)
    
    return col_sums_list

def plot_all_motif_overlaid_plots(df, m_infos, o_name, up, fasta_F, down):
    """
    Parameters
    ----------
    df : dataframe
        a full dataframe which considers all chromosomes in a given sample and class.
        each row: fixed representative cleavage site where it is tuple of (chromId, direction, representative_fc)  
        each col: motif
        
        each cell value = 1 if a motif is present in the range defined by that particular fixed representative cleavage site. 
        (does not matter where you have the motif as long as it is within the range)
        
        each cell value = 0 if a motif is not present in the range defined by that particular fixed representative cleavage site. 
        
    m_infos : a list of tuple
        a list of tuple that contains the motif infomation.
        1st element of a tuple: motif name. e.g. AAUAAA
        2nd element of a tuple: length of a motif. e.g. 6

    o_name : string
        motif frequency plot output file name.
            
    up : int (positive integer)
        how many number of base pairs do we go upstream of fixed representative cleavage site?
    
    fasta_F : a fasta flie
        contains the reference genome sequence.
    
    down : int (negative integer)
        how many number of base pairs do we go downstream of fixed representative cleavage site?    
         
    Returns
    -------
    returns nothing but draws a plot that overlays all 18 motive frequency plots at once.
    """                
    motives_only = [elem[0] for elem in m_infos]
    lengths_of_motives = [elem[1] for elem in m_infos]
    motif_color_dict = {}
    
    # if using canonical motives
    if len(motives_only) == 12:
        # setting colors for motif frequency plots
        colors_manual = ['#00FFFF', '#0000BB', '#DAA520', '#FF00FF', '#9A0EEA', '#9ACD32', '#808000', '#F0E68C', '#40E0D0', '#C0C0C0', '#A52A2A', '#008000']
    
    # if using all motives 
    elif len(motives_only) == 18:
        colors_manual = ['#00FFFF', '#0000BB', '#DAA520', '#FF00FF', '#9A0EEA', '#9ACD32', '#808000', '#F0E68C', '#40E0D0', '#C0C0C0', '#A52A2A', '#008000',\
                         '#01153E', '#BBF90F', '#9A0EEA', '#C79FEF', '#AAA662', '#A9561E']
    
    # assigning colors to motives.
    for i, elem in enumerate(motives_only):
        motif_color_dict[elem] = colors_manual[i]
    
    # width of the line. (gets thicker if you have higher frequency)
    width = 2
    plt.figure()
    while len(motives_only) != 0:
        col_sums = get_colsum_list(df, motives_only)
        idx = col_sums.index(max(col_sums))
        
        max_col_sum = col_sums[idx]
        assert(max(col_sums) == max_col_sum)
        # select a motif with the highest frequency
        max_motif = motives_only[idx]
        length_max_motif = int(lengths_of_motives[idx])
        
        # subset df only if a fixed cleavage site has count of the selected motif > 0
        # and return row_names (i.e fixed_cleavage_site)
        # df.index returns the row name (fixed cleavage site)
        # df[max_motif] > 0 returns you a set of indices 
        # where it returns true if a max_motif column has value > 0 at that row and else otherwise 
        
        # return row name (i.e cleavage sites) where count in max_motif is bigger than 0
        used_fixed_cleavage_sites = list(df.index[df[max_motif] > 0])
        frequency_dict = get_frequency_by_position(used_fixed_cleavage_sites, max_motif, up, length_max_motif, fasta_F, down)
                 
        # overlay all motif plots on a single figure
        plot_linegraph_overlay(frequency_dict, max_motif, up, length_max_motif, down, max_col_sum, motif_color_dict[max_motif], width)    
        
        # reduce df by removing rows that were used for plotting and the column (max_motif) used
        # reduce column as well. Otherwise, same motif can be used multiple times
        unused_df = df.loc[df[max_motif] == 0, ~df.columns.isin([max_motif])]

        print('succesfully reduced df')        
        motives_only.remove(max_motif)
        df = unused_df
        width -= 0.1
    
    # for an overlay plot
    plt.legend(bbox_to_anchor=(1.04, 1), loc="upper left")
    plt.xticks(fontsize='x-large')
    plt.yticks(fontsize='x-large')
    plt.xlabel('bp upstream', fontsize = 'x-large')
    plt.ylabel('frequency', fontsize = 'x-large')
    
    plt.rcParams['font.family'] = "Arial"
    
    overlay_output = o_name + "_overlaid.png"
    plt.savefig(overlay_output, bbox_inches='tight')    
    
def plot_all_motif_plots(df, m_infos, o_name, up, fasta_F, down, peaks_dict):
    """
    Parameters
    ----------
    df : dataframe
        a full dataframe which considers all chromosomes in a given sample and class.
        each row: fixed representative cleavage site where it is tuple of (chromId, direction, representative_fc)  
        each col: motif
        
        each cell value = 1 if a motif is present in the range defined by that particular fixed representative cleavage site. 
        (does not matter where you have the motif as long as it is within the range)
        
        each cell value = 0 if a motif is not present in the range defined by that particular fixed representative cleavage site. 
        
    m_infos : a list of tuple
        a list of tuple that contains the motif infomation.
        1st element of a tuple: motif name. e.g. AAUAAA
        2nd element of a tuple: length of a motif. e.g. 6

    o_name : string
        motif frequency plot output file name.
            
    up : int (positive integer)
        how many number of base pairs do we go upstream of fixed representative cleavage site?
    
    fasta_F : a fasta flie
        contains the reference genome sequence.
    
    down : int (negative integer)
        how many number of base pairs do we go downstream of fixed representative cleavage site?    

    peak_dictionary : dictionary of tuple
        key = motif
        value = (lower boundary, upper boundary)
         
    Returns
    -------
    ordered_motives : list
        order of the motives in that sample and class (from highest frequency to lowest frequency)
        
    total_score : int
        total motif score per sample (using 18 motives).
        
        For each motif : Whenever a peak is located within a certain range (produced by gtf), you increment score by 1.
        min score : 0
        max score : 18
    
    peaks : list of tuple
        a list of tuple where
        1st element = motif
        2nd element = peak position of that motif
        
    P.S. Additionally this function also plots all 18 motif frequency plots.       
    """                
    motives_only = [elem[0] for elem in m_infos]
    lengths_of_motives = [elem[1] for elem in m_infos]
    ordered_motives = []
    
    peaks = []
    scores = {}

    while len(motives_only) != 0:
        col_sums = get_colsum_list(df, motives_only)
        idx = col_sums.index(max(col_sums))
        
        max_col_sum = col_sums[idx]
        assert(max(col_sums) == max_col_sum)
        # select a motif with the highest frequency
        max_motif = motives_only[idx]
        length_max_motif = int(lengths_of_motives[idx])
        
        # subset df only if a fixed cleavage site has count of the selected motif > 0
        # and return row_names (i.e fixed_cleavage_site)
        # df.index returns the row name (fixed cleavage site)
        # df[max_motif] > 0 returns you a set of indices 
        # where it returns true if a max_motif column has value > 0 at that row and else otherwise 
        
        # return row name (i.e cleavage sites) where count in max_motif is bigger than 0
        used_fixed_cleavage_sites = list(df.index[df[max_motif] > 0])
        frequency_dict = get_frequency_by_position(used_fixed_cleavage_sites, max_motif, up, length_max_motif, fasta_F, down)
         
        (score, peak) = get_score(frequency_dict, up, length_max_motif, down, peaks_dict[max_motif])
        scores[max_motif] = score
        peaks.append((max_motif, peak))
           
        # plot 1 motif frequency plot at a time
        plot_linegraph(frequency_dict, max_motif, o_name, up, length_max_motif, down, max_col_sum)
        
        # reduce df by removing rows that were used for plotting and the column (max_motif) used
        # reduce column as well. Otherwise, same motif can be used multiple times
        unused_df = df.loc[df[max_motif] == 0, ~df.columns.isin([max_motif])]

        print('succesfully reduced df')
        
        motives_only.remove(max_motif)
        ordered_motives.append(max_motif)
        df = unused_df
    
    total_score = sum(scores.values())
        
    return ordered_motives, total_score, peaks

def read_peak_csvs(input_directory):
    """
    Parameters
    ----------
    input_directory : string
        directory toward input csv file (gtf_peaks.csv) which is constructed as follows:
        row[0] = motif
        row[1] = lower boundary for the peak of a motif.
        It is defined as amongst genomic positions smaller than the peak position, 
        the maximum position which has approximately 90% frequency of the peak.
        
        row[2] = upper boundary for the peak of a motif.
        It is defined as amongst genomic positions bigger than the peak position, 
        the minimum position which has approximately 90% frequency of the peak.
        
        In total 18 rows (i.e. 18 motives)
        
    Returns
    -------
    peak_dictionary : dictionary of tuple
        key = motif
        value = (lower boundary, upper boundary)
    """    
    peaks_dictionary = {}
    with open (input_directory) as csvfile:
        spamreader = csv.reader(csvfile)
        for row in spamreader:
            peaks_dictionary[row[0]] = (row[1], row[2])
    return peaks_dictionary

def change_dna_to_rna(sequence, direction):
    """
    Parameters
    ----------    
    sequence : string
        a current genome DNA sub-sequence with same length as motif 
        (always + strand because reference genome is always + strand).
        
        we need to convert this DNA into RNA so that we can decide
        whether this sub-sequence is identical to the motif or not.
        
    direction : character
        direction of DNA in which a read maps to.
        
    Returns
    -------
    corrected_string : string
        RNA version (5' -> 3') of the current genome DNA sub-sequence.
        Now you can directly compare it with the motif.
        
        i.e. change a subesequence of DNA into RNA so that it becomes compatible with motif (5' -> 3')
    """     
    corrected_sequence = []
    # if a read is mapping to - strand,
    # revert the DNA sequence and then make a complementary
    if direction == '-':
        reverted_subsequence = list(reversed(sequence))
        for elem in reverted_subsequence:
            if elem == 'A':
                corrected_sequence.append('U')
            elif elem == 'T':
                corrected_sequence.append('A')
            elif elem == 'G':
                corrected_sequence.append('C')
            elif elem == 'C':
                corrected_sequence.append('G')
    
    # if a read is mapping to + strand
    # only change T in the DNA -> U in RNA. other nucleotides stay the same
    elif direction == '+':
        for elem in sequence:
            if elem == 'A':
                corrected_sequence.append('A')
            elif elem == 'T':
                corrected_sequence.append('U')
            elif elem == 'G':
                corrected_sequence.append('G')
            elif elem == 'C':
                corrected_sequence.append('C')
                
    # convert a list of characters into a single string
    corrected_string = ''.join(corrected_sequence)
    return corrected_string

def search_motif_in_cs(representative_fc_concatenate, upstream, fasta, motif_tuple, down_bp):
    """
    Parameters
    ----------    
    representative_fc_concatenate : string
        a string made of 3 components (chromID, direction, representative_fc)
        which is structred as: chromID_direction_representative_fc
    
    upstream : int (positive integer)
        how many number of base pairs do we go upstream of a fixed representative cleavage site?
    
    fasta : a fasta flie
        contains the reference genome sequence.
    
    motif_tuple : tuple
        a tuple that contains current motif infomation.
        1st element of a tuple: motif name. e.g. AAUAAA
        2nd element of a tuple: length of a motif. e.g. 6
    
    down_bp : int (negative integer)
        how many number of base pairs do we go downstream of a fixed representative cleavage site?
    
    Returns
    -------
    True if motif is present within the range specified by upstream and down_bp.
    False if motif is not present within the range specified by upstream and down_bp.
    
    in other words this function checks whether a motif exist in a particular fixed representative cleavage sites.
    """     
    representative_fc_tuple = representative_fc_concatenate.split('_')
    chromosome = representative_fc_tuple[0]
    direction = representative_fc_tuple[1]
    cleavage_site = int(representative_fc_tuple[2])
    
    motif = motif_tuple[0]
    length_motif = int(motif_tuple[1])
    
    # if down_bp is 0, you dont search beyond cleavage site
    if down_bp == 0:
        # create a list from for example from 5(=6-1) to 40
        # because you dont want your motif to exceed the cleavage site
        # So we want to look at from 5 bp upstream to 40 bp upstream.
        upstream_basepairs = list(range(length_motif - 1, upstream + 1))
        
    # if down_bp is not 0, you search beyond cleavage site
    # search from down_bp ~ upstream
    else:
        upstream_basepairs = list(range(down_bp, upstream + 1))
    
    # search ranges from 5 to 40 (i.e. -40 ~ -5) if down_bp == 0
    # search ranges from -20 to 40 (i.e. -40 ~ +20) if down_bp != 0   
    for upstream_bp in upstream_basepairs:
        if direction == '-':
            start_point = cleavage_site + upstream_bp 
            end_point = cleavage_site + upstream_bp + length_motif - 1                
        
        elif direction == '+':
            start_point = cleavage_site - upstream_bp
            end_point = cleavage_site - upstream_bp + length_motif - 1
            
        subsequence = fasta.fetch(reference = chromosome, start = start_point, end = end_point + 1)
        # correct DNA subsequence into RNA so that it becomes compatible with motif
        corrected_subsequence = change_dna_to_rna(subsequence, direction)
        
        if corrected_subsequence == motif:
            return True
        else:
            continue      
    # by the time one comes here, it means one couldnt find that motif (e.g. 'AAUAAA') in the upstream (typically -40bp to -1bp) of this subsequence
    return False
    
def get_full_dataframe(total_representative_cs, until, Fasta_file, m_list, down):
    """
    Parameters
    ----------    
    total_representative_cs : list of string
        a list where each element is chromID_direction_representative_fc
    
    until : int (positive integer)
        how many number of base pairs do we go upstream of a fixed representative cleavage site?
    
    Fasta file : a fasta flie
        contains the reference genome sequence.
    
    m_list : a list of tuple
        a list of tuple that contains the motif infomation.
        1st element of a tuple: motif name. e.g. AAUAAA
        2nd element of a tuple: length of a motif. e.g. 6
    
    down : int (negative integer)
        how many number of base pairs do we go downstream of a fixed representative cleavage site?
    
    Returns
    -------
    full_dataframe : dataframe
        a full dataframe which considers all chromosomes in a given sample and class.
        each row: fixed representative cleavage site where it is tuple of (chromId, direction, representative_fc)  
        each col: motif
        
        each cell value = 1 if a motif is present in the range defined by that particular fixed representative cleavage site. 
        (does not matter where you have the motif as long as it is within the range)
        
        each cell value = 0 if a motif is not present in the range defined by that particular fixed representative cleavage site. 
    """       
    
    row_names = [elem for elem in total_representative_cs]
    col_names = [elem[0] for elem in m_list]
    
    # create an empty dataframe
    full_dataframe = pd.DataFrame(np.zeros((len(row_names), len(col_names))))
    
    # give column and row names to the empty dataframe
    full_dataframe.columns = col_names
    full_dataframe.index = row_names
    
    for representative_fc in full_dataframe.index:
        for motif_tuple in m_list:
            
            motif_present = search_motif_in_cs(representative_fc, until, Fasta_file, motif_tuple, down)
            
            if motif_present:
                
                # acessing the element by row_name and col_name
                # += in dataframe is not recommended.
                # If you want to update, copy the dataframe and then increment by 1
                copy_dataframe = full_dataframe.copy()
                full_dataframe.loc[representative_fc, motif_tuple[0]] = copy_dataframe.loc[representative_fc, motif_tuple[0]] + 1
                                
            elif not motif_present:
                continue

    return full_dataframe

def get_representative_cleavage_sites(tempfile):
    """
    Parameters
    ----------    
    tempfile : string
        bed file of specific class (e.g. intergenic). each row is a cluster of polyA reads with score.
        cluster_id in each row contains information about representative fixed cleavage site.
        
    Returns
    -------
    representative_FC_list : list of string
        a list where each element is chromID_direction_representative_fc
    """      
    representative_FC_list = []
    with open(tempfile, "r") as t:
        for line in t:          
            cluster_id = line.split()[3]
            chromId = cluster_id.split(':')[0]
            representative_fc = cluster_id.split(':')[1]
            direction = cluster_id.split(':')[2]
            
            representative_fc_concatenate = chromId + '_' + direction + '_' + representative_fc
            representative_FC_list.append(representative_fc_concatenate)            
    
    return representative_FC_list
       
def get_motif_info(input_directory):
    """
    Parameters
    ----------    
    input_directory : string
        input directory to csv file containing motif infomation.
       
    Returns
    -------
    motif_list : a list of tuple
        1st element of a tuple: a motif
        2nd element of a tuple: length of the motif 
    """      
    motif_list = []
    with open (input_directory) as csvfile:
        spamreader = csv.reader(csvfile, delimiter = ',')
        for row in spamreader:

            motif = row[0]
            length_motif = row[1]

            motif_list.append((motif, length_motif))
            
    return motif_list 

def get_args():        
    parser = argparse.ArgumentParser(description="searching for motives upstream of cleavage sites")

    parser.add_argument('--bed', dest = 'bed',
                        required = True,
                        help = 'bed file containing pA sites in a specific sample and class')
        
    parser.add_argument('--fasta', dest = 'fasta',
                        required = True,
                        help = 'fasta file dir')
   
    parser.add_argument('--out_name', dest = 'out_name',
                        required = True,
                        help = 'name template for motif frequency plots and motive orders.csv')

    parser.add_argument('--annotated', type = int, dest = 'annotated',
                        required = True,
                        help = 'which reads did you recieve?')    
    
    parser.add_argument('--window', dest = 'window',
                        required = True,
                        help = 'how many bp upstream for the searching motives')

    parser.add_argument('--motif_info_dir', dest = 'motif_info_dir',
                        required = True,
                        help = 'directory to a file which contains info to look for')

    parser.add_argument('--downstream', dest = 'downstream',
                        required = True,
                        help = 'whether you also look into downstream. and if yes, how far')    

    parser.add_argument('--peaks', dest = 'peaks',
                        required = True,
                        help = 'csv file that contains dictionary of key = motif, value = position of the peak in the frequency plot. based on gtf file')  
  
    args = parser.parse_args()
    
    bed = args.bed
   
    fasta_dir = args.fasta

    out_name = args.out_name
    annotated = args.annotated

    window = int(args.window)
    motif_info_dir = args.motif_info_dir
    downstream =  int(args.downstream)
    peaks = args.peaks
    
    return bed, fasta_dir, out_name, annotated, window, motif_info_dir, downstream, peaks

def run_process():

    bed, fasta_dir, out_name, annotated, window, motif_info_dir, downstream, peaks = get_args()
        
    fasta_file = pysam.FastaFile(fasta_dir)
    motif_infos = get_motif_info(motif_info_dir) 
    print('successfully got motif infos')
    
    print('motif info is: ' + str(motif_infos))
    total_represenatative_fc = []
    
    total_represenatative_fc = get_representative_cleavage_sites(bed)
    print('successfully got total representative cleavage_sites.')
    
    print('total_representative_fixed_cleavage sites: ' + str(total_represenatative_fc))
    
    full_dataframe_all = get_full_dataframe(total_represenatative_fc, window, fasta_file, motif_infos, downstream)
    print('successfully got full dataframe_all')
    
    print('full dataframe is: ' + str(full_dataframe_all))

    peaks_dictionary = read_peak_csvs(peaks)
    print("successfully read the peaks dictionary")
    
    ordered_motives, tot_score, peaks_in_data = plot_all_motif_plots(full_dataframe_all, motif_infos, out_name, window, fasta_file, downstream, peaks_dictionary)
    print('successfully drew all plots, got all ordered motives and scores for all motives')
    
    ordered_motif_csv_name = out_name + "_ordered_motives.csv"
    # e.g. 10X_P4_7
    sample_name = '_'.join(out_name.split('_')[0:3])
    # e.g. annotated
    sub_header_name = out_name.split('_')[3]
    # modify class name so that it is compatible with the terminology used in the paper.
    sub_header_name = modify_class_name(annotated)
    
    ordered_motives.insert(0, sample_name)
    ordered_motives.insert(1, sub_header_name)
    
    write_out(ordered_motives, ordered_motif_csv_name)
    print('successfully saved the ordered_motives')
    
    score_out_name = out_name + "_scores.csv"

    score_data = peaks_in_data
    score_data.insert(0, sample_name)
    score_data.insert(1, sub_header_name)
    score_data.insert(2, tot_score)
    write_out(score_data, score_out_name)
    print('successfully saved total score for the sample')
    
    plot_all_motif_overlaid_plots(full_dataframe_all, motif_infos, out_name, window, fasta_file, downstream)
    print('successfully drew an overlaid motif frequency plots')
    
if __name__ == "__main__":
    run_process()
    print("success")