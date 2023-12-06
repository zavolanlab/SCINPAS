
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
import os
from collections import Counter
import pandas as pd
import statistics
"""
Aim 1 : Using all terminal exons from gtf, draw all 18 motif frequency plots.
Aim 2 : Using all terminal exons from gtf, draw an overlaid motif frqeuncy plot which contains all 18 motif frequency plots.
Aim 3 : Using all terminal exons from gtf, save order of motives (from highest frequency to lowest frequency) in the csv.
Aim 4 : Using all terminal exons from gtf, save peak positions of all motives at the csv file.
      
        18 motives were used (12 + 6). 
"""
def write_peaks(dictionary, o_file):
    """
    Parameters
    ----------
    dictionary : dictionary of tuple
        key = motif
        value = (lower boundary, upper boundary)
            
    o_file : string
        output file name of peak positions.    
       
    Returns
    -------        
    returns nothing but writes peak positions of all motives at the csv file with output name: o_file.
    """       
    w = csv.writer(open(o_file, "w", newline=''))
    for key in dictionary.keys():
        
        val = dictionary[key]
        row_data = [key, val[0], val[1]]
        w.writerow(row_data)
        
def write_out(row_data, o_file):
    """
    Parameters
    ----------
    row_data : list
        a list that contains 20 elements. 
        1st element: sample name
        2nd element: class name. In this case: 'annotated(known)'
        3rd ~ 20th element: ordered motives. (from highest frequency to lowest frequency)
            
    o_file : string
        motif orders output file name    
       
    Returns
    -------        
    returns nothing but writes motif orders at the csv file with output name: o_file.
    """     
    w = csv.writer(open(o_file, "w"))
    w.writerow(row_data)
        
def find_nearest(frequency, positions):
    """
    Parameters
    ----------
    frequency : list
        a list of smoothen frequencies of the motif occurrence.        
    
    positions : a list
        a list of positions in bp.
        
    Returns
    ------- 
    peak_boundary_tuple : tuple
        a tuple which contains lower and upper boundary of the peak position of the motif.
        i.e. (lower boundary, upper boundary)
        
        lower boundary is defined as amongst genomic positions smaller than the peak position, 
        the maximum position which has approximately 90% frequency of the peak.
        
        Upper boundary is defined as amongst genomic positions bigger than the peak position, 
        the minimum position which has approximately 90% frequency of the peak.    
        
        Assuming unimodal peak because annotated gtf file has good profile in general (especially in the canonical motives)
    """      
    peak = max(frequency)
    peak_position_idx = frequency.index(max(frequency))
    # value = 90% of peak
    value = peak * 0.9

    array = np.asarray(frequency)

    # from index 0 to peak position find the index of the nearest value to 0.9*peak
    lower_boundary_idx = (np.abs(array - value))[0 : peak_position_idx].argmin()
    uncorrected_upper_boundary_idx = (np.abs(array - value))[peak_position_idx + 1 : ].argmin()
    
    # because index starts from 0, you have to correct it to get real index.
    upper_boundary_idx = uncorrected_upper_boundary_idx + peak_position_idx + 1
    
    peak_boundary_tuple = (positions[lower_boundary_idx], positions[upper_boundary_idx])       
    
    return peak_boundary_tuple

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

def plot_linegraph_overlay(frequency_dictionary, motif_of_interest, up_bp, length_motif_interest, down_bp, tot_freq, motif_color, width):
    """
    Parameters
    ----------
    frequency_dictionary : dictionary
        key = position in bp
        value = frequency. (considering all fixed representative cleavage sites that has the motif)
        
        e.g. if representativeCleavageSite_1 has the motif at -40, -30, -20 and if representativeCleavageSite_2 has motif at -30 and -20bp,
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
        
def plot_linegraph(frequency_dictionary, motif_of_interest, file, up_bp, length_motif_interest, down_bp, tot_freq):
    """
    Parameters
    ----------
    frequency_dictionary : dictionary
        key = position in bp
        value = frequency. (considering all representative cleavage sites that has the motif)
        
        e.g. if representativeCleavageSite_1 has the motif at -40, -30, -20 and if representativeCleavageSite_2 has motif at -30 and -20bp,
        frequency_dict: -40: 1/3, -30: 5/6, -20: 5/6
        
    motif_of_interest : string
        a motif of interest (the most frequently appearing motif at that moment)
        
    file : string
        ouput figure file name template.
        
    up_bp : int (positive integer)
        how many number of base pairs do we go upstream of representative cleavage site?
    
    length_motif_interest : int
        length of the most frequently appearing motif at that moment.
            
    down_bp : int (negative integer)
        how many number of base pairs do we go downstream of representative cleavage site?
    
    tot_freq : int
        the maximum column sum of full_dataframe at that moment.
        i.e. the number of occurrences of the most frequently appearing motif at that moment.       
        
    Returns
    ------- 
    peak_position_tuple : tuple
        a tuple which contains lower and upper boundary of the peak position of the motif.
        i.e. (lower boundary, upper boundary)
        
        lower boundary is defined as amongst genomic positions smaller than the peak position, 
        the maximum position which has approximately 90% frequency of the peak.
        
        Upper boundary is defined as amongst genomic positions bigger than the peak position, 
        the minimum position which has approximately 90% frequency of the peak.
    
    Additionally this function plots a motif frequency plot of interest.
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
    
    plt.legend(loc = 'upper right')
    
    plt.xticks(fontsize='x-large')
    plt.yticks(fontsize='x-large')
    plt.xlabel('bp upstream', fontsize = 'x-large')
    plt.ylabel('frequency', fontsize = 'x-large')
    
    plt.rcParams['font.family'] = "Arial"
    plt.savefig(full_out_name, bbox_inches='tight')
    
    peak_position_tuple = find_nearest(smoothen_list, base_pairs_upstream)
    
    return peak_position_tuple

def get_frequency_by_position(selected_row_names, selected_motif_only, up_until, length_motif, fasta_f, down_bp):
    """
    Parameters
    ----------
    selected_row_names : list
        a list of representative cleavage sites which has the motif of interest 
        (the most frequently appearing motif at that moment) within the range specified.
        
    selected_motif_only : string
        a motif of interest (the most frequently appearing motif at that moment)
        
    up_until : int (positive integer)
        how many number of base pairs do we go upstream of representative cleavage site?
    
    length_motif : int
        length of the motif
        
    fasta_f : a fasta flie
        contains the reference genome sequence.
    
    down_bp : int (negative integer)
        how many number of base pairs do we go downstream of representative cleavage site?
         
    Returns
    -------        
    frequency_dict : dictionary
        key = position in bp
        value = frequency. (considering all representative cleavage sites that has the motif)
        
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
        a full dataframe which considers all terminal exons from the terminal_exons.bed (a subset of gtf).
        each row: representative cleavage site where it is tuple of (chromId, direction, representative_fc)  
        each col: motif
        
        each cell value = 1 if a motif is present in the range defined by that particular representative cleavage site. 
        (does not matter where you have the motif as long as it is within the range)
        
        each cell value = 0 if a motif is not present in the range defined by that particular representative cleavage site. 
    
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

def plot_all_motif_overlaid_plots(df, m_infos, file_template, up, fasta_F, down):
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

    file_template : string
        output file template
            
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
        unused_df = df.loc[df[max_motif] == 0, ~df.columns.isin([max_motif])].copy()

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
    
    plt.rcParams['font.family'] = "DejaVu Sans"
    
    file1 = file_template + "_overlaid.png"
    file2 = file_template + "_overlaid.svg"
    
    plt.savefig(file1, bbox_inches='tight')
    plt.savefig(file2, bbox_inches='tight') 
    
def plot_all_motif_plots(df, m_infos, o_name, up, fasta_F, down):
    """
    Parameters
    ----------
    df : dataframe
        a full dataframe which considers all terminal exons from the terminal_exons.bed (a subset of gtf).
        each row: representative cleavage site where it is tuple of (chromId, direction, representative_fc)  
        each col: motif
        
        each cell value = 1 if a motif is present in the range defined by that particular representative cleavage site. 
        (does not matter where you have the motif as long as it is within the range)
        
        each cell value = 0 if a motif is not present in the range defined by that particular representative cleavage site.
        
        Note here we use representative cleavage site instead of fixed representative cleavage site
        because terminal exons come from gtf but not bam file.
        
    m_infos : a list of tuple
        a list of tuple that contains the motif infomation.
        1st element of a tuple: motif name. e.g. AAUAAA
        2nd element of a tuple: length of a motif. e.g. 6
    
    o_name : string
        motif frequency plot output file name.
        
    up : int (positive integer)
        how many number of base pairs do we go upstream of representative cleavage site?
    
    fasta_F : a fasta flie
        contains the reference genome sequence.
    
    down : int (negative integer)
        how many number of base pairs do we go downstream of representative cleavage site?    

    Returns
    -------
    ordered_motives : list
        order of the motives in the terminal_exons.bed (from highest frequency to lowest frequency)
    
    peaks : dictionary of tuple
        key = motif
        value = (lower boundary, upper boundary)
        
    P.S. Additionally this function also plots all 18 motif frequency plots based on terminal_exons.bed.
        
    """       
    motives_only = [elem[0] for elem in m_infos]
    lengths_of_motives = [elem[1] for elem in m_infos]
    ordered_motives = []
    peaks = {}
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
        
        peak = plot_linegraph(frequency_dict, max_motif, o_name, up, length_max_motif, down, max_col_sum)
        
        # reduce df by removing rows that were used for plotting and the column (max_motif) used
        # You should reduce column as well. Otherwise, same motif can be used multiple times
        # df.columns.isin([max_motif]) returns True only if you are in that column, else False.
        # ~df.columns.isin([max_motif]) reverses the result of df.columns.isin([max_motif])
        unused_df = df.loc[df[max_motif] == 0, ~df.columns.isin([max_motif])].copy()
                
        motives_only.remove(max_motif)
        
        ordered_motives.append(max_motif)
        peaks[max_motif] = peak
        df = unused_df
    
    return ordered_motives, peaks

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
        how many number of base pairs do we go upstream of representative cleavage site?
    
    fasta : a fasta flie
        contains the reference genome sequence.
    
    motif_tuple : tuple
        a tuple that contains current motif infomation.
        1st element of a tuple: motif name. e.g. AAUAAA
        2nd element of a tuple: length of a motif. e.g. 6
    
    down_bp : int (negative integer)
        how many number of base pairs do we go downstream of representative cleavage site?
    
    Returns
    -------
    True if motif is present within the range specified by upstream and down_bp.
    False if motif is not present within the range specified by upstream and down_bp.
    
    in other words this function checks whether a motif exist in a particular representative cleavage sites.
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
        how many number of base pairs do we go upstream of representative cleavage site?
    
    Fasta file : a fasta flie
        contains the reference genome sequence.
    
    m_list : a list of tuple
        a list of tuple that contains the motif infomation.
        1st element of a tuple: motif name. e.g. AAUAAA
        2nd element of a tuple: length of a motif. e.g. 6
    
    down : int (negative integer)
        how many number of base pairs do we go downstream of representative cleavage site?
    
    Returns
    -------
    full_dataframe : dataframe
        a full dataframe which considers all terminal exons from the terminal_exons.bed (a subset of gtf).
        each row: representative cleavage site where it is tuple of (chromId, direction, representative_fc)  
        each col: motif
        
        each cell value = 1 if a motif is present in the range defined by that particular representative cleavage site. 
        (does not matter where you have the motif as long as it is within the range)
        
        each cell value = 0 if a motif is not present in the range defined by that particular representative cleavage site.
        
        Note here we use representative cleavage site instead of fixed representative cleavage site
        because terminal exons come from gtf but not bam file.
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
        bed file containing all terminal exons of a gtf file. each row is a terminal exon.
        
    Returns
    -------
    representative_FC_list : list of string
        a list where each element is chromID_direction_representative_fc
    """      
    representative_FC_list = []
    with open(tempfile, "r") as t:
        for line in t:          
            chromId = line.split()[0]
            read_start = line.split()[1]
            read_end = line.split()[2]
            print('line: ' + str(line.split()))
            print('score: ' + str(line.split()[4]))
            direction = line.split()[5]
            
            if direction == '-':
                representative_fc = read_start
                
            elif direction == '+':
                representative_fc = read_end  
                      
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

    parser.add_argument('--terminal_exons', dest = 'terminal_exons',
                        required = True,
                        help = 'terminal exons.bed input file')
        
    parser.add_argument('--fasta', dest = 'fasta',
                        required = True,
                        help = 'fasta file dir')
    
    parser.add_argument('--out_name', dest = 'out_name',
                        required = True,
                        help = 'name template for motif frequency plots, peak_positions.csv and motive orders.csv')
    
    parser.add_argument('--window', dest = 'window',
                        required = True,
                        help = 'how many bp upstream for the searching motives')

    parser.add_argument('--motif_info_dir', dest = 'motif_info_dir',
                        required = True,
                        help = 'directory to a file which contains info to look for')

    parser.add_argument('--downstream', dest = 'downstream',
                        required = True,
                        help = 'whether you also look into downstream. and if yes, how far')    
    
    args = parser.parse_args()
    
    terminal_exons = args.terminal_exons
    fasta_dir = args.fasta
    out_name = args.out_name
    
    window = int(args.window)
    motif_info_dir = args.motif_info_dir
    downstream =  int(args.downstream)
    
    return terminal_exons, fasta_dir, out_name, window, motif_info_dir, downstream

def run_process():

    terminal_exons, fasta_dir, out_name, window, motif_info_dir, downstream = get_args()
        
    fasta_file = pysam.FastaFile(fasta_dir)
    motif_infos = get_motif_info(motif_info_dir)
    print('successfully got motif infos')
    
    total_represenatative_fc = []
        
    total_represenatative_fc = get_representative_cleavage_sites(terminal_exons)
    print('successfully got total representative cleavage_sites')
     
    full_dataframe = get_full_dataframe(total_represenatative_fc, window, fasta_file, motif_infos, downstream)
    print('successfully got full dataframe')
        
    ordered_motives, peaks = plot_all_motif_plots(full_dataframe, motif_infos, out_name, window, fasta_file, downstream)
    print('successfully drew all plots')

    ordered_motif_csv_name = out_name + "_ordered_motives.csv"
    
    # i.e. gtf
    header_name = out_name
    # i.e. annotated
    sub_header_name = 'annotated(known)'
    
    ordered_motives.insert(0, header_name)
    ordered_motives.insert(1, sub_header_name)
    
    write_out(ordered_motives, ordered_motif_csv_name)
    print('successfully saved ordered motives')
    
    peaks_name = out_name + "_peaks.csv"
    write_peaks(peaks, peaks_name)
    print('successfully written peaks position')
    
    plot_all_motif_overlaid_plots(full_dataframe, motif_infos, out_name, window, fasta_file, downstream)
    print('successfully drew an overlaid motif frequency plots')
    
if __name__ == "__main__":
    run_process()
    print("success")