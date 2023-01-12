# -*- coding: utf-8 -*-
"""

Created on Mon Jan 10 10:26:42 2022

@author: Youngbin Moon (y.moon@unibas.ch)
"""

import csv
import matplotlib.pyplot as plt
import numpy as np
import argparse

"""
Aim: get a histogram that shows distribution of 'absolute' distance between end of a read and end of terminal exon.
"""
def read_distances (input_directory):
    """
    Parameters
    ----------
    input_directory : string
        An input directory to the csv file that contains 
        a distribution of 'absolute' distance between end of a read and end of terminal exon.
    
    Returns
    -------
    min_distances : list of int
        contains 'absolute' distances between end of a read and end of terminal exon.
    
    are_minus : list of bool
        contains directions of a reads which correspond to min_distances.
    
    min_log_distances : list of float
        log10 value of min_distances
    """     
    min_distances = []
    are_minus = []
    min_log_distances = []
    with open (input_directory) as csvfile:
        spamreader = csv.reader(csvfile)
        for row in spamreader:
            min_distances.append(int(row[1]))
            are_minus.append(row[2])
            min_log_distances.append(float(row[3]))
    return min_distances, are_minus, min_log_distances

def get_mean(scores_list):
    """
    Parameters
    ----------
    scores_list : list of float
        A list that contains log10 'absolute' distances between end of a read and end of terminal exon.
        
    Returns
    -------
    mean_val : float
        mean value of log10 'absolute' distances between end of a read and end of terminal exon.    
    """
    mean_val = np.mean(scores_list)         
    return mean_val
        
def plot_histogram(file, log_distances, cntrl_log_distances, use_Fc):
    """
    Parameters
    ----------
    file : string
        output histogram file name
    
    log_distances : list of float
        A list that contains log10 'absolute' distances between end of a read of a sample and end of terminal exon.

    cntrl_log_distances : list of float
        A list that contains log10 'absolute' distances between end of a read of a negative control and end of terminal exon.
        
    use_fc : bool
        whether you use fixed cleavage site or not.
        
    Returns
    -------
    returns nothing but plots a histogram that shows 
    a distribution of 'absolute' distance between end of a read and end of terminal exon.
    """
    plt.figure()
    plt.xticks(fontsize='x-large')
    plt.yticks(fontsize='x-large')
    plt.xlabel('log10 distance', fontsize = 'x-large')
    plt.ylabel('frequency', fontsize = 'x-large')
    
    plt.rcParams['font.family'] = "Arial"
        
    bins = np.arange(min(min(log_distances), min(cntrl_log_distances)), max(max(log_distances), max(cntrl_log_distances)), 0.05)
                      
    plt.hist(log_distances, bins = bins, histtype = "step", cumulative = True, density = True, label = "sample",  ec = 'blue', fc = 'None') 
    plt.hist(cntrl_log_distances, bins = bins, histtype = "step", cumulative = True, density = True, label = "negative control", ec = 'purple', fc = 'None')                     

    plt.legend(bbox_to_anchor=(1.04, 1), loc="upper left")
    
    plt.savefig(file, bbox_inches='tight')

def get_inputs():
    
    parser = argparse.ArgumentParser(description="get distance histogram")
    parser.add_argument('--csv_input', dest = 'csv_input',
                        required = True,
                        help = 'distances csv file of sample')
    
    parser.add_argument('--cntrl_csv', dest = 'cntrl_csv',
                        required = True,
                        help = 'distances csv file of negative control')    

    parser.add_argument('--file_name', dest = 'file_name',
                        required = True,
                        help = 'output file_name')   
    
    parser.add_argument('--use_fc', dest = 'use_fc',
                        required = True,
                        help = 'use fixed cleavage sites')
    
    args = parser.parse_args()
    
    input_dir = args.csv_input
    control_input_dir = args.cntrl_csv
    file_name = args.file_name
    use_fc = args.use_fc
    
    return input_dir, control_input_dir, file_name, use_fc
                        
def run_process():
    
    input_dir, control_input_dir, file_name, use_fc = get_inputs()
    print('successfully got inputs')
    
    min_distances, are_minus, min_log_distances = read_distances(input_dir)
    print('successfully got distanes for the data')
    
    _, _, control_min_log_distances = read_distances(control_input_dir)
    print('successfully got distanes for the negative control')
    
    plot_histogram(file_name, min_log_distances, control_min_log_distances, use_fc)
    print('successfully got a histograms')

if __name__ == "__main__":
    
    run_process()
    print('success')
    
    
    