# -*- coding: utf-8 -*-
"""

Created on Mon Jan 10 10:26:42 2022

@author: Youngbin Moon (y.moon@unibas.ch)
"""

import csv
import matplotlib.pyplot as plt
import numpy as np
import argparse
import seaborn as sns
import pandas as pd




"""
Aim: get a histogram that shows 2 distributions of 'absolute' distance between a pA site(representative cleavage site) and an end of a terminal exon.
1 distribution for sample.
1 distribution for negative control.
"""
def read_distances (input_directory):
    """
    Parameters
    ----------
    input_directory : string
        An input directory to the csv file that contains 
        a distribution of 'absolute' distance between a pA site(representative cleavage site) and an end of a terminal exon.
    
    Returns
    -------
    min_distances : list of int
        contains 'absolute' distance between a pA site(representative cleavage site) and an end of a terminal exon.
    
    are_minus : list of bool
        contains direction of a PA site
    
    min_log_distances : list of float
        log10 value of min_distances
    """     
    min_distances = []
    are_minus = []
    min_log_distances = []
    with open (input_directory) as csvfile:
        spamreader = csv.reader(csvfile)
        for row in spamreader:
            min_distances.append(int(row[0]))
            are_minus.append(row[1])
            min_log_distances.append(float(row[2]))
    return min_distances, are_minus, min_log_distances

def get_mean(scores_list):
    """
    Parameters
    ----------
    scores_list : list of float
        A list that contains log10 'absolute' distance between a pA site(representative cleavage site) and an end of a terminal exon..
        
    Returns
    -------
    mean_val : float
        mean value of log10 'absolute' distance between a pA site(representative cleavage site) and an end of a terminal exon. 
    """
    mean_val = np.mean(scores_list)         
    return mean_val
        
def plot_histogram(file_template, log_distances, cntrl_log_distances):
    """
    Parameters
    ----------
    file_template : string
        output file template
    
    log_distances : list of float
        A list that contains log10 'absolute' distance between a pA site(representative cleavage site) from sample and an end of a terminal exon.

    cntrl_log_distances : list of float
        A list that contains log10 'absolute' distance between a pA site(representative cleavage site) from negative control and an end of a terminal exon.
                
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
            
    sample_names = ['sample'] * len(log_distances)
    control_names = ['control'] * len(cntrl_log_distances)
                       
    log_distances += cntrl_log_distances        
    sample_names += control_names
      
    distance_data = {'distance': log_distances, 'name': sample_names}
    distance_df = pd.DataFrame(distance_data)
    
    sns.ecdfplot(data = distance_df, x = 'distance', hue = 'name')                    

    plt.legend(bbox_to_anchor=(1.04, 1), loc="upper left")
    
    plt.rcParams['font.family'] = "DejaVu Sans"
    
    file1 = file_template + ".png"
    file2 = file_template + ".svg"
    
    plt.savefig(file1, bbox_inches='tight')
    plt.savefig(file2, bbox_inches='tight')

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
        
    args = parser.parse_args()
    
    input_dir = args.csv_input
    control_input_dir = args.cntrl_csv
    file_name = args.file_name
    
    return input_dir, control_input_dir, file_name
                        
def run_process():
    
    input_dir, control_input_dir, file_name = get_inputs()
    print('successfully got inputs')
    
    min_distances, are_minus, min_log_distances = read_distances(input_dir)
    print('successfully got distanes for the data')
    
    _, _, control_min_log_distances = read_distances(control_input_dir)
    print('successfully got distanes for the negative control')
    
    plot_histogram(file_name, min_log_distances, control_min_log_distances)
    print('successfully got a histograms')

if __name__ == "__main__":
    
    run_process()
    print('success')
    
    
    