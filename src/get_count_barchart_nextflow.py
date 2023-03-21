# -*- coding: utf-8 -*-
"""

Created on Mon Jan 10 10:26:42 2022

@author: Youngbin Moon (y.moon@unibas.ch)
"""
import matplotlib.pyplot as plt
import argparse
import pandas as pd
import numpy as np






"""
Aim1: draw bar plot of read counts in each class (for a given sample)
Aim2: draw bar plot of read counts in each class (for a negative control)
"""
def plot_bar(dataframe, file_name):
    """
    Parameters
    ----------
    dataframe : dataframe
        a merged dataframe from 
        1) dataframe containing the number of reads in each class of a sample
        2) dataframe containing the number of reads in each class of a negative control
          
    file_name : string
        output file name
    
    Returns
    -------
    returns nothing but draws barcharts of 
    1) the number of reads in each class of a sample
    2) the number of reads in each class of a negative control
    in 1 figure.
    """   
    labels = list(dataframe['class'])
    num_reads_samples = list(dataframe['#_reads_sample'])
    num_reads_control = list(dataframe['#_reads_control'])
    
    # the label locations                                               
    # shifted by 1 to the right because on the first position, you will put the total number of genes.
    # i.e. at index: 0
    x = np.array((range(0, len(labels))))
    # the width of the bars
    width = 0.35
    fig, axis = plt.subplots()
    
    plt.bar(x - width/2, num_reads_samples, width, label = "sample")
    plt.bar(x + width/2, num_reads_control, width, label = "control")                   
                             
    plt.xticks(
    rotation=45, 
    horizontalalignment='right',
    fontweight='light',
    fontsize='x-large'  
    )
    
    plt.tight_layout()
    plt.yticks(fontsize='x-large')
    plt.ylabel('Number of reads', fontsize = 'x-large')
    plt.yscale('log')
    plt.xticks(x, labels)                           
    
    plt.legend(bbox_to_anchor=(1.04, 1), loc="upper left")
    
    plt.rcParams['font.family'] = "Arial"
    plt.savefig(file_name, bbox_inches='tight')
        
def merge_count_dataframe(sample_dir, control_dir):
    """
    Parameters
    ----------
    sample_dir : string
        directory towards sample full_count.csv containing the number of reads in each class of a sample
    
    control_dir : string
        directory towards control full_count.csv containing the number of reads in each class of a negative control
        
    Returns
    -------
    merged_dataframe : dataframe
        a merged dataframe from 
        1) dataframe containing the number of reads in each class of a sample
        2) dataframe containing the number of reads in each class of a negative control
         
    """       
    sample_count_csv = pd.read_csv(sample_dir, header = None, names = ['sample', 'class', '#_reads_sample'])
    control_count_csv = pd.read_csv(control_dir, header = None, names = ['sample', 'class', '#_reads_control'])

    subset_sample_count_csv = sample_count_csv[['class', '#_reads_sample']]
    subset_control_count_csv = control_count_csv[['class', '#_reads_control']]
    print("subset_sample_count_csv")
    print(str(subset_sample_count_csv))
    print("subset_control_count_csv")
    print(str(subset_control_count_csv))    
    merged_dataframe = subset_sample_count_csv.merge(subset_control_count_csv, left_on = 'class', right_on = 'class', how = 'inner')
    print(str(merged_dataframe))
    return merged_dataframe

def get_inputs():
    
    parser = argparse.ArgumentParser(description="get barchart of read counts in each BAM files")
    parser.add_argument('--input_dir', dest = 'input_dir',
                        required = True,
                        help = 'input directory of sample counts_full.csv')

    parser.add_argument('--control_input_dir', dest = 'control_input_dir',
                        required = True,
                        help = 'input directory of negative control counts_full.csv')
    
    parser.add_argument('--out_file', dest = 'out_file',
                        required = True,
                        help = 'creates a bar chart')

    args = parser.parse_args()
    
    input_dir = args.input_dir
    control_input_dir = args.control_input_dir
    out_file = args.out_file

    return input_dir, control_input_dir, out_file
    
def run():
    input_dir, control_input_dir, out_file = get_inputs()
    print('successfully got inputs')
    
    merged_dataframe = merge_count_dataframe(input_dir, control_input_dir)
    print('successfully read the counts info')
    
    plot_bar(merged_dataframe, out_file)
    print('successfully plot the bar plot')

if __name__ == "__main__":
    run()
    print('success')
    
