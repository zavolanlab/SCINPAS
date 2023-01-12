# -*- coding: utf-8 -*-
"""

Created on Mon Jan 10 10:26:42 2022

@author: Youngbin Moon (y.moon@unibas.ch)
"""
import csv
import matplotlib.pyplot as plt
import argparse







"""
Aim : get a barchart of the number of polyA sites in each sample.
"""
def read_counts(input_directory):
    """
    Parameters
    ----------
    input_directory : string
        directory toward input csv file.

    Returns
    -------
    names_list : list
        a list that contains sample names used to run our pipeline.
        
    counts_list : list
        a list that contains the number of polyA sites in each sample.
    """    
    names_list = []
    counts_list = []
    with open (input_directory) as csvfile:
        spamreader = csv.reader(csvfile)
        for row in spamreader:
            names_list.append(row[0])
            counts_list.append(float(row[1]))
    
    return names_list, counts_list

def autolabel(rects, ax):
    """
    Parameters
    ----------
    rects : BarContainer object
        contains all bars
        
    ax : matplotlib.axes._subplots.AxesSubplot
        axes of the figure.
        
    Returns
    -------
    returns nothing but attaches a text label above each bar displaying its height
    """       
    for rect in rects:
        height = rect.get_height()
        ax.text(rect.get_x() + rect.get_width()/2., 1.05*height,
                '%d' % float(height),
                ha='center', va='bottom')
        
def plot_bar (names, counts, file_name):
    """
    Parameters
    ----------
    names : list
        a list that contains all sample names used to run our pipeline.
    
    counts : list
        a list that contains the number of polyA sites in each sample.
    
    file_name : string
        output barchart name.
        
    Returns
    -------
    returns nothing but draws a barchart of the number of polyA sites from different samples.
    """   
    fig, axis = plt.subplots()
    rectangles = axis.bar(names, counts) 
    autolabel(rectangles, axis)
    
    plt.xticks(
    rotation=45, 
    horizontalalignment='right',
    fontweight='light',
    fontsize='x-large'  
    )
    
    plt.tight_layout()
        
    plt.yticks(fontsize='x-large')
    plt.ylabel('number of polyA sites', fontsize = 'x-large')
    plt.rcParams['font.family'] = "Arial"
    plt.yscale("log") 
    plt.savefig(file_name, bbox_inches='tight')

def get_inputs():
    
    parser = argparse.ArgumentParser(description="get a barchart of number of polyA sites for all samples")
    parser.add_argument('--input_dir', dest = 'input_dir',
                        required = True,
                        help = 'input csv file')
    
    parser.add_argument('--out_file', dest = 'out_file',
                        required = True,
                        help = 'creates a bar chart')
    
    args = parser.parse_args()
    
    input_dir = args.input_dir
    out_file = args.out_file

    return input_dir, out_file
    
def run_process():
    input_dir, file = get_inputs()
    print('successfully got inputs')
    
    names_list, counts_list = read_counts(input_dir)
    print('successfully read the data')
    
    plot_bar(names_list, counts_list, file)
    print('successfully plot barchart')

if __name__ == "__main__":
    run_process()
    print('success')