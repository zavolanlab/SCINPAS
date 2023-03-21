# -*- coding: utf-8 -*-
"""

Created on Mon Jan 10 10:26:42 2022

@author: Youngbin Moon (y.moon@unibas.ch)
"""
import csv
import matplotlib.pyplot as plt
import argparse
import numpy as np





"""
Aim: draw bar plot of read counts in each class (for a given sample)
"""

def read_counts(input_directory):
    """
    Parameters
    ----------
    input_directory : string
        input directory to the counts.csv that contains read counts in each class for a given sample.

    Returns
    -------
    names_list : list
        sorted list of each class name (type) of reads according to our need.
        we want in this order: 'original', 'deduplicated', 'softclipped', 'nonPolyA', 'allPolyA', 'terminal+A', 'exonic_A'
        
    counts_list : list
        sorted list of the number of reads in each class according to our need.
        we want in this order: 'original', 'deduplicated', 'softclipped', 'nonPolyA', 'allPolyA', 'terminal+A', 'exonic_A'
    """
    names_list = []
    counts_list = []
    with open (input_directory) as csvfile:
        spamreader = csv.reader(csvfile)
        for row in spamreader:
            names_list.append(row[1])
            counts_list.append(float(row[2]))
    return names_list, counts_list

def autolabel(rects, ax):
    """
    Parameters
    ----------
    rects : matplotlib.container.BarContainer
        Bars representing the read counts in each class.
        
    ax : matplotlib.axes._subplots.AxesSubplot
        The figure panel.

    Returns
    -------
    returns nothing but it attaches a text label above each bar displaying its height. (to the figure panel)
    """
    count = 0
    variation = 1.05
    for rect in rects:
        height = rect.get_height()
        # instead of fixing to 1.05, fiddle a little bit where text is likely to overlap
        if count == 2:
            variation += 1.2
        elif count > 2:
            variation = 1.05    
        # %s and %d are placeholders for a string and a number respectively. 
        ax.text(rect.get_x() + rect.get_width()/2., variation*height,
                '%s' % '{:.2E}'.format(float(height)),
                ha='center', va='bottom')
        count += 1
        
def plot_bar (names, counts, file_name):
    """
    Parameters
    ----------
    names : list
        sorted list of each class name (type) of reads according to our need.
        
    counts : list
        sorted list of the number of reads in each class according to our need.
        
    file_name : string
        output bar plot name.

    Returns
    -------
    returns nothing but draws a bar plot.
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
    plt.yscale('log')
    
    plt.yticks(fontsize='x-large')
    plt.ylabel('Number of reads', fontsize = 'x-large')
    plt.rcParams['font.family'] = "Arial"
    
    plt.savefig(file_name, bbox_inches='tight')

def get_inputs():
    
    parser = argparse.ArgumentParser(description="get barchart of read counts in each BAM files")
    parser.add_argument('--input_dir', dest = 'input_dir',
                        required = True,
                        help = 'input file for counts_full.csv')
    
    parser.add_argument('--out_file', dest = 'out_file',
                        required = True,
                        help = 'creates a bar chart')
    args = parser.parse_args()
    
    input_dir = args.input_dir
    out_file = args.out_file

    return input_dir, out_file
    
def run():
    input_dir, file = get_inputs()
    print('successfully got inputs')
    
    names_list, counts_list = read_counts(input_dir)
    print('successfully read the counts info')
    
    plot_bar(names_list, counts_list, file)
    print('successfully plot the bar plot')

if __name__ == "__main__":
    run()
    print('success')
    
