# -*- coding: utf-8 -*-
"""

Created on Mon Jan 10 10:26:42 2022

@author: Youngbin Moon (y.moon@unibas.ch)
"""
import csv
import matplotlib.pyplot as plt
import argparse

"""
Aim : Draw a barchart of motif scores from different samples.

18 motives were used (12 + 6) for scoring. 
For each motif : Whenever a peak is located within a certain range (produced by gtf), you increment score by 1.
min score : 0
max score : 18

Score was computed using all "polyA" reads.

"""
def read_counts(input_directory):
    """
    Parameters
    ----------
    input_directory : string
        directory toward input csv file (polyA_combined_motif_scores.csv).

    Returns
    -------
    corrected_names_list : list
        a list that contains all sample names used to run our pipeline.
        "UmiDedup" is a negative control. Hence name is converted into "negative_control"
        
    counts_list : list
        a list that contains motif scores of each sample.
    """    
    names_list = []
    counts_list = []
    with open (input_directory) as csvfile:
        spamreader = csv.reader(csvfile)
        for row in spamreader:
            names_list.append(row[0])
            counts_list.append(float(row[1]))
    
    corrected_names_list = [elem if "UmiDedup" not in elem else "negative_control" for elem in names_list]
    return corrected_names_list, counts_list

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
        "UmiDedup" is a negative control. Hence name is converted into "negative_control".        
    
    counts : list
        a list that contains motif scores of each sample.
    
    file_name : string
        output barchart name.
    Returns
    -------
    returns nothing but draws a barchart of motif scores from different samples.
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
    plt.ylabel('motif scores', fontsize = 'x-large')
    plt.rcParams['font.family'] = "Arial"
    plt.savefig(file_name, bbox_inches='tight')

def get_inputs():
    
    parser = argparse.ArgumentParser(description="get a barchart of motif scores from different samples. (Using all polyA reads)")
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
