# -*- coding: utf-8 -*-
"""
Created on Sat Feb 19 18:49:24 2022

@author: Youngbin Moon (y.moon@unibas.ch)
"""
import matplotlib.pyplot as plt
import argparse
import csv
import numpy as np

"""
Aim: Plot a distribution of span of a cluster with same CB + UB tags or CB + YB tags.
"""
def draw_histogram(logspans, file_template): 
    """
    Parameters
    ----------
    logspans : a list
        A list of spans in log10 scale.
    
    file_template : string
        output file template

    Returns
    -------
    returns nothing but plots a histogram.

    """
    bins = np.arange(min(logspans), max(logspans), 0.1)
    plt.figure()
    histogram = plt.hist(logspans, bins=bins, histtype='bar', facecolor = 'blue')                  
    # set y above 0 because log100 is not defined
    plt.ylim(1, 100*max(histogram[0]))
    plt.yscale("log")
    
    plt.xticks(fontsize='x-large')
    plt.yticks(fontsize='x-large')
    plt.xlabel('log spans', fontsize = 'x-large')
    plt.ylabel('frequency', fontsize = 'x-large')
        
    plt.rcParams['font.family'] = "DejaVu Sans"
    
    file1 = file_template + ".png"
    file2 = file_template + ".svg"
    
    plt.savefig(file1, bbox_inches='tight')
    plt.savefig(file2, bbox_inches='tight')
    
def read_spans (input_directory):
    """
    Parameters
    ----------
    input_directory : string
        input directory to the csv file that contains log10 span of the cluster with same CB + UB tags or CB + YB tags.

    Returns
    -------
    logspans : a list
        A list of spans in log10 scale. 

    """
    logspans = []
    with open (input_directory) as csvfile:
        spamreader = csv.reader(csvfile)
        for row in spamreader:
            logspans.append(float(row[0]))
    return logspans

def get_input():
    parser = argparse.ArgumentParser(description="get span histogram")
    parser.add_argument('--input_csv', dest = 'input_csv',
                        required = True,
                        help = 'input csv file')
    
    parser.add_argument('--out_hist', dest = 'out_hist',
                        required = True,
                        help = 'output histogram')
        
    args = parser.parse_args()
    input_csv = args.input_csv
    out_hist = args.out_hist
    
    return input_csv, out_hist

def run_process():
    input_csv, out_hist = get_input()
    ("successfully got input")
    
    log_spans = read_spans(input_csv)
    ("successfully got log_spans")
    
    draw_histogram(log_spans, out_hist)
    ("successfully got spans histogram")

if __name__ == "__main__":
    run_process()
    print("success")