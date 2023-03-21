# -*- coding: utf-8 -*-
"""

Created on Mon Jan 10 10:26:42 2022

@author: Youngbin Moon (y.moon@unibas.ch)
"""
import csv
import matplotlib.pyplot as plt
import argparse
import pandas as pd
import numpy as np




"""
Aim : to draw barcharts of 
    1) total number of expressed genes from the species of interest
    2) the number of genes expressed by each sample from dedup bam file
    3) the number of genes covered by each sample from polyA reads mapping to terminal exons
    in 1 figure.
"""

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
        
def plot_barchart(dataframe, file_name, tot_num_genes):
    """
    Parameters
    ----------
    dataframe : dataframe
        a merged dataframe from 
        1) dataframe containing the number of genes expressed by each sample from dedup bam file.
        2) dataframe containing the number of genes covered by each sample from polyA reads mapping to terminal exons.
          
    file_name : string
        output file name
    
    tot_num_genes : int
        total number of expressed genes from the species of interest
        
    Returns
    -------
    returns nothing but draws barcharts of 
    1) total number of expressed genes from the species of interest
    2) the number of genes expressed by each sample from dedup bam file
    3) the number of genes covered by each sample from polyA reads mapping to terminal exons
    in 1 figure.
    """   
    labels = list(dataframe['sample'])
    num_genes_expressed = list(dataframe['#_total_genes_expressed'])
    num_gene_covered = list(dataframe['gene_covered'])
    
    # the label locations                                               
    # shifted by 1 to the right because on the first position, you will put the total number of genes.
    # i.e. at index: 0
    x = np.array((range(1, len(labels)+1)))
    # the width of the bars
    width = 0.35
    fig, axis = plt.subplots()
    
    rectangles1 = plt.bar(x - width/2, num_genes_expressed, width, label = "# genes expressed")
    rectangles2 = plt.bar(x + width/2, num_gene_covered, width, label = "# of genes covered")
                   
    # a single bar of total number of genes 
    rectangles3 = plt.bar(0, tot_num_genes, width, label = "total # genes")
                              
    plt.xticks(
    rotation=45, 
    horizontalalignment='right',
    fontweight='light',
    fontsize='x-large'  
    )
    
    # Add some text for labels, title and custom x-axis tick labels, etc.
    # plt.ylim((1, 50000))
    plt.yticks(fontsize='x-large')
    plt.ylabel('number of genes or number of reads', fontsize = 'x-large')
    plt.yscale('log')
    plt.xticks(x, labels, fontsize='x-large')                           
#    autolabel(rectangles1, axis)
#    autolabel(rectangles2, axis)
#    autolabel(rectangles3, axis)
    
    plt.legend(bbox_to_anchor=(1.04, 1), loc="upper left")
    
    plt.rcParams['font.family'] = "Arial"
    plt.savefig(file_name, bbox_inches='tight')
  
def merge_dataframe(Dedup_csv_dir, PolyA_csv_dir, tot_gene_dir):
    """
    Parameters
    ----------
    Dedup_csv_dir : string
        a directory towards csv file containing number of genes expressed by each sample from dedup bam file.
    
    PolyA_csv_dir : string
        a directory towards csv file containing number of genes covered by each sample from polyA reads mapping to terminal exons
    
    tot_gene_dir : string
        a directory towards csv file containing total number of expressed genes from the species of interest
        
    Returns
    -------
    merged_dataframe : dataframe
        a merged dataframe from 
        1) dataframe containing the number of genes expressed by each sample from dedup bam file.
        2) dataframe containing the number of genes covered by each sample from polyA reads mapping to terminal exons.
        
    total_num_genes : int
        total number of expressed genes from the species of interest    
    """       
    dedup_csv = pd.read_csv(Dedup_csv_dir, header = None, names = ['sample', '#_total_genes_expressed'])
    polyA_csv = pd.read_csv(PolyA_csv_dir, header = None, names = ['sample', 'gene_covered'])
    total_gene_csv = pd.read_csv(tot_gene_dir, header = None, names = ['gtf', 'total_#_genes'])
        
    merged_dataframe = dedup_csv.merge(polyA_csv, left_on = 'sample', right_on = 'sample', how = 'inner')
    
    total_num_genes = int(total_gene_csv['total_#_genes'])
                             
    return merged_dataframe, total_num_genes

def get_inputs():
    
    parser = argparse.ArgumentParser(description="get a barchart of gene coverage for all samples")
    parser.add_argument('--dedup_csv_dir', dest = 'dedup_csv_dir',
                        required = True,
                        help = 'csv file containing # genes expressed generated from dedup bam file')

    parser.add_argument('--polyA_csv_dir', dest = 'polyA_csv_dir',
                        required = True,
                        help = 'csv file containing gene coverage generated from polyA + terminal exons reads')

    parser.add_argument('--total_gene_dir', dest = 'total_gene_dir',
                        required = True,
                        help = 'csv file containing total number of genes')
   
    parser.add_argument('--coverage_out_file', dest = 'coverage_out_file',
                        required = True,
                        help = 'creates a bar chart of gene coverage by each sample')

    args = parser.parse_args()
    
    dedup_csv_dir = args.dedup_csv_dir
    polyA_csv_dir = args.polyA_csv_dir
    total_gene_dir = args.total_gene_dir
        
    coverage_out_file = args.coverage_out_file
    
    
    return dedup_csv_dir, polyA_csv_dir, total_gene_dir, coverage_out_file
    
def run_process():
    dedup_csv_dir, polyA_csv_dir, total_gene_dir, coverage_out_file = get_inputs()
    print('successfully got inputs')
    
    merged_dataframe, total_num_genes = merge_dataframe(dedup_csv_dir, polyA_csv_dir, total_gene_dir)
    print("successfully got merged_dataframe")
    print("merged data frame is: " + str(merged_dataframe))
    
    plot_barchart(merged_dataframe, coverage_out_file, total_num_genes)
    print("successfully got the plot")
    
if __name__ == "__main__":
    run_process()
    print('success')