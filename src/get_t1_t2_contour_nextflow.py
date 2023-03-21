# -*- coding: utf-8 -*-
"""

Created on Mon Jan 10 10:26:42 2022

@author: Youngbin Moon (y.moon@unibas.ch)
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import argparse
import seaborn as sns
import math





"""
Aim1 : To plot sns KDE and REG of avg distances between start of terminal exon to representative cleavage site of pA site in type1 vs type2.
Aim2 : To plot histogram of avg distances between start of terminal exon to representative cleavage site of pA site in type1 vs type2.
Aim3 : To plot cumulative histogram of ratio of avg distances between type1 and type2. (log2(avg_d1/avg_d2))
where distance is a distance between start of terminal exon to representative cleavage site of pA site.
"""
def plot_cumul_histogram(distances_df, file, t1, t2):
    """
    Parameters
    ----------
    distances_df : dataframe
        1st column : a list of distances between start of terminal exon to representative cleavage site of pA site in type 1 cell.
        2nd coluimn : a list of distances between start of terminal exon to representative cleavage site of pA site in type 2 cell.
                
    file : string
        histogram plot file name
        
    t1 : string
        type1 cell name
        
    t2 : string
        type2 cell name   

    Returns
    -------
    returns nothing but plots a cumulative histogram of distances between terminal exon start to representative cleavage site of pA site in type1 vs type2.
    """   
    plt.figure()
    pre_x = distances_df[t1]
    pre_y = distances_df[t2]
    
    x = [elem for elem in pre_x]
    y = [elem for elem in pre_y]
    assert(len(x)==len(y))

    log_ratio_x_y = [math.log2((x[i]+1)/(y[i]+1)) for i in range(len(x))]
    bins = np.linspace(min(log_ratio_x_y), max(log_ratio_x_y), 100)
    
    ratio_label = t1 + "/" + t2
    plt.hist(log_ratio_x_y, bins = bins, histtype = "step", cumulative = True, density = True, label = ratio_label, ec = 'blue', fc = 'None')
    plt.xticks(fontsize='x-large')
    plt.yticks(fontsize='x-large')
    plt.xlabel('log2 fold change', fontsize = 'x-large')
    plt.ylabel('frequency', fontsize = 'x-large')
    
    plt.rcParams['font.family'] = "Arial"
    plt.legend(bbox_to_anchor=(1.04, 1), loc="upper left")    
    
    plt.savefig(file, bbox_inches='tight')
    
def plot_histogram(distances_df, file, t1, t2):
    """
    Parameters
    ----------
    distances_df : dataframe
        1st column : a list of distances between start of terminal exon to representative cleavage site of pA site in type 1 cell.
        2nd coluimn : a list of distances between start of terminal exon to representative cleavage site of pA site in type 2 cell.
                
    file : string
        histogram plot file name
        
    t1 : string
        type1 cell name
        
    t2 : string
        type2 cell name   

    Returns
    -------
    returns nothing but plots a histogram of distances between terminal exon start to representative cleavage site of pA site in type1 vs type2.
    """    
    plt.figure()
    pre_x = distances_df[t1]
    pre_y = distances_df[t2]
    
    x = [math.log10(elem) if elem != 0 else math.log10(1) for elem in pre_x]
    y = [math.log10(elem) if elem != 0 else math.log10(1) for elem in pre_y]
    bins = np.linspace(min(min(x), min(y)), max(max(x), max(y)), 100)
    
    plt.hist(x, bins, alpha=0.5, label = t1)
    plt.hist(y, bins, alpha=0.5, label = t2)
    
    plt.xticks(fontsize='x-large')
    plt.yticks(fontsize='x-large')
    plt.xlabel('log10 distance', fontsize = 'x-large')
    plt.ylabel('frequency', fontsize = 'x-large')
    
    plt.rcParams['font.family'] = "Arial"
    plt.legend(bbox_to_anchor=(1.04, 1), loc="upper left")    
    
    plt.savefig(file, bbox_inches='tight')
    
def plot_Contour(distances_df, file, t1, t2, plot_type):
    """
    Parameters
    ----------
    distances_df : dataframe
        1st column : a list of distances between start of terminal exon to representative cleavage site of pA site in type 1 cell.
        2nd coluimn : a list of distances between start of terminal exon to representative cleavage site of pA site in type 2 cell.
                
    file : string
        tri-contour plot file name
        
    t1 : string
        type1 cell name
        
    t2 : string
        type2 cell name   

    Returns
    -------
    returns nothing but plots a contour plot of distances between terminal exon start to representative cleavage site of pA site in type1 vs type2.
    """
    plt.figure()
    
    if plot_type == "kde":     
        sns.kdeplot(data = distances_df, x = t1, y = t2, fill = True, log_scale = True)
    
    else:
        sns.jointplot(data = distances_df, x = t1, y = t2, kind = plot_type)  
        
    plt.xlabel(t1, fontsize = 'x-large')                         
    plt.ylabel(t2, fontsize = 'x-large')
    
    x_data = list(distances_df[t1])
    x_min = min(x_data)
    x_max = max(x_data)
    
    plt.plot([x_min, x_max], [x_min, x_max], linewidth=2, color = 'brown')
    plt.rcParams['font.family'] = "Arial"
    
    plt.savefig(file, bbox_inches='tight')
    
def run ():    
    parser = argparse.ArgumentParser(description="draw terminal exon average distance plot (cell type1 vs cell type2)")
    parser.add_argument('--distance_csv_dir', dest = 'distance_csv_dir',
                        required = True,
                        help = 'distance csv containing tuple of distance to terminal exon (cell type1, cell type2)')

    parser.add_argument('--output_name', dest = 'output_name',
                        required = True,
                        help = 'output name of the scatter plot')                     

    parser.add_argument('--type1', dest = 'type1',
                        required = True,
                        help = 'first cell type')

    parser.add_argument('--type2', dest = 'type2',
                        required = True,
                        help = 'Second cell type')
    
    args = parser.parse_args()
    
    distance_csv_dir = args.distance_csv_dir 
    output_template = args.output_name
    
    type1 = args.type1 + "_distance(bp)"
    type2 = args.type2 + "_distance(bp)"
    
    distance_csv = pd.read_csv(distance_csv_dir, delimiter = '\t', header = None)
    distance_csv.columns = [type1, type2]
    
    return distance_csv, output_template, type1, type2
    
def run_process():
    distance_csv, output_template, type1, type2 = run()
    print('successfully initialized')
    
    # kde plot
    out_name1 = output_template + "KDE.png"
    p_type1 = "kde"
    plot_Contour(distance_csv, out_name1, type1, type2, p_type1)
    print('successfully plotted contour plot1')

    # regression plot
    out_name2 = output_template + "REG.png"
    p_type2 = "reg"
    plot_Contour(distance_csv, out_name2, type1, type2, p_type2)
    print('successfully plotted contour plot2')
    
    # distribution plot (histogram)
    out_name3 = output_template + "HIST.png"
    plot_histogram(distance_csv, out_name3, type1, type2)
    print('successfully plotted histogram')
    
    # cumulative distribution plot (histogram)
    out_name4 = output_template + "cumul_HIST.png"
    plot_cumul_histogram(distance_csv, out_name4, type1, type2)
    print('successfully plotted cumulative histogram')    
    
if __name__ == "__main__":
    run_process()
    print('success')