# -*- coding: utf-8 -*-
"""

Created on Sun Dec 19 22:55:32 2021

@author: Youngbin Moon (y.moon@unibas.ch)
"""










import argparse
import csv
import pandas as pd
import multiprocessing as mp
"""
Aim1 : Count the number of polyA sites in a given sample and in a given class 
Aim2 : and save it in csv file.
"""

def change_class_name(class_type):
    """
    Parameters
    ----------        
    class_type : string
        class of reads. (either annotated(4), unannotated(5), intronic(7), intergenic(8), exonic(10) or all_polyA(12))
           
    Returns
    -------
    modified_class_name : string
        modified class name of reads in order to use same terminology as paper. 
        (either ATE, UTE, I, IG, NTE, PATR)
    """  
    if class_type == 4:
        modified_class_name = "ATE"

    elif class_type == 5:
        modified_class_name = "UTE"
 
    elif class_type == 7:
        modified_class_name = "I"

    elif class_type == 8:
        modified_class_name = "IG"
        
    elif class_type == 10:
        modified_class_name = "NTE"        

    elif class_type == 12:
        modified_class_name = "PATR"   

    elif class_type == 14:
        modified_class_name = "TE"   
        
    return modified_class_name
   
def write_out(row_data, o_file):
    """
    Parameters
    ----------
    row_data : list
        a list of tuple where each element is
        [sample_name, class_name, num_polyA_sites]
            
    o_file : string
        output file name    
       
    Returns
    -------        
    returns nothing but writes row_data at the csv file with output name: o_file.
    """        
    w = csv.writer(open(o_file, "w"))
    w.writerow(row_data)

def num_polyA(bed, threshold):
    """
    Parameters
    ----------
    bed : bed file
        a bed file which contains pA sites for a particular strand and chromosome.
        
    threshold : int
        threshold for scores of a pA site (cluster) in order to be considered as a true pA site.
       
    Returns
    -------        
    num_polyA : int
        the number of pA sites within a bed file (specific sample, specific class)
    """            
    print('bed: ' + str(bed))
    filtered_bed = bed[bed['score'] >= threshold].copy()

    print('filtered_bed: ' + str(filtered_bed))
    num_polyA = len(filtered_bed)
    
    return num_polyA

def get_args():        
    parser = argparse.ArgumentParser(description="get number of polyA sites")

    parser.add_argument('--bed', dest = 'bed',
                        required = True,
                        help = 'bed file containing pA sites of a specific class')
            
    parser.add_argument('--csv_out_name', dest = 'csv_out_name',
                        required = True,
                        help = 'number of polyA sites csv file name')
    
    parser.add_argument('--annotated', type = int, dest = 'annotated',
                        required = True,
                        help = 'which reads did you recieve?')  

    parser.add_argument('--cluster_threshold', type = int, dest = 'cluster_threshold',
                        required = True,
                        help = 'a threshold on cluster size in order to be considered as a pA site') 

    parser.add_argument('--n_cores', type = int, dest = 'n_cores',
                        required = True,
                        help = 'the number of cores available') 
    
    args = parser.parse_args()
    
    bed_dir = args.bed    
    csv_out_name = args.csv_out_name
    annotated = args.annotated 
    cluster_threshold = args.cluster_threshold
    
    n_cores = args.n_cores
    
    return bed_dir, csv_out_name, annotated, cluster_threshold, n_cores

def run_process():

    bed_dir, csv_out_name, annotated, cluster_threshold, n_cores = get_args()
    print('successfully got arguments')
        
    bedfile = pd.read_csv(bed_dir, delimiter = '\t', header = None, names = ['seqid', 'start', 'end', 'id', 'score', 'strand'])
        
    groupedByChrStrand = bedfile.groupby(['seqid', 'strand'])
    splitted_dataframes = [group for name, group in groupedByChrStrand]
    # set the number of cores available
    pool = mp.Pool(n_cores)
    
    # call apply_async() without callback
    result_objects = [pool.apply_async(num_polyA, args = (splitted_df, cluster_threshold)) for splitted_df in splitted_dataframes]

    # result_objects is a list of pool.ApplyResult objects
    overlap_results = [r.get() for r in result_objects]
    
    print("successfully got num_pA" )
    
    pool.close()
    pool.join()

    num_polyA_sites = sum(overlap_results)
    
    print('successfully got total num_pA')
    # e.g. 10X_P4_7
    sample_name = '_'.join(csv_out_name.split('_')[0:3])
    if "UmiDedup" in sample_name:
        # sample name at this point is : "L10X_P7_14_UmiDedup" where L can be any alphabetic letter. you want original sample name which is 10X_P7_14.
        class_name = sample_name[1:]
        sample_name = "negative_control"
        
    else:
        class_name = change_class_name(annotated)
        
    num_polyA_data = [sample_name, class_name, num_polyA_sites]
        
    write_out(num_polyA_data, csv_out_name)
    print('successfully saved the result')
    
if __name__ == "__main__":
    run_process()
    print("success")