# -*- coding: utf-8 -*-
"""

Created on Sun Dec 19 22:55:32 2021

@author: Youngbin Moon (y.moon@unibas.ch)
"""

import csv
import argparse






"""
Aim: fuse 3 csv files together so that
each row is sample name, class name, motif score and # of PAS. (including negative control)
"""
def write_fused_csv (fused_csv, o_file):
    """
    Parameters
    ----------
    fused_csv : dictionary
    returns a fused dictionary where
        1st key : sample name
        2nd key : class name (e.g. annotated, unannotated, intronic, intergenic, exonic, control)
        value : [motif_score, #_PAS]
    
    o_file : string
        output csv file name
        
    Returns
    -------
    returns nothing but writes output (fused_csv) with o_file name.
    """    
    w = csv.writer(open(o_file, "w", newline=''))
   
    for sample in fused_csv.keys():
        for class_name in fused_csv[sample].keys():
            row_data = fused_csv[sample][class_name]
            row_data.insert(0, sample)
            row_data.insert(1, class_name)
            w.writerow(row_data)
            
def fuse_csv(m_score_dict, polyA_m_score_dict, num_polyA_csv):
    """
    Parameters
    ----------
    m_score_dict : dictionary
        1st key = sample name (does not include negative control)
        2nd key = class name
        value = motif score at each class
            
    polyA_m_score_dict : dictionary
        1st key = sample name (includes negative control)
        value = motif score at polyA level 
    
    num_polyA_csv : dictionary
        1st key = sample name (includes negative control)
        2nd key = class name
        value = # of PAS
        
    Returns
    -------
    new_dict : dictionary
    returns a fused dictionary where
        1st key : sample name
        2nd key : class name (e.g. annotated, unannotated, intronic, intergenic, exonic, control)
        value : [motif_score, #_PAS]
    """
    new_dict = m_score_dict.copy()
    for sample_name in m_score_dict.keys():
        for class_name in m_score_dict[sample_name].keys():
            new_dict[sample_name][class_name].append(num_polyA_csv[sample_name][class_name])
    
    sample_name = list(num_polyA_csv['negative_control'].keys())[0]
    new_dict[sample_name] = {}
                
    if 'control' not in new_dict[sample_name].keys():
        new_dict[sample_name]['control'] = []
        
    new_dict[sample_name]['control'].append(polyA_m_score_dict['negative_control'])
    new_dict[sample_name]['control'].append(num_polyA_csv['negative_control'][sample_name])
    
    return new_dict

def read_counts(input_directory, type_csv):
    """
    Parameters
    ----------
    input_directory : string
        input directory to the csv file that contains read counts in each class for a given sample.
        input directory can be one of three types.
        1) directory of total_motif_scores_csv: contains motif scores for each sample and for each class.
        2) directory of polyA_csv: contains motif scores of all polyA reads for each sample. (including negative control)
        3) directory of total_num_poly_csv: contains number of PAS for each sample and for each class (including negative control)
    
    Returns
    -------
    dictionary : dictionary
    returns a dictionary where it can be of 3 types
    1) total_motif_scores_csv:
        1st key = sample name (does not include negative control)
        2nd key = class name
        value = motif score at each class
            
    2) polyA_csv:
        1st key = sample name (includes negative control)
        value = motif score at polyA level
    
    3) total_num_polyA_csv:
        1st key = sample name (includes negative control)
        2nd key = class name
        value = # of PAS
    """
    dictionary = {}
    with open (input_directory) as csvfile:
        spamreader = csv.reader(csvfile)
        for row in spamreader:
            # for polyA_csv
            # row[0] =  sample name, row[1] = score
            if len(row) == 2:
                dictionary[row[0]] = row[1]
            
            # for total_motif_scores_csv and total_num_polyA_csv
            # row[0] = sample name, row[1] = class name, row[2] = motif score at polyA level or # of PAS.
            elif len(row) > 2:
                if row[0] not in dictionary.keys():
                    dictionary[row[0]] = {}
                
                if type_csv == "motif_score":
                    if row[1] not in dictionary[row[0]].keys():
                        dictionary[row[0]][row[1]] = []
                    dictionary[row[0]][row[1]].append(row[2])
                
                elif type_csv == "num_polyA":
                    dictionary[row[0]][row[1]] = row[2]

    return dictionary

def get_args():        
    parser = argparse.ArgumentParser(description="fuse 3 csv files together")

    parser.add_argument('--total_motif_scores_csv', dest = 'total_motif_scores_csv',
                        required = True,
                        help = 'csv file containing all motif scores in each class and sample')
        
    parser.add_argument('--polyA_csv', dest = 'polyA_csv',
                        required = True,
                        help = 'csv file containing all motif scores at all polyA reads in each sample')

    parser.add_argument('--total_num_polyA_csv', dest = 'total_num_polyA_csv',
                        required = True,
                        help = 'csv file containing all number of PAS in each sample')  
    
    parser.add_argument('--csv_out_name', dest = 'csv_out_name',
                        required = True,
                        help = 'a new csv file containing motif score and number of PAS for each sample and each class')
    
    args = parser.parse_args()
    
    total_motif_scores_csv = args.total_motif_scores_csv
    polyA_csv = args.polyA_csv
    total_num_polyA_csv = args.total_num_polyA_csv
    csv_out_name = args.csv_out_name
    
    return total_motif_scores_csv, polyA_csv, total_num_polyA_csv, csv_out_name

def run_process():

    total_motif_scores_csv, polyA_csv, total_num_polyA_csv, csv_out_name = get_args()
    print('successfully got arguments')

    motif_score_dict = read_counts(total_motif_scores_csv, "motif_score")
    print('successfully read motif_score_csv')
    
    polyA_motif_score_dict = read_counts(polyA_csv, "polyA_motif_score")
    print('successfully read polyA_motif_score_csv')
    
    num_polyA_csv = read_counts(total_num_polyA_csv, "num_polyA")
    print('successfully read num_polyA_csv')
    
    updated_m_score_dict = fuse_csv(motif_score_dict, polyA_motif_score_dict, num_polyA_csv)
    print('successfully fused csv file')
    
    write_fused_csv(updated_m_score_dict, csv_out_name)    
    print('successfully wrote output')
    
if __name__ == "__main__":
    run_process()
    print("success")