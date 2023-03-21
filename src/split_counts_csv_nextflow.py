# -*- coding: utf-8 -*-
"""

Created on Mon May  9 16:28:57 2022

@author: Youngbin Moon (y.moon@unibas.ch)
"""

import csv
import argparse






"""
Aim: From a count csv file which contains the number of counts in each class and in each sample,
split according to the sample such that you have csv file containing the number of counts in all classes.
"""

def find_index(n_list, element):
    """
    Parameters
    ----------
    n_list : list
        list of each class name (type) of reads
    
    element : string
        current class you are searching for
        
    Returns
    -------
    idx = int
        index of where your current class is in the n_list
    """
    
    idx = n_list.index(element)
    return idx
    
def sort_lists(names, counts):
    """
    Parameters
    ----------
    names : list
        list of each class name (type) of reads.
    
    counts : list
        list of the number of reads in each class.
        It matches 1:1 with "names" list
        
    Returns
    -------
    organized_names : list
        sorted list of each class name (type) of reads according to our need.
        we want in this order: 'original', 'deduplicated', 'non-PATR', 'softclipped', 'PATR', 'TE', 'ATE', 'UTE', 'NTE', 'I', 'IG'
    
    organized_counts : list
        sorted list of the number of reads in each class according to our need.
        we want in this order: 'original', 'deduplicated', 'non-PATR', 'softclipped', 'PATR', 'TE', 'ATE', 'UTE', 'NTE', 'I', 'IG'
    """
    
    organized_names = ['original', 'deduplicated', 'non-PATR', 'softclipped', 'PATR', 'TE', 'ATE', 'UTE', 'NTE', 'I', 'IG']
    
    print("names: " + str(names))
    print("counts: " + str(counts))
    indices = [find_index(names, elem) for elem in organized_names]
    organized_counts = [counts[elem] for elem in indices]
    
    return organized_names, organized_counts
    
def write_out(samples, names_list, counts_list):
    """
    Parameters
    ----------
    samples : string
        class of the reads
        
    names_list : list 
        sorted list of each class name (type) of reads according to our need.
        
    counts_list : list
        sorted list of the number of reads in each class according to our need.

    Returns
    -------
    returns nothing but write the csv file such that it contains the number of counts in all classes for a given sample. 
    """
    
    o_file = samples[0] + '.csv'
    # newline is used to avoid having blank lines
    w = csv.writer(open(o_file, "w", newline=''))
    # sort according to the workflow
    sorted_names, sorted_counts = sort_lists(names_list, counts_list)
    
    for i in range(len(sorted_names)):
        w.writerow([samples[i], sorted_names[i], sorted_counts[i]])

def split_csv(dict_count):
    """
    Parameters
    ----------
    counts_dict : dictionary
        a dictionary which contains the number of counts in each class and in each sample.

    Returns
    -------
    returns nothing but it calls write_out function
    to split according to the sample such that you have csv file containing the number of counts in all classes.
    """
    
    for key, val in dict_count.items():
        sample_name = key
        names = val[0]
        counts = val[1]
        
        if "negative_control" in sample_name:
            # recover original file name so that later you can put csv and figure in the original folder.
            # sample_name = negative_control(A10X_P7_14UmiDedup)
            sample_name = sample_name.split('(')[1].split(')')[0]
                          
        if len(names) == len(counts):
            sample_names = [sample_name for i in range(len(names))]
            write_out(sample_names, names, counts)
            print('successfully generated split csv file for 1 sample')
            
        else:
            print('lengths do not match')

def read_counts_csv(input_directory):
    """
    Parameters
    ----------
    input_directory : string
        input count csv file which contains the number of counts in each class and in each sample.

    Returns
    -------
    counts_dict : dictionary
        a dictionary which contains the number of counts in each class and in each sample.
    """
    
    counts_dict = {}
    with open (input_directory) as csvfile:
        spamreader = csv.reader(csvfile)
        # key (1st column) is "sample name"
        for row in spamreader:
            key = row[0]
            # 2nd column is "type"
            # 3rd column is "counts"
            if key not in counts_dict.keys():
                counts_dict[key] = ([], [])
            # UMI-tools deduplicated data does not have intergenic or intronic reads.
            # Hence do not plot bar plot for intergenic or intronic reads
            if "negative_control" in row[0]:
                if int(row[2]) == 0 or int(row[2]) == 0:
                    print("Passed")
                    # continue makes you go next iteration at this point
                    continue
            counts_dict[key][0].append(row[1])
            counts_dict[key][1].append(int(row[2]))
    return counts_dict

def get_inputs():
    
    parser = argparse.ArgumentParser(description="get barchart of read counts in each BAM files")
    parser.add_argument('--total_csv_input', dest = 'total_csv_input',
                        required = True,
                        help = 'total_csv_input')
    
    args = parser.parse_args()
    
    total_csv_input = args.total_csv_input
    return total_csv_input
            
def run():
    
    input_csv = get_inputs()
    print('successfully got the input')
    
    counts_dict = read_counts_csv(input_csv)
    print('successfully generated counts_dict')
    
    split_csv(counts_dict)
    print('successfully split all csvs according to sample names')

if __name__ == "__main__":    
    run()
    print('success')