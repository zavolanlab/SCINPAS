# -*- coding: utf-8 -*-
"""

Created on Mon May  9 16:28:57 2022

@author: Youngbin Moon (y.moon@unibas.ch)
"""
import csv
import argparse

"""
Aim : From total_merged_motif_orders.csv, separate motif orders by sample.

total_merged_motif_orders.csv:
    row[0] = sample name
    row[1] = read class
    row[2] = most frequent motif
    row[3] = 2nd most frequent motif
    row[4] = 3rd most frequent motif
    row[5] = 4th most frequent motif
    ......
    row[19] = the least frequent motif
"""

def read_counts_csv(input_directory):
    """
    Parameters
    ----------
    input_directory : string
        directory toward input csv file (total_merged_motif_orders.csv)
        which contains motif orders for all samples and for all classes
        
    Returns
    -------
    counts_dict : dictionary of tuple where each element is a list.
    
        key = sample name
        1st value of the tuple = read type (class)
        2nd value of the tuple = motif orders (1st-18th) in the form of list of list
    """    
    counts_dict = {}
    with open (input_directory) as csvfile:
        spamreader = csv.reader(csvfile)
        # key (1st column) is "sample name"
        for row in spamreader:
            key = row[0]
            # 2nd column is "type"
            # 3rd to 20th column are "motif orders"
            if key not in counts_dict.keys():
                counts_dict[key] = ([], [])
            counts_dict[key][0].append(row[1])
            counts_dict[key][1].append(row[2:])
    print(counts_dict)
    return counts_dict

def find_index(n_list, element):
    """
    Parameters
    ----------
    n_list : list
        a list containing all classes which are ordered as:
        ["annotated", "unannotated", "intronic", "intergenic","exonic"].
        
    element : string
        One of the read classes. It can be either "annotated", "unannotated", 
        "intronic", "intergenic" or "exonic".
        
    Returns
    -------
    idx : int
        index of the location of element in the n_list
    """
    idx = n_list.index(element)
    return idx
    
def sort_lists(names, counts):
    """
    Parameters
    ----------    
    names : list of strings
        a list containing all classes which are not ordered.
    
    counts : list of list
        It contains motif orders (1st-18th) in the form of list of list
    
    Returns
    -------
    organized_names : list of strings
        a list containing all classes which are ordered as :
        ["annotated", "unannotated", "intronic", "intergenic","exonic"].
        
    organized_counts : list of list
        It contains motif orders (1st-18th) in the form of list of list
        but ordered such that it corresponds to organized_names.
        
    """   
    organized_names = ['annotated', 'unannotated', 'intergenic', 'intronic', 'exonic']
    indices = [find_index(names, elem) for elem in organized_names]
    organized_counts = [counts[elem] for elem in indices]
    return organized_names, organized_counts
    
def write_out(samples, names_list, motives_list):
    """
    Parameters
    ----------
    samples : list
        a list which contains a specific sample name len(names_list) times (or len(motives_list) times).
        e.g. ['intergenic', 'intergenic', 'intergenic', 'intergenic', intergenic' ......]
    
    names_list : list of strings
        a list containing all classes which are not ordered.
        
    motives_list : list of list
        It contains motif orders (1st-18th) in the form of list of list
    
    Returns
    -------
    returns nothing but writes motif orders for each class in a "specific" sample.
    """       
    o_file = samples[0] + '.csv'
    # newline is used to avoid having blank lines
    w = csv.writer(open(o_file, "w", newline=''))
    
    # sort according to the workflow
    # sorted_motives is a list of list
    sorted_names, sorted_motives = sort_lists(names_list, motives_list)
    
    for i in range(len(sorted_names)):
        sample = samples[i]
        sorted_name = sorted_names[i]
        sorted_motives_per_sample = sorted_motives[i]
        
        print("sorted_motives: " + str(sorted_motives_per_sample))
        sorted_motives_per_sample.insert(0, sample)
        sorted_motives_per_sample.insert(1, sorted_name)
        
        w.writerow(sorted_motives_per_sample)

def split_csv(dict_count):
    """
    Parameters
    ----------
    dict_count : dictionary of tuple where each element is a list.
    
        key = sample name
        1st value of the tuple = read type (class)
        2nd value of the tuple = motif orders (1st-18th) in the form of list of list
        
    Returns
    -------
    returns nothing but it calls the writing function for each sample. (i.e. write_out)
    """       
    for key, val in dict_count.items():
        sample_name = key
        # class of reads (list)
        names = val[0]
        # motives is a list of list where each element is a list of ordered motives.
        motives = val[1]
        if len(names) == len(motives):
            # creating a specific sample name list with same length 
            sample_names = [sample_name for i in range(len(names))]
            write_out(sample_names, names, motives)
            print('successfully generated split csv file for 1 sample')
        else:
            print('lengths do not match')

def get_inputs():
    
    parser = argparse.ArgumentParser(description="separate motif orders by sample")
    parser.add_argument('--total_csv_input', dest = 'total_csv_input',
                        required = True,
                        help = 'total_merged_motif_orders.csv')
    
    args = parser.parse_args()
    
    total_csv_input = args.total_csv_input
    return total_csv_input
            
def run_process():
    
    input_csv = get_inputs()
    print('successfully got the input')
    
    counts_dict = read_counts_csv(input_csv)
    print('successfully generated counts_dict')
    
    split_csv(counts_dict)
    print('successfully split all csvs according to sample names')

if __name__ == "__main__":
    run_process() 
    print('success')

