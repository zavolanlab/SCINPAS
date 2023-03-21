# -*- coding: utf-8 -*-
"""

Created on Mon Jan 10 10:26:42 2022

@author: Youngbin Moon (y.moon@unibas.ch)
"""

import pysam
import csv
import argparse






"""
Aim: Get csv output file that contains number of read having softclipped region(direction does not matter) in a specific sample and a specific bam type.
Only used for deduplicated bam file. because we want total number of softclipped reads.

"""
def write_out(row_data, o_file):
    """
    Parameters
    ----------
    row_data : list
        A list that contains sample, bam type and count. This will be 1 row in the csv.
        Later, this will be merged across all samples and all bam types.
    
    o_file : string
        output csv file name    

    Returns
    -------
    returns nothing but saves the output (sample, bam type and count) in csv format.
    """
    
    w = csv.writer(open(o_file, "w"))
    w.writerow(row_data)

def count_softclipped_reads(samFile):
    """
    Parameters
    ----------
    samFile : bam file
        input bam file of a specific class of reads

    Returns
    -------
    softclipped_count : int
        The number of reads that have softclipped reads in that bam file.
        Does not matter if softclipped direction does not match with direction of a read.
    """
    softclipped_count = 0
    for read in samFile.fetch():
        tuples = read.cigartuples
        left_end = tuples[0]
        right_end = tuples[-1]
        
        if left_end[0] == 4 or right_end[0] == 4:
            softclipped_count += 1
        
        else:
            print('not softclipped read')
            continue
        
    return softclipped_count

def get_count_data(sam, sample, bamType):
    """
    Parameters
    ----------
    sam : bam file
        input bam file of a specific class of reads
    
    sample : string
        sample name.
        
    bamType : string
        specific class of reads.
        bam type can be original, dedup, nonpolyA, all polyA, 
        polyA mapping to terminal exon, exonic, intronic and intergenic.
            
    Returns
    -------
    row_data : list
        A list that contains sample, bam type and count. This will be 1 row in the csv.
        Later, this will be merged across all samples and all bam types. 
    """
    count = count_softclipped_reads(sam)
        
    row_data = [sample, bamType, count]
    return row_data
    
def run ():    
    parser = argparse.ArgumentParser(description="get counts for all types of reads")
    parser.add_argument('--bam_input', dest = 'bam_input',
                        required = True,
                        help = 'bam file')
    
    parser.add_argument('--csv_output', dest = 'csv_output',
                        required = True,
                        help = 'output csv name')
  
    parser.add_argument('--sample_name', dest = 'sample_name',
                        required = True,
                        help = 'sample_name')
    
    parser.add_argument('--bam_type', dest = 'bam_type',
                        required = True,
                        help = 'type of the bam file') 
                 
    args = parser.parse_args()
    
    bam_dir = args.bam_input
    bam = pysam.AlignmentFile(bam_dir, "rb")
    out_dir = args.csv_output
    sample_name = args.sample_name
    bam_type = args.bam_type

    return (bam, out_dir, sample_name, bam_type)
    
def run_process():
    (bam, out_dir, sample_name, bam_type) = run()
    print('successfully initialized')

    modified_sample_name =  sample_name    
    # UmiRaw is the raw data used for negative control.
    # UmiDedup is the UMI-tools deduplicated data that is used for negative control.
    # because these are negative controls and they are bam file joined at different stages of the pipeline, the name has to be unified.     
    if "UmiDedup" in sample_name:
        modified_sample_name = "negative_control" + "(" + sample_name + ")"
        
    if "UmiRaw" in sample_name:
        modified_sample_name = "negative_control" + "(" + sample_name.replace("UmiRaw", "UmiDedup") + ")"
        
    row = get_count_data(bam, modified_sample_name, bam_type)
    print('successfully got dictionary')
    
    write_out(row, out_dir)
    print('successfully saved result')

if __name__ == "__main__":
    run_process()
    print('success')