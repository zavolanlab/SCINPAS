# -*- coding: utf-8 -*-
"""

Created on Mon Jan 10 10:26:42 2022

@author: Youngbin Moon (y.moon@unibas.ch)
"""

import pysam
import csv
import argparse

"""
Aim: Get csv output file that contains number of read counts in a specific sample and a specific bam type.
bam type can be original, dedup, nonpolyA, all polyA, polyA mapping to terminal exon, exonic, intronic and intergenic.
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
    
def count_total_reads (samFile):
    """
    Parameters
    ----------
    samFile : bam file
        input bam file of a specific class of reads

    Returns
    -------
    count : int
        The number of reads in that bam file.
    """
    count = 0
    for read in samFile.fetch():
        count += 1
    return count

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
    count = count_total_reads(sam)
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
    
    row = get_count_data(bam, sample_name, bam_type)
    print('successfully got dictionary')
    
    write_out(row, out_dir)
    print('successfully saved result')

if __name__ == "__main__":
    run_process()
    print('success')