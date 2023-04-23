# -*- coding: utf-8 -*-
"""

Created on Wed Dec 15 14:41:02 2021

@author: Youngbin Moon (y.moon@unibas.ch)
"""

import pysam
import matplotlib.pyplot as plt
import numpy as np
import math
import argparse
import csv

"""
Aim: get a span of a cluster with same CB + UB tags + same chromosome or CB + YB tags + same chromosome (log10 scale).

YB tag: when you split a cluster of CB and UB tags into several sub-clusters because the span of the cluster is too large, 
you assign different UB tags to each sub-clusters and name it as YB tag.

YB = UB + k where k means kth subcluster.
 
"""

def write_out_to_csv(data, file_name):
    """
    Parameters
    ----------
    data : list
        a list of cluster span of reads with same CB, UB and chromosome (log10 scale).
    
    file_name : string
        output span csv file name.

    Returns
    -------
    returns nothing but saves computed log10 spans in the csv file.
    """
    w = csv.writer(open(file_name, 'a'))
    for row in data:
        w.writerow([row])
        
def get_distal_pos (cigar, rev, list_c):
    """
    Parameters
    ----------
    cigar : string
        a string that explains how a read is mapped.
        e.g. 4S76M151N21M, 38M1D48M795N12M and 73M25S.
        
        M: 0, alignment match (can be a sequence match or mismatch)
        I: 1, insertion to the reference
        D: 2, deletion from the reference
        N: 3, skipped region from the reference
        S: 4, soft clipping (clipped sequences present in SEQ)
        H: 5, hard clipping (clipped sequences NOT present in SEQ)
        P: 6, padding (silent deletion from padded reference)
        =: 7, sequence match
        X: 8, sequence mismatch
    
    rev : bool
        True if a read maps to (-) strand of DNA.
        False if a read maps to (+) strand of DNA.
    
    list_c : list
        a list which is ['I', 'D', 'N', 'S', 'H', 'P', 'X']
    
    Returns
    -------
    distal : int
        length of distal part of the read that is mapped.
    """    
    distal = ''
    intermediate = []
    # it is checked that all reads end with either 'M' or 'S' in both direction.
    if rev == True:
        # cigar example: 4S76M151N21M
        # In this case: distal is 76. (4 is reset because distal = '' at S)
        for c in cigar:
            if c not in ['M', 'S']:
                distal += c
            elif c == 'S':
                distal = ''
            elif c == 'M':
                break
    if rev == False:
        # cigar example: 38M1D48M795N12M
        if cigar[-1] == 'M':
            # remove 'M' from the end of cigar string. 38M1D48M795N12
            cigar_m = cigar[:-1]
            # reverse cigar_m. i.e. 21N597M84D1M83
            for c in cigar_m[::-1]:
                # list_c = ['I', 'D', 'N', 'S', 'H', 'P', 'X']
                if c not in list_c:
                    # intermediate: 21
                    intermediate.append(c)
                else:
                    break
        # if you have soft-clipped region, you will have mapped region in the next vicinity
        # cigar example: 73M25S
        elif cigar[-1] == 'S':
            # remove 'S' from the end of cigar string
            # cigar_m: 73M25
            cigar_m = cigar[:-1]
            # if after soft clipped region other element in list_c appeared, it will cause an error.
            # but it didnt cause an error. Hence no problem.
            # cigar_m[::-1] = 52M37
            # in this case intermediate will be: 37 (52 is reset when you see 'M'. This is because 25bp are softclipped. i.e. 25S)
            for c in cigar_m[::-1]:
                if c != 'M' and c not in list_c:
                    intermediate.append(c)
                elif c == 'M':
                    intermediate = []
                elif c in list_c:
                    break
        # if cigar: 38M1D48M795N12M, then interval: 21, distal: 12.
        # if cigar: 73M25S, then interval: 37, distal: 73.
        for elem in intermediate[::-1]:
            distal += elem

    distal = int(distal)
    return distal

def get_distribution(unique_pairs, reads_per_ub):
    """
    Parameters
    ----------
    unique_pairs : list of tuple
        It is an unique pair of a cell barcode (CB), UMI tag (UB or YB) and a chromosome
       
    reads_per_ub : dictionary of dictionary.
        first key: CB, UB and chromosome
        second key: 'reads', 'starts', 'ends' ('starts and 'ends' are extended from original input)
        value: 
            if second key is 'reads': reads that have same CB, UB and chromosome
            if second key is 'starts': start positions of reads that have same CB, UB and chromosome
            if second key is 'ends': end positions of reads that have same CB, UB and chromosome
    Returns
    -------
    spans : list
        a list of cluster span of reads with same CB and UMI tags (log10 scale).

    """
    spans = []
    list_c = ['I', 'D', 'N', 'S', 'H', 'P', 'X']
    for unique_pair in unique_pairs:
        # get a list of reads with the same CB and (UB or YB) tags and same chromosome.
        read_list = reads_per_ub[unique_pair]['reads']
        start_positions = []
        end_positions = []
        reads_per_ub[unique_pair]["start"] = []
        reads_per_ub[unique_pair]["end"] = []
        
        for read in read_list:
            cigar = read.cigarstring
            
            # if a read does not span introns, save start and end positions as it is.
            if 'N' not in cigar:
                
                alignedRefPositions = read.get_reference_positions()
                refStart = alignedRefPositions[0]
                refEnd = alignedRefPositions[-1]
                start_positions.append(refStart)
                end_positions.append(refEnd)
            
            # if a read spans introns:
            else:
                
                alignedRefPositions = read.get_reference_positions()
                refStart = alignedRefPositions[0]
                refEnd = alignedRefPositions[-1]
                rev = read.is_reverse
                
                # if a read spans introns, use distal part of read only
                distal = get_distal_pos(cigar, rev, list_c)
                
                if rev == True:
                    start = refStart
                    end = refStart + distal -1
                    start_positions.append(start)
                    end_positions.append(end)
                    
                elif rev == False:
                    start = refEnd - distal + 1
                    end = refEnd
                    start_positions.append(start)
                    end_positions.append(end)    
                    
        if len(start_positions) != 0 and len(end_positions) != 0:
            reads_per_ub[unique_pair]["start"] = start_positions
            reads_per_ub[unique_pair]["end"] = end_positions
            span = max(end_positions) - min(start_positions)
                        
            spans.append(math.log10(span + 1))
            
    return spans

def get_all_reads_per_umi (sam, type_input):
    """
    Parameters
    ----------
    sam : BAM file
        input bam file.
    
    type_input : string
        type of input which can be either raw bam file (raw), intermediate bam file of our deduplication process(split)
        or final output of our deduplication process (dedup).
        
        raw: raw bam file
        split: bam file where where it splits a cluster of CB and UB tags into several sub-clusters in case span is too large.
        dedup: splitted & deduplicated bam file.

    Returns
    -------
    unique_pairs : list of tuple
        It is an unique pair of a cell barcode (CB), UMI tag (UB or YB) and a chromosome
        UB = UMI tag before splitting a cluster into several sub-clusters. (for raw)
        YB = UMI tag after splitting a cluster into several sub-clusters. (for splitted and dedup)
       
    reads_per_ub : dictionary of dictionary.
        first key: CB, UMI and chromosome.
        second key: 'reads'
        value: reads that have same CB, UMI and chromosome.

    """
    reads_per_ub = {}
    unique_pairs=set()
    
    for read in sam.fetch():
        chrom = read.reference_name
        if type_input == "dedup" or type_input == "split":
            if read.has_tag('YB') and read.has_tag('CB'):
                umi = read.get_tag('YB')
                cb = read.get_tag('CB')
                
                if (cb, umi, chrom) not in reads_per_ub.keys():
                    reads_per_ub[(cb, umi, chrom)] = {}
                    reads_per_ub[(cb, umi, chrom)]['reads'] = []
                reads_per_ub[(cb, umi, chrom)]['reads'].append(read)
                # unique umis
                unique_pairs.add((cb, umi, chrom))
                
        elif type_input == "raw":
            if read.has_tag('UR') and read.has_tag('CB'):
                umi = read.get_tag('UR')
                cb = read.get_tag('CB')
                
                if (cb, umi, chrom) not in reads_per_ub.keys():
                    reads_per_ub[(cb, umi, chrom)] = {}
                    reads_per_ub[(cb, umi, chrom)]['reads'] = []
                reads_per_ub[(cb, umi, chrom)]['reads'].append(read)
                #unique umis
                unique_pairs.add((cb, umi, chrom))
                
    unique_pairs = list(unique_pairs)
    
    return unique_pairs, reads_per_ub
        
def get_input():

    parser = argparse.ArgumentParser(description="compute span")
    parser.add_argument('--input', dest = 'input',
                        required = True,
                        help = 'bam file')

    parser.add_argument('--output', dest = 'output',
                        required = True,
                        help = 'output span csv')   

    parser.add_argument('--action_type', dest = 'action_type',
                        required = True,
                        help = 'raw or dedup or split')
       
    args = parser.parse_args()
    
    input_dir = args.input
    input_bam = pysam.AlignmentFile(input_dir, "rb")    
    output_csv = args.output
    action_type = args.action_type
    
    return input_bam, output_csv, action_type

def run_process():
    
    input_bam, output_csv, action_type = get_input()
    print("successfully initialized")
    
    # action_type = 'raw', 'split' or 'dedup'    
    Unique_pairs, Reads_per_ub = get_all_reads_per_umi (input_bam, action_type)
    print("successfully done get_all_reads_per_umi")
        
    logspans = get_distribution(Unique_pairs, Reads_per_ub)
    print("successfully got final spans")
    
    write_out_to_csv(logspans, output_csv)
    print("successfully got output csv file")    
    
if __name__ == "__main__":
    run_process()
    print("success")
