# -*- coding: utf-8 -*-
"""

Created on Sun Jan 23 07:02:13 2022

@author: Youngbin Moon (y.moon@unibas.ch)

modified on Apr 20 2026
@author: Aleksei Mironov (aleksei.mironov@unibas.ch)
"""

import pysam
import argparse
import re
import statistics
import numpy as np
import pandas as pd
import time
import tracemalloc
"""
Aim : From deduplicated bam file of a given sample, collect all polyA reads
and save it to a new bam file.
"""    
def write_output(final_reads, o_name, o_mode, bam):
    """
    Parameters
    ----------
    final_reads : list
        a list of reads that that only have polyA reads.
    
    out_name : string  
        output file name.
    
    out_mode : string
        output format. 'wb' refers to wrtie in bam file.
    
    bam : bam file
        the original bam file (input)
        This is to use original bam file as a template.
        
    Returns
    -------        
    returns nothing but saves polyA reads in the bam format.
    """              
    outfile = pysam.AlignmentFile(o_name, o_mode, template=bam)
    for read in final_reads:
        outfile.write(read)

def get_median_phred_polyA(r, is_fc, tag_orig_cs="XO", tag_fixed_cs="XF"):
    rev = r.is_reverse
    tuples = r.cigartuples
    left_end = tuples[0]
    right_end = tuples[-1]
    
    if rev == True:
        if not is_fc:
            qualities = r.get_forward_qualities()
            len_pA = left_end[1]
            len_read = len(qualities)
            mapped_qualities = qualities[0 : len_read - len_pA]
            softclipped_qualities = qualities[len_read - len_pA : len_read]
        elif is_fc:
            qualities = r.get_forward_qualities()
            OCS = r.get_tag(tag_orig_cs)
            FCS = r.get_tag(tag_fixed_cs)
            difference = OCS - FCS
            assert(difference >= 0)
            len_pA = left_end[1] - difference
            len_read = len(qualities)
            mapped_qualities = qualities[0 : len_read - len_pA]
            softclipped_qualities = qualities[len_read - len_pA : len_read]
    else:
        if not is_fc:
            qualities = r.get_forward_qualities()
            len_pA = right_end[1]
            len_read = len(qualities)
            mapped_qualities = qualities[0 : len_read - len_pA]
            softclipped_qualities = qualities[len_read - len_pA : len_read]  
        elif is_fc:
            qualities = r.get_forward_qualities()
            OCS = r.get_tag(tag_orig_cs)
            FCS = r.get_tag(tag_fixed_cs)
            difference = FCS - OCS
            assert(difference >= 0)
            len_pA = right_end[1] - difference
            len_read = len(qualities)
            mapped_qualities = qualities[0 : len_read - len_pA]
            softclipped_qualities = qualities[len_read - len_pA : len_read]
    
    phred_median_mapped = statistics.median(mapped_qualities)
    phred_median_softclipped = statistics.median(softclipped_qualities)
    return phred_median_mapped, phred_median_softclipped
        
def count_A(sub_sequence):
    """
    Parameters
    ----------
    sub_sequence : string
        a soft clipped region of a read from 5' -> 3' (left to right).
        This is original sequence. (regardless of the read direction)
    
    Returns
    -------        
    number_A : int
        The number of "A"s in the sub_sequence.
    """       
    list_count = [1 if elem == 'A' else 0 for elem in sub_sequence]
    number_A = sum(list_count)
    return number_A

"""
get_forward_sequence(): 

This function is used to get real transcript sequence
if it is a read mapping to (-) strand, you have to reverse complement a read to get "original real transcript sequence"
which is handled by the get_forward_sequence function. This is because BAM file always saves reads in 5' -> 3' direction.

query_sequence():

This function is used when you want to compare a read sequence (+ or -) to reference genome(+)
because you have to use sequence that is saved in BAM file.
(You can still compare reads mapping to - strand, because a read is always saved in 5->3')

Summary:
If you want original sequence in 5'->3', use read.get_forward_sequence()
If you want reverse complemented sequence so that you can compare it to genome 5'-> 3', use read.query_sequence
Need to consdier direction as well
"""
def check_polyA(read, left_end, right_end, percentage_threshold, length_threshold, use_fc, tag_orig_cs="XO", tag_fixed_cs="XF"):
    """
    Parameters
    ----------
    read : pysam object
        a deduplicated read of interest.

    left_end : tuple
        contains whether a read has a soft clipped region and if yes, how long?
        This is for a read mapping to (-) strand of the genome.
        
        if left_end[0] == 4 -> there is a soft clipped region in the left side of a read.
        left_end[1] -> gives you length of the soft clipped region in the left side of a read.
        
    right_end : tuple
        contains whether a read has a soft clipped region and if yes, how long?
        This is for a read mapping to (+) strand of the genome.   
        
        if right_end[0] == 4 -> there is a soft clipped region in the right side of a read.
        right_end[1] -> gives you length of the soft clipped region in the right side of a read.
        
    percentage_threshold : int
        a percentage threshold for a deduplicated read to be considered as polyA read.
        percentage of "A" nucleotide in the softclipped region has to be over this threshold
        in order to be considered as polyA reads.
        
    length_threshold : int
        a length threshold for a deduplicated read to be considered as polyA read.
        The number of "A" nucleotide in the softclipped region has to be over this threshold
        in order to be considered as polyA reads.
        (Note: it does not have to be consecutive number of "A"s)
    
    use_FC : bool
        whether you use fixed softclipped region or original soft clipped region.
        True if you want to use fixed softclipped region.
        False if you do not want to use fixed softclipped region.
    tag_orig_cs : string
        the name of the custom tag that contains the original cleavage site position.
    tag_fixed_cs : string
        the name of the custom tag that contains the fixed cleavage site position.        
    
    Returns
    -------        
    True if a read has polyA tail.
    False if a read does not have a polyA tail.
    """      
    full_sequence = read.get_forward_sequence()
    rev = read.is_reverse

    if rev == True and left_end[0] == 4:
        if not use_fc:
            potential_polyA = full_sequence[len(full_sequence) - left_end[1] : len(full_sequence)]
            len_pA = left_end[1]
        elif use_fc:
            OCS = read.get_tag(tag_orig_cs)
            FCS = read.get_tag(tag_fixed_cs)
            difference = OCS - FCS
            assert(difference >= 0)
            potential_polyA = full_sequence[len(full_sequence) - left_end[1] + difference : len(full_sequence)]
            len_pA = left_end[1] - difference
            
        num_A = count_A(potential_polyA)
        percentage_A = (num_A/len_pA)*100
        if len_pA <= 5:
            percentage_threshold = 100
        if len_pA >= length_threshold and percentage_A >= percentage_threshold:
            return True
        else:
            return False
            
    elif rev == False and right_end[0] == 4:
        if not use_fc:
            potential_polyA = full_sequence[len(full_sequence) - right_end[1] : len(full_sequence)]
            len_pA = right_end[1]
        elif use_fc:
            OCS = read.get_tag(tag_orig_cs)
            FCS = read.get_tag(tag_fixed_cs)
            difference = FCS - OCS
            assert(difference >= 0)
            potential_polyA = full_sequence[len(full_sequence) - right_end[1] + difference : len(full_sequence)]
            len_pA = right_end[1] - difference
            
        num_A = count_A(potential_polyA)
        percentage_A = (num_A/len_pA)*100
        if len_pA <= 5:
            percentage_threshold = 100
        if len_pA >= length_threshold and percentage_A >= percentage_threshold:
            return True
        else:
            return False 
    else:
        return False

def find_polyA_seq(sam, percentage_threshold, length_threshold, fasta, use_FC, min_phred=30.0, 
                   tag_phred_mapped="", tag_phred_softclipped="", tag_orig_cs="XO", tag_fixed_cs="XF"):
    """
    Parameters
    ----------
    sam : bam file
        a bam file that contains deduplicated reads in which their alignment is fixed
        and saved in 'FC' tag.
    
    percentage_threshold : int
        a percentage threshold for a deduplicated read to be considered as polyA read.
        percentage of "A" nucleotide in the softclipped region has to be over this threshold
        in order to be considered as polyA reads.
        
    length_threshold : int
        a length threshold for a deduplicated read to be considered as polyA read.
        The number of "A" nucleotide in the softclipped region has to be over this threshold
        in order to be considered as polyA reads.
        (Note: it does not have to be consecutive number of "A"s)
    
    fasta : a fasta flie
        contains the reference genome sequence.
        
    use_FC : bool
        whether you use fixed softclipped region or original soft clipped region.
        True if you want to use fixed softclipped region.
        False if you do not want to use fixed softclipped region.

    min_phred : int
        the minimum PHRED quality threshold for a deduplicated read to be considered as high-quality polyA read.

    tag_phred_mapped : string or None
        if not None, the name of the custom tag to store the median PHRED quality of the mapped part of a read.

    tag_phred_softclipped : string or None
        if not None, the name of the custom tag to store the median PHRED quality of the softclipped part of a read.
    tag_orig_cs : string
        the name of the custom tag that contains the original cleavage site position.
    tag_fixed_cs : string
        the name of the custom tag that contains the fixed cleavage site position.
    Returns
    -------        
    polyA_reads : list
        a list of reads that have polyA tail.
    
    non_polyA_reads : list
        a list of reads that do not have polyA tail.
    
    low_quality_pA_reads : list
        a list of reads that have polyA tail but do not meet the minimum PHRED quality threshold.
    """       
    polyA_reads = []
    non_polyA_reads = []
    low_quality_pA_reads = []
    
    for read in sam.fetch():
        tuples = read.cigartuples
        left_end = tuples[0]
        right_end = tuples[-1]
        
        if left_end[0] == 4 or right_end[0] == 4:           
            is_polyA = check_polyA(read, left_end, right_end, percentage_threshold, length_threshold, use_FC, tag_orig_cs, tag_fixed_cs)
            if is_polyA:
                phred_median_mapped, phred_median_softclipped = get_median_phred_polyA(read, use_FC, tag_orig_cs, tag_fixed_cs)
                
                if tag_phred_mapped:
                    read.set_tag(tag_phred_mapped, float(phred_median_mapped), value_type='f')
                if tag_phred_softclipped:
                    read.set_tag(tag_phred_softclipped, float(phred_median_softclipped), value_type='f')

                if phred_median_mapped > min_phred and phred_median_softclipped > min_phred:    
                    polyA_reads.append(read)
                else:
                    low_quality_pA_reads.append(read)
                    
            elif not is_polyA:
                non_polyA_reads.append(read)
        elif left_end[0] != 4 and right_end[0] != 4:
            non_polyA_reads.append(read)
            
    return polyA_reads, non_polyA_reads, low_quality_pA_reads

def get_all_polyA_input():
    parser = argparse.ArgumentParser(description="get filtered polyA reads")
    parser.add_argument('--bam_input', dest='bam_input', required=True)
    parser.add_argument('--o_polyA', dest='o_polyA', required=True)
    parser.add_argument('--o_nonpolyA', dest='o_nonpolyA', required=True)      
    parser.add_argument('--fasta', dest='fasta', required=True)    
    parser.add_argument('--percentage_threshold', type=int, dest='percentage_threshold', required=True)
    parser.add_argument('--length_threshold', type=int, dest='length_threshold', required=True)  
    parser.add_argument('--use_fc', type=int, dest='use_fc', required=True)   
    parser.add_argument('--o_low_q_polyA', dest='o_low_q_polyA', required=False)
    parser.add_argument('--exact_out', action='store_true')
    parser.add_argument('--stats_tsv', dest='stats_tsv', required=False, help="Output TSV with polyA statistics")
    parser.add_argument('--sample_id', dest='sample_id', required=False, help="Sample ID for stats tracking")
    parser.add_argument('--min_phred', type=float, default=30.0, help="Minimum median PHRED score to pass")
    parser.add_argument('--tag_phred_mapped', type=str, default="", help="Custom SAM tag for mapped PHRED (e.g., ZM)")
    parser.add_argument('--tag_phred_softclipped', type=str, default="", help="Custom SAM tag for softclipped PHRED (e.g., ZC)")
    parser.add_argument('--tag_orig_cs', dest='tag_orig_cs', default="XO")
    parser.add_argument('--tag_fixed_cs', dest='tag_fixed_cs', default="XF")
    
    args = parser.parse_args()
    
    bamFile = args.bam_input
    sam = pysam.AlignmentFile(bamFile, "rb")
    fasta_file = pysam.FastaFile(args.fasta)
    if args.exact_out:
        number = ""
    else:
        number = re.split('_chr', bamFile)[1].split('_')[0]
        
    return sam, "wb", fasta_file, args.o_polyA, args.o_nonpolyA, args.o_low_q_polyA, \
            args.percentage_threshold, args.length_threshold, bool(args.use_fc), number, \
            args.exact_out, args.stats_tsv, args.sample_id, \
            args.min_phred, args.tag_phred_mapped, args.tag_phred_softclipped, args.tag_orig_cs, args.tag_fixed_cs

def run_process():
    start = time.time()
    tracemalloc.start()
    
    sam, out_mode, fasta_file, out_polyA, out_non_polyA, o_low_q_polyA, \
    percentage_threshold, length_threshold, use_fc, number, exact_out, \
    stats_tsv, sample_id, min_phred, tag_phred_mapped, tag_phred_softclipped, \
    tag_orig_cs, tag_fixed_cs = get_all_polyA_input()
    
    polyA_reads, non_polyA_reads, low_quality_pA_reads = find_polyA_seq(
        sam, percentage_threshold, length_threshold, fasta_file, use_fc, 
        min_phred, tag_phred_mapped, tag_phred_softclipped, tag_orig_cs, tag_fixed_cs
    )
    print('successfully got all polyA reads')
    
    if exact_out:
        corrected_out_polyA = out_polyA
        corrected_out_non_polyA = out_non_polyA
        corrected_out_low_q_polyA = o_low_q_polyA if o_low_q_polyA else "low_q.bam"
    else:
        corrected_out_polyA = out_polyA.split('.')[0] + '_chr' + str(number) + '.bam'
        corrected_out_non_polyA = out_non_polyA.split('.')[0] + '_chr' + str(number) + '.bam'
        corrected_out_low_q_polyA = out_polyA.split('.')[0] + '_lowQualityChrom' + str(number) + '.bam'
    
    write_output(polyA_reads, corrected_out_polyA, out_mode, sam)
    write_output(non_polyA_reads, corrected_out_non_polyA, out_mode, sam)
    write_output(low_quality_pA_reads, corrected_out_low_q_polyA, out_mode, sam)
    
    if stats_tsv:
        import csv
        import os
        with open(stats_tsv, 'w', newline='') as tsv_file:
            writer = csv.writer(tsv_file, delimiter='\t')
            writer.writerow(["sample_id", "chunk_filename", "polyA_alignments", "non_polyA_alignments", "low_q_polyA_alignments"])
            chunk_name = os.path.basename(sam.filename.decode() if isinstance(sam.filename, bytes) else sam.filename)
            s_id = sample_id if sample_id else "unknown"
            writer.writerow([s_id, chunk_name, len(polyA_reads), len(non_polyA_reads), len(low_quality_pA_reads)])
        
    print('elapsed time: ' + str(time.time() - start))
    print('memory usage is: ' + str(tracemalloc.get_traced_memory()))
    tracemalloc.stop()
    
if __name__ == "__main__":    
    run_process()
    print('success')