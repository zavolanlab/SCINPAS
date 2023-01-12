# -*- coding: utf-8 -*-
"""

Created on Wed Jul 20 17:26:26 2022

@author: Youngbin Moon (y.moon@unibas.ch)
"""

import pysam
import argparse
from collections import Counter
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import math
import statistics
"""
Aim : From deduplicated bam file of a sample, create a new deduplicated bam file
in which their alignments are fixed.

i.e. fixed cleavage site is saved in a new 'FC' tag per each read.
"""
def write_output(final_reads, out_name, out_mode, sam):
    """
    Parameters
    ----------
    final_reads : list
        a list of deduplicated reads in which their alignments are fixed.
        (cleavage sites are fixed and saved in the 'FC' tag)
    
    out_name : string  
        output file name.
    
    out_mode : string
        output format. 'wb' refers to wrtie in bam file.
    
    bam : bam file
        the original bam file (input)
        This is to use original bam file as a template.
        
    Returns
    -------        
    returns nothing but saves fixed deduplicated reads in the bam format.
    """       
    outfile = pysam.AlignmentFile(out_name, out_mode, template=sam)
    for read in final_reads:
        outfile.write(read)
        
def find_positions_to_fix (genome_seq, read_Softclipped, tolerance_threshold, rev):
    """
    Parameters
    ----------    
    genome_seq : string
        a genome sequence that corresponds to the soft-clipped region of a read
        
    read_Softclipped : string
        a soft clipped sequence of a read.
    
    tolerance_threshold : int
        a maximum number of allowed mismatches between a read and genome during amendment of the softclipped region.
        As you go through the soft-clipped region, if the number of mismatches between a read and genome exceeds this
        threshold, we stop fixing the soft-clipped region.
        
        Threshold changes depending on the length of the soft-clipped region.
    
    rev : bool
        True if a read maps to (-) strand of the genome.
        False if a read maps to (+) strand of the genome.
        
    Returns
    -------        
    n_proceed : int
        From old cleavage site, how many base pairs do we have to go? in order to get a fixed cleavage site.
    """      
    assert len(genome_seq) == len(read_Softclipped)
    n_mismatch = 0
    previous_move = 'match'
    # last position will have at least 1 mismatch.
    # Hence consecutive_mismatch >= 1
    consecutive_mismatch = 1
    end_pos = 0
    
    # if a read maps to - strand, you have to start editing from the right end
    # Hence you need to reverse the genomic sequence and softclipped region
    # if you use reverse(), original list is reversed
    if rev == True:
        genome_seq = [elem for elem in genome_seq]
        read_Softclipped = [elem for elem in read_Softclipped]
        genome_seq.reverse()
        read_Softclipped.reverse()

    # for each read find a new start of the softclipped region
    for i in range(len(genome_seq)):
        genome_nt = genome_seq[i]
        read_nt = read_Softclipped[i]
        
        if read_nt != genome_nt:
            n_mismatch += 1
            if previous_move == 'mismatch':
                consecutive_mismatch += 1
            previous_move = 'mismatch'
            
            if n_mismatch > tolerance_threshold or i == len(genome_seq) - 1:
                end_pos = i
                break
            
            continue
        
        elif read_nt == genome_nt:
            # consecutive_mismatch == 1 means you have mismatch = 1. You incremented in advance.
            # Because when you have mismatch, consecutive mismatch only increases when you had mismatch just before that position.
            # It also function as reseting consecutive_mismatch in case you have a match.
            consecutive_mismatch = 1
            previous_move = 'match'
            
            if n_mismatch > tolerance_threshold or i == len(genome_seq) - 1:
                end_pos = i
                break
            
            continue
        
    n_proceed = end_pos - consecutive_mismatch

    return n_proceed  

def extract_sequences(Read, fasta):
    """
    Parameters
    ----------    
    Read : a pysam object
        current read from deduplicated bam file.
        
    fasta : a fasta flie
        contains the reference genome sequence.
        
    Returns
    -------        
    genome_sequence : string
        a genome sequence that corresponds to the soft-clipped region of a read
        
    read_softclipped : string
        a soft clipped sequence of a read.
    
    threshold : int
        a maximum number of allowed mismatches between a read and genome during amendment of the softclipped region.
        As you go through the soft-clipped region, if the number of mismatches between a read and genome exceeds this
        threshold, we stop fixing the soft-clipped region.
        
        Threshold changes depending on the length of the soft-clipped region.
    """     
    rev = Read.is_reverse
    chromosome = Read.reference_name
    alignedRefPositions = Read.get_reference_positions()
    refStart = alignedRefPositions[0]
    refEnd = alignedRefPositions[-1]
    tuples = Read.cigartuples
    left_end = tuples[0]
    right_end = tuples[-1]
        
    # reads that are checked should have polyA tail
    if rev == True and left_end[0] == 4:
        # if a read is mapped reversly,  reads should be reverse complemented
        # so that you can compare read vs genome
        # for that, you need to use query_sequence (but U is converted into T)
        full_sequence = Read.query_sequence
        
        # start_point should be smaller than the end_point to extract genomic sequence
        end_point = int(refStart) - 1
        start_point = int(refStart)  - left_end[1]
        
        genome_sequence = fasta.fetch(reference = chromosome, start = start_point, end = end_point + 1)
        read_softclipped = full_sequence[0 : left_end[1]]
        print("read softclipped (-) is: " + str(read_softclipped))
        print("genomic_sequence is: " + str(genome_sequence))
        threshold = max(left_end[1]/10, 2)
        
    elif rev == False and right_end[0] == 4:
        # get original transcript sequence (but U is converted into T)
        full_sequence = Read.get_forward_sequence()
        start_point = int(refEnd) + 1
        end_point = int(refEnd) + right_end[1]
        
        genome_sequence = fasta.fetch(reference = chromosome, start = start_point, end = end_point + 1)
        read_softclipped = full_sequence[len(full_sequence) - right_end[1] : len(full_sequence)]
        print("read softclipped is (+) : " + str(read_softclipped))
        print("genomic_sequence is: " + str(genome_sequence))
        threshold = max(right_end[1]/10, 2)
        
    return genome_sequence, read_softclipped, threshold

def fix_soft_clipped(sam, fasta_file):
    """
    Parameters
    ----------    
    sam : bam file
        a bam file containing all deduplicated reads from a sample.
    
    fasta_file : a fasta flie
        contains the reference genome sequence.
        
    Returns
    -------        
    changed_reads : list
        a list of deduplicated reads in which their alignments are fixed.
        (cleavage sites are fixed and saved in the 'FC' tag)
    """       
    changed_reads = []
    for read in sam.fetch():
        rev = read.is_reverse
        alignedRefPositions = read.get_reference_positions()
        refStart = alignedRefPositions[0]
        refEnd = alignedRefPositions[-1]
        tuples = read.cigartuples
        left_end = tuples[0]
        right_end = tuples[-1]
                        
        if rev == True and left_end[0] == 4:
            genome_sequence, read_softclipped, threshold = extract_sequences(read, fasta_file)          
            n_proceed = find_positions_to_fix(genome_sequence, read_softclipped, threshold, rev)
            
            cleavage_site = int(refStart)
            fixed_cleavage_site = cleavage_site - n_proceed - 1
        
        # if there is no softclipped region, use original cleavage site.
        # because there isnt a thing to extend mapped region.
        elif rev == True and left_end[0] != 4:            
            cleavage_site = int(refStart)
            fixed_cleavage_site = cleavage_site
            
        elif rev == False and right_end[0] == 4:
            genome_sequence, read_softclipped, threshold = extract_sequences(read, fasta_file)          
            n_proceed = find_positions_to_fix(genome_sequence, read_softclipped, threshold, rev)
            
            cleavage_site = int(refEnd)
            fixed_cleavage_site = cleavage_site + n_proceed + 1
            
        # if there is no softclipped region, use original cleavage site.
        # because there isnt a thing to extend mapped region.        
        elif rev == False and right_end[0] != 4:          
            cleavage_site = int(refEnd)
            fixed_cleavage_site = cleavage_site 
            
        # save the original and fixed cleavage site into OC and FC tag respectively.     
        read.set_tag("OC", cleavage_site)
        read.set_tag("FC", fixed_cleavage_site)

        changed_reads.append(read)
        
    return changed_reads
            
def get_inputs():
    parser = argparse.ArgumentParser(description = "fix softclipped regions to get better cleavage sites" )
    
    parser.add_argument('--bam_file', dest = 'bam_file',
                        required = True,
                        help = 'bam_file that can be dedup, internal priming, annotated polyA reads')    
        
    parser.add_argument('--fasta', dest = 'fasta',
                        required = True,
                        help = 'fasta file')  
             
    parser.add_argument('--bam_out', dest = 'bam_out',
                        required = True,
                        help = 'output bam file')
     
    args = parser.parse_args()

    bamFile = args.bam_file  
    bam = pysam.AlignmentFile(bamFile, "rb")     
       
    fasta_dir = args.fasta
    fasta = pysam.FastaFile(fasta_dir)
            
    bam_out = args.bam_out
    
    return bam, fasta, bam_out
        
def run_process():
    
    bam, fasta, bam_out = get_inputs()
    print('successfully got inputs')
    
    changed_reads = fix_soft_clipped(bam, fasta)
    print('successfully added new tags')
    print('successfully got dictionary of new cleavage_sites')
    
    write_output(changed_reads, bam_out, "wb", bam)
    print('successfully wrote a new bam file')
    
if __name__ == "__main__":
    run_process()
    print('success')
    
    