# -*- coding: utf-8 -*-
"""

Created on Sun Jan  9 07:02:13 2022

@author: Youngbin Moon (y.moon@unibas.ch)
"""

import pysam
import argparse






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
def check_polyA(read, left_end, right_end, percentage_threshold, length_threshold, use_fc):
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
    
    Returns
    -------        
    True if a read has polyA tail.
    False if a read does not have a polyA tail.
    """      
    full_sequence = read.get_forward_sequence()
    rev = read.is_reverse

    # if read is mapped to negative strand, polyA tail should be on the left end
    if rev == True and left_end[0] == 4:
        if not use_fc:
            potential_polyA = full_sequence[len(full_sequence) - left_end[1] : len(full_sequence)]
            len_pA = left_end[1]
        
        elif use_fc:
            OCS = read.get_tag('XO')
            FCS = read.get_tag('XF')
            difference = OCS - FCS
            potential_polyA = full_sequence[len(full_sequence) - left_end[1] + difference : len(full_sequence)]
            # length of a softclipped region
            len_pA = left_end[1] - difference
            
        # you should use num_A not num_T because you use full_sequence rather than fasta.fetch()
        num_A = count_A(potential_polyA)
        percentage_A = (num_A/len_pA)*100
        
        # for softclipped length of <=6, you want to have 100% "A"        
        if len_pA <= 5:
            percentage_threshold = 100
                
        # decide whether a read is polyA or not
        if len_pA >= length_threshold and percentage_A >= percentage_threshold:
            return True
        
        else:
            return False
        
    # if right_end[0] == 4, it means soft clipped on the right side of a read.
    # right_end[1] gives how many bases are soft clipped on the right side.
    elif rev == False and right_end[0] == 4:
        if not use_fc:
            potential_polyA = full_sequence[len(full_sequence) - right_end[1] : len(full_sequence)]
            len_pA = right_end[1]
        
        elif use_fc:
            OCS = read.get_tag('XO')
            FCS = read.get_tag('XF')
            difference = FCS - OCS
            potential_polyA = full_sequence[len(full_sequence) - right_end[1] + difference : len(full_sequence)]
            # length of a softclipped region
            len_pA = right_end[1] - difference
            
        num_A = count_A(potential_polyA)
        percentage_A = (num_A/len_pA)*100
        
        # for softclipped length of <=6, you want to have 100% "A"  
        if len_pA <= 5:
            percentage_threshold = 100
        
        # decide whether a read is polyA or not
        if len_pA >= length_threshold and percentage_A >= percentage_threshold:
            return True
        
        else:
            return False 
        
    # a read does not have softclipped region
    # or it has soft clipped region but not at the right direction.
    else:
        return False

def find_polyA_seq(sam, percentage_threshold, length_threshold, fasta, use_FC):
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
    
    Returns
    -------        
    polyA_reads : list
        a list of reads that have polyA tail.
    
    non_polyA_reads : list
        a list of reads that do not have polyA tail.
    """       
    polyA_reads = []
    non_polyA_reads = []
    for read in sam.fetch():
        # list of tuples where each element is tuple.
        # assume soft clipp happens either left end, right end or both ends. it cant be softclipped in the middle.
        # first element of each tuple = type of cigar block
        # 2nd element of each tuple = length of each cigar block
        tuples = read.cigartuples
        left_end = tuples[0]
        right_end = tuples[-1]
        
        if left_end[0] == 4 or right_end[0] == 4:           
            is_polyA = check_polyA(read, left_end, right_end, percentage_threshold, length_threshold, use_FC)
            if is_polyA:
                polyA_reads.append(read)               
            # 3 possiblities not to have a polyA tail
            # 1) dont have soft clipped part at all
            # 2) soft clipped part in the wrong direction
            # 3) you have soft clipped part in the right direction but do not exceed threshold to become a polyA reads
            # case 2) and 3)
            elif not is_polyA:
                non_polyA_reads.append(read)
        # case 1)                
        elif left_end[0] != 4 and right_end[0] != 4:
            non_polyA_reads.append(read)
            
    return polyA_reads, non_polyA_reads

def get_all_polyA_input():
    
    parser = argparse.ArgumentParser(description="get full polyA reads")
    parser.add_argument('--bam_input', dest = 'bam_input',
                        required = True,
                        help = 'full deduplicated bam file in which the alignment is fixed')
    
    parser.add_argument('--o_polyA', dest = 'o_polyA',
                        required = True,
                        help = 'full poly A reads')

    parser.add_argument('--o_nonpolyA', dest = 'o_nonpolyA',
                        required = True,
                        help = 'full nonpoly A reads')      
    
    parser.add_argument('--fasta', dest = 'fasta',
                        required = True,
                        help = 'fasta file directory')    
    
    parser.add_argument('--percentage_threshold', type = int, dest = 'percentage_threshold',
                        required = True,
                        help = 'minimal percentage to become polyA')

    parser.add_argument('--length_threshold', type = int, dest = 'length_threshold',
                        required = True,
                        help = 'minimum number of A to become polyA')  
   
    parser.add_argument('--use_fc', type = int, dest = 'use_fc',
                        required = True,
                        help = 'whether use fixed cleavage site or not')   
      
    args = parser.parse_args()
    
    bamFile = args.bam_input
    sam = pysam.AlignmentFile(bamFile, "rb")
    
    out_mode = "wb"
    
    fasta_dir = args.fasta
    fasta_file = pysam.FastaFile(fasta_dir)
    
    out_polyA = args.o_polyA
    out_non_polyA = args.o_nonpolyA
    
    percentage_threshold = args.percentage_threshold
    length_threshold = args.length_threshold
    
    use_fc = bool(int(args.use_fc))
    
    return sam, out_mode, fasta_file, out_polyA, out_non_polyA,\
            percentage_threshold, length_threshold, use_fc

def run_process():
    
    sam, out_mode, fasta_file, out_polyA, out_non_polyA,\
    percentage_threshold, length_threshold, use_fc = get_all_polyA_input()
    print('successfully got inputs')
    
    polyA_reads, non_polyA_reads = find_polyA_seq(sam, percentage_threshold, length_threshold, 
                                                   fasta_file, use_fc)
    print('successfully got all polyA reads')
    
    write_output(polyA_reads, out_polyA, out_mode, sam)
    print('successfully got all polyA reads bamfile')
    
    write_output(non_polyA_reads, out_non_polyA, out_mode, sam)
    print('successfully got all nonpolyA reads bamfile')
        
if __name__ == "__main__":    
    run_process()
    print('success')