# -*- coding: utf-8 -*-
"""

Created on Sun Jan  9 07:02:13 2022

@author: Youngbin Moon (y.moon@unibas.ch)
"""

import pysam
import argparse








"""
Aim : collect all reads that do not have softclipped region from deduplicated bam file.
"""

def write_non_softclipped_output(final_reads, o_name, o_mode, bam):
    """
    Parameters
    ----------
    final_reads : list
        a list of reads that does not softclipped region.
    
    out_name : string  
        output file name.
    
    out_mode : string
        output format. 'wb' refers to wrtie in bam file.
    
    bam : bam file
        the original bam file (input)
        This is to use original bam file as a template.
        
    Returns
    -------        
    returns nothing but saves reads that do not have softclipped regions in the bam format.
    """              
    outfile = pysam.AlignmentFile(o_name, o_mode, template=bam)
    for read in final_reads:
        outfile.write(read)
        
def get_no_softclipped(bam):
    """
    Parameters
    ----------    
    bam : bam file
        a bam file containing all deduplicated reads.
        
    Returns
    -------        
    non_softclipped_reads : list
        a list of reads that does not softclipped region.
    """   
    
    non_softclipped_reads = []
    
    for read in bam.fetch():
        tuples = read.cigartuples
        left_end = tuples[0]
        right_end = tuples[-1]
        
        if left_end[0] != 4 and right_end[0] != 4:
            print('no softclipped reads')
            non_softclipped_reads.append(read)
    
    return non_softclipped_reads
            
def get_dedup_input():
    
    parser = argparse.ArgumentParser(description="get reads that do not have softclipped reads")
    parser.add_argument('--bam_input', dest = 'bam_input',
                        required = True,
                        help = 'full deduplicated bam file in which the alignment is fixed')
    
    parser.add_argument('--o_file', dest = 'o_file',
                        required = True,
                        help = 'output bam file')


    args = parser.parse_args()
    
    bamFile = args.bam_input
    sam = pysam.AlignmentFile(bamFile, "rb")
    
    out_mode = "wb"
    
    o_file = args.o_file
    
    return sam, out_mode, o_file

def run_process():
    
    sam, out_mode, o_file = get_dedup_input()
    print('successfully got inputs')
    
    non_softclipped_reads = get_no_softclipped(sam)
    print('successfully got non softclipped reads')
    
    write_non_softclipped_output(non_softclipped_reads, o_file, out_mode, sam)
    print('successfully wrote output')
    
if __name__ == "__main__":    
    run_process()
    print('success')