
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 28 15:45:14 2022

@author: Youngbin Moon (y.moon@unibas.ch)
"""

import pysam
import argparse

"""
Aim : To add a pseudo tag to the negative control.

Our pipeline does deduplication and fixes soft clipped region.
Hence we use our own fixed cleavage sites saved in 'FC' tag.

But this negative control is deduplicated by UMI tools and does not fix soft clipped region.
Hence, they do not have fixed cleavage site saved in 'FC' tag.
Therefore, here we add 'FC' tag to these negative control reads. (but their content is original cleavage sites not fixed cleavage sites.)
In addition we add 'OC' tag (old cleavage site) to these negative control reads (to be compatible with the rest of the workflow.)

For negative control: FC == OC
"""

def write_output(final_reads, out_name, o_mode, bam):
    """
    Parameters
    ----------
    final_reads : list of reads
        Save the reads that have the 'FC' and 'OC' tag added.
    
    out_name : string
        output bam file name
    
    o_mode : string
        output mode. we use 'wb' which means writing in bam format.
    
    bam : pysam object 
        contains all the reads in the bam file. This is used as a template when writing a new bam file oputput.

    Returns
    -------
    returns nothing but writes output reads (that have 'FC' and 'OC' tag) in the bam format.
    """
    
    outfile = pysam.AlignmentFile(out_name, o_mode, template = bam)
    for read in final_reads:
        outfile.write(read)
        
def trim(bam):
    """
    Parameters
    ----------
    bam : pysam object
        contains all the reads in the negative control bam file.

    Returns
    -------
    changed_reads : a list of reads
        contains all the reads that have 'FC','OC' and 'YB' tag added.
        For negative control, FC == OC and YB == UB.
    """
    changed_reads = []
    for read in bam.fetch():  
        rev = read.is_reverse
        chromosome = read.reference_name
        alignedRefPositions = read.get_reference_positions()
        refStart = alignedRefPositions[0]
        refEnd = alignedRefPositions[-1]
        tuples = read.cigartuples
        
        if rev == True:
            cleavage_site = int(refStart)
            
        elif rev == False:  
            cleavage_site = int(refEnd)
        
        cleavage_site += 1
        read.set_tag("XO", cleavage_site)
        read.set_tag("XF", cleavage_site)
        
        if read.has_tag("UB"):
            UB = read.get_tag("UB")
            read.set_tag("YB", UB)
            
        changed_reads.append(read)
    
    return changed_reads

def get_all_polyA_input():
    
    parser = argparse.ArgumentParser(description="add pseudo tags to the negative control bam file")
    parser.add_argument('--bam_input', dest = 'bam_input',
                        required = True,
                        help = 'full deduplicated bam file')
    
    parser.add_argument('--bam_out', dest = 'bam_out',
                        required = True,
                        help = 'full poly A reads')

    args = parser.parse_args()
    
    bamFile = args.bam_input
    sam = pysam.AlignmentFile(bamFile, "rb")
    
    out_mode = "wb"
        
    bam_out = args.bam_out
    
    return sam, out_mode, bam_out

def run_process():
    
    sam, out_mode, bam_out = get_all_polyA_input()
    print('successfully got inputs')
    
    modified_reads = trim(sam)
    print('successfully changed the reads')
    
    write_output(modified_reads, bam_out, out_mode, sam)
    print('successfully wrote the bam file')
    
if __name__ == "__main__":
    run_process()
    print('success')
    