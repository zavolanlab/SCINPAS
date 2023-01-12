# -*- coding: utf-8 -*-
"""

Created on Wed Jan 26 01:00:23 2022

@author: Youngbin Moon (y.moon@unibas.ch)
"""
import pysam
import csv
import argparse






"""
Aim : collect polyA reads that do not map to exon part of the genes (polyA reads with no GX tag)
and save it as a bam file
"""
def convert_into_dict(list_of_tuples):
    """
    Parameters
    ----------
    list_of_tuples : list
        an object that contains infomation about 1 read from BAM file.
        1st element of a tuple: tag name. For example, 'GX' tag
        2nd element of a tuple: value of the tag. if it is 'GX' tag,
        then value is gene name(s) that a read maps to. For example, ENSMUSG00000022995.
        
    Returns
    -------        
    tag_dict : dictionary
        key = tag name. e.g. 'GX' tag.
        value = value that the tag has. e.g. if key is 'GX', value is gene name.
    """      
    tag_dict = dict()
    for key, tag in list_of_tuples:
        tag_dict.setdefault(key, []).append(tag)
    return tag_dict           

def get_noGX_polyA(bam):
    """
    Parameters
    ----------
    bam : bam file    
       a bam file containing all polyA reads.
       
    Returns
    -------        
    noGX_polyA_reads : a list of reads
        a list of reads that do not have GX tag. (polyA reads that do not map to exon part of the genes)
    """       
    noGX_polyA_reads = []    
    for read in bam.fetch():
        list_tuples = read.tags
        tag_dict = convert_into_dict(list_tuples)
        # those reads that have GX tag we covered them already in polyA reads mapping to terminal exon.
        if 'GX' not in tag_dict.keys():
            noGX_polyA_reads.append(read)
            
    return noGX_polyA_reads

def save_output(final_reads, out_name, out_mode, bam):
    """
    Parameters
    ----------
    final_reads : list
        a list of reads that do not have GX tag. (reads that do not map to the genes)
    
    out_name : string  
        output file name
    
    out_mode : string
        output format. 'wb' refers to wrtie in bam file.
    
    bam : bam file
        the original bam file (input)
        This is to use original bam file as a template.
        
    Returns
    -------        
    returns nothing but saves a list of polyA reads that do not have GX tag in the bam format. 
    (polyA reads that do not map to exon part of the genes)
    """           
    outfile = pysam.AlignmentFile(out_name, out_mode, template = bam)
    for read in final_reads:
        outfile.write(read)
        
def get_input_noGX():
    
    parser = argparse.ArgumentParser(description="get noGX reads (pro intronic intergenic reads)")
    parser.add_argument('--bam_input', dest = 'bam_input',
                        required = True,
                        help = 'full all polyA reads bam file')
      
    parser.add_argument('--bam_out', dest = 'bam_out',
                        required = True,
                        help = 'output bam file')
        
    args = parser.parse_args()
    
    bamFile = args.bam_input
    sam = pysam.AlignmentFile(bamFile, "rb")
    out_mode = "wb"
        
    bam_out = args.bam_out

    return sam, out_mode, bam_out
        
def run_process():
    
    sam, out_mode, bam_out = get_input_noGX()
    print('successfully got inputs')
                
    noGX_polyA_reads = get_noGX_polyA(sam)
    print('successfully got resulting reads')
    
    save_output(noGX_polyA_reads, bam_out, out_mode, sam)
    print('successfully saved into bam file')
        
if __name__ == "__main__":
    run_process()
    print('success')
       
    