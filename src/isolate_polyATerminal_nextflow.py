# -*- coding: utf-8 -*-
"""

Created on Sat Jan 22 14:46:41 2022

@author: Youngbin Moon (y.moon@unibas.ch)
"""
import pysam
import pandas as pd
import argparse





"""
Aim : split polA reads mapping to terminal exon by a certain distance: e.g. 200 bp
i.e. annotated or unannotated polyA reads mapping to terminal exons
"""
def write_output(final_reads, out_name, o_mode, sam):
    """
    Parameters
    ----------
    final_reads : list
        a list of reads that that are either annotated or uannotated polyA reads mapping to terminal exons.
    
    out_name : string  
        output file name (annotated or unannotated)
    
    out_mode : string
        output format. 'wb' refers to wrtie in bam file.
    
    bam : bam file
        the original bam file (input)
        This is to use original bam file as a template.
        
    Returns
    -------        
    returns nothing but saves annotated or unannotated polyA reads mapping to terminal exons
    in the bam format.
    """         
    outfile = pysam.AlignmentFile(out_name, o_mode, template=sam)
    for read in final_reads:
        outfile.write(read)
        
def convert_to_dict(list_of_tuples):
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

def get_distance(current_read, current_terminal_exon, fixed):
    """
    Parameters
    ----------
    current_read : pysam object
        a current object that contains infomation about 1 read from BAM file.
    
    current_terminal_exon : dataframe
        terminal exon in which current read is thought to map to.
        It is 1 row (subset) of the terminal_exons.bed file.
    
    fixed : bool
        whether you use fixed softclipped region or original soft clipped region.
        True if you want to use fixed softclipped region.
        False if you do not want to use fixed softclipped region.
        
    Returns
    -------
    distance : int
        An 'absolute' distance between end of a read and end of terminal exon.
    
    """
    alignedRefPositions = current_read.get_reference_positions()
    refStart = alignedRefPositions[0]
    refEnd = alignedRefPositions[-1]   
    rev = current_read.is_reverse
    # read_end = cleavage_site
    # if soft-clipped region is fixed, you just need to retrieve the FC tag (fixed_cleavage_site)
    if fixed == True:
        read_end = current_read.get_tag("FC")
                
    elif fixed == False:
        if rev == True:
            read_end = int(refStart)
            assert(read_end == current_read.get_tag("OC"))
        
        elif rev == False:
            read_end = int(refEnd)
            assert(read_end == current_read.get_tag("OC"))
    
    # need to consider the direction of terminal exon as well.
    if rev == True:
        terminal_exon_end = int(current_terminal_exon['start'])
    
    elif rev == False:
        terminal_exon_end = int(current_terminal_exon['end'])
                
    distance = abs(terminal_exon_end - read_end)
    
    return distance

def split_polyA_reads(sam, bed, threshold, is_fixed):
    """
    Parameters
    ----------
    sam : bam file.
        A bam file that contains polyA reads mapping to terminal exons.
    
    bed : bed file
        A bed file that contains terminal exons only.
    
    threshold : int
        a distance threshold which is used to split polyA reads mapping to terminal exons.
        
    is_fixed : bool
        whether you use fixed softclipped region or original soft clipped region.
        True if you want to use fixed softclipped region.
        False if you do not want to use fixed softclipped region.  
        
    Returns
    -------
    unannotated_reads : list of reads
        a list of polyA reads mapping to terminal exons with distance >= threshold
    
    annotated_reads : list of reads
        a list of polyA reads mapping to terminal exons with distance < threshold
        
    """    
    bed.columns = ['seqid', 'start', 'end', 'id', 'score', 'strand']
    unannotated_reads = set()
    annotated_reads = set()
    for read in sam.fetch():
        list_tuples = read.tags
        tag_dict = convert_to_dict(list_tuples)
        
        min_distance_candidates = []
     
        # there could be reads that do not have GX tag despite they map to terminal exon and have poly A tail.
        if 'GX' in tag_dict.keys():
            # GX_tags is a list of tags
            GX_tags = tag_dict['GX']
            # GX_tags = [ENSMUSG00000022995;ENSMUSG00000025779]            
            if GX_tags[0] != GX_tags[-1] or len(GX_tags) > 1:
                print('you have more than 1 element')

            # GX_tags = [ENSMUSG00000022995, ENSMUSG00000025779]
            GX_tags = GX_tags[0].split(';')
            # length of GX_tags in general is 1 (sometimes you can have more than 2 gene tags for a given read)
            # for each gene
            for GX_tag in GX_tags:
                # For Tabula muris, GX tag is written ENSMUSG00000023143.10
                # but algorithm uses ENSMUSG00000023143
                # if a tag doesnt have ".", it will just output the same thing.
                GX_tag = GX_tag.split(".")[0]
                distances = []
                terminal_exons = bed[bed['id'] == GX_tag]
                filtered_terminal_exons = terminal_exons[terminal_exons['score'] <= 3]

                # If one of the genes (gene A) do not exist in gtf file or if it is filtered out because of high score
                # then there will be no terminal exons for that gene A.(filtered_terminal_exons list is empty)
                # In such case you should skip this gene A. Hence you need to check filtered_terminal_exons.
                # e.g. ENSMUSG000000116048
                if len(filtered_terminal_exons) >= 1:
                    # for each terminal exon
                    # this 'for' phrase is needed because a read can be mapped to more than 1 terminal exons of the same gene.
                    # this can happen even if we only keep reads that are primarily aligned.
                    for index, terminal_exon in filtered_terminal_exons.iterrows():
                        distance = get_distance(read, terminal_exon, is_fixed)
                        distances.append(distance)
                        
                    # find minimum distance in case a read maps to multuple terminal exons of the same gene.
                    # Assume that the terminal exon that has minimum distance to the read is
                    # the true terminal exon that a read maps to.                        
                    min_distance_candidate = min(distances)
                    min_distance_candidates.append(min_distance_candidate)
                    
                else:
                    continue
            # pick min of min (in case there are multiple GX tags per read. Assume only the gene that has minimum of the minimum distance with the read is true mapping.) 
            # if phrase needed in case a read does not have a GX tag or does not have a GX_tag that satisfies the condition.
            # bedtools only check if the regions overlap between a read and terminal exon.
            if len(min_distance_candidates)!= 0:    
                min_index = min_distance_candidates.index(min(min_distance_candidates))
                min_distance = min_distance_candidates[min_index]
            
            # go to next iteration in case the condition not met.
            else:
                continue
            
            # split polyA reads (mapping to terminal exons) according to the distance threshold.
            if min_distance >= threshold:
                # have to add a read once and only once
                unannotated_reads.add(read)
                
            elif min_distance < threshold:
                annotated_reads.add(read)            
        # in case a read does not have GX tag, just skip that read.    
        else:
            continue
    unannotated_reads = list(unannotated_reads)
    annotated_reads = list(annotated_reads)

    return unannotated_reads, annotated_reads

def generate_input_isolate(args):
    """
    Parameters
    ----------
    args : list
        a list of input arguments
    
    Returns
    -------
    bam : a bam file
        polyA reads mapping to terminal exons in a specific chromosome.
    
    bedFile : a bed file
        a bed file containing terminal exons in gtf.
    
    unannotated_out_name : string
        output bam file name for unannoated polyA reads (mapping to terminal exon)

    annotated_out_name : string
        output bam file name for annoated polyA reads (mapping to terminal exon)
    
    thres : int
        distance threshold used for splitting polyA reads mapping to terminal exons into annotated and unannotated.
    
    out_mode : string
        output format. 'wb' refers to wrtie in bam file.
    
    use_fc : bool
      whether use fixed cleavage site or not.
      
    """       
    bamFile = args.bam_input
    number = bamFile.split('_')[-1].split('.')[0]
    bam = pysam.AlignmentFile(bamFile, "rb")    
    
    bed_dir = args.bed_dir
    bedFile = pd.read_csv(bed_dir, delimiter = '\t', header = None)  
    
    unannotated_out_template = args.unannotated_template
    unannotated_out_name = unannotated_out_template +  number + '.bam'

    annotated_out_template = args.annotated_template
    annotated_out_name = annotated_out_template +  number + '.bam' 
    
    thres = args.distance_threshold
    out_mode = "wb"
    
    use_fc = bool(args.use_fc)
    
    return bam, bedFile, unannotated_out_name, annotated_out_name, thres, out_mode, use_fc
        
def get_isolate_args():
    parser = argparse.ArgumentParser(description="isolate polyA reads mapping to terminal exon")
    
    parser.add_argument('--bam_input', dest = 'bam_input',
                        required = True,
                        help = 'polyA reads mapping to terminal exons in a specific chromosome')
    
    parser.add_argument('--bed_dir', dest = 'bed_dir',
                        required = True,
                        help = 'directory towards terminal_exons bed file')    
    
    parser.add_argument('--unannotated_template', dest = 'unannotated_template',
                        required = True,
                        help = 'unannotated reads output name template')
    
    parser.add_argument('--annotated_template', dest = 'annotated_template',
                        required = True,
                        help = 'annotated reads output name template')
    
    parser.add_argument('--distance_threshold', type = int, dest = 'distance_threshold',
                        required = True,
                        help = 'distance threshold for splitting polyA reads mapping to terminal exon into 2 categories')           

    parser.add_argument('--use_fc', dest = 'use_fc',
                        required = True,
                        help = 'use fixed cleavage sites')
    
    args = parser.parse_args()
    return args
    
def run_process():
    
    arguments = get_isolate_args()
    
    bam, bedFile, unannotated_out_name, annotated_out_name, thres, out_mode, use_fc = generate_input_isolate(arguments)
    print("successfully initialized")
    
    un_reads, a_reads = split_polyA_reads(bam, bedFile, thres, use_fc)
    print("successfully got final reads for chr")
    
    write_output(un_reads, unannotated_out_name, out_mode, bam)
    print("successful generation of unannotated polyA sites bam file")
    
    write_output(a_reads, annotated_out_name, out_mode, bam)
    print("successful generation of annotated polyA sites bam file")
        
if __name__ == "__main__":
    run_process()
    print("success")


