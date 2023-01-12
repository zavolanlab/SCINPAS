# -*- coding: utf-8 -*-
"""

Created on Mon Jan 10 10:26:42 2022

@author: Youngbin Moon (y.moon@unibas.ch)
"""
import pandas as pd
import pysam
import csv
import argparse
import math

"""
Aim: to get distance between end of the terminal exon and end of the read 
and save the output in csv file.
"""
def get_read_properties(each_read):
    """
    Parameters
    ----------
    each_read : pysam object
        an object that contains infomation about 1 read from BAM file.

    Returns
    -------
    rev : bool
        a direction in which a read maps to.
        if a read maps to a negative direction -> True.
        if a read maps to a positive direction -> False.
    
    chromosome : string
        e.g. chr6, chrY etc......
    
    refStart : int 
        start position of a mapped part of a read.
    
    refEnd : int 
        end position of a mapped part of a read.
    
    CB : string
        Cell bar code. Used to distinguish cell identity.
    
    UMI : string
        Unique molecular identifier. Used to distinguish same transcript.
    """
    rev = each_read.is_reverse
    chromosome = each_read.reference_name
    alignedRefPositions = each_read.get_reference_positions()
    refStart = alignedRefPositions[0]
    refEnd = alignedRefPositions[-1]
    CB = each_read.get_tag('CB')
    UMI = each_read.get_tag('UB')
    
    return rev, chromosome, refStart, refEnd, CB, UMI

def convert_dict(list_of_tuples):
    """
    Parameters
    ----------
    list_of_tuples : a list of tuples
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
    
    is_minus : bool
        whether distance between end of a read to end of terminal exon is - or +.
        if distance is -, is_minus is True.
        if distance is +, is_minus is False.
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
        
    if terminal_exon_end >= read_end:
        is_minus = False
    else:
        is_minus = True
        
    distance = abs(terminal_exon_end - read_end)
    return distance, is_minus 

def get_d_for_polyA_reads(sam, bed, is_fixed):
    """
    Parameters
    ----------
    sam : bam file.
        A bam file that contains polyA reads mapping to terminal exons.
    
    bed : bed file
        A bed file that contains terminal exons only.
    
    is_fixed : bool
        whether you use fixed softclipped region or original soft clipped region.
        True if you want to use fixed softclipped region.
        False if you do not want to use fixed softclipped region.  
        
    Returns
    -------
    min_distances : a list of tuple.
        A list of tuples in which each tuple is (GX, min_distance, min_is_minus, min_log_distance).
        
        GX : GX tag in which distance between a read of interest and several genes is minimal.
        
        min_distance : If a read has 2 GX tags(rarely it happens), we assume a gene that is closest to the read
        is the true gene that a read atually maps to. min_distance is the minimum distance 
        between a read of interest and several genes (terminal exons). 
        
        min_is_minus : whether minimum distance between end of a read to end of terminal exon is - or +
        
        min_log_distance : log10 value of min_distance
    
    subset_reads : a list of tuple.
        A list of tuples which contains info about reads that have distance > 5 (outliers).
        each tuple is (rev, chromosome, refStart, refEnd, CB, UMI).

        rev : bool
            a direction in which a read maps to.
            if a read maps to a negative direction -> True.
            if a read maps to a positive direction -> False.
        
        chromosome : string
            e.g. chr6, chrY etc......
        
        refStart : int 
            start position of a mapped part of a read.
        
        refEnd : int 
            end position of a mapped part of a read.
        
        CB : string
            Cell bar code. Used to distinguish cell identity.
        
        UMI : string
            Unique molecular identifier. Used to distinguish same transcript.        
    """
    min_distances = []
    bed.columns = ['seqid', 'start', 'end', 'id', 'score', 'strand']
    subset_reads = []
    for read in sam.fetch():
        list_tuples = read.tags
        tag_dict = convert_dict(list_tuples)
        
        min_distance_candidates = []
        min_is_minus_candidates = []
        candidates_GX = []
        # there could be reads that do not have GX tag despite they map to terminal exon and have poly A tail.
        if 'GX' in tag_dict.keys():
            # GX_tags is a list of tags
            GX_tags = tag_dict['GX']
            # GX_tags = ['ENSMUSG00000022995;ENSMUSG00000025779']            
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
                are_minus = []
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
                        distance, is_minus = get_distance(read, terminal_exon, is_fixed)
                        distances.append(distance)
                        are_minus.append(is_minus)
                    
                    # find minimum distance in case a read maps to multuple terminal exons of the same gene.
                    # Assume that the terminal exon that has minimum distance to the read is
                    # the true terminal exon that a read maps to.
                    min_distance_candidate = min(distances)
                    min_d_candidate_index = distances.index(min_distance_candidate)
                    min_is_minus_candidate = are_minus[min_d_candidate_index]
                    
                    min_distance_candidates.append(min_distance_candidate)
                    min_is_minus_candidates.append(min_is_minus_candidate)
                    candidates_GX.append(GX_tag)

                else:
                    continue    
            # pick min of min (in case there are multiple GX tags per read. Assume only the gene that has minimum of the minimum distance with the read is true mapping.) 
            # if phrase needed in case a read does not have a GX tag or does not have a GX_tag that satisfies the condition.
            # bedtools only check if the regions overlap between a read and terminal exon.
            if len(min_distance_candidates)!= 0:    
                min_index = min_distance_candidates.index(min(min_distance_candidates))
                min_distance = min_distance_candidates[min_index]
                min_is_minus = min_is_minus_candidates[min_index]
                GX = candidates_GX[min_index]
            
            # go to next iteration in case condition not met.
            else:
                continue
            
            # min_log_distance always >= 0
            if min_distance != 0:
                min_log_distance = math.log10(min_distance)
            elif min_distance == 0:
                min_log_distance = math.log10(1)
            
            min_distances.append((GX, min_distance, min_is_minus, min_log_distance))
            
            # for debugging purpose.
            if min_log_distance > 5:
                if read.has_tag('CB') and read.has_tag('UB'):
                    rev, chromosome, refStart, refEnd, CB, UMI = get_read_properties(read)
                    subset_reads.append((rev, chromosome, refStart, refEnd, CB, UMI))
        
        # in case a read does not have GX tag, just skip that read.
        else:
            continue
    
    return min_distances, subset_reads

def write_out_to_csv(data, file_name):
    """
    Parameters
    ----------
    data : a list of tuples
        A list of tuples in which each tuple is (GX, min_distance, min_is_minus, min_log_distance).
    
    file_name : string 
        output csv file that contains a distribution of 'absolute' distance between end of a read and end of terminal exon.
    
    Returns
    -------
    returns nothing but it writes output to the csv file.
    """    
    w = csv.writer(open(file_name, 'w'))
    for row in data:
        w.writerow(list(row))

def get_inputs():
    
    parser = argparse.ArgumentParser(description="get distance between polyA read and its closest terminal exon")
    parser.add_argument('--bam_input', dest = 'bam_input',
                        required = True,
                        help = 'FULL polyA reads mapping to terminal exons')
    
    parser.add_argument('--bed_input', dest = 'bed_input',
                        required = True,
                        help = 'bed file containing terminal exons')   
    
    parser.add_argument('--distance_out', dest = 'distance_out',
                        required = True,
                        help = 'distance distribution output csv file name')

    parser.add_argument('--use_fc', dest = 'use_fc',
                        required = True,
                        help = 'use fixed cleavage sites')

    parser.add_argument('--sample_or_control', dest = 'sample_or_control',
                        required = True,
                        help = 'is it sample or control')
    
    args = parser.parse_args()
    
    bamFile = args.bam_input
    bed_dir = args.bed_input
    distance_out = args.distance_out
    use_fc = bool(args.use_fc)
    
    sample_or_control = args.sample_or_control
    
    # Umi-dedup negative control is too huge, hence needs to do chromosome wise.
    if sample_or_control == "control":
        number = bamFile.split('_')[7].split('.')[0]
        out_prefix = distance_out.split('.')[0]
        out_suffix = distance_out.split('.')[1]
        distance_out = out_prefix + number + '.' + out_suffix
        outlier_out = 'outlier' + number + '.' + 'csv'
    
    elif sample_or_control == "sample":
        outlier_out = "outlier.csv"
        
    samFile = pysam.AlignmentFile(bamFile, "rb")
    bedFile = pd.read_csv(bed_dir, delimiter = '\t', header = None)
    
    return (samFile, bedFile, distance_out, use_fc, outlier_out)

def run_process():

    (samFile, bedFile, distance_out, use_fc, outlier_out) = get_inputs()
    print('successfully got inputs')
    
    min_distances_distribution, subset_reads = get_d_for_polyA_reads(samFile, bedFile, use_fc)
    print("successfully got distance distribution from read end to end of terminal exon")
    
    write_out_to_csv(min_distances_distribution, distance_out)
    print('successfully wrote distribution to csv')
    
    write_out_to_csv(subset_reads, outlier_out)
    print('successfully wrote outlier')
    
if __name__ == "__main__":
    
    run_process()
    print("success")
    
    