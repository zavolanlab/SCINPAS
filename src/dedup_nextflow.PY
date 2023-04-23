# -*- coding: utf-8 -*-
"""

Created on Sun Dec 19 22:55:32 2021

@author: Youngbin Moon (y.moon@unibas.ch)
"""
import pysam
import numpy as np
import math
import argparse

"""
Aim 1 : Within a cluster of reads with same CB, UR tag and chromosome if span is too large, split a cluster into several sub-clusters
and assign different UR tags (YB tag). (Assuming they had same UR by chance)

Aim 2 : For a given CB, merge 2 clusters of reads that are nearby and have hamming distance = 1 in their UR tags.
the new UR tag is the one with the majority.

After merging step, UR tags are saved in the UB tag.
Note: we merge 2 clusters at a time.

Aim 3 : deduplicate reads (remove PCR artefact) by choosing the most distal read amongst reads with same CB and UB tag same chromosome. (choosing a read that is closest to the terminal exon)
This is to capture true polyA read.

Aim 4 : save all reads from original bam file but a new tag YB tag attached.
Later we will use this 'YB' tag to compute span in span_alter_nextflow.py.
"""
    
def flatten_reads(splitted_dict):
    """
    Parameters
    ----------
    splitted_dict: dictionary of a list
        key = UR, start, end position and chromosome of a new cluster
        value =  a list where it contains all reads in a new sub-cluster (with same CB and UR tag).
        e.g. splitted_dict = { (UR1.1, (1, 5, chr1)):[r1, r2, r3, r4, r5], (UR1.2, (10^5, 10^5 + 500, chr1)):[r6, r7, r8], (UR6, (0, 6, chr20)):[r9, r10, r11, r12, r13, r14, r15, r16, r17, r18] }
        
        NOTE: This splitted_dict is a dictionary for a given cell (CB).   

    Returns
    -------
    partial_split_reads : list
        a list of all reads in a given cell (CB) but a new tag 'YB' attached.
        This 'YB' tag is a result of splitting a cluster with same CB and UR tag into sub-clusters.
        We use this 'YB' tag to compute span in span_alter_nextflow.py.
        
    """
    partial_split_reads = []
    # list_of_clusters = [[r1, r2, r3, r4, r5], [r6, r7, r8], [r9, r10, r11, r12, r13, r14, r15, r16, r17, r18]]
    list_of_clusters = list(splitted_dict.values())
    for cluster in list_of_clusters:
        for read in cluster:
            partial_split_reads.append(read)

    return partial_split_reads

def get_splitted_reads(split_dicts, Cbs):
    """
    Parameters
    ----------
    split_dicts: dictionary of a list10^5, 10^5 + 500
        1st key = CB
        2nd key = UR, start, end position and chromosome of a new cluster
        value =  a list where it contains all reads in a new sub-cluster (with same CB and UR tag, start, end and chromosome).
        e.g. {'cb1':{ (UR1.1, (1, 5, chr1)):[r1, r2, r3, r4, r5], (UR1.2, (10^5, 10^5 + 500, chr1)):[r6, r7, r8], (UR6, (0, 6, chr20)):[r9, r10, r11, r12, r13, r14, r15, r16, r17, r18] }}
        
    Cbs : list
        a list that contains all unique cell-barcodes. 

    Returns
    -------
    total_splitted_reads : list
        a list of all reads in a given input bam file but a new tag 'YB' attached.
        This 'YB' tag is a result of splitting a cluster with same CB and UR tag into sub-clusters.
        We use this 'YB' tag to compute span in span_alter_nextflow.py.
    """
    total_splitted_reads = []
    for cb in Cbs:
        partial_splitted_reads = flatten_reads(split_dicts[cb])
        total_splitted_reads += partial_splitted_reads
    
    return total_splitted_reads

def write_output(final_reads, out_name, out_mode, sam):
    """
    Parameters
    ----------
    final_reads : list
       
        1) either a list of reads which contains all deduplicated reads in all cells (CBs).
        We choose a single read per CB and UB that is closest to the end of the terminal exon.
        
        2) or a list of all reads in a given input bam file but a new tag 'YB' attached.
    
    out_name : string
        1) either output bamfile name of deduplication.
        
        2) or output bamfile name which contains all reads in the original bam file but a new tag 'YB' attached.
    
    out_mode : string
        output format. 'wb' refers to write in bam format.
        
    sam : bam file
        the original raw input bam file. We use this as a template when saving outputs in the bam format.

    Returns
    -------
    returns nothing but 
    
    1) save all deduplicated reads in a new bam file. (deduplicated bam file)
    2) save all reads in a given input bam file but a new tag 'YB' attached. (splitted bam file)
    """     
    outfile = pysam.AlignmentFile(out_name, out_mode, template=sam)
    for read in final_reads:
        outfile.write(read)

def extract_spans(UR_and_spans):
    """
    Parameters
    ----------
    UR_and_spans : list
        a list of tuple that contains URs, start positions and end positions of all sub-clusters.

    Returns
    -------
    starts : list
        a list of start positions of all sub-clusters.
        
    ends : list
        a list of end positions of all sub-clusters.
        
    URs : list
        a list of URs of all sub-clusters.
    
    chroms : list
        a list of choromosomes of all sub-clusters.

    """
    starts = [UR_and_span[1][0] for UR_and_span in UR_and_spans]
    ends = [UR_and_span[1][1] for UR_and_span in UR_and_spans]
    URs = [UR_and_span[0] for UR_and_span in UR_and_spans]
    chroms = [UR_and_span[1][2] for UR_and_span in UR_and_spans]
    
    return starts, ends, URs, chroms


def deduplication(new_dict):
    """
    Parameters
    ----------
    new_dict : dictionary of a list
        same structure as splitted_dict but some sub-clusters are fused with a new key.                
        key = UR, start, end position and chromosome of a new cluster
        value =  a list where it contains all reads in a new sub-cluster (with same CB and UR tag, start, end and chromosome).    
        Note: Since span is also part of the key, even if you have same CB and UR, if you have different spans, you consider as a different cluster.
        Hence if you have a cluster oF CB and UR that has too much distance, they are subdivided into different clusters and treated as different clusters (Same effect as using YB tag).
        
        e.g. splitted_dict = { (UR1.1, (1, 5, chr1)):[r1, r2, r3, r4, r5], (UR2, (4, 25, chr1)):[r6, r7, r8], (UR6, (0, 6, chr20)):[r9, r10, r11, r12, r13, r14, r15, r16, r17, r18] } 
        new_dict = {(UR1.1, (1, 25, chr1)):[r1, r2, r3, r4, r5, r6, r7, r8], (UR6, (0, 6, chr20)):[r9, r10, r11, r12, r13, r14, r15, r16, r17, r18]}
        Here, UR1.1 and UR2 overlap and assuming their hamming distance = 1 they are fused to a single cluster.
        Since UR1.1 is the majority, the new cluster UR is UR1.1.
        Note that span is also changed because of fusion.
        
    NOTE: This splitted_dict and new_dict are dictionary for a given cell (CB).    
    
    Returns
    -------
    partial_final_reads : list
        a list of reads which contains all deduplicated reads in a given cell (CB).
        We choose a single read per CB and UB that is closest to the end of the terminal exon.
    
    Also, the UB tag is changed to the majority UR as a result of fusing sub-clusters.
    """
    partial_final_reads = []
    # list_of_clusters_r = [[r1, r2, r3, r4, r5, r6, r7, r8], [r9, r10, r11, r12, r13, r14, r15, r16, r17, r18]]
    list_of_clusters_r = list(new_dict.values())
    UR_and_spans = list(new_dict.keys())
    print("UR_and_spans: " + str(UR_and_spans))
    starts, ends, URs, chroms = extract_spans(UR_and_spans)
    
    for i, cluster in enumerate(list_of_clusters_r):
        start_positions = []
        end_positions = []
        are_negatives = []
        reads = []
        # this UR is the majority UR that you decided in the 'merge' step.
        UR = URs[i]
        
        for read in cluster:
            alignedRefPositions = read.get_reference_positions()
            refStart = alignedRefPositions[0]
            refEnd = alignedRefPositions[-1]
            rev = read.is_reverse
            
            start_positions.append(refStart)
            end_positions.append(refEnd)
            are_negatives.append(rev)
            
            # change the UMI of the read
            read.set_tag("UB", UR)
            reads.append(read)
            
        num_neg = sum([1 if elem == True else 0 for elem in are_negatives])
        num_pos = sum([1 if elem == False else 0 for elem in are_negatives])
        
        # if all reads are mapped to negative strand, choose the read that has minimum start position (closest to the end of a terminal exon)
        if num_neg == len(are_negatives) and num_pos == 0:
            read_idx = start_positions.index(min(start_positions))
            partial_final_reads.append(reads[read_idx])
        
        # if all reads are mapped to positive strand, choose the read that has maximuim end position (closest to the end of a terminal exon)
        elif num_pos == len(are_negatives) and num_neg == 0:
            read_idx = end_positions.index(max(end_positions))
            partial_final_reads.append(reads[read_idx])
            
        # if the majority of reads are mapped to negative strand (but not all), choose the read that has minimum start position amongst those reads mapping to - strand (closest to the end of a terminal exon).
        elif num_neg >= num_pos:
            neg_indices = [i for i in range(len(are_negatives)) if are_negatives[i] == True]
            reads_subset = [reads[index] for index in neg_indices]
            starts_subset = [start_positions[index] for index in neg_indices]
            
            read_idx = starts_subset.index(min(starts_subset))
            partial_final_reads.append(reads_subset[read_idx])
        
        # if the majority of reads are mapped to positive strand (but not all), choose the read that has maxium end position amongst those reads mapping to + strand (closest to the end of a terminal exon).
        elif num_pos > num_neg:
            pos_indices = [i for i in range(len(are_negatives)) if are_negatives[i] == False]
            reads_subset = [reads[index] for index in pos_indices]
            ends_subset = [end_positions[index] for index in pos_indices]
            
            read_idx = ends_subset.index(max(ends_subset))
            partial_final_reads.append(reads_subset[read_idx])   
            
    return partial_final_reads

def hamming_distance(UR1, UR2):
    """
    Parameters
    ----------
    UR1 : string
        a UR tag of a current sub-cluster.
        
    UR2 : string
        a UR tag of a next sub-cluster.

    Returns
    -------
    ham_dist : int
        hamming distance between UR1 and UR2 tags.
    """
    assert(len(UR1)==len(UR2))
    count = [1 if UR1[i] != UR2[i] else 0 for i in range(len(UR1))]
    ham_dist = sum(count)
    
    return ham_dist

def check_overlap (start1, end1, start2, end2, chromosome1, chromosome2):
    """
    Parameters
    ----------
    start1 : int
        start position of a sub-cluster 1.
        
    end1 : int
        end position of a sub-cluster 1.
        
    start2 : int
        start position of a sub-cluster 2.
    
    end2 : int
        end position of a sub-cluster 2.
    
    chromosome1: str
        chromosome of a sub-cluster 1.
        
    chromosome2: str
        chromosome of a sub-cluster 2.
        
    Returns
    -------
    is_overlap : bool
        tells you whether a sub-cluster 1 overlaps with a sub-cluster 2 or not.
        True if a sub-cluster 1 overlaps with a sub-cluster 2.
        False if a sub-cluster 1 does not overlap with a sub-cluster 2.
    """
    is_overlap = (start2 <= end1) and (end2 >= start1) and (chromosome1 == chromosome2)
    return is_overlap

def extract(UR_and_span):
    """
    Parameters
    ----------
    UR_and_span : a tuple
        a tuple where it is made of UR, start, end position and chromosome of a sub-cluster.
        e.g. (UR6, (0, 6, chr1))

    Returns
    -------
    UR : string
        a UR tag of the current sub-cluster.
        
    start : int
        start position of a sub-cluster.
        
    end : int
        end position of a sub-cluster.
    
    chrom : str
        chromosome in which a sub-cluster belongs to.
    """      
    UR = UR_and_span[0]
    start = UR_and_span[1][0]
    end = UR_and_span[1][1]
    chrom = UR_and_span[1][2]
        
    return UR, start, end, chrom

def merge(splitted_dict):
    """
    Parameters
    ----------
    splitted_dict: dictionary of a list
        key = UR, start, end position and chromosome of a new cluster
        value =  a list where it contains all reads in a new sub-cluster (with same CB and UR tag).
        e.g. splitted_dict = { (UR1.1, (1, 5, chr1)):[r1, r2, r3, r4, r5], (UR2, (4, 25, chr1)):[r6, r7, r8], (UR6, (0, 6, chr20)):[r9, r10, r11, r12, r13, r14, r15, r16, r17, r18] }

    Returns
    -------
    new_dict : dictionary of a list
        same structure as splitted_dict but some sub-clusters are fused with a new key.
        We fuse 2 sub-clusters at a time but fusing only happens if you satisfy 2 conditions.
        1) hamming distance between 2 URs = 1
        2) 2 sub-clusters overlap.
    
        e.g. splitted_dict = { (UR1.1, (1, 5, chr1)):[r1, r2, r3, r4, r5], (UR2, (4, 25, chr1)):[r6, r7, r8], (UR6, (0, 6, chr20)):[r9, r10, r11, r12, r13, r14, r15, r16, r17, r18] } 
        new_dict = {(UR1.1, (1, 25, chr1)):[r1, r2, r3, r4, r5, r6, r7, r8], (UR6, (0, 6, chr20)):[r9, r10, r11, r12, r13, r14, r15, r16, r17, r18]}
        Here, UR1.1 and UR2 overlap and assuming their hamming distance = 1 they are fused to a single cluster.
        Since UR1.1 is the majority, the new cluster UR is UR1.1.
        Note that span is also changed because of fusion.
        
    
    NOTE: This splitted_dict and new_dict are dictionary for a given cell (CB).
    """
    new_dict={}
    # splitted_dicts = {'cb1':{ (UR1.1, (1, 5, chr1)):[r1, r2, r3, r4], (UR1.2, (4, 25, chr1)):[r5, r6, r7, r8], (UR6, (0, 6, chr20)):[r9, r10, r11, r12, r13, r14, r15, r16, r17, r18] }}
    # splitted_dict = { (UR1.1, (1, 5, chr1)):[r1, r2, r3, r4], (UR1.2, (4, 25, chr1)):[r5, r6, r7, r8], (UR6, (0, 6, chr20)):[r9, r10, r11, r12, r13, r14, r15, r16, r17, r18] }
    UR_and_spans = list(splitted_dict.keys())
    
    # For a given CB if you have only 1 UR cluster, then dont merge. just return the same thing.
    if len(UR_and_spans) == 1:
        new_dict[UR_and_spans[0]] = splitted_dict[UR_and_spans[0]]
        
    else:
        # use stable sort to get rid of randomness
        # need to sort first by chromosome and then start and then end
        # 1 element: (UR1.1, (1, 5, chr1)
        sorted_UR_and_spans = sorted(UR_and_spans, key = lambda x: (x[1][2], x[1][0], x[1][1]))
        
        # reserve == 'not changed' means previous 2 clusters are not fused in the previous iteration.
        reserve = 'not changed'
        
        for i in range(len(sorted_UR_and_spans)-1):
            # This means previous 2 clusters are fused into 1 cluster. Hence you need to use the new key.
            # We do this iteratively one by one, so all keys are considered. It is just a matter of using old or new key.
            if reserve != 'not changed':
                UR_and_span1 = reserve
                
            else:
                UR_and_span1 = sorted_UR_and_spans[i]
                
            UR_and_span2 = sorted_UR_and_spans[i+1]
            print("UR_and_span1 is: " + str(UR_and_span1))
            UR1, start1, end1, chrom1 = extract(UR_and_span1)
            UR2, start2, end2, chrom2 = extract(UR_and_span2)
                
            is_overlap = check_overlap(start1, end1, start2, end2, chrom1, chrom2)
            ham_d = hamming_distance(UR1, UR2)
            
            # if reserve is changed to a new key it means previous 2 clusters are fused into 1 cluster(changed). 
            # So you need to get reads from new dictionary with a new key.
            # We do this iteratively one by one, so all reads are considered.
            if reserve != 'not changed':
                reads_C1 = new_dict[reserve]
                
            # if reserve is not changed, it means it is not a fused cluster. so you get reads from old original dictionary
            else:
                reads_C1 = splitted_dict[UR_and_span1]
                
            # reads_C2 are new reads. they are not fused yet. so you get reads from the old original dictionary
            reads_C2 = splitted_dict[UR_and_span2]        
            
            if is_overlap and ham_d == 1 :
    
                if len(reads_C1) >= len(reads_C2):
                    assert(chrom1 == chrom2)
                    new_span = (min(start1, start2), max(end1, end2), chrom1)
                    # need to change to a key with UR1 because UR1 is majority. 
                    # also you have a new 'merged' span.
                    new_key = (UR1, new_span)
    
                    # delete reads in the old key if reads exist with old key(new_dict[UR_and_span] not empty)
                    # delete key and values
                    if UR_and_span1 in new_dict.keys():
                        del new_dict[UR_and_span1]
                        
                    if UR_and_span2 in new_dict.keys():
                        del new_dict[UR_and_span2]
                        
                    # 2 cases: 
                    # 1) didnt have new key in the new_dict (e.g.because 2 clusters are fused and its span is changed) -> delete old key values (because fused) -> add key values
                    # 2) have same key in the new_dict (because 2 clusters are fused but span is not changed) -> delete key values anyway. -> add key values

                    # your span might not change. you might have same key.
                    # but it doesnt matter because in such case, you deleted the old key values beforehand.
                    
                    if new_key not in new_dict.keys():
                        new_dict[new_key] = []                    
                     
                    new_dict[new_key] += reads_C1                
                    new_dict[new_key] += reads_C2
    
                    reserve = new_key
                    
                else:
                    
                    # need to change key with UR2 becuase UR2 is majority.
                    # also you have a new 'merged' span.
                    assert(chrom1 == chrom2)
                    new_span = (min(start1, start2), max(end1, end2), chrom2)
                    new_key = (UR2, new_span)
                    
                    # delete reads in the old key if reads exist with old key (i.e. new_dict[UR_and_span] not empty)
                    # delete key and values
                    if UR_and_span1 in new_dict.keys():
                        del new_dict[UR_and_span1]
                        
                    if UR_and_span2 in new_dict.keys():
                        del new_dict[UR_and_span2]
                        
                    if new_key not in new_dict.keys():
                        new_dict[new_key] = []
                        
                    new_dict[new_key] += reads_C1                
                    new_dict[new_key] += reads_C2
    
                    reserve = new_key
                    
            # else put them separately (do not merge 2 clusters)
            else:     
                
                if UR_and_span1 not in new_dict.keys():
                    new_dict[UR_and_span1] = []
                    
                if UR_and_span2 not in new_dict.keys():
                    new_dict[UR_and_span2] = []
                    
                # Add reads to this new_dict with the key if and only if new_dict is empty with this key. (Else, do nothing because reads are already there)
                # new_dict[UR_and_span1] can already have a key value if in the previous iteration, clusters are fused.
                # In such case, you should not add reads.
                # UR_and_span1 should already exist in the new dict unless it is really first element.
                # UR_and_span2 is always new key and hence new_dict should not have this key beforehand.
                if not new_dict[UR_and_span1]:
                    new_dict[UR_and_span1] += reads_C1
                    
                if not new_dict[UR_and_span2]:
                    new_dict[UR_and_span2] += reads_C2
                    
                reserve = 'not changed'
    # splitted_dict = { (UR1.1, (1, 5, chr1)):[r1, r2, r3, r4, r5], (UR1.2, (4, 25, chr1)):[r6, r7, r8], (UR6, (0, 6, chr20)):[r9, r10, r11, r12, r13, r14, r15, r16, r17, r18] }            
    # new_dict = {(UR1.1, (1, 25, chr1)):[r1, r2, r3, r4, r5, r6, r7, r8], (UR6, (0, 6, chr20)):[r9, r10, r11, r12, r13, r14, r15, r16, r17, r18]}
    return new_dict

def merge_two_dicts(x, y):
    """

    Parameters
    ----------
    x : dictionary of a list
        key = UR, start, end position and chromosome of a new cluster
        value =  a list where it contains all reads in a new sub-cluster (with same CB, UR tag, start, end and chromosome).
        e.g. { (UR1, (1, 5, chr1)):[r1, r2, r3, r4], (UR1, (10^5, 10^5 + 500, chr1)):[r5, r6, r7, r8] }
   
    y : dictionary of a list
        key = UR, start, end position and chromosome of a new cluster
        value =  a list where it contains all reads in a new sub-cluster (with same CB, UR tag, start, end and chromosome).
        e.g. {(UR6, (0, 6, chr20)):[r9, r10, r11, r12, r13, r14, r15, r16, r17, r18] }

    Returns
    -------
    z : dictionary of a list
        key = UR, start, end position and chromosome of a new cluster
        value =  a list where it contains all reads in a new sub-cluster (with same CB, UR tag, start, end and chromosome).
        
        z refers to a new merged dictionary between x and y.
        e.g. { (UR1, (1, 5, chr1)):[r1, r2, r3, r4], (UR1, (10^5, 10^5 + 500, chr1)):[r5, r6, r7, r8], (UR6, (0, 6, chr20)):[r9, r10, r11, r12, r13, r14, r15, r16, r17, r18] }        
    """
    z = x.copy()
    z.update(y)
    return z

def get_span(cluster_p):
    """
    Parameters
    ----------
    cluster_p : a list of tuple
        a list of tuple where each element is start and end position of reads within a new sub-cluster (with same CB and UR tag and chromosome)

    Returns
    -------
    start : int
        start of the span of a new sub-cluster.
        
    end : int
        end of the span of a new sub-cluster.

    """
    starts = [tupl[0] for tupl in cluster_p]
    ends = [tupl[1] for tupl in cluster_p]
    start = min(starts)
    end = max(ends)
    return (start, end)

def split(sorted_tuple_positions, start_positions, end_positions, threshold, sorted_reads, chromosome):
    """
    Parameters
    ----------
    sorted_tuple_positions : a list
        a list of tuple where each element is: 
        (start, end) of the read within a cluster of the same CB and UR tag.
    
    start_positions :a list
        a list of start positions of a read within a cluster of the same CB and UR tag.
            
    end_positions : a list
        a list of end positions of a read within a cluster of the same CB and UR tag.
    
    threshold : int
        Span threshold for a cluster of reads (with same CB and UR tag) to be considered as 1 cluster.
    
    sorted_reads :a list
        a list of reads with same CB and UR tags that are sorted by start and then end positions.
    
    chromosome: str
        chromosome number
        
    Returns
    -------
    dict_of_clusters_r : dictionary of list.
    
        first key :  tuple of UR, start, end and chromosome of the new cluster. Each key represents a new sub-cluster. 
        i.e. (UR, (new_cluster_start, new_cluster_end, new_cluster_chromosome))
        
        value : a list where it contains all reads in a new sub-cluster (with same CB, UR tag, start, end and chromosome).
        
        dict_of_clusters_r[('UR', (new_cluster_start, new_cluster_end, chromosome))] = [r19, r20, r21, r22, r23]
    """
    first_tuple = sorted_tuple_positions[0]
    min_start = first_tuple[0]
    assert(min_start == min(start_positions))
    max_end = max(end_positions)
    dict_of_clusters_r = {}
    cluster_p = []
    cluster_r = []
    span = max_end - min_start
    
    logspan = math.log10(span + 1)    
    # if logspan >= threshold, you need to split a cluster into subclusters at some point.
    if logspan >= threshold:
        k = 1
        for i, sorted_tuple in enumerate(sorted_tuple_positions):
            end = sorted_tuple[1]
            # if span is larger than threshold, no need to consider overlapping region
            # if span is larger than threshold, no need to consider if 2 reads have similar UR
            mini_span = end - min_start
                        
            mini_logspan = math.log10(mini_span + 1)
            if mini_logspan > threshold:
                # Need to use UR of the previous read not current read because you are splitting until previous read and current read onward. 
                # i-1 works because first element will not have mini_logspan > threshold.
                UR = sorted_reads[i-1].get_tag('UR')
                (new_cluster_start, new_cluster_end) = get_span(cluster_p)

                # dict_of_clusters_r[('UR', (new_cluster_start, new_cluster_end, chromosome))] = [r19, r20, r21, r22, r23]
                # e.g. cluster_r =[r19, r20, r21, r22, r23]
                if (UR, (new_cluster_start, new_cluster_end, chromosome)) not in dict_of_clusters_r.keys():
                    dict_of_clusters_r[(UR, (new_cluster_start, new_cluster_end, chromosome))] = []
                dict_of_clusters_r[(UR, (new_cluster_start, new_cluster_end, chromosome))] += cluster_r
                
                YB = UR + str(k)
                k += 1
                
                # assign a new tag within a new splitted cluster. e.g. cluster_r = [r19, r20, r21, r22, r23]
                [read.set_tag("YB", YB) for read in cluster_r]
                
                #reset a new cluster
                cluster_p =[]
                cluster_r =[]
                
                # reset the first read start position
                min_start = sorted_tuple[0]
                
                # you have to add this read to the new cluster as a seed
                # seed should be always included in the new cluster
                cluster_p.append(sorted_tuple)
                cluster_r.append(sorted_reads[i])     
                
            else:
                # Since it is sorted, the first element will always enter here and first element will be added to the first cluster as a seed.
                # (because if it is first element, min_logpsan can't be > threshold)
                # no need to check overlap or ham distance because they all have same UR.
                cluster_p.append(sorted_tuple)
                cluster_r.append(sorted_reads[i])
                 
        # you do this because always all reads in the last cluster are not saved in the dict_of_clusters_r in the above statments.
        # Hence you add all reads in the last cluster here.
        UR = sorted_reads[-1].get_tag('UR')
        (new_cluster_start, new_cluster_end) = get_span(cluster_p)

        if (UR, (new_cluster_start, new_cluster_end, chromosome)) not in dict_of_clusters_r.keys():
            dict_of_clusters_r[(UR, (new_cluster_start, new_cluster_end, chromosome))] = []
        dict_of_clusters_r[(UR, (new_cluster_start, new_cluster_end, chromosome))] += cluster_r        
        
        YB = UR + str(k)
        [read.set_tag("YB", YB) for read in cluster_r]
                           
    else:
        k = 1
        # nothing to split
        # they have all same UR
        UR = sorted_reads[0].get_tag('UR')
        (new_cluster_start, new_cluster_end) = get_span(sorted_tuple_positions)
        
        if (UR, (new_cluster_start, new_cluster_end, chromosome)) not in dict_of_clusters_r.keys():
            dict_of_clusters_r[(UR, (new_cluster_start, new_cluster_end, chromosome))] = []
        
        # nothing to split. put every sorted reads in here    
        dict_of_clusters_r[(UR, (new_cluster_start, new_cluster_end, chromosome))] += sorted_reads
        
        YB = UR + str(k)
        [read.set_tag("YB", YB) for read in sorted_reads]
        
    return dict_of_clusters_r

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
        
def purify_UR(reads_per_ub, unique_pairs, threshold):
    """
    Parameters
    ----------
    reads_per_ub : a dictionary of dictionary
    1st key : a tuple of CB, UR and chromosome.
    2nd key : 'reads'
    value : a list of reads with the same CB and UR tag and a chromosome
        
    unique_pairs : a list of tuple
        a list of tuple where each element is unique tuple of CB, UR and chromosome.
    
    threshold : int
        Span threshold for a cluster of reads (with same CB and UR tag and chromosome) to be considered as 1 cluster.
          
    Returns
    -------
    final_reads : list
        a list of reads which contains all deduplicated reads in all cells (CBs)
        We choose a single read per CB and UB that is closest to the end of the terminal exon.
    
    splitted_dicts: dictionary of a list
        1st key = CB
        2nd key = UR, start, end position and chromosome of a new cluster
        value =  a list where it contains all reads in a new sub-cluster (with same CB and UR tag and chromosome).
        e.g. {'cb1':{ (UR1, (1, 5, chr1)):[r1, r2, r3, r4], (UR1, (10^5, 10^5 + 500, chr1)):[r5, r6, r7, r8], (UR6, (0, 6, chr20)):[r9, r10, r11, r12, r13, r14, r15, r16, r17, r18] }}
        
    CBs : list
        a list that contains all unique cell-barcodes. 
    """       
    final_reads = []
    list_c = ['I', 'D', 'N', 'S', 'H', 'P', 'X']
    splitted_dicts ={}
    # to avoid randomness, fixate it by stable sort
    unique_pairs = sorted(unique_pairs, key = lambda x: (x[0], x[1], x[2]))
    
    for unique_pair in unique_pairs:
        read_list = reads_per_ub[unique_pair]['reads']
        start_positions = []
        end_positions = []
        positions_reads = []
        
        for i in range(len(read_list)):     
            read = read_list[i]
            cigar = read.cigarstring
            # N means Skipped region from the reference.
            # checking if there is a skipped region (span to intronic region)
            if 'N' not in cigar:
                alignedRefPositions = read.get_reference_positions()
                refStart = alignedRefPositions[0]
                refEnd = alignedRefPositions[-1]
                start_positions.append(refStart)
                end_positions.append(refEnd)
                positions_reads.append((refStart, refEnd, read))
                
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
                    positions_reads.append((start, end, read))
                    
                elif rev == False:
                    start = refEnd - distal + 1
                    end = refEnd
                    start_positions.append(start)
                    end_positions.append(end)
                    positions_reads.append((start, end, read))
                            
        # sort by start position first and then by end positions   
        sorted_positions_reads = sorted(positions_reads, key = lambda x: (x[0], x[1]))
        sorted_tuple_positions = [(sorted_tuple[0], sorted_tuple[1]) for sorted_tuple in sorted_positions_reads]
        # get sorted reads
        sorted_reads = [sorted_tuple[2] for sorted_tuple in sorted_positions_reads]
        chrom = unique_pair[2]
        
        # within same CB, UR tag and a chromosome. you have to add dictionary
        splitted_per_ur = split(sorted_tuple_positions, start_positions, end_positions, threshold, sorted_reads, chrom)
        
        # put dictionary to the splitted_dicts under the same cb outermost key
        # (UR1, (1,5, chr1)) and (UR1, (10^5, 10^5 + 500, chr1)) are different subclusters despite having same CB and UR. because they have too large span.
        # Their tags are differently saved in 'YB' tag.
        cb = unique_pair[0]
        
        if cb not in splitted_dicts.keys():
            splitted_dicts[cb] = splitted_per_ur
            
        else:
            # merge two dictionaries under same cb. so that you use this in the merge step
            # This is within the same chromosome:
            # splitted_dicts[cb] = { (UR1, (1, 5, chr1)):[r1, r2, r3, r4], (UR1, (10^5, 10^5 + 500, chr1)):[r5, r6, r7, r8] }
            # splitted_per_ur = {(UR6, (0, 6, chr20)):[r9, r10, r11, r12, r13, r14, r15, r16, r17, r18] }
            splitted_dicts[cb] = merge_two_dicts(splitted_dicts[cb], splitted_per_ur)
    
    # after merginx two dictionaries:
    # splitted_dicts = {'cb1':{ (UR1, (1, 5, chr1)):[r1, r2, r3, r4], (UR1, (10^5, 10^5 + 500, chr1)):[r5, r6, r7, r8], (UR6, (0, 6, chr20)):[r9, r10, r11, r12, r13, r14, r15, r16, r17, r18] }}
    # This new_dict should be within same CB tag
    cbs = splitted_dicts.keys()
    for cb in cbs:
        new_dict = merge(splitted_dicts[cb])
        partial_final_reads = deduplication(new_dict)
        final_reads += partial_final_reads
        
    return final_reads, splitted_dicts, cbs

def get_all_reads_per_umi (sam):
    """
    Parameters
    ----------
    sam : bam file
        an input raw bam file.
    
    Returns
    -------
    reads_per_ub : a dictionary of dictionary
    1st key : a tuple of CB, UMI and a chromosome
    2nd key : 'reads'
    value : a list of reads with the same CB and UR tag and a chromosome.
        
    unique_pairs : a list of tuple
        a list of tuple where each element is unique tuple of CB, UR tag and a chromosome.
    """        
    reads_per_ub = {}
    unique_pairs = set()
    
    for read in sam.fetch(until_eof = True):
        chrom = read.reference_name
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
    return (reads_per_ub, unique_pairs)
        
def generate_partial_input(number, bam_template, out_template, split_read_template):
    """
    Parameters
    ----------
    number : character
        chromosome number.
    
    bam_template : string 
        input bam template. You have to concatenate chromosome number and make a full input directory. 
    
    out_template : string
        output bam template that is splitted and deduplicated. 
        You have to concatenate chromosome number and make a full splitted & deduplicated output directory. 
    
    split_read_template : string
        output bam template that is splitted only.
        You have to concatenate chromosome number and make a full splitted output directory.
    
    Returns
    -------
    sam : bam file
        an input raw bam file.
    
    out_name : string
        a full splitted & deduplicated output directory. 
    
    split_out_name : string
        a full splitted output directory.
    """    

    sam = pysam.AlignmentFile(bam_template, "rb")
    out_name = out_template + number + '_output.bam'        
    split_out_name = split_read_template + number + '_output.bam'
    return sam, out_name, split_out_name

def get_dedup_args():
    
    parser = argparse.ArgumentParser(description="deduplication")
    parser.add_argument('--bam_template', dest = 'bam_template',
                        required = True,
                        help = 'BAM input template')
    
    parser.add_argument('--out_template', dest = 'out_template',
                        required = True,
                        help = 'deduplication output template')
    
    parser.add_argument('--span_threshold', type = int, dest = 'span_threshold',
                        required = True,
                        help = 'span_threshold for splitting a cluster')
    
    parser.add_argument('--split_read_template', dest = 'split_read_template',
                        required = True,
                        help = 'split_read_template')    
    args = parser.parse_args()
    return args

def run_process():
    
    args = get_dedup_args()
    
    b_template = args.bam_template
    o_template = args.out_template
    number = b_template.split('.')[0].split('_')[-1]
    thres = args.span_threshold
    s_read_template = args.split_read_template
    
    o_mode = "wb"
    
    bam, o_name, split_out_name = generate_partial_input(number, b_template, o_template, s_read_template)
    print("successfully initialized for chr")
    
    r_per_ub, u_pairs = get_all_reads_per_umi (bam)
    print("successfully done get_all_reads_per_umi for chr" + str(number))
    
    resulting_reads, splitted_dicts, cbs = purify_UR(r_per_ub, u_pairs, thres)
    print("successfully got splitted and deduped reads for chr" + str(number))
    
    write_output(resulting_reads, o_name, o_mode, bam)
    print("successfully saved splitted and deduped reads for chr" + str(number))
    
    total_splitted_reads = get_splitted_reads(splitted_dicts, cbs)
    print("successfully got splitted reads for chr" + str(number))
    
    write_output(total_splitted_reads, split_out_name, o_mode, bam)
    print("successfully saved splitted reads for chr" + str(number))
    
if __name__ == "__main__":
    run_process()
    print("success")
    
