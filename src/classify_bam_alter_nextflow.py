# -*- coding: utf-8 -*-
"""
Created on Fri Feb 24 17:10:42 2023


@author: Youngbin Moon (y.moon@unibas.ch)
"""










import pandas as pd
import argparse
import pysam
import multiprocessing.pool
"""
Aim : Using a set of bed files (terminal exons.bed, exons.bed, intronic.bed or intergenic.bed from sample) from a sample,
separate all polyA reads from a sample into specific classes (TE, exonic, intronic and intergenic) and save them as bam files.
"""
def write_output(final_reads, o_name, o_mode, bam):
    """
    Parameters
    ----------
    final_reads : list
        a list of reads that belongs to a specific class (TE, exonic, intronic, intergenic)
    
    out_name : string  
        output file name.
    
    out_mode : string
        output format. 'wb' refers to wrtie in bam file.
    
    bam : bam file
        the original bam file (input)
        This is to use original bam file as a template.
        
    Returns
    -------        
    returns nothing but saves the output in the bam format.
    """              
    outfile = pysam.AlignmentFile(o_name, o_mode, template=bam)
    for read in final_reads:
        outfile.write(read)

def get_full_reads(filtered, remaining):
    """
    Parameters
    ----------
    filtered: a list of list
        list where each element is a list of classified reads
    
    remaining: a list of list
        list where each element is a list of un-classified reads
        
    Returns
    -------        
    classified_reads: a list
        a list containing all classified reads of an input bam file

    remaining_reads: a list
        a list containing all un-classified reads of an input bam file
        
    """    
    classified_reads = []
    remaining_reads = []
    for filtered_list in filtered:
        classified_reads += filtered_list
    
    for remaining_list in remaining:
        remaining_reads += remaining_list
    
    return classified_reads, remaining_reads

def get_cs(current_read, use_FC):
    """
    Parameters
    ----------    
    current_read : pysam object
        a current pysam object
    
    use_FC : bool
        whether use a fixed cleavage site or not.
        
    Returns
    -------        
    chrom : string
        chromosome where a read maps to. (e.g. 'chr19')
    
    rev : string
        a direction of a strand in which a read maps to (e.g. '+')
    
    end : int
        cleavage site of a read.
    """        
    chrom = current_read.reference_name
    rev = current_read.is_reverse
    alignedRefPositions = current_read.get_reference_positions()
    refStart = alignedRefPositions[0]
    refEnd = alignedRefPositions[-1]    
    
    if rev == True:
        rev = '-'
        if not use_FC:
            # assert(refStart == int(current_read.get_tag('OC')))
            assert(refStart + 1 == int(current_read.get_tag('XO')))
            end = int(current_read.get_tag('XO'))
        
        elif use_FC:
            end = int(current_read.get_tag('XF'))      
    else:
        rev = '+'
        if not use_FC:
            # assert(refEnd == int(current_read.get_tag('XO')))
            assert(refEnd + 1 == int(current_read.get_tag('XO')))
            end = int(current_read.get_tag('XO'))
        
        elif use_FC:
            end = int(current_read.get_tag('XF'))
    
    return chrom, rev, end

def classify_bam(polyA_reads, bed, use_fixed_cs, polyA_cluster_id):
    """
    Parameters
    ----------    
    polyA_reads : list
        a list which contains all/partial polyA reads in pysam object format.
    
    bed : bed file
        a bed file that is either terminal exons.bed, exons.bed, intronic.bed or intergenic.bed from sample.
    
    use_fixed_cs : bool
        whether use a fixed cleavage site or not.
    
    polyA_cluster_id : string
        a string used when labeling a read with a tag which indicates 
        which class a read belongs to. ("TE", "exonic", "intronic", "intergenic")
        
    Returns
    -------        
    final reads : list
        a list which contains reads that belong to a specific class (polyA_cluster_id)
    
    leftover_reads : list
        a list which contains reads that do not belong to a specific class (polyA_cluster_id)
    """       
    final_reads = []
    leftover_reads = []
    overlap_count = 0
    for read in polyA_reads:
        chrom, rev, end = get_cs(read, use_fixed_cs)
        filtered_bed = bed[(bed['seqid'] == chrom) & (bed['start'] <= end) & (bed['end'] >= end) & (bed['strand'] == rev)].copy()
        
        if len(filtered_bed) == 1:
            # e.g. intergenic_chr19:4833275:+
            sequence_id = filtered_bed['seqid'].values[0]
            start_pos = filtered_bed['start'].values[0]
            end_pos = filtered_bed['end'].values[0]
            direction = filtered_bed['strand'].values[0]
            c_id = filtered_bed['id'].values[0]
            score = filtered_bed['score'].values[0]
            
            cluster_id = str(sequence_id) + str(':') + str(start_pos) + str(':') + str(end_pos) + str(':') + str(direction)\
            + str(':') + str(c_id)
            
            pA_id = polyA_cluster_id + '_' + cluster_id
            # set cluster ID
            read.set_tag("ZI", pA_id)
            # set score of the cluster
            read.set_tag("ZS", int(score))
            # set boolean mark that specifies whether a read is located in both clusters or not. 
            # 0 -> a read is located in 1 cluster.
            # 1 -> a read is located in both clusters.
            read.set_tag("ZD", 0)
            
            final_reads.append(read)
            
        elif len(filtered_bed) == 2:
            # a read is mapping between 2 clusters exactly
            print('overlapping cluster......')
            overlap_count += 1
            ind = list(filtered_bed['score']).index(max(list(filtered_bed['score'])))
            further_filtered_bed = filtered_bed.iloc[ind]
            
            sequence_id = further_filtered_bed['seqid']
            start_pos = further_filtered_bed['start']
            end_pos = further_filtered_bed['end']
            direction = further_filtered_bed['strand']
            c_id = further_filtered_bed['id']
            score = further_filtered_bed['score']
            
            cluster_id = str(sequence_id) + str(':') + str(start_pos) + str(':') + str(end_pos)\
            + str(':') + str(direction) + str(':') + str(c_id)
            
            pA_id = polyA_cluster_id + '_' + cluster_id
            
            read.set_tag("ZI", pA_id)
            
            read.set_tag("ZS", int(score))
            
            read.set_tag("ZD", 1)
            
            final_reads.append(read)
        
        elif len(filtered_bed) > 2:
            print('really wierd, not expected to happen......')
        
        else:
            leftover_reads.append(read)
    
    print('overlapping count: ' + str(overlap_count))
    return final_reads, leftover_reads

def convert_bam_to_list(sam):
    """
    Parameters
    ----------    
    sam : bam file
        the original bam file (input)
        which contains all polyA reads.
        
    Returns
    -------        
    final reads : list
        a list which contains all polyA reads in pysam object format.
    """      
    final_reads = []
    for read in sam.fetch():
        final_reads.append(read)

    return final_reads
        
def get_inputs():
    parser = argparse.ArgumentParser(description="From polyA bam file, split into bam file with different classes")
    parser.add_argument('--bam', dest = 'bam',
                        required = True,
                        help = 'input bam file containing all polyA reads')
    
    parser.add_argument('--input_bed', dest = 'input_bed',
                        required = True,
                        help = 'bed file containing clusters of pA sites in a specific sample and class')    

    parser.add_argument('--use_fc', type = int, dest = 'use_fc',
                        required = True,
                        help = 'whether use fixed cleavage site or not')  
    
    parser.add_argument('--out_bam', dest = 'out_bam',
                        required = True,
                        help = 'output bam file name template')

    parser.add_argument('--remaining_bam', dest = 'remaining_bam',
                        required = True,
                        help = 'remainder bam file name template')
    
    parser.add_argument('--class_reads', dest = 'class_reads',
                        required = True,
                        help = 'class of reads that you want to separate')
        
    args = parser.parse_args()
    
    bam_dir = args.bam
    bam = pysam.AlignmentFile(bam_dir, "rb")
    
    input_bed_dir = args.input_bed
        
    try:
        input_bed = pd.read_csv(input_bed_dir, delimiter = '\t', header = None)
        input_bed.columns = ['seqid', 'start', 'end', 'id', 'score', 'strand']
    
    except pd.errors.EmptyDataError:
        input_bed = pd.DataFrame()
        
    use_fc = bool(args.use_fc)
    class_reads = args.class_reads
    
    if class_reads == 'TE':
        number = bam_dir.split('.')[0].split('_')[-1]
    else:
        number = bam_dir.split('.')[0].split('_')[-2]
        
    out_bam = args.out_bam + '_' + number + '.bam'
    remaining_bam = args.remaining_bam + '_' + number + '.bam'

    return bam, input_bed, use_fc, out_bam, remaining_bam, class_reads

def run_process():
    
    bam, input_bed, use_fc, out_bam, remaining_bam, class_reads = get_inputs()
    print('successfully got inputs')
            
    all_polyA_reads = convert_bam_to_list(bam)
    print('successfully converted bam to a list')
    
    len_all = len(all_polyA_reads)
    print("input bed: " + str(len(input_bed)))
    print("all_polyA_readsL " + str(len_all))
   
    if len(input_bed) == 0:
        classified_reads = []
        remaining_reads = all_polyA_reads
    
    else:
        # create a default thread pool
        pool = multiprocessing.pool.ThreadPool()
        n_cores = pool._processes
        
        sublist_length = max(int(len_all/n_cores), 1000)
        print('cores: ' + str(n_cores))
        print('sublist_length: ' + str(sublist_length))
        sublists = [all_polyA_reads[i:i + sublist_length] for i in range(0, len(all_polyA_reads), sublist_length)]
        
        # call apply_async() without callback
        result_objects = [pool.apply_async(classify_bam, args = (sublist, input_bed, use_fc, class_reads)) for sublist in sublists]       

        final_reads_list = [r.get()[0] for r in result_objects]
        leftover_reads_list = [r.get()[1] for r in result_objects]
        pool.close()
        pool.join()
        
        classified_reads, remaining_reads = get_full_reads(final_reads_list, leftover_reads_list)
    
    assert(len_all == len(classified_reads) + len(remaining_reads))
    print('successfully separated reads')
        
    out_mode = "wb"
    
    write_output(classified_reads, out_bam, out_mode, bam)
    print('successfully wrote classified reads')
    
    write_output(remaining_reads, remaining_bam, out_mode, bam)
    print('successfully wrote remaining reads')
            
if __name__ == "__main__":
    run_process()
    print('success')