# -*- coding: utf-8 -*-
"""
Created on Fri Feb 24 17:10:42 2023


@author: Youngbin Moon (y.moon@unibas.ch)
"""
import pandas as pd
import argparse
import pysam








"""
Aim : Using a set of bed files (annotated.bed, unannotated.bed) from a sample,
1. separate all polyA reads mapping to terminal exons from a sample into sub-classes (annotated, unannotated)
2. and save them as bam files.
"""

def write_bam(final_reads, o_name, o_mode, bam):
    """
    Parameters
    ----------
    final_reads : list
        a list of reads that belongs to a specific class (annotated, unannotated) or remaining reads.
    
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
        
def get_cleavage_site(current_read, use_FC):
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
            # assert(refEnd == int(current_read.get_tag('OC')))
            assert(refEnd + 1 == int(current_read.get_tag('XO')))
            end = int(current_read.get_tag('XO'))
        
        elif use_FC:
            end = int(current_read.get_tag('XF'))
    
    return chrom, rev, end

def separate_polyAT(polyA_reads, bed, use_fixed_cs, polyA_cluster_id):
    """
    Parameters
    ----------    
    polyA_reads : list
        a list which contains all/partial polyA reads mapping to terminal exons in pysam object format.
    
    bed : bed file
        a bed file that is either annotated.bed or unannotated.bed from sample.
    
    use_fixed_cs : bool
        whether use a fixed cleavage site or not.
    
    polyA_cluster_id : string
        a string used when labeling a read with a tag which indicates 
        which class a read belongs to. ("annotated", "unannotated")
        
    Returns
    -------        
    final reads : list
        a list which contains reads that belong to a specific class (annotated. unannotated).
    
    leftover_reads : list
        a list which contains reads that do not belong to a specific class.
    """       
    final_reads = []
    leftover_reads = []
    overlap_count = 0
    for read in polyA_reads:
        chrom, rev, end = get_cleavage_site(read, use_fixed_cs)
        # use copy() to avoid any misbehaviour
        filtered_bed = bed[(bed['seqid'] == chrom) & (bed['start'] <= end) & (bed['end'] >= end) & (bed['strand'] == rev)].copy()
        
        if len(filtered_bed) == 1:
            sequence_id = filtered_bed['seqid'].values[0]
            start_pos = filtered_bed['start'].values[0]
            end_pos = filtered_bed['end'].values[0]
            direction = filtered_bed['strand'].values[0]
            c_id = filtered_bed['id'].values[0]
            
            # cluster_id = str(sequence_id) + str(':') + str(start_pos) + str(':') + str(end_pos) + str(':') + str(direction)\
            # + str(':') + str(c_id)
            # pA_id = polyA_cluster_id + '_' + cluster_id
            
            pA_id = polyA_cluster_id
            # set sub-cluster ID
            read.set_tag("Zi", pA_id)
            # set boolean mark that specifies whether a read is located in both clusters or not. 
            # 0 -> a read is located in 1 cluster.
            # 1 -> a read is located in both clusters.
            read.set_tag("Zd", str(0))
            
            final_reads.append(read)
            
        elif len(filtered_bed) == 2:
            # a read is mapping between 2 clusters exactly
            print('overlapping cluster......')
            print(str(filtered_bed))
            overlap_count += 1
            ind = list(filtered_bed['score']).index(max(list(filtered_bed['score'])))
            further_filtered_bed = filtered_bed.iloc[ind]

            sequence_id = further_filtered_bed['seqid']
            start_pos = further_filtered_bed['start']
            end_pos = further_filtered_bed['end']
            direction = further_filtered_bed['strand']
            c_id = further_filtered_bed['id']
            
            # cluster_id = str(sequence_id) + str(':') + str(start_pos) + str(':') + str(end_pos)\
            # + str(':') + str(direction) + str(':') + str(c_id)
            # pA_id = polyA_cluster_id + '_' + cluster_id
            
            pA_id = polyA_cluster_id
            read.set_tag("Zi", pA_id)
            
            read.set_tag("Zd", str(1))
            
            final_reads.append(read)
        
        elif len(filtered_bed) > 2:
            print('really wierd, not expected to happen......')
        
        else:
            leftover_reads.append(read)
    
    print('overlapping count: ' + str(overlap_count))
    return final_reads, leftover_reads

def bam_to_list(sam):
    """
    Parameters
    ----------    
    sam : bam file
        the original bam file (input)
        which contains all polyA reads mapping to terminal exons.
        
    Returns
    -------        
    final reads : list
        a list which contains all polyA reads mapping to terminal exons in pysam object format.
    """          
    final_reads = []
    for read in sam.fetch():
        final_reads.append(read)

    return final_reads

def get_inputs():
    parser = argparse.ArgumentParser(description="From polyA + T bam file, split into annotated.bam or unannotated.bam")
    parser.add_argument('--bam', dest = 'bam',
                        required = True,
                        help = 'input bam file containing all polyA reads mapping to terminal exons')
    
    parser.add_argument('--input_bed', dest = 'input_bed',
                        required = True,
                        help = 'annotated bed file or unannotated bed file')    
    
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
    input_bed = pd.read_csv(input_bed_dir, delimiter = '\t', header = None)
    input_bed.columns = ['seqid', 'start', 'end', 'id', 'score', 'strand']
    
    use_fc = bool(args.use_fc)
    class_reads = args.class_reads

    number = bam_dir.split('.')[0].split('_')[-2]    
    out_bam = args.out_bam + '_' + number + '.bam'
    remaining_bam = args.remaining_bam + '_' + number + '.bam'
    
    return bam, input_bed, use_fc, out_bam, remaining_bam, class_reads

def run_process():
    
    bam, input_bed, use_fc, out_bam, remaining_bam, class_reads = get_inputs()
    print('successfully got inputs')
    
    all_polyAT_reads = bam_to_list(bam)
    print('successfully converted bam to a list')
    
    classified_reads, remaining_reads = separate_polyAT(all_polyAT_reads, input_bed, use_fc, class_reads)
    print('successfully separated reads')

    out_mode = "wb"
    
    write_bam(classified_reads, out_bam, out_mode, bam)
    print('successfully wrote classified reads')

    write_bam(remaining_reads, remaining_bam, out_mode, bam)
    print('successfully wrote remaining reads')
        
if __name__ == "__main__":
    run_process()
    print('success')