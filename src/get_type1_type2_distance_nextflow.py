# -*- coding: utf-8 -*-
"""

Created on Mon Jan 10 10:26:42 2022

@author: Youngbin Moon (y.moon@unibas.ch)
"""
import pandas as pd
import argparse
import pysam
import numpy as np







"""
Aim : write a csv file where:

1st column = average distance from start of longest terminal exon to representative cleavage site of a current pA site.
2nd column = average distance from start of longest terminal exon to representative cleavage site of a current pA site.

NOTE: type1 average distance and type 2 average distance should be matched with the same gene.
"""
def write_out(final_df, out_file):
    """
    Parameters
    ----------
    final_df : dataframe
        dataframe where 
        1st column = average distance from start of longest terminal exon to representative cleavage site of a current pA site.
        2nd column = average distance from start of longest terminal exon to representative cleavage site of a current pA site.
        
        NOTE: type1 average distance and type 2 average distance should be matched with the same gene. 
        
    out_file : string
        output file name.

    Returns
    -------
    returns nothing but writes the output in the csv format.
    """
    final_df.to_csv(out_file, sep = '\t', header = False, index = False,
    columns = ['type1_distance', 'type2_distance'])
    
def make_distances_tuple(dictionary1, dictionary2):
    """
    Parameters
    ----------    
    dictionary1 : dictionary
        a dictionary for "type1 cell" where
        key = GX tag
        value = average distance from start of the longest terminal exon (of that gene) to representative cleavage site of a current pA site.    

    dictionary2 : dictionary
        a dictionary for "type2 cell" where
        key = GX tag
        value = average distance from start of the longest terminal exon (of that gene) to representative cleavage site of a current pA site. 
        
    Returns
    -------
    final_df : dataframe
        dataframe where 
        1st column = average distance from start of longest terminal exon to representative cleavage site of a current pA site.
        2nd column = average distance from start of longest terminal exon to representative cleavage site of a current pA site.
        
        NOTE: type1 average distance and type 2 average distance should be matched with the same gene.
    """       
    distance_tuples_list = []
    if len(list(dictionary1.keys())) >= len(list(dictionary2.keys())):
        for GX_tag in dictionary1.keys():
            if GX_tag in dictionary2.keys():
                distance1 = dictionary1[GX_tag]
                distance2 = dictionary2[GX_tag]
                distance_tuples_list.append((distance1, distance2))
            
    elif len(list(dictionary2.keys())) >= len(list(dictionary1.keys())):
        for GX_tag in dictionary2.keys():
            if GX_tag in dictionary1.keys():
                distance1 = dictionary1[GX_tag]
                distance2 = dictionary2[GX_tag]                
                distance_tuples_list.append((distance1, distance2))
    
    final_df = pd.DataFrame(distance_tuples_list, columns = ['type1_distance', 'type2_distance'])
    return final_df
    
def start_to_end_distance(current_cleavage_site, start, end, rev):
    """
    Parameters
    ----------    
    current_cleavage_site : int
        representative cleavage site of a current pA site (cluster). 
    
    start : int
        start position of the longest terminal exon that the current pA site maps to.        

    end : int
        end position of the longest terminal exon that the current pA site maps to. 
        
    use_FC : bool
        whether use a fixed cleavage site or not.  
        
    Returns
    -------
    distance : int
        distance from start of longest terminal exon to representative cleavage site of a current pA site.
    """       
    if rev == '-':
        terminal_exon_start = end
    
    elif rev == '+':
        terminal_exon_start = start
                
    distance = abs(current_cleavage_site - terminal_exon_start)
    
    return distance
    
def get_longest_te_positions(filtered_terminal_Exons):
    """
    Parameters
    ----------    
    filtered_terminal_Exons : bed file
        bed file that contains 1 or more terminal exons of the same gene where current pA site maps to (from gtf).
        
    Returns
    -------
    start_longest_te: start position of the longest terminal exon where current pA site maps to.
    
    end_longest_te: end position of the longest terminal exon where current pA site maps to.
    """      
    start_positions = np.asarray(filtered_terminal_Exons['start'].copy())
    end_positions = np.asarray(filtered_terminal_Exons['end'].copy())
    assert(len(start_positions) == len(end_positions))
    
    spans = (end_positions - start_positions)
    assert(all(val > 0 for val in spans))
    longest_terminal_exon_idx = np.argmax(spans)
    
    # get the longest terminal exon span
    start_longest_te = start_positions[longest_terminal_exon_idx]
    end_longest_te = end_positions[longest_terminal_exon_idx]
    
    return start_longest_te, end_longest_te

def get_GX_distance_dict(dictionary, terminal_exons):
    """
    Parameters
    ----------    
    dictionary : dictionary
        a dictionary of tuples where each elements is (cluster_id, cluster_start, cluster_end, 1/n_genes).
        cluster_id : id of a pA site
        cluster_start : start position of a pA site
        cluster_end : end position of a pA site
        1/n_genes : probability that a pA cluster is mapping to this gene. (in case this pA cluster is mapping to multiple genes)
                    max value == 1.
                    n_genes is the number of genes that a current pA cluster is mapping to.

    terminal_exons : dataframe
        a dataframe containing all terminal exons from gtf.
        
    Returns
    -------
    GX_distance_dict : dictionary
        a dictionary where key = gene and value = average distance between pA site and start of the terminal exon.
        The average distance is computed because for a given terminal exon of a gene, multiple pA sites can map to different locations of the same terminal exon.
    """       
    GX_distance_dict = {}
    for gene in dictionary.keys():
        read_tuples = dictionary[gene]
        frequencies = [elem[3] for elem in read_tuples]
        tot_frequency = sum(frequencies)
        avg_distance = 0
        
        # each read tuple represents a single pA cluster
        for read_tuple in read_tuples:
            cluster_id = read_tuple[0]
            
            chrom = cluster_id.split(':')[0]
            cs = int(cluster_id.split(':')[1])
            direction = cluster_id.split(':')[2]
            
            cluster_start = int(read_tuple[1])
            cluster_end = int(read_tuple[2])
            
            # 0 < frequency <=1
            frequency = float(read_tuple[3])
            # convert frequency into probability so that we compute average easily.
            probability = frequency/tot_frequency
            
            # find all terminal exons of a particular gene only that overlaps with the pA cluster.
            filtered_terminal_exons = terminal_exons[(terminal_exons['seqid'] == chrom) & (terminal_exons['start'] <= cluster_end)\
                                                     & (terminal_exons['end'] >= cluster_start) & (terminal_exons['strand'] == direction)\
                                                         & (terminal_exons['id'] == gene)].copy()   
            
            print(filtered_terminal_exons)    
            # whther a current pA site maps to single or multiple terminal exons, always choose the longest terminal exon.
            start_longest_te, end_longest_te = get_longest_te_positions(filtered_terminal_exons)
            start_to_end_d = start_to_end_distance(cs, start_longest_te, end_longest_te, direction)
            avg_distance += start_to_end_d * probability
            
        GX_distance_dict[gene] = avg_distance
    
    return GX_distance_dict
            
def group_by_gene(cell_type_bed_input, terminal_exons):
    """
    Parameters
    ----------    
    cell_type_bed_input : string
        directory of clusters of a polyA reads mapping to terminal exons in a specific cell type (either cell type 1 or 2). 
    
    terminal_exons : dataframe
        a dataframe containing all terminal exons from gtf.
        
    Returns
    -------
    GX_read_dict : dictionary
        a dictionary of tuples where each elements is (cluster_id, cluster_start, cluster_end, 1/n_genes).
        cluster_id : id of a pA site
        cluster_start : start position of a pA site
        cluster_end : end position of a pA site
        1/n_genes : probability that a pA cluster is mapping to this gene. (in case this pA cluster is mapping to multiple genes)
                    max value == 1.
                    n_genes is the number of genes that a current pA cluster is mapping to.
    """      
    GX_read_dict = {}
    with open(cell_type_bed_input, "r") as t:
        for line in t:
            # use cluster span to find which terminal exons that a cluster intersects with.
            chrom = line.split()[0]
            cluster_start = int(line.split()[1])
            cluster_end = int(line.split()[2])
            cluster_id = line.split()[3]
            score = int(line.split()[4])
            direction = line.split()[5]
            
            # subset terminal exon bed file (gtf) with overlapping span.
            filtered_terminal_exons = terminal_exons[(terminal_exons['seqid'] == chrom) & (terminal_exons['start'] <= cluster_end)\
                                                     & (terminal_exons['end'] >= cluster_start) & (terminal_exons['strand'] == direction)].copy()
            
            genes = list(set(filtered_terminal_exons['id']))
            n_genes = len(genes)
            # frequency = 1/n_genes
            frequency = score/n_genes 
            # if a pA cluster maps to multiple genes, normalize it.
            for gene in genes:
                if gene not in GX_read_dict.keys():
                    GX_read_dict[gene] = []
                
                GX_read_dict[gene].append((cluster_id, cluster_start, cluster_end, frequency))
                
    return GX_read_dict

def run ():    
    parser = argparse.ArgumentParser(description="compute average distance between terminal exon and pA site in a specific cell type")
    parser.add_argument('--type1_polyA_bed', dest = 'type1_polyA_bed',
                        required = True,
                        help = 'directory of clusters of a polyA reads mapping to terminal exon in type 1 cells')
    
    parser.add_argument('--type2_polyA_bed', dest = 'type2_polyA_bed',
                        required = True,
                        help = 'directory of clusters of a polyA reads mapping to terminal exon in type 2 cells')
  
    parser.add_argument('--bed_input', dest = 'bed_input',
                        required = True,
                        help = 'bed file containing terminal exons')   
    
    parser.add_argument('--distance_out', dest = 'distance_out',
                        required = True,
                        help = 'distance distribution output csv file name')
                     
    args = parser.parse_args()
    
    type1_polyA_bed = args.type1_polyA_bed

    type2_polyA_bed = args.type2_polyA_bed
    
    bed_input = args.bed_input
    bedFile = pd.read_csv(bed_input, delimiter = '\t', header = None)
    bedFile.columns = ['seqid', 'start', 'end', 'id', 'score', 'strand']
    
    distance_out = args.distance_out

    return type1_polyA_bed, type2_polyA_bed, bedFile, distance_out
    
def run_process():
    type1_polyA_bed, type2_polyA_bed, bedFile, distance_out = run()
    print('successfully initialized')
    
    GX_read_dict_type1 = group_by_gene(type1_polyA_bed, bedFile)
    print('successfully generated GX read dict for type 1')
    
    GX_read_dict_type2 = group_by_gene(type2_polyA_bed, bedFile)
    print('successfully generated GX read dict for type 2')
    
    GX_distance_dict_type1 = get_GX_distance_dict(GX_read_dict_type1, bedFile)
    print('successfully dictionary for cell type1 where key = GX, value = start to end distance in terminal exon')
    print(str(GX_distance_dict_type1))

    GX_distance_dict_type2 = get_GX_distance_dict(GX_read_dict_type2, bedFile)
    print('successfully dictionary for cell type2 where key = GX, value = start to end distance in terminal exon')
    print(str(GX_distance_dict_type2))
    
    output_dataframe = make_distances_tuple(GX_distance_dict_type1, GX_distance_dict_type2)
    print('successfully converted into a list of tuple')
    print(str(output_dataframe))
    
    write_out(output_dataframe, distance_out)
    print('successfully wrote the output')
    
if __name__ == "__main__":
    run_process()
    print('success')