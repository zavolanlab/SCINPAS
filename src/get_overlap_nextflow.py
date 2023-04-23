# -*- coding: utf-8 -*-
"""

Created on Sun Dec 19 22:55:32 2021

@author: Youngbin Moon (y.moon@unibas.ch)
"""
import pandas as pd
import matplotlib.pyplot as plt
import argparse
import multiprocessing as mp








"""
Aim1 : compare how many overlaps found between 
all PAS clusters with score >= k of all available samples from the dataset (e.g. Tabula Muris) and clusters of polyA catalog.

Aim2 : visualize it with a barchart.
"""

def autolabel(rects, ax):
    """
    Parameters
    ----------
    rects : matplotlib.container.BarContainer
        Bars representing the read counts in each class.
        
    ax : matplotlib.axes._subplots.AxesSubplot
        The figure panel.

    Returns
    -------
    returns nothing but it attaches a text label above each bar displaying its height. (to the figure panel)
    """
    variation = 1.05
    for rect in rects:
        height = rect.get_height()
        # %s and %d are placeholders for a string and a number respectively. 
        ax.text(rect.get_x() + rect.get_width()/2., variation*height,
                '%d' % float(height),
                ha='center', va='bottom')

def plot_bar(names, counts, file_template):
    """
    Parameters
    ----------
    names : list
        sorted list of each class name (total, overlap, non_overlap)
        
    counts : list
        sorted list of the number of PAS in each class (total, overlap, non_overlap)
        
    file_template : string
        output file template

    Returns
    -------
    returns nothing but draws a bar plot.
    """
    fig, axis = plt.subplots()
    rectangles = axis.bar(names, counts)
    autolabel(rectangles, axis)
    
    plt.xticks(
    rotation=45, 
    horizontalalignment='right',
    fontweight='light',
    fontsize='x-large'  
    )
    
    plt.tight_layout()
    
    plt.yticks(fontsize='x-large')
    plt.ylabel('Number of PAS', fontsize = 'x-large')
    plt.rcParams['font.family'] = "DejaVu Sans"
    
    file1 = file_template + ".png"
    file2 = file_template + ".svg"
    
    plt.savefig(file1, bbox_inches='tight')
    plt.savefig(file2, bbox_inches='tight')
    
def find_overlap_PAS(sample_bed, catalog_bed):
    """
    Parameters
    ----------
    sample_bed : a bed file
        a partial bed file containing PAS clusters of all available samples from the dataset (e.g. Tabula Muris).
        (partial means: 1. PAS filtered by score >= k; 2. PAS from a specific chromosomes)
        
    catalog_bed : a bed file
        a bed file containing all PAS from the polyA catalog
                
    Returns
    -------
    overlap : int
        the number of PAS clusters that overlap with the catalog
    
    no_overlap : int
        the number of PAS clusters that do not overlap with the catalog.
    """    
    no_overlap = 0
    overlap = 0
    
    if len(sample_bed) == 0:
        overlap = 0
        no_overlap = 0
        
    elif len(sample_bed) != 0:
        for index, row in sample_bed.iterrows():
            # use cluster span to find which terminal exons that a cluster intersects with.
            chrom = row[0]
            cluster_start = int(row[1])
            cluster_end = int(row[2])
            cluster_id = row[3]
            score = int(row[4])
            direction = row[5]
            
            overlapped_PAS = catalog_bed[(catalog_bed['seqid'] == chrom) & (catalog_bed['start'] <= cluster_end)\
                                                     & (catalog_bed['end'] >= cluster_start) & (catalog_bed['strand'] == direction)].copy()
            
            if len(overlapped_PAS) == 0:    
                no_overlap += 1
             
            else:
                if len(overlapped_PAS) == 1 :
                    overlap += 1
                    print(overlapped_PAS)
                    
                elif len(overlapped_PAS) > 1 :
                    overlap += 1
        
    return (overlap, no_overlap)

def split_bed(subset_sample, chromosomes):
    """
    Parameters
    ----------
    subset_sample : a bed file
        partial bed file containing PAS clusters of all available samples from the dataset. (e.g. Tabula Muris)
        partial means: the bed file is filtered by score >= k. (i.e. PAS has to be supported by more than k reads)
    
    chromosomes : a list
        a list of chromosomes of the species.
        
    Returns
    -------
    splitted_bed: a bed file
        partial bed file containing PAS clusters of all available samples from the dataset. (e.g. Tabula Muris)
        (partial means: 1. PAS filtered by score >= k; 2. PAS from a specific chromosomes)
    """       
    splitted_beds = []
    
    for chrom in chromosomes:
        sample_per_chrom = subset_sample[subset_sample['seqid'] == chrom].copy()
        print(sample_per_chrom)
        splitted_beds.append(sample_per_chrom)
        
    return splitted_beds

def get_args():
        
    parser = argparse.ArgumentParser(description="get how many overlaps are found between PAS clusters and clusters of polyA catalog")

    parser.add_argument('--pas_bed_dir', dest = 'pas_bed_dir',
                        required = True,
                        help = 'directory towards PAS bedfile of our pipeline')
            
    parser.add_argument('--catalog_dir', dest = 'catalog_dir',
                        required = True,
                        help = 'directory towards catalog')

    parser.add_argument('--overlap_barchart', dest = 'overlap_barchart',
                        required = True,
                        help = 'output file name')

    parser.add_argument('--species', dest = 'species',
                        required = True,
                        help = 'species of the sample')
 
    parser.add_argument('--n_cores', dest = 'n_cores',
                        required = True,
                        help = 'number of cores')

    parser.add_argument('--threshold', dest = 'threshold',
                        required = True,
                        help = 'minimum number of reads that should support a cluster in order to consider as a true PAS')
    
    args = parser.parse_args()
    
    pas_bed_dir = args.pas_bed_dir
    catalog_dir = args.catalog_dir
    overlap_barchart = args.overlap_barchart
    species = args.species
    n_cores = int(args.n_cores)
    threshold = int(args.threshold)
    
    sample = pd.read_csv(pas_bed_dir, delimiter = '\t', header = None, low_memory=False)
    sample.columns = ['seqid', 'start', 'end', 'id', 'score', 'direction'] 
    
    catalog = pd.read_csv(catalog_dir, delimiter = '\t', header = None, low_memory=False)
    catalog.columns = ['seqid', 'start', 'end', 'id', 'score', 'strand', 'percentage', 'num_protocol', 'avg_expression', 'annotation', 'info_polyA']

    catalog['seqid'] = catalog['seqid'].apply(lambda x: 'chr' + str(x))
    
    return sample, catalog, overlap_barchart, species, n_cores, threshold

def run_process():
    
    sample, catalog, overlap_barchart, species, n_cores, threshold = get_args()
    print('successfully got inputs')

    # filter PAS that is supported more than k reads
    subset_sample = sample[sample['score'] >= threshold].copy()
    
    if species == 'mouse':
        chromoSomes = ['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','X','Y']
        chromoSomes = ['chr' + str(elem) for elem in chromoSomes]

    elif species == 'human':
        chromoSomes = ['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19', '20', '21', '22', 'X','Y']
        chromoSomes = ['chr' + str(elem) for elem in chromoSomes]   
        
    splitted_beds = split_bed(subset_sample, chromoSomes)
    print('successfully got splitted sample bed')
    
    # set the number of cores available
    pool = mp.Pool(n_cores)
    
    # call apply_async() without callback
    result_objects = [pool.apply_async(find_overlap_PAS, args = (split_bed, catalog)) for split_bed in splitted_beds]

    # result_objects is a list of pool.ApplyResult objects
    overlap_results = [r.get()[0] for r in result_objects]
    no_overlap_results = [r.get()[1] for r in result_objects]
    
    pool.close()
    pool.join()
    
    total_overlap = sum(overlap_results)
    total_non_overlap = sum(no_overlap_results)
    print('successfully got results')
    
    total = len(subset_sample)
    total_PAS = total_overlap + total_non_overlap
    assert(total == total_PAS)
    print('parallelization results are correct and successful!')
    
    names_list = ['total', 'overlap', 'no_overlap']
    counts_list = [total, total_overlap, total_non_overlap]
    
    plot_bar(names_list, counts_list, overlap_barchart)
    print('successfully saved results')

if __name__ == "__main__":
    run_process()    
    print('success')
