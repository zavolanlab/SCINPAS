# -*- coding: utf-8 -*-
"""

Created on Sun Dec 19 22:55:32 2021

@author: Youngbin Moon (y.moon@unibas.ch)
"""
import pandas as pd
import argparse
import matplotlib.pyplot as plt
import multiprocessing as mp







"""
Aim1: compute % overlap between partial PAS clusters of all available samples from the dataset (e.g. Tabula Muris) and clusters of polyA catalog
(partial means: 1. PAS with top 100, top 500 etc scores; 2. PAS filtered by score >= k;)

Aim2 : visualize them with a barchart.
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

def plot_barchart(names, counts, file_template):
    """
    Parameters
    ----------
    names : list
        sorted list of each class name e.g.(top100, top500, top1000,.......... all)
        
    counts : list
        sorted list of the number of PAS in each class (% overlap of top100, % overlap of top500, % overlap of top1000,.......... % overlap of all PAS)
        
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
    plt.ylabel('percentage overlap', fontsize = 'x-large')
    plt.rcParams['font.family'] = "DejaVu Sans"
    
    file1 = file_template + ".png"
    file2 = file_template + ".svg"
    
    plt.savefig(file1, bbox_inches='tight')
    plt.savefig(file2, bbox_inches='tight')

def find_PAS_overlap(sample_bed, catalog_bed):
    """
    Parameters
    ----------
    sample_bed : a bed file
        a partial or full bed file containing PAS clusters of all available samples from the dataset. (e.g. Tabula Muris)
        (partial means: 1. PAS with top 100, top 500 etc scores; 2. PAS filtered by score >= k; 3. PAS from a specific chromosomes)
        
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
            
            # overlapped_PAS = catalog_bed[(catalog_bed['seqid'] == chrom) & (catalog_bed['start'] <= cluster_end)\
            #                                          & (catalog_bed['end'] >= cluster_start) & (catalog_bed['strand'] == direction)].copy()

            overlapped_PAS = catalog_bed[(catalog_bed['seqid'] == chrom) & (catalog_bed['start'] <= cluster_end + 1)\
                                                     & (catalog_bed['end'] >= cluster_start - 1) & (catalog_bed['strand'] == direction)].copy()                 
            if len(overlapped_PAS) == 0:    
                no_overlap += 1
             
            else:
                if len(overlapped_PAS) == 1 :
                    overlap += 1
                    print(overlapped_PAS)
                    
                elif len(overlapped_PAS) > 1 :
                    overlap += 1
    
    return (overlap, no_overlap)

def split_by_chromosomes(sample_bed, chromosomes):
    """
    Parameters
    ----------
    sample_bed : a bed file
        partial bed file containing PAS clusters of all available samples from the dataset. (e.g. Tabula Muris)
        (partial means: 1. PAS with top 100, top 500 etc scores; 2. PAS filtered by score >= k)
        
    chromosomes : a list
        a list of chromosomes of the species.
        
    Returns
    -------
    splitted_beds_by_chromosomes : a list
        a list of bed files where each element is a partial bed file containing PAS clusters of all available samples from the dataset.
        (partial means: 1. PAS with top 100, top 500 etc scores; 2. PAS filtered by score >= k; 3. PAS from a specific chromosomes)
    """        
    splitted_beds_by_chromosomes = []
    
    for chrom in chromosomes:
        sample_per_chrom = sample_bed[sample_bed['seqid'] == chrom].copy()
        print(sample_per_chrom)
        splitted_beds_by_chromosomes.append(sample_per_chrom)
        
    return splitted_beds_by_chromosomes

def make_subset(subset_sample):
    """
    Parameters
    ----------
    subset_sample : a bed file
        a partial bed file containing PAS clusters of all available samples from the dataset. (e.g. Tabula Muris)
        This bed file is filtered by score >= k. (i.e. PAS has to be supported by more than k reads)
                    
    Returns
    -------
    splitted_samples : a list
        a list of bed files where each element is a partial bed file containing PAS clusters of all available samples from the dataset.
        (partial means: 1. PAS with top 100, top 500 etc scores; 2. PAS filtered by score >= k)
    
    """      
    splitted_samples = []
    num_clusters = len(subset_sample)

    if num_clusters >= 20000:
        names_list = ["top100", "top500", "top1000", "top5000", "top10000", "top20000", "all"]
        sample_sizes = [100, 500, 1000, 5000, 10000, 20000]
        
    elif num_clusters >= 10000:
        names_list = ["top100", "top500", "top1000", "top5000", "top10000", "all"]
        sample_sizes = [100, 500, 1000, 5000, 10000]
        
    elif num_clusters >= 5000:
        names_list = ["top100", "top500", "top1000", "top5000", "all"]
        sample_sizes = [100, 500, 1000, 5000]
    
    elif num_clusters >= 1000:
        names_list = ["top100", "top500", "top1000", "all"]
        sample_sizes = [100, 500, 1000]

    elif num_clusters >= 500:
        names_list = ["top100", "top500", "all"]
        sample_sizes = [100, 500]
    
    else:
        names_list = ["top100", "all"]
        sample_sizes = [100]    
            
    for sample_size in sample_sizes:
        splitted_samples.append(subset_sample[0:sample_size].copy())
                
    splitted_samples.append(subset_sample.copy())
    
    return splitted_samples, names_list

def get_args():
        
    parser = argparse.ArgumentParser(description="get how many overlaps are found between top-scored PAS clusters and clusters of polyA catalog")

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
    sample.sort_values(by=['score'], ascending = False, inplace=True, ignore_index=True)
    
    return sample, catalog, overlap_barchart, species, n_cores, threshold

def run_process():
    
    sample, catalog, overlap_barchart, species, n_cores, threshold = get_args()
    print('successfully got inputs')

    if species == 'mouse':
        chromoSomes = ['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','X','Y']
        chromoSomes = ['chr' + str(elem) for elem in chromoSomes]

    elif species == 'human':
        chromoSomes = ['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19', '20', '21', '22', 'X','Y']
        chromoSomes = ['chr' + str(elem) for elem in chromoSomes] 

    # filter PAS that is supported more than k reads
    subset_sample = sample[sample['score'] >= threshold].copy()
        
    top_splitted_samples, names_list = make_subset(subset_sample)
    print('successfully got top 100, top 200 samples......etc')
    
    percentages =[]
    for top_split_sample in top_splitted_samples:
        top_split_sample_by_chrs = split_by_chromosomes(top_split_sample, chromoSomes)
        # set the number of cores available
        pool = mp.Pool(n_cores)
        
        # call apply_async() without callback
        result_objects = [pool.apply_async(find_PAS_overlap, args = (top_split_bed_by_chr, catalog)) for top_split_bed_by_chr in top_split_sample_by_chrs]
           
        # result_objects is a list of pool.ApplyResult objects
        # the number of PAS that overlaps with catalog (1st element: chr1.........., 19th element = chr19..........)
        overlap_results = [r.get()[0] for r in result_objects]
        no_overlap_results = [r.get()[1] for r in result_objects]
        
        pool.close()
        pool.join()
        
        partial_overlap = sum(overlap_results)
        partial_no_overlap = sum(no_overlap_results)
        
        partial_total = partial_overlap + partial_no_overlap
        assert(partial_total == len(top_split_sample))
        
        percentage_overlap = (partial_overlap*100)/partial_total
        percentages.append(percentage_overlap)
    
    assert(len(names_list) == len(percentages))
    
    plot_barchart(names_list, percentages, overlap_barchart)
    print('successfully generated bar chart')

if __name__ == "__main__":
    run_process()    
    print('success')