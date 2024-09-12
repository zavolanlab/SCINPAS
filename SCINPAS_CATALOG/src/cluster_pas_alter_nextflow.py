# -*- coding: utf-8 -*-
"""

Created on Sun Oct 22 22:55:32 2023

@author: Youngbin Moon (y.moon@unibas.ch)
"""











import argparse
import itertools
from pathlib import Path
from typing import Tuple
import numpy as np
import pandas as pd
from multiprocessing import Pool

def read_BED(file: Path) -> pd.DataFrame:
    """Read a BED file.
    
    Args:
        posix path to file

    Returns:
        pandas Dataframe
    """
    df = pd.read_csv(file, sep='\t', header = 0)
    return(df)

def write_BED(out_path: Path, df: pd.DataFrame) -> None:
    """Write a dataframe to BED.
    
    Args:
        dataframe containing the defined columns
    """
    df.to_csv(out_path, sep="\t", header=True, index=False)

def rank_PAS(df: pd.DataFrame) -> pd.DataFrame:
    """Rank the polyA sites by supp and score which are:
    
    supp = the number of samples (experiments) supporting that cleavage site.
    score = RPM
    
    Inplace sorting.

    Args:
        df: BED formatted file with column "score" and "supp".

    Returns:
        sorted df (by score and supp)
    """
    # df = df.sort_values(by=['score', 'supp'], ascending=False, inplace=False, ignore_index=True)
    df = df.sort_values(by=['supp', 'score'], ascending=False, inplace=False, ignore_index=True)
    return(df)

def _traverse_PAS(zipped: Tuple[Tuple[str, str], pd.DataFrame],
    dist_upstream: int, dist_downstream: int) -> pd.DataFrame:
    """Traverse df of PAS and cluster.
    """
    name = zipped[0]
    group = zipped[1]
    out = pd.DataFrame(columns=group.columns)
    
    if name == '+':
        dist_start = dist_upstream
        dist_end = dist_downstream
    else:
        # on minus strand: up- and downstream is switched
        dist_start = dist_downstream
        dist_end = dist_upstream

    # Sort PAS
    group = rank_PAS(group)
    
    total_df = pd.DataFrame(columns = ['seqid', 'start', 'end', 'id', 'score', 'strand', 'supp'])
    while len(group) > 0:
        # find representative site as first entry
        next_site = group.iloc[0]
        # group = group.iloc[1:]
        # get PAS in vicinity of next_site
        ind = (group.loc[:,'start'] >= next_site.start - dist_start) & (
            group.loc[:,'end'] <= next_site.end + dist_end)
        # construct one BED entry for the coordinates
        # chr, name, strand already given.
        # print('next_site' + str(next_site.to_frame().T))
                
        entry = next_site.copy()
        # print('entry:' + str(entry.to_frame().T))
        
        # entry positions and score updated here.
        # entry.to_frame().T converts entry to dataframe
        # you dont change cluster id because this cleavage site (with the highest score) 
        # will be representative cleavage site of a PAS.
        entry.start = np.min(group.loc[ind,'start'])
        entry.end = np.max(group.loc[ind,'end'])
        entry.score = np.sum(group.loc[ind, 'score'])
        # print('corrected entry:' + str(entry.to_frame().T))
        entry['id'] = str(entry['id']) + ':' + str(entry.start) + ':' + str(entry.end) + ':' + str(entry.score) + ':' + str(entry.supp)
        
        # selected group includes entry
        selected_grp = group.loc[ind, :].copy()
        
        # print('selected_grp: ' + str(selected_grp))
        # print('name of entry: ' + str(entry['id']))
        
        # update cleavage site id -> cluster id
        selected_grp['id'] = entry['id']
        # print('corrected_selected_grp: ' + str(selected_grp))
        
        total_df = pd.concat([total_df, selected_grp])
        # print('total_df: ' + str(total_df))
        
        out = pd.concat([out, entry.to_frame().T], ignore_index=True)
        # remove the current entries from the dataframe
        group = group.loc[-ind,:]
        
    return(out, total_df)

def cluster_PAS(df: pd.DataFrame, dist_upstream: 
    int, dist_downstream: int, n_cores: int) -> pd.DataFrame:
    """Obtain distance between adjacent PAS.

    Sort PAS by position, and compute distance
    to next PAS

    Args:
        df: BED file containing individual polyA sites
        dist_upstream: distance to look for PAS upstream.
        dist_downstream: distance to look for PAS downstream.
        n_cores: number of cores to use for multiprocessing.

    Returns:
        Pandas dataframe with PAS clusters
    """
    groupedByChrStrand = df.groupby(['seqid', 'strand'])
    
    # manual subsetting:
    # group = df[(df.loc[:,"chr"] == "chr19") & (df.loc[:,"strand"] == "+")]
    with Pool(n_cores) as pool:
        result = pool.starmap(_traverse_PAS ,zip(groupedByChrStrand, 
          itertools.repeat(dist_upstream),
          itertools.repeat(dist_downstream)))
    
    outs = [elem[0] for elem in result]
    modified_originals = [elem[1] for elem in result]
    
    outs = [elem[0] for elem in result]
    modified_originals = [elem[1] for elem in result]
    # print('outs:' + str(outs))
    # print('modified_originals: ' + str(modified_originals))
    
    total_outs = pd.concat(outs)
    # print('total_outs: ' + str(total_outs))
    
    total_modified_originals = pd.concat(modified_originals)
    print('total_modified_originals: ' + str(total_modified_originals))
    
    return (total_outs, total_modified_originals)

def arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Cluster PAS as in Gruber et al. 2016 for the polyAsite atlas.")

    parser.add_argument('--in_bed', dest='in_bed', 
      required = True,
      help='Input PAS in BED format.')
    
    parser.add_argument('--out', dest='out_file', 
      default = 'cluster_pas',
      help='Output PAS clusters in BED format. Default cluster_pas.')
    
    parser.add_argument('--du', dest='dist_upstream',
      default = 25, type=int,
      help='check up to INT distance upstream for other sites. Default: 25.')
    
    parser.add_argument('--dd', dest='dist_downstream',
      default = 25, type=int,
      help='check up to INT distance downstream for other sites. Default: 25.')
    
    parser.add_argument('--c', dest='n_cores',
      default=4, type=int,
      help='Number of cores for multiprocessing. Default: 4.')

    parser.add_argument('--num', dest='num', required = True,
      help='Chromosome number')

    parser.add_argument('--strand', dest='strand', required = True,
      help='direction of PAS')

    parser.add_argument('--cs_out', dest='cs_out', 
      default = 'modified_unique_cs',
      help='modified unique cleavage sites in BED format containing PAS cluster id. Default modified_unique_cs.')
    
    args = parser.parse_args()
    
    return(args)

def main():
    
    args = arguments()
    
    df = read_BED(args.in_bed)   
    
    clusters, modified_unique_cs = cluster_PAS(df, args.dist_upstream, args.dist_downstream, args.n_cores)
    clusters.sort_values(by=['seqid', 'start', 'end'], inplace=True, ignore_index=True)
    modified_unique_cs.sort_values(by=['seqid', 'start', 'end'], inplace=True, ignore_index=True)
    
    print('successfully generated clusters')
    
    chromosome = args.num
    direction = args.strand
    out_name = args.out_file + '_' + str(chromosome) + '_' + str(direction) + '.bed'
    
    write_BED(out_name, clusters)
    print('successfully saved cluster result for chromosome: ' + str(chromosome) + str(direction))
    
    out_name2 = args.cs_out + '_' + str(chromosome) + '_' + str(direction) + '.bed'

    write_BED(out_name2, modified_unique_cs)
    print('successfully saved modified unique cleavage sites for chromosome: ' + str(chromosome) + str(direction))
    
if __name__ == '__main__':
    main()
