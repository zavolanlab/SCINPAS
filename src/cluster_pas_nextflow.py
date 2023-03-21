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
    df = pd.read_csv(file, sep='\t',
                     names=['chr', 'start', 'end', 'name', 'score', 'strand'])
    return(df)


def write_BED(out_path: Path, df: pd.DataFrame) -> None:
    """Write a dataframe to BED.
    
    Args:
        dataframe containing the defined 6 columns
    """
    df.to_csv(out_path, sep="\t", header=False, index=False)


def rank_PAS(df: pd.DataFrame) -> pd.DataFrame:
    """Rank the polyA sites by score.

    That is, by the number of reads supporting the PAS.
    Inplace sorting.

    Args:
        df: BED formatted file with column "score".

    Returns:
        sorted df (by score)
    """
    df = df.sort_values(by=['score'], ascending=False, inplace=False, ignore_index=True)
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
    
    while len(group) > 0:
        # find representative site as first entry
        next_site = group.iloc[0]
        # group = group.iloc[1:]
        # get PAS in vicinity of next_site
        ind = (group.loc[:,'start'] >= next_site.start - dist_start) & (
            group.loc[:,'end'] <= next_site.end + dist_end)
        # construct one BED entry for the coordinates
        # chr, name, strand already given.
        entry = next_site.copy()
        entry.start = np.min(group.loc[ind,'start'])
        entry.end = np.max(group.loc[ind,'end'])
        entry.score = np.sum(group.loc[ind, 'score'])
        out = pd.concat([out, entry.to_frame().T], ignore_index=True)
        # remove the current entries from the dataframe
        group = group.loc[-ind,:]
    return(out)


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
    groupedByChrStrand = df.groupby(['chr', 'strand'])
    
    # manual subsetting:
    # group = df[(df.loc[:,"chr"] == "chr19") & (df.loc[:,"strand"] == "+")]
    with Pool(n_cores) as pool:
        result = pool.starmap(_traverse_PAS ,zip(groupedByChrStrand, 
          itertools.repeat(dist_upstream),
          itertools.repeat(dist_downstream)))

    return(pd.concat(result))

def arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Cluster PAS as in Gruber et al. 2016 for the polyAsite atlas.")
    parser.add_argument('in_bed',
      nargs='*', help='Input PAS in BED format.')
    parser.add_argument('--out', dest='out_file', 
      default = 'cluster_pas.bed',
      help='Output PAS clusters in BED format. Default cluster_pas.bed.')
    parser.add_argument('--du', dest='dist_upstream',
      default = 25, type=int,
      help='check up to INT distance upstream for other sites. Default: 25.')
    parser.add_argument('--dd', dest='dist_downstream',
      default = 25, type=int,
      help='check up to INT distance downstream for other sites. Default: 25.')
    parser.add_argument('--c', dest='n_cores',
      default=4, type=int,
      help='Number of cores for multiprocessing. Default: 4.')
    args = parser.parse_args()
    return(args)

def main():
    args = arguments()
    if len(args.in_bed) == 2:
        df = read_BED(args.in_bed[1])
    elif len(args.in_bed) > 2:
        # multiple BED files provided
        print(f"Concatenate BED files")
        df = pd.DataFrame()
        for bed in args.in_bed:
            print(f"file: {bed}")
            tmpdf = read_BED(bed)
            df = pd.concat([df, tmpdf])
        
        # alternative: include as separate samples

    clusters = cluster_PAS(df, args.dist_upstream, args.dist_downstream, args.n_cores)
    clusters.sort_values(by=['chr', 'start', 'end'], inplace=True, ignore_index=True)
    write_BED(args.out_file, clusters)

if __name__ == '__main__':
    main()
