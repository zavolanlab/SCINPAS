# -*- coding: utf-8 -*-
"""
Created on Mon Jul 18 19:25:43 2024

@author: Youngbin Moon (y.moon@unibas.ch)
"""







import pandas as pd
import pybedtools
import argparse
from multiprocessing import Pool
import os
import tempfile
import shutil
import gc

def write_result(df, out_name):
    df.to_csv(out_name, sep = '\t', header = True, index = False) 
    
def change_format(df):
    intermediate_df = df[['o_seqid', 'o_start', 'o_end', 'o_id', 'o_score', 'o_strand']]
    final_df  = intermediate_df.rename(columns = {'o_seqid': 'seqid', 'o_start': 'start', 'o_end': 'end', 'o_id': 'id', 'o_score': 'score', 'o_strand': 'strand'})
    return final_df

def non_overlapping_region(df, base_temp_dir):

    print('subset: ' + str(df))
    genes = pybedtools.BedTool.from_dataframe(df)
    results = []
    for gene in genes:
        # Create a unique temporary directory within the base_temp_dir
        temp_dir = tempfile.mkdtemp(dir=base_temp_dir)
        pybedtools.helpers.set_tempdir(temp_dir)      
        # taking a single genomic interval (gene), converting it to a string that represents the interval in standard BED format, 
        # and then creating a new BedTool object from that string. 
        # The new BedTool object (gene_bed) is now ready to be used in subsequent operations (intersections, merging, and subtractions) with other genomic intervals.
        # gene_bed is the current gene
        gene_bed = pybedtools.BedTool(str(gene), from_string = True)
        # find genes that overlap with current gene
        overlaps = genes.intersect(b=gene_bed, wa=True, wb=True, s=True)
        
        # Handling ovelaps (intersect leads to duplicate columns. need to drop those columns)
        columns = ['o_seqid', 'o_start', 'o_end', 'o_id', 'o_score', 'o_strand', 
                   'c_seqid', 'c_start', 'c_end', 'c_id', 'c_score', 'c_strand']
        
        overlap_df = pd.read_table(overlaps.fn, header = None, names = columns)
        changed_overlap_df = change_format(overlap_df)
        
        # remove the current gene from the overlaps because this is not technically overlap
        overlaps_filtered = changed_overlap_df[changed_overlap_df['id'] != gene.name]
        
        if not overlaps_filtered.empty:
            overlaps_filtered_bed = pybedtools.BedTool.from_dataframe(overlaps_filtered)
            # merge all gene regions that overlap with a current gene
            # you need to sort before merging
            # also by default, merging just gives chrom, start and end. you need to add id, score, direction with collapse, distinct, distinct
            merged_overlaps = overlaps_filtered_bed.sort().merge(s=True, c='4,5,6', o='collapse,distinct,distinct')
            
            # Subtraction with strand consideration after merging: this will give non_overlapping regions
            # gene_bed is the current gene
            non_overlapping = gene_bed.subtract(merged_overlaps, s=True)

            # Debug outputs
            print(f"Current gene: {gene}")
            print(f"Merged overlaps: {merged_overlaps}")
            print(f"Non-overlapping: {non_overlapping}")
            
            # non_overlapping region of a current gene
            # only use the relevant region. the rest is the same as current gene
            # .chrom, .name are default names by pybedtool
            for non_overlapping_region in non_overlapping:
                results.append({
                    'seqid': non_overlapping_region.chrom,
                    'start': non_overlapping_region.start,
                    'end': non_overlapping_region.end,
                    'id': gene.name,
                    'score': gene.score,
                    'strand': gene.strand
                    })
                        
        else:
            # if there is no overlaping, the entire current gene is non-overlapping region
            results.append({
                'seqid': gene.chrom,
                'start': gene.start,
                'end': gene.end,
                'id': gene.name,
                'score': gene.score,
                'strand': gene.strand
                })
        
        shutil.rmtree(temp_dir)

    print('successfully finished getting non-overlapping regions')
    return pd.DataFrame(results)

def process_genes(df, chrom, direction, temp_dir):
    copy_df = df[['seqid', 'start', 'end', 'id', 'score', 'strand']]
    subset_df = copy_df[(copy_df['seqid'] == chrom) & (copy_df['strand'] == direction)]
    
    final_df = non_overlapping_region(subset_df, temp_dir)
    return final_df

def get_args():        
    parser = argparse.ArgumentParser(description="get regions of genes that are not overlapping with other genes")

    parser.add_argument('--genes_dir', dest = 'genes_dir',
                        required = True,
                        help = 'genes.bed directory')

    parser.add_argument('--out', dest = 'out',
                        required = True,
                        help = 'output name')
    
    parser.add_argument('--chrom', dest = 'chrom',
                        required = True,
                        help = 'chromosome')

    parser.add_argument('--direction', dest = 'direction',
                        required = True,
                        help = 'strand')

    parser.add_argument('--species', dest = 'species',
                        required = True,
                        help = 'species')
    
    args = parser.parse_args()
    
    genes_dir = args.genes_dir
    out = args.out

    number = args.chrom
    species = args.species
    if species == 'worm':
        chrom = number
    else:    
        chrom = 'chr' + number
    
    direction = args.direction
    
    genes = pd.read_csv(genes_dir, delimiter = '\t', names = ['seqid', 'start', 'end', 'id', 'score', 'strand', 'gene_type'])
    print('genes: ' + str(genes))
    
    return genes, out, chrom, direction, number

def run_process():

    genes, out, chrom, direction, number = get_args()
    print('successfully got arguments')

    # Determine the base temporary directory from SLURM_TMPDIR or fallback to /scratch or /tmp
    base_temp_dir = os.getenv('SLURM_TMPDIR', '/scratch')
    
    # Ensure the fallback directory exists and is writable
    if not os.path.exists(base_temp_dir) or not os.access(base_temp_dir, os.W_OK):
        base_temp_dir = os.path.join('/tmp', os.getenv('USER', 'default_user'))
        os.makedirs(base_temp_dir, exist_ok=True)
        
    final_df = process_genes(genes, chrom, direction, base_temp_dir)
    print('successfully got all non overlapping regions of genes')
    
    print(final_df)
    
    out_name = out + '_' + number + '_' + direction + '.bed'
    write_result(final_df, out_name)
    print('successfully saved the output')
    
    pybedtools.cleanup(remove_all = True)
    
if __name__ == "__main__":
    run_process()
    print("success")