# -*- coding: utf-8 -*-
"""
Created on Sun Jul 21 11:47:38 2024

@author: Youngbin Moon (y.moon@unibas.ch)
"""











import argparse
import pandas as pd
from multiprocessing import Pool
import itertools

def write_to_bed(final_df, out_file):
    final_df.to_csv(out_file, sep = '\t', header = True, index = False)
    
def filter_te(te_df, chrom, direction):
    subset_te = te_df[(te_df['seqid'] == chrom) & (te_df['strand'] == direction)]
    return subset_te

def get_pas_info(pas_row):
    chrom = pas_row['seqid']
    rcs = int(pas_row['id'].split(':')[1])
    direction = pas_row['strand']
    
    return chrom, rcs, direction

def find_closest_gene_id(pas_row, t_exons):
    chrom, rcs, direction = get_pas_info(pas_row)
    
    subset_te = filter_te(t_exons, chrom, direction)
    
    if direction == '+':
        # valid_exons = subset_te[subset_te['end'] <= rcs]
        valid_exons = subset_te.reset_index(drop=True)
        if not valid_exons.empty:
            distances = abs(rcs - valid_exons['end'])
            # closest_exon is the entire row 
            # valid_exons.iloc[distances.idxmin()] gives the same result (accessing the relevant row)
            closest_exon = valid_exons.loc[distances.idxmin()]
            
            print('distances:' + str(distances))
            print('idx: ' + str(distances.idxmin()))
            print('valid_exons: ' +str(valid_exons))
            print('closest_exon: ' + str(closest_exon))
            return closest_exon['id']
        
        else:
            print('no_gene')
            return 'no_gene'
    
    elif direction == '-':
        # valid_exons = subset_te[subset_te['start'] >= rcs]
        valid_exons = subset_te.reset_index(drop=True)
        if not valid_exons.empty:
            distances = abs(valid_exons['start'] - rcs)
            # closest_exon is the entire row 
            # valid_exons.iloc[distances.idxmin()] gives the same result (accessing the relevant row)
            closest_exon = valid_exons.loc[distances.idxmin()]
            print('distances:' + str(distances))
            print('idx: ' + str(distances.idxmin()))
            print('valid_exons: ' +str(valid_exons))
            print('closest_exon: ' + str(closest_exon))
            return closest_exon['id']
        
        else:
            print('no_gene')
            return 'no_gene'
        
def assign_nearest_gene(df, te):
    copy_df = df.copy()
    copy_df['reassigned_g'] = df.apply(find_closest_gene_id, t_exons = te, axis = 1)
    return copy_df

def get_all_nearest_gene(df, te_bed, n_cores):
    groupedBy_chrom = df.groupby(['seqid', 'strand'])
    # if you dont do this you need to use [1] to get the content because each group is a tuple of seqid and df
    contents = [contents_of_group for name_of_group, contents_of_group in groupedBy_chrom]
    
    with Pool(n_cores) as pool:
        result = pool.starmap(assign_nearest_gene, zip(contents, itertools.repeat(te_bed)))
    
    final_df = pd.concat(result)
    # final_df['re_class'] = 'intergenic'
    final_df['re_class'] = final_df['class']
    return final_df

def left_join(pas, original):
    print('length of pas: ' + str(len(pas)))
    print('length of original: ' + str(len(original)))
    assert(len(pas) == len(original))
    # Define the keys for merging
    keys = ['seqid', 'start', 'end', 'id', 'score', 'strand']
    merged_df = pd.merge(left = original, right = pas[['seqid', 'start', 'end', 'id', 'score', 'strand', 'overlap']], how = 'left', on = keys)
    print('merged_df: ' + str(merged_df))
    print('merged_df length:' + str(len(merged_df)))
    assert(len(pas) == len(merged_df))
    return merged_df

def get_args():        
    parser = argparse.ArgumentParser(description="assign PAS gene iD to the nearest gene")
           
    parser.add_argument('--pas', dest = 'pas',
                        required = True,
                        help = 'total pas bed')

    parser.add_argument('--original', dest = 'original',
                        required = True,
                        help = 'original full pas bed')
    
    parser.add_argument('--te', dest = 'te',
                        required = True,
                        help = 'terminal_exons.bed')
    
    parser.add_argument('--out', dest = 'out',
                        required = True,
                        help = 'output name')

    parser.add_argument('--n', dest = 'n',
                        required = True,
                        help = 'number of cores')

    args = parser.parse_args()
    
    pas_dir = args.pas
    original_dir = args.original
    te_dir = args.te
    out = args.out
    n = int(args.n)
    
    pas = pd.read_csv(pas_dir, delimiter = '\t', header = 0)
    original = pd.read_csv(original_dir, delimiter = '\t', header = 0)
    te = pd.read_csv(te_dir, delimiter = '\t', names = ['seqid', 'start', 'end', 'id', 'score', 'strand'])
    
    intergenic_pas = pas[pas['overlap'] == 'intergenic']
    intergenic_original = original[(original['class'] == 'true_intergenic') | (original['class'] == 'antisense_TE') | (original['class'] == 'antisense_intronic') | (original['class'] == 'antisense_exonic')]
    return intergenic_pas, intergenic_original, te, out, n

def run_process_alter():
    intergenic_pas, intergenic_original, te, out, n = get_args()
    print('successfully got inputs')

    left_joined_intergenic_pas = left_join(intergenic_pas, intergenic_original)
    print('successfully got old gene id')
    
    final_df = get_all_nearest_gene(left_joined_intergenic_pas, te, n)
    print('successfully got nearest gene id')
        
    write_to_bed(final_df, out)
    print('successfully saved the output')
    
if __name__ == "__main__":
    run_process_alter()
    print('success')
