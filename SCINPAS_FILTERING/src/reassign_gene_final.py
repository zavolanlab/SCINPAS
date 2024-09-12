# -*- coding: utf-8 -*-
"""
Created on Mon Jul 22 13:44:47 2024

@author: Youngbin Moon (y.moon@unibas.ch)
"""











import argparse
import pandas as pd
from multiprocessing import Pool
from scipy.stats import kendalltau

def write_to_bed(final_df, out_file):
    final_df.to_csv(out_file, sep = '\t', header = True, index = False)

def find_overlapping_gene_random(pas_row, non_overlapping_sites, genes, organ_columns, which_method, species):
    overlaps = genes[(genes['seqid'] == pas_row['seqid']) & (genes['start'] <= pas_row['end'] + 1) & (genes['end'] >= pas_row['start'] - 1) & (genes['strand'] == pas_row['strand'])]
    assert(len(overlaps) > 1)
    # Randomly select one overlapping gene
    random_gene = overlaps.sample(n = 1).iloc[0]['id']
    return (pas_row['id'], random_gene)

def get_pas_info(pas_row):
    chrom = pas_row['seqid']
    rcs = int(pas_row['id'].split(':')[1])
    direction = pas_row['strand']
    
    return chrom, rcs, direction

def compute_distance(pas_row, te_row):
    chrom, rcs, direction = get_pas_info(pas_row)
    if direction == '+':
        distance = abs(rcs - te_row['end'])
    
    elif direction == '-':
        distance = abs(rcs - te_row['start'])
    
    return distance

def find_overlapping_closest(pas_row, te, genes, organ_columns, which_method, species):
    overlaps = genes[(genes['seqid'] == pas_row['seqid']) & (genes['start'] <= pas_row['end'] + 1) & (genes['end'] >= pas_row['start'] - 1) & (genes['strand'] == pas_row['strand'])]
    assert(len(overlaps) > 1)
    
    if species == 'human':
        long_noncoding_rna = 'lncRNA'
        
    elif species == 'mouse':
        long_noncoding_rna = 'lincRNA'
        
    elif species == 'worm':
        long_noncoding_rna = 'lincRNA'
    
    # Prioritize by type
    overlaps.loc[:, 'priority'] = overlaps['gene_type'].apply(lambda x: 1 if x == 'protein_coding' else 2 if x == long_noncoding_rna else 3)
    print('overlaps: ' + str(overlaps['priority']))
    # Check for the highest priority group with exactly one overlap
    highest_priority = overlaps['priority'].min()
    # reset index true in order to this properly highest_priority_overlaps.iloc[largest_span_index]['id']
    highest_priority_overlaps = overlaps[overlaps['priority'] == highest_priority].reset_index(drop=True)
    # If there is exactly one highest priority overlap and highest priority is protein_coding or lncRNA, return the gene name
    if len(highest_priority_overlaps) == 1 and highest_priority < 3:
        # df.iloc[0] accesses the entire first row of the DataFrame.
        # df.iloc[0, 'gene'] is incorrect because iloc expects integer positional indexing for both row and column.
        print('only 1 highest_priority: ' + str(highest_priority_overlaps.iloc[0]['id']))
        return (pas_row['id'], highest_priority_overlaps.iloc[0]['id'])
    
    highest_priority_genes = list(highest_priority_overlaps['id'])
    overlapping_te = te[(te['id'].isin(highest_priority_genes)) & (te['seqid'] == pas_row['seqid']) & (te['strand'] == pas_row['strand'])]
    print('overlapping_te: ' + str(overlapping_te))
    assert(len(overlapping_te) > 1)
    
    distances = []
    for _, te_row in overlapping_te.iterrows():
        distance = compute_distance(pas_row, te_row)
        distances.append((distance, te_row['id']))
    
    distances_df = pd.DataFrame(distances, columns = ['distance', 'gene_id']).reset_index(drop=True)
    print('distances_df: ' + str(distances_df))
    min_idx = distances_df['distance'].idxmin()
    print('min_idx: ' + str(min_idx))
    print('selected_gene: ' + str(distances_df.iloc[min_idx]['gene_id']))
    return (pas_row['id'], distances_df.iloc[min_idx]['gene_id'])

    raise ValueError("Unable to determine the reassigned gene id. This should never happen.") 

# pas_row = pas row in overlapping region of genes
def find_overlapping_gene(pas_row, non_overlapping_sites, genes, organ_columns, which_method, species):
    overlaps = genes[(genes['seqid'] == pas_row['seqid']) & (genes['start'] <= pas_row['end'] + 1) & (genes['end'] >= pas_row['start'] - 1) & (genes['strand'] == pas_row['strand'])]
    assert(len(overlaps) > 1)

    if species == 'human':
        long_noncoding_rna = 'lncRNA'
        
    elif species == 'mouse':
        long_noncoding_rna = 'lincRNA'
        
    elif species == 'worm':
        long_noncoding_rna = 'lincRNA'
        
    # Multiple overlaps: prioritize by type and then by Kendall coefficient
    
    # Prioritize by type
    overlaps.loc[:, 'priority'] = overlaps['gene_type'].apply(lambda x: 1 if x == 'protein_coding' else 2 if x == long_noncoding_rna else 3)
    print('overlaps: ' + str(overlaps['priority']))
    # Check for the highest priority group with exactly one overlap
    highest_priority = overlaps['priority'].min()
    # reset index true in order to this properly highest_priority_overlaps.iloc[largest_span_index]['id']
    highest_priority_overlaps = overlaps[overlaps['priority'] == highest_priority].reset_index(drop=True)
    print(f'pas_row: {pas_row}\noverlaps: {overlaps}\nhighest_priority_overlaps: {highest_priority_overlaps}')
    print('overlap genes: ' + str(overlaps['id']))
    print('highest_priority_overlaps genes: ' + str(highest_priority_overlaps['id']))
    # If there is exactly one highest priority overlap and highest priority is protein_coding or lncRNA, return the gene name
    if len(highest_priority_overlaps) == 1 and highest_priority < 3:
        # df.iloc[0] accesses the entire first row of the DataFrame.
        # df.iloc[0, 'gene'] is incorrect because iloc expects integer positional indexing for both row and column.
        print('only 1 highest_priority: ' + str(highest_priority_overlaps.iloc[0]['id']))
        return (pas_row['id'], highest_priority_overlaps.iloc[0]['id'])
    
    if which_method == 'kendall':
        # If there are multiple highest priority overlaps or the highest priority is 3, compare Kendall's tau correlations
        kendall_similarities = []
        for _, gene_row in highest_priority_overlaps.iterrows():
            # current PASes that are in non_overlapping regions of gene = gene_row
            current_pas = non_overlapping_sites[non_overlapping_sites['gene_id'] == gene_row['id']]
            print('current_pas: ' + str(current_pas))
            # It could be that there is no pas in the non overlapping region for this gene
            if not current_pas.empty:
                similarities = []
                for _, current_pas_row in current_pas.iterrows():
                    # comparing organ scores of a current pas in "overlapping" region to that of a particular current PAS of "non overlapping" region of this gene
                    similarity, _ = kendalltau(current_pas_row[organ_columns], pas_row[organ_columns])
                    print('similarity: ' + str(similarity))
                    print('pd.notna: ' + str(pd.notna(similarity)))
                    
                    if pd.notna(similarity):
                        similarities.append(similarity)
                print('similarities: ' + str(similarities))
                if similarities:
                    median_similarity = pd.Series(similarities).median()
                    print('median_similarity: ' + str(median_similarity))
                    kendall_similarities.append((median_similarity, gene_row['id']))
            
        # Check if there are valid Kendall's tau correlations
        if kendall_similarities:
            kendall_similarities = pd.DataFrame(kendall_similarities, columns = ['similarity', 'gene_id']).reset_index(drop=True)
            # should only consider positive correlation
            kendall_similarities = kendall_similarities[kendall_similarities['similarity'] > 0].reset_index(drop=True)
            if not kendall_similarities.empty:
                max_index = kendall_similarities['similarity'].idxmax()
                print('kendall_similarities: ' + str(kendall_similarities))
                print('max_index: ' + str(max_index))
                # If you are sure max_index is an integer index, both methods will yield the same result.
                # The idxmax() function in pandas returns the index label of the first occurrence of the maximum value in the Series. 
                # If the index of the Series is the default integer index, idxmax() will return an integer. 
                # However, if the index is a custom label, it will return that label.
                # In the context of the function, since kendall_similarities is a DataFrame created from a list, 
                # it will have the default integer index. Therefore, idxmax() will return an integer position.
                print('chosen geneid based on kendall: ' + str(kendall_similarities.iloc[max_index]['gene_id']))
                return (pas_row['id'], kendall_similarities.iloc[max_index]['gene_id']) # or kendall_similarities.loc[max_index, 'gene']
            
    # Fallback if no valid Kendall's tau correlations can be computed
    print('reached here: ..... still ambigous')
    # case1: multiple genes with highest_priority and none of genes have pas at non_overlapping regions
    # case2: all genes had highest_priority = 3 and none of genes have pas at non_overlapping regions
    # case3: when used method == 'longest' (only choose by highest priority and longest gene)
    # In such case choose the longest gene
    highest_priority_overlaps['span'] = highest_priority_overlaps['end'] - highest_priority_overlaps['start']
    largest_span_index = highest_priority_overlaps['span'].idxmax()
    print('highest_priority_overlaps span: ' + str(highest_priority_overlaps['span']))
    print('largest_span_index: ' + str(largest_span_index))
    print('highest_priority_overlaps: ' + str(highest_priority_overlaps))
    print('highest_priority_overlaps id: ' + str(highest_priority_overlaps['id']))
    
    print('chosen gene id based on span: ' + str(highest_priority_overlaps.iloc[largest_span_index]['id']))
    return (pas_row['id'], highest_priority_overlaps.iloc[largest_span_index]['id'])
    raise ValueError("Unable to determine the reassigned gene id. This should never happen.") 

def worker(args):
    pas_row, non_overlapping_pas, genes_bed, organ_cols, method, te_bed, species = args
    if method == 'random':
        return find_overlapping_gene_random(pas_row, non_overlapping_pas, genes_bed, organ_cols, method, species)
    
    elif method == 'closest':
        return find_overlapping_closest(pas_row, te_bed, genes_bed, organ_cols, method, species)
    
    else:
        return find_overlapping_gene(pas_row, non_overlapping_pas, genes_bed, organ_cols, method, species)
    
def find_all_overlapping_regions(df, genes_bed, organ_cols, method, te_bed, species):
    non_overlapping_pas = df[df['overlap'] == 'non_overlapping']
    overlapping_pas = df[df['overlap'] == 'overlapping']
    print('non_overlapping_pas: ' + str(non_overlapping_pas))
    print('non_overlapping_pas geneid: ' + str(non_overlapping_pas['gene_id']))
    
    overlapping_copy = overlapping_pas.copy()
    non_overlapping_copy = non_overlapping_pas.copy()
    
    if not overlapping_pas.empty:
        tasks = [(row, non_overlapping_pas, genes_bed, organ_cols, method, te_bed, species) for _, row in overlapping_pas.iterrows()]
        with Pool() as pool:
            results = pool.map(worker, tasks)
        
        result_df = pd.DataFrame(results, columns = ['id', 'reassigned_g'])
        print('result_df: ' + str(result_df))
        overlapping_copy = pd.merge(left = overlapping_copy, right = result_df, on = 'id', how = 'left')
    else:
        print("No overlapping regions to process.")
    non_overlapping_copy['reassigned_g'] = non_overlapping_pas['gene_id']
    
    reassigned_df = pd.concat([overlapping_copy, non_overlapping_copy])
    return reassigned_df
    
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
    parser = argparse.ArgumentParser(description="assign PAS gene ID according to 2 criterion")
           
    parser.add_argument('--pas', dest = 'pas',
                        required = True,
                        help = 'total pas bed')

    parser.add_argument('--original', dest = 'original',
                        required = True,
                        help = 'original full pas bed')
    
    parser.add_argument('--genes', dest = 'genes',
                        required = True,
                        help = 'filtered_genes.bed')
    
    parser.add_argument('--out', dest = 'out',
                        required = True,
                        help = 'output name')
    
    parser.add_argument('--chrom', dest = 'chrom',
                        required = True,
                        help = 'chromosome')

    parser.add_argument('--strand', dest = 'strand',
                        required = True,
                        help = 'direction of PAS')

    parser.add_argument('--species', dest = 'species',
                        required = True,
                        help = 'species for example human, mouse etc')

    parser.add_argument('--method', dest = 'method',
                        required = True,
                        help = 'which method to use')

    parser.add_argument('--te', dest = 'te',
                        required = True,
                        help = 'terminal_exons.bed')
    
    args = parser.parse_args()
    
    pas_dir = args.pas
    original_dir = args.original
    genes_dir = args.genes
    out = args.out
    chrom = args.chrom
    strand = args.strand
    te_dir = args.te
    
    species = args.species
    if species == 'worm':
        chromosome = chrom
    else:    
        chromosome = 'chr' + chrom
        
    pas = pd.read_csv(pas_dir, delimiter = '\t', header = 0)
    original = pd.read_csv(original_dir, delimiter = '\t', header = 0)
    genes = pd.read_csv(genes_dir, delimiter = '\t', names = ['seqid', 'start', 'end', 'id', 'score', 'strand', 'gene_type'])
    te = pd.read_csv(te_dir, delimiter = '\t', names = ['seqid', 'start', 'end', 'id', 'score', 'strand'])
    
    subset_pas = pas[(pas['seqid'] == chromosome) & (pas['strand'] == strand)]
    subset_original = original[(original['seqid'] == chromosome) & (original['strand'] == strand)]
    subset_genes = genes[(genes['seqid'] == chromosome) & (genes['strand'] == strand)]
    subset_te = te[(te['seqid'] == chromosome) & (te['strand'] == strand)]
    
    method = args.method
    
    out_name = method + '_' + out + '_' + chrom + '_' + strand + '.bed'
    
    organ_columns = {
        'human': ['nose', 'trachea', 'heart', 'intestine', 'breast', 'bone', 'pancreas', 'eye', 'kidney', 'penis', 'ureter', 'lung', 'liver', 'skin', 'prostate', 'uterus', 'bloodImmune', 'brain'],
        'mouse': ['Bladder', 'Tongue', 'unknown', 'Kidney', 'Spleen', 'Fat', 'Marrow', 'Lung', 'Aorta', 'Heart', 'LimbMuscle', 'MammaryGland', 'Liver', 'Skin', 'Pancreas', 'Thymus', 'LargeIntestine', 'Trachea'],
        'worm': ['EmbryoVC2010', 'EmbryoN2']
    }
    
    selected_organ_cols = organ_columns[species]
    
    return subset_pas, subset_original, subset_genes, out_name, selected_organ_cols, method, subset_te, species
    
def run_process_alter():
    subset_pas, subset_original, subset_genes, out_name, selected_organ_cols, method, subset_te, species = get_args()
    print('successfully got inputs')
    
    left_joined_pas = left_join(subset_pas, subset_original)
    print('successfully got old gene id')
    
    geneid_reassigned_genic_pas = find_all_overlapping_regions(left_joined_pas, subset_genes, selected_organ_cols, method, subset_te, species)
    print('successfully reassigned gene id')
    
    write_to_bed(geneid_reassigned_genic_pas, out_name)
    print('successfully saved the output')
    
if __name__ == "__main__":
    run_process_alter()
    print('success')