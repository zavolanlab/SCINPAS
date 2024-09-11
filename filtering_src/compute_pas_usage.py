# -*- coding: utf-8 -*-
"""
Created on Mon Jul 22 13:44:47 2024

@author: Youngbin Moon (y.moon@unibas.ch)
"""











import argparse
import pandas as pd
from multiprocessing import Process, Manager
# compute pas usage per each organ
def write_to_bed(final_df, out_file):
    final_df.to_csv(out_file, sep = '\t', header = True, index = False)
    
def compute_pas_usage_for_organ(df, organ, return_list):
    organ_pas = df[df[organ] > 0]
    
    organ_specific_pas = organ_pas.copy()
    # Compute the number of polyA sites for each gene and append as a new column 'num_pas'
    organ_specific_pas['num_pas'] = organ_specific_pas.groupby('reassigned_g')['id'].transform('count')   
    
    # Group by 'gene_id' and compute the sum of scores for each gene
    organ_specific_pas['total_score'] = organ_specific_pas.groupby('reassigned_g')[organ].transform('sum')
    
    # Rank the total scores and convert to quantile-like categories
    organ_specific_pas['rank'] = organ_specific_pas['total_score'].rank(method='first')
    
    # Create 20 quantiles based on 'total_score'
    quantile_labels = [f'q{i+1}' for i in range(20)]
    organ_specific_pas['quantile'] = pd.qcut(organ_specific_pas['rank'], 20, labels = quantile_labels)
    
    # Compute the ratio of each row's score to the total scores for its respective gene
    organ_specific_pas['usage'] = organ_specific_pas[organ]/organ_specific_pas['total_score']
    print(f"Computed pas usage for {organ}:")
    print('computed pas usage: ' + str(organ_specific_pas[['id', 'reassigned_g', 'num_pas', 'total_score', 'quantile', 'usage']]))
    
    return_list.append((organ, organ_specific_pas))

def compute_pas_usage(df, organs):
    with Manager() as manager:
        return_list = manager.list()
        processes = []
        
        for organ in organs:
            p = Process(target = compute_pas_usage_for_organ, args = (df, organ, return_list))
            p.start()
            processes.append(p)
        
        for p in processes:
            p.join()
        
        organ_dfs = list(return_list)
    
    return organ_dfs
        
def get_args():        
    parser = argparse.ArgumentParser(description="compute PAS usage per organ")
           
    parser.add_argument('--pas', dest = 'pas',
                        required = True,
                        help = 'total pas bed')

    parser.add_argument('--species', dest = 'species',
                        required = True,
                        help = 'species for example human, mouse etc')
    
    parser.add_argument('--out', dest = 'out',
                        required = True,
                        help = 'output name')
    
    args = parser.parse_args()
    
    pas_dir = args.pas
    species = args.species
    out = args.out  
        
    pas = pd.read_csv(pas_dir, delimiter = '\t', header = 0)
        
    organ_columns = {
        'human': ['nose', 'trachea', 'heart', 'intestine', 'breast', 'bone', 'pancreas', 'eye', 'kidney', 'penis', 'ureter', 'lung', 'liver', 'skin', 'prostate', 'uterus', 'bloodImmune', 'brain'],
        'mouse': ['Bladder', 'Tongue', 'unknown', 'Kidney', 'Spleen', 'Fat', 'Marrow', 'Lung', 'Aorta', 'Heart', 'LimbMuscle', 'MammaryGland', 'Liver', 'Skin', 'Pancreas', 'Thymus', 'LargeIntestine', 'Trachea'],
        'worm': ['EmbryoVC2010', 'EmbryoN2']
    }
    
    selected_organ_cols = organ_columns[species]
    
    return pas, out, selected_organ_cols

def run_process_alter():
    pas, out, selected_organ_cols = get_args()
    print('successfully got inputs')
    
    organ_dfs = compute_pas_usage(pas, selected_organ_cols)
    print('successfully computed organ_dfs')
    
    for elem in organ_dfs:
        organ = elem[0]
        organ_df = elem[1]
        out_name = organ + '_' + out + '.bed'
        
        write_to_bed(organ_df, out_name)
        print(f'successfully got the {organ} dataframe')
        
if __name__ == "__main__":
    run_process_alter()
    print('success')