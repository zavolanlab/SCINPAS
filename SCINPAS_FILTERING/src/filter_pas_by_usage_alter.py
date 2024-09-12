# -*- coding: utf-8 -*-
"""
Created on Mon Jul 22 13:44:47 2024

@author: Youngbin Moon (y.moon@unibas.ch)
"""











import argparse
import pandas as pd
from multiprocessing import Pool
import itertools

def write_to_bed(final_df, out_file):
    final_df.to_csv(out_file, sep = '\t', header = True, index = False)
    
def filter_by_usage(group, threshold, organ):
    # Sort each group by 'usage' from highest to lowest
    group = group.sort_values( by = ['usage', organ, 'id'], ascending = [False, False, True]).reset_index(drop = True)
    
    # Compute cumulative motif_presence and cumulative number of PAS within each group
    group['cumul_motif'] = group['all_motif'].cumsum()
    group['cumul_pas'] = group.index + 1
    
    # Compute motif percentage
    group['motif_pct'] = group['cumul_motif']/group['cumul_pas'] * 100
    print('organ: ' + str(organ))
    print('group: ' + str(group[['id', 'num_pas', 'quantile', 'cumul_motif', 'cumul_pas', 'motif_pct']]))
    # Find the last index where motif_pct exceeds the threshold        
    exceeding_indices = group[group['motif_pct'] >= threshold].index
    if not exceeding_indices.empty:
        last_exceeding_idx = exceeding_indices[-1]
        print('last exceeding idx: ' + str(last_exceeding_idx))
        # Filter to include rows up to and including the last exceeding index
        filtered_group = group.iloc[:last_exceeding_idx+1]
        print('filtered_group: ' + str(filtered_group[['id', 'num_pas', 'quantile', 'cumul_motif', 'cumul_pas', 'motif_pct']]))
    
    else:
        filtered_group = pd.DataFrame()
        print('empty df')

    return filtered_group

def process_organ(args):
    organ, df, threshold = args
    filtered_dict = {}
    results = []    
    # group_keys = False is used to control whether the keys from the groupby operation are included in the index of the result. 
    # Setting group_keys=False ensures that the resulting DataFrame does not have the grouping columns as part of its index, 
    # which can make subsequent operations simpler.
    grouped = df.groupby(['num_pas', 'quantile'], group_keys = False)
    
    for name, group in grouped:
        filtered_group = filter_by_usage(group, threshold, organ)
        print(f"Data for {organ}, group {name} after computing motif_pct:")
        print(filtered_group)
        print("\n")
        
        if not filtered_group.empty:
            # need to do this because you constantly update filtered_dict
            if organ not in filtered_dict:
                filtered_dict[organ] = pd.DataFrame()
            filtered_dict[organ] = pd.concat([filtered_dict[organ], filtered_group])
            
            # Gather results for the last PAS that satisfied the criteria (access the last row using iloc)
            last_good_pas = filtered_group.iloc[-1]
            results.append({'organ': organ, 'quantile': last_good_pas['quantile'], 'pas_usage': last_good_pas['usage'], 'num_good_pas': len(filtered_group)})    
    
    return filtered_dict, results

def process_polyA_sites(organ_dfs, threshold = 70):

    
    tasks = [(organ, df, threshold) for organ, df in organ_dfs.items()]
    
    with Pool() as pool:
        all_results = pool.map(process_organ, tasks)

    combined_filtered_dict = {}
    combined_results = []
    for filtered_dict, results in all_results:
        # if same key exist, update can rewrite the value. but in this case, since you do it for each different organ it is safe
        combined_filtered_dict.update(filtered_dict)
        combined_results.extend(results)
        
    # Combine all filtered PAS from different organs
    combined_filtered_df = pd.concat(combined_filtered_dict.values()).drop_duplicates(subset = ['id'])
    results_df = pd.DataFrame(combined_results)
    
    return results_df, combined_filtered_df

def load_organ_pas(beds_list):
    organ_dfs = {}
    for file in beds_list:
        organ = file.split('/')[-1].split('_')[0]
        print(file.split('/')[-1].split('_')[0])
        organ_df = pd.read_csv(file, delimiter = '\t', header = 0)
        organ_dfs[organ] = organ_df
        
    return organ_dfs

def get_args():        
    parser = argparse.ArgumentParser(description="filter PAS by the usage per organ")
           
    parser.add_argument('in_bed',
      nargs='*', help='a list of organ specific pas with pas usage')

    parser.add_argument('--species', dest = 'species',
                        required = True,
                        help = 'species for example human, mouse etc')
    
    parser.add_argument('--out', dest = 'out',
                        required = True,
                        help = 'output name')

    parser.add_argument('--thres', dest = 'thres',
                        required = True,
                        help = 'threshold')
    
    args = parser.parse_args()
    
    bed_in = args.in_bed
    species = args.species
    out = args.out  
    threshold = int(args.thres)
        
    first_elem = bed_in.pop(0)
    assert(first_elem == 'in_bed')
    
    organ_dfs = load_organ_pas(bed_in)
    out_name = str(threshold) + '_' + species + '_' + out + '.bed'
        
    organ_columns = {
        'human': ['nose', 'trachea', 'heart', 'intestine', 'breast', 'bone', 'pancreas', 'eye', 'kidney', 'penis', 'ureter', 'lung', 'liver', 'skin', 'prostate', 'uterus', 'bloodImmune', 'brain'],
        'mouse': ['Bladder', 'Tongue', 'unknown', 'Kidney', 'Spleen', 'Fat', 'Marrow', 'Lung', 'Aorta', 'Heart', 'LimbMuscle', 'MammaryGland', 'Liver', 'Skin', 'Pancreas', 'Thymus', 'LargeIntestine', 'Trachea'],
        'worm': ['EmbryoVC2010', 'EmbryoN2']
    }
    
    selected_organ_cols = organ_columns[species]
    
    assert(len(organ_dfs) == len(selected_organ_cols))
    return organ_dfs, out_name, selected_organ_cols, threshold

def run_process_alter():
    organ_dfs, out_name, selected_organ_cols, threshold = get_args()
    print('successfully got inputs')
        
    results_df, combined_filtered_df = process_polyA_sites(organ_dfs, threshold)
    print('successfully got filtered pas')
    
    combined_filtered_df.drop(columns=['num_pas', 'total_score', 'rank', 'quantile', 'usage', 'cumul_motif', 'cumul_pas', 'motif_pct'], inplace=True)
    
    write_to_bed(combined_filtered_df, out_name)
    print('successfully saved the output')
    
    out_name2 = str(threshold) + '_filtering_summary_results.bed'
    write_to_bed(results_df, out_name2)
    print('successfully saved the output')    
        
if __name__ == "__main__":
    run_process_alter()
    print('success')