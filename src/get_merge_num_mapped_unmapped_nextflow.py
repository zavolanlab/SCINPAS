# -*- coding: utf-8 -*-
"""

Created on Sun Nov 26 22:55:32 2023

@author: Youngbin Moon (y.moon@unibas.ch)
"""











import argparse
import pandas as pd
import os
import numpy as np

def get_filtered_input(df, original_samples_df):
    subset = df[['scinpas', 'organ', 'total', 'unique']]
    subset.columns = ['sample', 'organ', 'total', 'unique']
    filtered_input = pd.merge(subset, original_samples_df, on = ['sample', 'organ'], how = 'left')
    filtered_input['percentage'] = filtered_input['unique']/filtered_input['total']
    return filtered_input[['percentage', 'dir', 'sample', 'organ']]

def filter_csv(df):
    projects = list(set((df['project'])))
    # print(projects)
    final_df = pd.DataFrame(columns= ['total', 'mapped', 'unmapped', 'unique', 'scinpas', 'organ', 'project'])
    filtered_out_df = pd.DataFrame(columns= ['total', 'mapped', 'unmapped', 'unique', 'scinpas', 'organ', 'project'])
    # print(final_df)
    
    for project in projects:
        subset_df = df[df['project'] == project]
        # print(len(subset_df))
       

        unique_mapped_ratio = np.divide(subset_df['unique'], subset_df['total'])*100
    
        print(len(unique_mapped_ratio))
        # print(unique_mapped_ratio)
        
        # calculate interquartile range 
        q3, q1 = np.percentile(sorted(unique_mapped_ratio), [75, 25])
        # quartiles = unique_mapped_ratio.quantile([0.25, 0.75])
        # iqr = quartiles[0.75] - quartiles[0.25]
        # q3 = quartiles[0.75]
        # q1 = quartiles[0.25]
        print(sorted(unique_mapped_ratio))
        print(project)
        print('q3 is: ' + str(q3))
        print('q1 is: ' + str(q1))
        
        iqr = q3 - q1
        outlier_threshold = q1 - 1*iqr
        print('iqr is: ' + str(iqr))
        print('outlier is: ' + str(outlier_threshold))
        less_df = subset_df[(unique_mapped_ratio < outlier_threshold) | (unique_mapped_ratio <= 70)]
        
        good_df = subset_df[(unique_mapped_ratio >= outlier_threshold) & (unique_mapped_ratio > 70)]
        
        final_df = pd.concat([final_df, good_df]).copy()
        filtered_out_df = pd.concat([filtered_out_df, less_df]).copy()
    
    return final_df, filtered_out_df
        
def merge_csv(in_csvs):
    
    first_csv_dir = in_csvs.pop(0)
    # print(first_csv_dir)
    first_df = pd.read_csv(first_csv_dir, header = None)
    first_df.columns = ['total', 'mapped', 'unmapped', 'unique']
    basename = os.path.basename(first_csv_dir)
    print(basename)
    sample = ('_'.join(basename.split('_')[0:3])).split('-')[0]
    print(sample)
    organ = ('_'.join(basename.split('_')[0:3])).split('-')[1]
    print(organ)
    project = sample.split('_')[1]
    
    if organ == 'brain':
        project = 'BICCN'

    
    first_df['scinpas'] = sample
    first_df['organ'] = organ
    first_df['project'] = project
    
    # print(first_df)
    total_df = first_df.copy()
    
    for file in in_csvs:
        df = pd.read_csv(file, header = None)
        df.columns = ['total', 'mapped', 'unmapped', 'unique']
        next_basename = os.path.basename(file)
        
        sample = ('_'.join(next_basename.split('_')[0:3])).split('-')[0]
        organ = ('_'.join(next_basename.split('_')[0:3])).split('-')[1]
        project = sample.split('_')[1]   
        
        if organ == 'brain':
            project = 'BICCN'  
            
        df['scinpas'] = sample
        df['organ'] = organ
        df['project'] = project  
      
        total_df = pd.concat([total_df, df]).copy()
    
    print(total_df) 
    return total_df
       
def get_args():        
    parser = argparse.ArgumentParser(description="get merged csv file of the number of mapped, unmapped and unique reads with sample info")

    parser.add_argument('csv_inputs',
      nargs='*', help='csv input files')
       
    parser.add_argument('--out_csv', dest = 'out_csv',
                        required = True,
                        help = 'csv output')    
 
    parser.add_argument('--original_input', dest = 'original_input',
                        required = True,
                        help = 'original input directory')  

    parser.add_argument('--modified_input_out', dest = 'modified_input_out',
                        required = True,
                        help = 'filtered_input out name')  
       
    args = parser.parse_args()
    
    csv_inputs = args.csv_inputs
    out_csv = args.out_csv
    original_input_dir = args.original_input
    modified_input_out = args.modified_input_out
    
    first_elem = csv_inputs.pop(0)
    assert(first_elem == 'csv_inputs')
    print('csv_inputs: ' + str(csv_inputs))
    print('length: ' + str(len(csv_inputs)))
    
    # os.chdir(r'C:/Users/Geniu/Downloads/num_mapped_folder/')
    # files_list = os.listdir(r"./")
    # csv_inputs = ["./" + elem for elem in files_list if elem.endswith(".csv")]
    # csv_inputs.insert(0, 'csv_inputs')
    
    # out_csv = r'C:/Users/Geniu/Downloads/num_mapped_unmapped_reads_all_samples_modifed.csv'    
    # assert(first_elem == 'csv_inputs')
    
    original_input = pd.read_csv(original_input_dir, header = 0)
    
    return csv_inputs, out_csv, original_input, modified_input_out

def run_process():

    csv_inputs, out_csv, original_input, modified_input_out = get_args()
    print('successfully got inputs')
    
    total_df = merge_csv(csv_inputs)
    print('successfully merged dataframe')
    
    out_csv1 = out_csv + '.csv'
    total_df.to_csv(out_csv1, header=True, index=False)
    print('successfully wrote the merged csv output')
        
    final_df, filtered_out_df = filter_csv(total_df)
    print('successfully filtered bad samples')
    
    out_csv2 = out_csv + '_good.csv'
    out_csv3 = out_csv + '_bad.csv'
    
    final_df.to_csv(out_csv2, header=True, index=False)
    filtered_out_df.to_csv(out_csv3, header=True, index=False)
    print('successfully saved bad and good samples')
    
    modified_input = get_filtered_input(final_df, original_input)
    print('successfully got modified input df')
    
    modified_input.to_csv(modified_input_out, header=True, index=False)
    print('successfully saved modified input df')
    
if __name__ == "__main__":
    run_process()
    print("success")