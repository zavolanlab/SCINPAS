# -*- coding: utf-8 -*-
"""

Created on Sun Nov 26 22:55:32 2023

@author: Youngbin Moon (y.moon@unibas.ch)
"""











from fastqcparser import FastQCParser
import statistics
import argparse
import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

def write_BED(out_path, df):
    df.to_csv(out_path, header=True, index=False)

def draw_swarm(df, out):
    
    plt.xticks(fontsize='x-large')
    plt.yticks(fontsize='x-large')
        
    plt.rcParams['font.family'] = "DejaVu Sans"    
    
    sns.swarmplot(data=df, x="sequencing quality", y="organ", hue="organ")
    
    plt.savefig(out, bbox_inches='tight')
    
def get_median_of_median(sample_fastqc_list):
    median_of_median_list = []
    n = 0
    for sample_dir in sample_fastqc_list:
        n += 1
        sample_fastqc = FastQCParser(sample_dir)
        per_base_dict = sample_fastqc.modules['Per base sequence quality']

        # print(per_base_dict)
        # print(per_base_dict['data'])
        data = per_base_dict['data']
        sample_median_values = [elem[2] for elem in data]
        # print(sample_median_values)
        sample_median_of_median = statistics.median(sample_median_values)
        # print(sample_median_of_median)

        file_name = os.path.basename(sample_dir)
        organ = file_name.split('_')[2].split('-')[1]
        scinpas_organ = file_name
        
        print(organ)
        
        median_of_median_list.append((organ, sample_median_of_median, scinpas_organ))
        
    final_df = pd.DataFrame(median_of_median_list, columns = ['organ', 'sequencing quality', 'scinpas_organ'])
    print(final_df)
    return final_df

def get_args():
    parser = argparse.ArgumentParser(description="draw swarmplot and make a csv output containing median of median sequence quality in each sample.")

    parser.add_argument('in_fastqc',
      nargs='*', help='input fastqc zip files (multiple).')
    
    parser.add_argument('--swarm_out', dest='swarm_out', 
      help='swarm plot output file')
    
    parser.add_argument('--csv_out', dest='csv_file', 
      help='csv output file')
    
    args = parser.parse_args()
    
    fastqc_list = args.in_fastqc
    swarm_out = args.swarm_out
    csv_out = args.csv_out
        
    return fastqc_list, swarm_out, csv_out
    
def run_process():    
    
    fastqc_list, swarm_out, csv_out = get_args()
    print('successfully got input data')

    first_elem = fastqc_list.pop(0)
    assert(first_elem == 'in_fastqc')

    final_df = get_median_of_median(fastqc_list)
    print('successfully got median-organ final_df')
    
    draw_swarm(final_df, swarm_out)
    print('successfully got swarm plot')
    
    write_BED(csv_out, final_df)
    print('successfully saved csv')
    
if __name__ == '__main__':
    run_process()
    print('success')
