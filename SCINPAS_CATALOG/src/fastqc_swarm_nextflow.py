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

# columns of first_filtered_df: 'percentage', 'dir', 'sample', 'organ'
# columns of df: organ, sequencing_quality, scinpas_organ
# left join first_filtered_df with df
# fill na = 0
# change sequencing quality > 0 as 2nd_filtering = 'good' else 'bad'
# intermediate_df: [[sample, organ, percentage, dir, sequencing_quality, 2nd_filtering]]
# 1st out_df: [[sample, organ, percentage, sequencing_quality, 2nd_filtering]]
# 2nd out_df: filter with 2nd_filtering = 'good' and then  [[dir, sample, organ]]
def make_final_input(df, first_filtered_df):
    df['sample'] = df['scinpas_organ'].apply(lambda x: x.split('-')[0])
    df = df.drop('scinpas_organ', axis = 1)
    
    intermediate_df = pd.merge(left = first_filtered_df, right = df, how = 'left', on = ['sample', 'organ'])
    
    intermediate_df.fillna(0, inplace = True)
    intermediate_df['second_filtering'] = intermediate_df['sequencing_quality'].apply(lambda x: 'good' if x > 0 else 'bad')
    
    final_df1 = intermediate_df[['sample', 'organ', 'percentage', 'sequencing_quality', 'second_filtering']]
    
    final_df2 = intermediate_df[intermediate_df['second_filtering'] == 'good']
    final_df2 = final_df2[['dir', 'sample', 'organ']]
    
    return final_df1, final_df2
    
def draw_swarm(df, out):
    
    plt.xticks(fontsize='x-large')
    plt.yticks(fontsize='x-large')
        
    plt.rcParams['font.family'] = "DejaVu Sans"    
    
    ax = sns.swarmplot(data=df, x="sequencing_quality", y="organ", hue="organ")
    sns.move_legend(ax, "upper left", bbox_to_anchor=(1, 1))
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
        scinpas_organ = '_'.join(file_name.split('_')[0:3])
        
        print(organ)
        
        median_of_median_list.append((organ, sample_median_of_median, scinpas_organ))
        
    final_df = pd.DataFrame(median_of_median_list, columns = ['organ', 'sequencing_quality', 'scinpas_organ'])
    print(final_df)
    return final_df

def get_args():
    parser = argparse.ArgumentParser(description="draw swarmplot and make a csv output containing median of median sequence quality in each sample.")

    parser.add_argument('in_fastqc',
      nargs='*', help='input fastqc zip files (multiple).')
    
    parser.add_argument('--swarm_out', dest='swarm_out', 
      help='swarm plot output file')
    
    parser.add_argument('--csv_out', dest='csv_out', 
      help='csv output file')

    parser.add_argument('--threshold', dest='threshold', 
      help='threshold for median of median phred score')

    parser.add_argument('--first_filtered_dir', dest='first_filtered_dir', 
      help='first filtered samples csv file')

    parser.add_argument('--second_filtered_out', dest='second_filtered_out', 
      help='second filtered samples output csv file name')    
    args = parser.parse_args()
    
    fastqc_list = args.in_fastqc
    swarm_out = args.swarm_out
    csv_out = args.csv_out
    threshold = int(args.threshold)
    
    first_filtered_dir = args.first_filtered_dir
    second_filtered_out = args.second_filtered_out
    
    first_filt_df = pd.read_csv(first_filtered_dir, delimiter = ',', header = 0)
    
    return fastqc_list, swarm_out, csv_out, threshold, first_filt_df, second_filtered_out
    
def run_process():    
    
    fastqc_list, swarm_out, csv_out, threshold, first_filt_df, second_filtered_out = get_args()
    print('successfully got input data')

    first_elem = fastqc_list.pop(0)
    assert(first_elem == 'in_fastqc')

    seq_quality_df = get_median_of_median(fastqc_list)
    print('successfully got median-organ seq_quality_df')
    
    draw_swarm(seq_quality_df, swarm_out)
    print('successfully got swarm plot')
    
    csv_out1 = csv_out + '_full.csv'
    write_BED(csv_out1, seq_quality_df)
    print('successfully saved csv')
    
    filtered_seq_quality_df = seq_quality_df[seq_quality_df['sequencing_quality'] >= threshold]
    print('successfully filtered bad samples')
    
    final_df1, final_df2 = make_final_input(filtered_seq_quality_df, first_filt_df)
    print('successfully got 2 final dtaframes')
    
    csv_out2 = csv_out + '_filtered.csv'
    write_BED(csv_out2, final_df1)
    print('successfully saved filtered csv')    
    
    write_BED(second_filtered_out, final_df2)
    print('successfully saved final filtered input samples.csv')
    
if __name__ == '__main__':
    run_process()
    print('success')
