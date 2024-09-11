# -*- coding: utf-8 -*-
"""
Created on Fri May 31 16:58:34 2024

@author: Youngbin Moon (y.moon@unibas.ch)
"""

import numpy as np
import pandas as pd
import argparse
from multiprocessing import Pool
import itertools

def write_to_bed(df, out_file):
    df.to_csv(out_file, sep = '\t', header = True, index = False)

# until 6th columns are normal bed format: seqid, start, end, id, score, strand.
# step 1): The sum of the sample columns (from the 7th column onward) calculated in row-wise. 
# This gives a Series where each value is the sum of sample scores for that row(cleavage site) 
# -> grouped_cs_sum
 
# step 2): This Series is then grouped by 'id' (pas id) followed by sum within each pas id group. 
# This gives a Series where each value is pas score of this particular organ
# -> grouped_pas_su,

# step 3): count number of samples in this organ
# -> num_samples
 
# step 4) calcaulate average scores by dividing the total number of samples in this organ
# Total scores for each pas id are divided by the number of samples in that organ
# The result is converted into dataframe. index (pas_id) is resetted as 'id'
# -> average_scores

# When you use the groupby function in pandas, the grouping key (in this case, 'id') becomes the index of the resulting DataFrame or Series.
# However, you can reset the index to make 'id' a regular column again.
def get_avg_score_per_organ(bed_dir, organ):
    bed = pd.read_csv(bed_dir, sep='\t', header=0)
    
    # the first sum gives you total score of a given cleavage site
    grouped_cs_sum = bed.iloc[:, 6:].sum(axis = 1)
    print('sample_sum: ' + str(grouped_cs_sum))

    # The grouped_cs_sum Series is aligned with the 'id' column from the original DataFrame.
    # The groupby method groups the values in "grouped_cs_sum" according to the corresponding 'id' values from bed df (pas id).
    # The sum() method is then applied to these groups, giving the total score for each unique 'id' (pas) in this organ.
    # When performing grouped_cs_sum.groupby(bed_df['id']).sum(), pandas aligns the grouped_cs_sum Series with the 'id' column from the original DataFrame.
    # Each value in grouped_cs_sum is grouped according to the corresponding 'id' from bed_df.
    
    grouped_pas_sum = grouped_cs_sum.groupby(bed['id']).sum()
    print('group_sum: ' + str(grouped_pas_sum))
    
    # Calculate the number of sample columns from the 7th cols onwards (sample columns)
    num_samples = len(bed.columns[6:])
    print('num_samples: ' + str(num_samples))
    print('organ: ' + str(organ))
    
    # Calculate average scores by dividing by total number of samples in that organ
    # The total scores for each 'id' (from grouped_pas_sum) are divided by the total number of samples in that organ
    # The result is converted to a DataFrame, and columns are renamed to 'id' (pas id) and 'average_score'.
    average_scores = (grouped_pas_sum/num_samples).reset_index()
    # Rename columns to 'id' and the specified organ name
    average_scores.columns = ['id', organ]
    print('average_scores: ' + str(average_scores))
    return average_scores
    
def get_all_average_scores_per_organ(bed_list, num_cores, organ):
    # zip(bed_list, itertools.repeat(organ)) creates an iterator that pairs each element of bed_list with the organ value. 
    # The itertools.repeat(organ) part creates an iterator that repeats the organ value indefinitely.
    with Pool(num_cores) as pool:
        results = pool.starmap(get_avg_score_per_organ, zip(bed_list, itertools.repeat(organ)))
    
    final_df = pd.concat(results)
    
    return final_df
    
def get_args():        
    
    parser = argparse.ArgumentParser(description="compute average score of a PAS per organ")

    parser.add_argument('in_bed',
      nargs='*', help='Input modified cleavage sites bed file of a specific organ')
   
    parser.add_argument('--out_name', dest = 'out_name',
                        required = True,
                        help = 'bed file out name')

    parser.add_argument('--organ', dest = 'organ',
                        required = True,
                        help = 'organ name')
    
    parser.add_argument('--n', dest = 'n',
                        required = True,
                        help = 'number of cores')
    args = parser.parse_args()
    
    in_bed_dir = args.in_bed
    out_name = args.out_name
    organ = args.organ
    n = int(args.n)
    
    print('in_bed_dir: ' + str(in_bed_dir))
    in_bed_dir.pop(0)
    print('in_bed_dir: ' + str(in_bed_dir))
    print('length of in_bed_dir: ' + str(len(in_bed_dir)))
            
    return in_bed_dir, out_name, organ, n

def run_process():
    in_bed_dir, out_name, organ, n = get_args()
    print('successfully got inputs')
    
    final_df = get_all_average_scores_per_organ(in_bed_dir, n, organ)
    print('successfully got final df')
    
    write_to_bed(final_df, out_name)
    print('successfully saved the result')
    
if __name__ == "__main__":
    run_process()
    print("success")



