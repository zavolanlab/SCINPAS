# -*- coding: utf-8 -*-
"""
Created on Wed Feb  7 13:54:32 2024

@author: Youngbin Moon (y.moon@unibas.ch)
"""











import argparse
import pandas as pd

def write_to_bed(final_df, out_file):
    
    final_df.to_csv(out_file, sep = '\t', header = False, index = False)

def convert_col_value(name):
    if name == 1:
        value = 'true_intergenic'
    
    elif name == 2:
        value = 'antisense_TE'
    
    elif name == 3:
        value = 'antisense_intronic'
        
    elif name == 4:
        value = 'antisense_exonic'
    
    elif name == 5:
        value = 'TE'
    
    elif name == 6:
        value = 'intronic'
    
    elif name == 7:
        value = 'exonic'
    
    return value

def convert_back_df(df):
    class_col = df['class']
    converted_class_col = [convert_col_value(elem) for elem in class_col]
    
    df['class'] = converted_class_col
    print('final final df: ' + str(df))
    
    return df
    
def get_column_value(name):
    if name == 'true_intergenic':
        value = 1
    
    elif name == 'antisense_TE':
        value = 2
    
    elif name == 'antisense_intronic':
        value = 3
        
    elif name == 'antisense_exonic':
        value = 4
    
    elif name == 'TE':
        value = 5
    
    elif name == 'intronic':
        value = 6
    
    elif name == 'exonic':
        value = 7
    
    return value

def add_column(original_df, df_list):
    
    for file in df_list:
        name_template = file.split('_')[1]
        # print(name_template)
        if name_template == 'antisense' or name_template == 'true':
            file_name = name_template + '_' + file.split('_')[2]
        else:
            file_name = name_template 
            
        df = pd.read_csv(file, delimiter = '\t', names = ['seqid', 'start', 'end', 'id', 'score', 'strand'])
        
        subset_df = df.iloc[:, [0, 1, 2, 3, 5]]
        # print(file_name)
        
        col_value = get_column_value(file_name)
        col_values = [col_value] * len(subset_df)
        
        subset_df[file_name] = col_values
        print('df: ' + str(subset_df))
        
        # left join means it is based on the left df. 
        # entries of the merged df that both left and right df has will have the column from the right df
        original_df = pd.merge(left = original_df, right = subset_df, how = 'left', on = ['seqid', 'start', 'end', 'id', 'strand'])
    
    intermediate_df = original_df.fillna(0, inplace = False)
    print('length of original df: ' + str(len(original_df)))
    print('length of added df(before filtering): ' + str(len(intermediate_df)))    
    print('intermediate df: ' + str(intermediate_df))
    
    semi_final_df = intermediate_df.iloc[:, [0, 1, 2, 3, 4, 5]]
    class_info_df = intermediate_df.iloc[:, [6, 7, 8, 9, 10, 11, 12]].copy()
    # Since you filled NA with 0 and since you have 1 class per PAS, max value will be the class of that PAS
    class_column = list(class_info_df.max(axis = 1))
    
    final_df = semi_final_df.copy()
    final_df.loc[:, 'class'] = class_column
    
    print('final df: ' + str(final_df))
    
    return final_df
    
def get_args():        
    parser = argparse.ArgumentParser(description="Add class columns to ALL PAS")

    parser.add_argument('in_bed',
      nargs='*', help='a list of bed files of TE_PAS, E_PAS, I_PAS, True_IG_PAS, antisense_TE_PAS, antisense_E_PAS, antisense_I_PAS and all_PAS')
    
    parser.add_argument('--bed_out', dest = 'bed_out',
                        required = True,
                        help = 'bed filename that contains all samples PAS but with additional column  (class column)')               
    
    args = parser.parse_args()
    
    in_bed = args.in_bed
    first_elem = in_bed.pop(0)
    assert(first_elem == 'in_bed')

    class_pas = []
    for elem in in_bed:
        if 'modified' in elem:
            all_pas_dir = elem
        
        else:
            class_pas.append(elem)
    
    bed_out = args.bed_out
    
    # all_pas_dir = r"C:/Users/Geniu/Downloads/modifiedAllsamples_polyA_cluster_out_21_+.bed"
    
    # dir1 = r"C:/Users/Geniu/Downloads/polyA_true_intergenic_all_samples_+_21.bed"
    # dir2 = r"C:/Users/Geniu/Downloads/polyA_antisense_TE_all_samples_+_21.bed"
    # dir3 = r"C:/Users/Geniu/Downloads/polyA_antisense_intronic_all_samples_+_21.bed"
    # dir4 = r"C:/Users/Geniu/Downloads/polyA_antisense_exonic_all_samples_+_21.bed"
    # dir5 = r"C:/Users/Geniu/Downloads/polyA_TE_all_samples_+_21.bed"
    # dir6 = r"C:/Users/Geniu/Downloads/polyA_intronic_all_samples_+_21.bed"
    # dir7 = r"C:/Users/Geniu/Downloads/polyA_exonic_all_samples_+_21.bed"
    
    # class_pas = [dir1, dir2, dir3, dir4, dir5, dir6, dir7]
    # bed_out = r"C:/Users/Geniu/Downloads/AdditionalColAllsamples_polyA_cluster_out_21_+.bed"
    all_pas = pd.read_csv(all_pas_dir, delimiter = '\t', names = ['seqid', 'start', 'end', 'id', 'score', 'strand'])
    
    print(bed_out)
    return all_pas, class_pas, bed_out

def run_process_alter():
    all_pas, class_pas, bed_out = get_args()
    print('successfully got inputs')
    
    final_df = add_column(all_pas, class_pas)
    print('successfully added the column')
    
    converted_final_df = convert_back_df(final_df)
    print('successfully converted class column into class name column')
    
    write_to_bed(converted_final_df, bed_out)
    print('successfully saved the result')
    
if __name__ == "__main__":
    run_process_alter()
    print('success!')