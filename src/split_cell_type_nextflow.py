# -*- coding: utf-8 -*-
"""
Created on Wed Feb  1 16:18:48 2023

@author: Youngbin Moon (y.moon@unibas.ch)
"""

import pandas as pd
import argparse









"""
Aim : subset a meta dataframe that contains all CB tags of a specific sample and cell type.
"""
def split_cell_type(df, cell_type, new_sample_name):
    """
    Parameters
    ----------
    df : dataframe
        a dataframe that contains all cell annotations across all samples and cell types
    
    cell_type : string
        cell type 
        
    new_sample_name : string
        modified sample name that is used in the meta data csv (df).

    Returns
    -------
    final_df : dataframe
        a dataframe containing CB tags of a specific sample and cell type.

    """    
    # subset dataframe which has a specific cell type and a specific sample name
    subset_df = df.loc[(df['cell_type'] == cell_type) & (df['sample'] == new_sample_name)]
    # collect 'CB' column only from the dataframe
    final_df = subset_df[['CB']]
    
    return final_df
    
def change_sample_name(sample, type_sample):
    """
    Parameters
    ----------
    sample : string
        sample name
        
    type_sample : string
        species of the sample. it should be either immune_cells or spermatocytes.

    Returns
    -------
    new_sample_name : string
        modified sample name that is used in the meta data csv.
        
    """    
    if type_sample == "spermatocyte":
        number = sample.split('_')[2]
        new_sample_name = 'mouse' + str(number)
    
    elif type_sample == "immune_cells":
        new_sample_name = sample.split('_')[1] + sample.split('_')[2]
    
    return new_sample_name

def run ():    
    parser = argparse.ArgumentParser(description="make 2 text files containing CB of cell type 1 and cell type 2 respecitvely")
    parser.add_argument('--cell_type_dir', dest = 'cell_type_dir',
                        required = True,
                        help = 'directory of csv containing cell type annotation')
    
    parser.add_argument('--txt_output_template', dest = 'txt_output_template',
                        required = True,
                        help = 'output txt file template')
  
    parser.add_argument('--type1', dest = 'type1',
                        required = True,
                        help = 'first cell type that you want to split')

    parser.add_argument('--type2', dest = 'type2',
                        required = True,
                        help = 'Second cell type that you want to split')

    parser.add_argument('--sample_name', dest = 'sample_name',
                        required = True,
                        help = 'Sample name that you want to split')

    parser.add_argument('--sample_type', dest = 'sample_type',
                        required = True,
                        help = 'Sample type: It should be either immune_cells or spermatocytes')
                     
    args = parser.parse_args()
    
    cell_type_dir = args.cell_type_dir
    cell_type_csv = pd.read_csv(cell_type_dir, header = None, names = ['sample', 'CB', 'cell_type', 'count1', 'count2'])
    
    txt_output_template = args.txt_output_template
    type1 = args.type1
    type2 = args.type2
    sample_name = args.sample_name
    sample_type = args.sample_type
    
    return cell_type_csv, txt_output_template, type1, type2, sample_name, sample_type
    
def run_process():
    cell_type_csv, txt_output_template, type1, type2, sample_name, sample_type = run()
    print('successfully initialized')
    
    new_sample = change_sample_name(sample_name, sample_type)
    print('successfully changed sample type')
    
    type1_df = split_cell_type(cell_type_csv, type1, new_sample)
    print('successfully got cell type 1 specific dataframe')

    type2_df = split_cell_type(cell_type_csv, type2, new_sample)
    print('successfully got cell type 2 specific dataframe')
    
    out_name1 = txt_output_template + '_' + new_sample + '_' + type1 + '.txt'
    out_name2 = txt_output_template + '_' + new_sample + '_' + type2 + '.txt'
    
    # without header and index save it as a txt file
    type1_df.to_csv(out_name1, header = False, index = False)
    print('successfully got cell type 1 specific dataframe')
 
    type2_df.to_csv(out_name2, header = False, index = False)
    print('successfully got cell type 1 specific dataframe')
    
if __name__ == "__main__":
    run_process()
    print('success')