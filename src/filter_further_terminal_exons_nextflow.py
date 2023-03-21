# -*- coding: utf-8 -*-
"""

Created on Fri May  6 19:05:00 2022

@author: Youngbin Moon (y.moon@unibas.ch)
"""
import pandas as pd
import argparse







"""
Aim : From terminal_exons.bed (generated from gtf file), filter further terminal exons with transcript support level <= 3. 
"""

def write_bed(final_df, out_file):
    """
    Parameters
    ----------
    final_df : dataframe
        dataframe containing terminal_exons with transcript support level <= 3.
    
    out_file : string
        Output file name. This output will be in bed format.

    Returns
    -------
    Returns nothing but writes the output in the bed format.
    """
    # chr, start, end are compulsory for bed file. others are accessory. adjust according to your needs
    final_df.to_csv(out_file, sep = '\t', header = False, index = False,
    columns = ['seqid', 'start', 'end', 'id', 'score', 'strand'])
    
def filter_further_terminal_exons(df):
    """
    Parameters
    ----------
    df : dataframe
        dataframe that only contains terminal exons from gtf.

    Returns
    -------
    result_df : dataframe
        dataframe that only contains terminal exons with transcript support level <= 3.
    """
    result_df = df[df['score'] <= 3]
    print(str(result_df))
    return result_df

def get_inputs():
    
    parser = argparse.ArgumentParser(description="filter terminal exons that have transcript support level <= 3")
    
    parser.add_argument('--terminal_exons_bed_input', dest = 'terminal_exons_bed_input',
                        required = True,
                        help = 'terminal exons.bed input generated from gtf file.')
      
    parser.add_argument('--bed_out', dest = 'bed_out',
                        required = True,
                        help = 'bed file output')  
        
    args = parser.parse_args()
        
    terminal_exons_bed_input = args.terminal_exons_bed_input
    terminal_exons_bed = pd.read_csv(terminal_exons_bed_input, delimiter = '\t', header = None)
    terminal_exons_bed.columns = ['seqid', 'start', 'end', 'id', 'score', 'strand']
    
    bed_out = args.bed_out

    return terminal_exons_bed, bed_out

def run_process():
    
    terminal_exons_bed, bed_out = get_inputs()
    print('successfully got inputs')
    
    filtered_terminal_exons = filter_further_terminal_exons(terminal_exons_bed)
    print('successfully filtered terminal exons')
    
    write_bed(filtered_terminal_exons, bed_out)
    print('successfully saved the output in the bed file')
    
if __name__ == "__main__":
    run_process()
    print('success')