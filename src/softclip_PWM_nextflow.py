# -*- coding: utf-8 -*-
"""

Created on Mon Jan 10 10:26:42 2022

@author: Youngbin Moon (y.moon@unibas.ch)
"""
import matplotlib.pyplot as plt
import logomaker as lm
import argparse
import numpy as np
import pysam
import pandas as pd




"""
Aim : plot a nucleotide logo of PWM of a fixed length of softclipped region.
"""
def plot_logo(dataframe, file_template):
    """
    Parameters
    ----------
    new_df : dataframe
        a normalized PWM dataframe where
        row : position from = 0 ~ length_softclip - 1
        col : A/T/G/C
        
        new_df.loc[i, j] means normalized frequency of nucleotide j happening at position i.
        summing new_df.loc[i, j] over all nucleotides j at a given position returns 1. 

    file_template : string
        output file template
        
    Returns
    -------        
    returns nothing but plots a nucleotide logo of PWM of a fixed length of softclipped region.
    """       
    plt.figure()
    
    lm.Logo(df = dataframe)
    plt.xlabel('position')
    plt.ylabel('frequency')
    
    plt.rcParams['font.family'] = "DejaVu Sans"
    
    file1 = file_template + ".png"
    file2 = file_template + ".svg"
    
    plt.savefig(file1, bbox_inches='tight')
    plt.savefig(file2, bbox_inches='tight')

def normalize_PWM(df):
    """
    Parameters
    ----------
    df : dataframe
        PWM dataframe where
        row : position from = 0 ~ length_softclip - 1
        col : A/T/G/C
        
        df.loc[i, j] means frequency of nucleotide j happening at position i.
        
    Returns
    -------        
    new_df : dataframe
        a normalized PWM dataframe where
        row : position from = 0 ~ length_softclip - 1
        col : A/T/G/C
        
        new_df.loc[i, j] means normalized frequency of nucleotide j happening at position i.
        summing new_df.loc[i, j] over all nucleotides j at a given position returns 1.        
    """         
    new_df = df.apply(lambda x: x/df.sum(axis = 1))
    return new_df
    
def update_frequency_df(sequence):
    """
    Parameters
    ----------
    sequence : string
        sequence of the softclipped region
       
    Returns
    -------        
    temp_df : dataframe
        a partial dataframe where
        row : position from = 0 ~ length_softclip - 1
        col : A/T/G/C
        
        temp_df.loc[i,j] = 1 if nucleotide j exist in position i.
    """           
    col_names = ["A", "T", "G", "C"]
    row_names = list(range(0, len(sequence)))    
    temp_df = pd.DataFrame(np.zeros((len(row_names), len(col_names))))
    temp_df.columns = col_names
    temp_df.index = row_names
    
    for i in range(len(sequence)):
        nucleotide = sequence[i]
        copy_df = temp_df.copy()
        
        if nucleotide == 'A':
            temp_df.loc[i, "A"] += copy_df.loc[i, "A"] + 1
            
        elif nucleotide == 'T':
            temp_df.loc[i, "T"] += copy_df.loc[i, "T"] + 1
            
        elif nucleotide == 'G':
            temp_df.loc[i, "G"] += copy_df.loc[i, "G"] + 1
            
        elif nucleotide =='C':
            temp_df.loc[i, "C"] += copy_df.loc[i, "C"] + 1
            
        else:
            print('not ATGC in making PWM')
    
    return temp_df

def compute_pwm(sam, length_softclip, use_fc):
    """
    Parameters
    ----------
    sam : bam file
        a deduplicated bam file
    
    length_softclip : int
        length of a softclipped region
    
    use_fc : bool
        whether use fixed cleavage site or not
        
    Returns
    -------        
    total_frequency_df : dataframe
        PWM dataframe where
        row : position from = 0 ~ length_softclip - 1
        col : A/T/G/C
        
        total_frequency_df.loc[i, j] means frequency of nucleotide j happening at position i.
    """           
    col_names = ["A", "T", "G", "C"]
    row_names = list(range(0, length_softclip))    
    total_frequency_df = pd.DataFrame(np.zeros((len(row_names), len(col_names))))
    total_frequency_df.columns = col_names
    total_frequency_df.index = row_names
    
    for read in sam.fetch():
        tuples = read.cigartuples
        left_end = tuples[0]
        right_end = tuples[-1]
        rev = read.is_reverse
        
        full_sequence = read.get_forward_sequence()
        
        # if read is mapped to negative strand, polyA tail should be on the left end
        if rev == True and left_end[0] == 4:
            if not use_fc:
                if left_end[1] == length_softclip:
                    potential_polyA = full_sequence[len(full_sequence) - left_end[1] : len(full_sequence)]
                
                else:
                    continue
            
            elif use_fc:
                OCS = read.get_tag('XO')
                FCS = read.get_tag('XF')
                difference = OCS - FCS
                if left_end[1] - difference == length_softclip:
                    potential_polyA = full_sequence[len(full_sequence) - left_end[1] + difference : len(full_sequence)]
                
                else:
                    continue
                
        # if right_end[0] == 4, it means soft clipped on the right side of a read.
        # right_end[1] gives how many bases are soft clipped on the right side.
        elif rev == False and right_end[0] == 4:
            if not use_fc:
                if right_end[1] == length_softclip:
                    potential_polyA = full_sequence[len(full_sequence) - right_end[1] : len(full_sequence)]
                
                else:
                    continue
            
            elif use_fc:
                OCS = read.get_tag('XO')
                FCS = read.get_tag('XF')
                difference = FCS - OCS
                if right_end[1] - difference == length_softclip:
                    potential_polyA = full_sequence[len(full_sequence) - right_end[1] + difference : len(full_sequence)]
                
                else:
                    continue
        else:
            continue
        
        temp_result_df = update_frequency_df(potential_polyA)
        total_frequency_df += temp_result_df
    
    return total_frequency_df
    
def get_inputs():
    parser = argparse.ArgumentParser(description="compute PWM of softclipped region of length = 5, 6, 7, 8, 9 or 10 separately and then draw signature plots" )
    
    parser.add_argument('--bam_file', dest = 'bam_file',
                        required = True,
                        help = 'deduplicated bam file')    
    
    parser.add_argument('--pwm_signature', dest = 'pwm_signature',
                        required = True,
                        help = 'output figure name')
    
    parser.add_argument('--length_softclipped', type = int, dest = 'length_softclipped',
                        required = True,
                        help = 'length of softclipped reads to make PWM') 

    parser.add_argument('--use_FC', dest = 'use_FC',
                        required = True,
                        help = 'use fixed cleavage site or not')      
    args = parser.parse_args()

    bamFile = args.bam_file  
    bam = pysam.AlignmentFile(bamFile, "rb")     
    
    pwm_signature = args.pwm_signature
    length_softclipped = args.length_softclipped
    
    use_FC = args.use_FC
    
    return bam, pwm_signature, length_softclipped, use_FC

def run_process():
    bam, pwm_signature, length_softclipped, use_FC = get_inputs()
    print("successfully got inputs")
    
    total_frequency_df = compute_pwm(bam, length_softclipped, use_FC)
    print("successfully got PWM")
    
    normalized_df = normalize_PWM(total_frequency_df)
    print("successfully normalized the PWM")
    
    out_name = pwm_signature + str(length_softclipped)
    
    plot_logo(normalized_df, out_name)
    print("successfully plot the PWM logo")
    
if __name__ == "__main__":    
    run_process()
    print("success")
    