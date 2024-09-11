#!/usr/bin/env python
# coding: utf-8

##########################
##### Import modules #####
##########################
import os
import time
import argparse
import numpy as np
import pandas as pd
from alive_progress import alive_bar

####################
##### argument #####
####################
parser = argparse.ArgumentParser(description='Counting mutation in read')

parser.add_argument('--sampleID')
parser.add_argument('--input_path')
parser.add_argument('--output_path', default=os.getcwd())
parser.add_argument('--R1_KnownSeq')
parser.add_argument('--R2_KnownSeq')
parser.add_argument('--R1_KnownSeq_start', default=1)
parser.add_argument('--R2_KnownSeq_start', default=1)

args = parser.parse_args()

sampleID = args.sampleID
input_path = args.input_path
output_path = args.output_path
R1_KnownSeq = args.R1_KnownSeq
R2_KnownSeq = args.R2_KnownSeq
R1_KnownSeq_start = int(args.R1_KnownSeq_start)
R2_KnownSeq_start = int(args.R2_KnownSeq_start)

##################################
##### Define Custom Function #####
##################################
###################################################################################
def time_conversion(duration_time):
    hour, remainder = divmod(duration_time, 3600)
    minute, second = divmod(remainder, 60)
    
    return f'{int(hour)}h {int(minute)}m {int(second)}s'
###################################################################################
def countMut_process(ref, compare, start_pos):
    length = len(ref)
    start = (start_pos - 1)
    end = (start + length)
    
    ref_seq = ref
    compare_seq = compare[start:end]
    
    mut_count = sum(1 for r, c in zip(ref_seq, compare_seq) if r != 'N' and r != c)
    return mut_count

def countMut(Seq_df, known_seq, start_pos):
    from alive_progress import alive_bar
       
    mutCount_dic = {}
    
    with alive_bar(Seq_df.shape[0], force_tty=True) as bar:
        for row in Seq_df.itertuples(index=True):
            i = row.Index
            seq = row.Sequence

            mut_count = countMut_process(known_seq, seq, start_pos)
            mutCount_dic[i] = mut_count

            bar()
    
    return mutCount_dic
###################################################################################

print(f'==================================================================================================================================================================')
print(f"Command :")
print(f"python 3-1_Counting_Mutation_inRead.py --prefix {sampleID} --input {input_path} --output {output_path} --R1_knownSeq {R1_KnownSeq} --R2_knownSeq {R2_KnownSeq} --R1_knownSeq_start {R1_KnownSeq_start} --R2_knownSeq_start {R2_KnownSeq_start}")
print(f'==================================================================================================================================================================')
print()

try: 
    if not os.path.exists(f'{output_path}'): 
        os.makedirs(f'{output_path}')
except OSError: 
    print("Error: Failed to create the directory.")

print(f'-----------------------------------------------------------------------------------------------------')
print(f"Start:\t{time.ctime()}"); start_time = time.time()
print(f'-----------------------------------------------------------------------------------------------------')
for read_N in [1, 2]:
    if read_N == 1: 
        polyBase = 'T'
        KnownSeq, KnownSeq_start = R1_KnownSeq, R1_KnownSeq_start
    elif read_N == 2: 
        polyBase = 'A'
        KnownSeq, KnownSeq_start = R2_KnownSeq, R2_KnownSeq_start
    print()
    print("######################################################")
    print(f"READ {read_N}")  
    print("######################################################")

    ### 1. Import data
    print("### 1. Importing CSV (Sequence / Split Sequence)")
    Seq_df = pd.read_csv(f"{input_path}/{sampleID}_Sequence.R{read_N}.txt.gz", header=None, sep="\t", dtype='category', compression='gzip')
    Seq_split_df = pd.read_csv(f"{input_path}/{sampleID}_Sequence_split.R{read_N}.txt.gz", header=None, sep="\t", dtype='category', compression='gzip')
    Seq_df.columns = ["Sequence"]
    Seq_split_df.columns = [i for i in range(1,Seq_split_df.shape[1]+1)]
    print()
    
    ### 2. Counting mutation in known sequence
    print("### 2. Counting mutation in known sequence")
    mutCount_dic = countMut(Seq_df, KnownSeq, KnownSeq_start)
    mutCount_df = pd.DataFrame.from_dict(mutCount_dic, orient='index')
    mutCount_df.columns = ["Mut_Count"]

    mutCount_df.to_csv(f'{output_path}/{sampleID}_knownSeq_mutCount_R{read_N}_df.csv', sep="\t")
    print("---------------------------- Save mutation count table ----------------------------")
    print()
    
    ### 3. Starting homopolymer as correct base
    print("### 3. Starting homopolymer as correct base")
    Homopolymer_correctStart = Seq_split_df.loc[Seq_split_df[len(KnownSeq) + 1] == polyBase]
    Homopolymer_correctStart_idx = Homopolymer_correctStart.index

    idx_file = f'{output_path}/{sampleID}_homopolymer_correctStart_R{read_N}_idx.csv'
    np.savetxt(idx_file, list(Homopolymer_correctStart_idx), fmt='%d', delimiter=',')
    print("------------- Save read index of starting homopolymer as correct base ------------- ")
    print()
    

print(f'-----------------------------------------------------------------------------------------------------')
print(f"End:\t{time.ctime()}"); end_time = time.time()
print(f'-----------------------------------------------------------------------------------------------------')
print()

duration_time = end_time - start_time
conversion_time = time_conversion(duration_time)
print(f"Time:\t{duration_time} ({conversion_time})")