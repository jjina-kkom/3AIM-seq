#!/usr/bin/env python
# coding: utf-8

##########################
##### Import modules #####
##########################
import os
import math
import time
import argparse
import numpy as np
import pandas as pd

####################
##### argument #####
####################
parser = argparse.ArgumentParser(description='Selection survival reads')

parser.add_argument('--sampleID')
parser.add_argument('--read_count')
parser.add_argument('--input_path_1')
parser.add_argument('--input_path_2')
parser.add_argument('--output_path', default=os.getcwd())
parser.add_argument('--R1_KnownSeq')
parser.add_argument('--R2_KnownSeq')
parser.add_argument('--NGS_ErrorRate', default=0.01)
parser.add_argument('--STdev_filter_whis', default=1.5)
parser.add_argument('--mutCount_filter_whis', default=1.5)

args = parser.parse_args()

sampleID = args.sampleID
read_count = int(args.read_count)
input_path_1 = args.input_path_1
input_path_2 = args.input_path_2
output_path = args.output_path
R1_KnownSeq = args.R1_KnownSeq
R2_KnownSeq = args.R2_KnownSeq
ErrorRate = float(args.NGS_ErrorRate)
STdev_filter_whis = float(args.STdev_filter_whis)
mutCount_filter_whis = float(args.mutCount_filter_whis)

##################################
##### Define Custom Function #####
##################################
###################################################################################
def time_conversion(duration_time):
    hour, remainder = divmod(duration_time, 3600)
    minute, second = divmod(remainder, 60)
    
    return f'{int(hour)}h {int(minute)}m {int(second)}s'
###################################################################################
def cutoff_outliers(values, whis):
    value_list = list(values)
    Q1, Q3 = np.percentile(value_list, [25,75])
    IQR = Q3 - Q1
    lower_bound = round(Q1 - (IQR * whis), 5)
    upper_bound = round(Q3 + (IQR * whis), 5) 
    
    dic = {"lower":lower_bound, 
           "upper":upper_bound}
    return dic

def cutoff_mut_ByErrorRate(ErrorRate, seqLen):
    mutCount_cutoff = math.ceil(seqLen * ErrorRate)
    return mutCount_cutoff
###################################################################################

print(f'==================================================================================================================================================================')
print(f"Command :")
print(f"python 3-2_Selection_SurvivalReads.py --prefix {sampleID} --input_1 {input_path_1} --input_2 {input_path_2} --output {output_path} --R1_knownSeq {R1_KnownSeq} --R2_knownSeq {R2_KnownSeq} --errorRate {ErrorRate} --STdev_filter_whis {STdev_filter_whis} --mutCount_filter_whis {mutCount_filter_whis}")
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
        KnownSeq_len = len(R1_KnownSeq)
    elif read_N == 2: 
        KnownSeq_len = len(R2_KnownSeq)
    print()
    
    print("######################################################")
    print(f"READ {read_N}")  
    print("######################################################")

    ### 1. Filtered by standard deviation of baseQ in read
    print("### 1. Filtered by standard deviation of baseQ in read")
    print("###### (1) Importing CSV (BaseQ standard)")
    BaseQ_stat_df = pd.read_csv(f"{input_path_1}/{sampleID}_BaseQ_stat.R{read_N}.txt", header=None, sep="\t", dtype='float16')[[0]]
    BaseQ_stat_df.columns = ['STdev']
    
    print("###### (2) Read filtering")
    Q_stdev_cutoff = cutoff_outliers(BaseQ_stat_df['STdev'], STdev_filter_whis)
    print(f'> Standard Deviation of BaseQ Cutoff (Lower):\t{Q_stdev_cutoff["lower"]}')
    print(f'> Standard Deviation of BaseQ Cutoff (Upper):\t{Q_stdev_cutoff["upper"]}')
    Survive_idx_1 = BaseQ_stat_df.loc[(BaseQ_stat_df.STdev >= Q_stdev_cutoff['lower']) & \
                                      (BaseQ_stat_df.STdev <= Q_stdev_cutoff['upper']), ].index
    Survive_idx_1 = pd.DataFrame({'Read':Survive_idx_1})
    
    print("###### (3) Save survival read")
    Survive_idx_1_file = f'{output_path}/{sampleID}_R{read_N}_Survive_idx_1_bySTdev.csv'
    np.savetxt(Survive_idx_1_file, Survive_idx_1, fmt='%d', delimiter=',')
    print(f'> 1st Survival Reads (By Standard Deviation):\t{len(Survive_idx_1)}\t({round(len(Survive_idx_1)/read_count*100, 5)} %)')
    print("------------------ Save 1st survived reads by standard deviation ------------------ ")
    print()
    
    ### 2. Filtered by base mutation in read
    print("### 2. Filtered by base mutation in read")
    print("###### (1) Importing CSV (Mutation count in known sequence & Starting correct homopolymer)")
    mutCount_df = pd.read_csv(f'{input_path_2}/{sampleID}_knownSeq_mutCount_R{read_N}_df.csv', sep="\t", index_col=0)
    correctStart_idx_file = f'{input_path_2}/{sampleID}_homopolymer_correctStart_R{read_N}_idx.csv'
    correctStart_idx = list(pd.read_csv(correctStart_idx_file, index_col=0, dtype='int').index)
    
    print("###### (2) Read filtering")
    mutCount_cutoff_byError = int(cutoff_mut_ByErrorRate(ErrorRate, KnownSeq_len))
    mutCount_cutoff_byOutlier = int(cutoff_outliers(mutCount_df['Mut_Count'], mutCount_filter_whis)['upper'])
    mutCount_cutoff = max([mutCount_cutoff_byError, mutCount_cutoff_byOutlier])
    print(f'> Mutation Count Cutoff (By Illumina Error Rate):\t{mutCount_cutoff_byError}')
    print(f'> Mutation Count Cutoff (By Mutation Count Outlier):\t{mutCount_cutoff_byOutlier}')
    Survive_idx_mutCount = mutCount_df.loc[(mutCount_df.Mut_Count <= mutCount_cutoff), ].index
    Survive_idx_2 = Survive_idx_mutCount.intersection(correctStart_idx)
    
    print("###### (3) Save survival read")
    Survive_idx_2_file = f'{output_path}/{sampleID}_R{read_N}_Survive_idx_2_byMut.csv'
    np.savetxt(Survive_idx_2_file, Survive_idx_2, fmt='%d', delimiter=',')
    print(f'> 2nd Survival Reads (By Mutation in Read):\t{len(Survive_idx_2)}\t({round(len(Survive_idx_2)/read_count*100, 5)} %)')
    print("----------------- Save 2nd survived reads by base mutation in read ---------------- ")
    print()

    ##### 3) UMI
    if read_N == 2:
        print("### 3. UMI Count ")
        print("###### (1) Importing CSV (UMI sequence & summary)")
        UMI_df = pd.read_csv(f"{input_path_1}/{sampleID}_UMI.txt.gz", header=None, compression='gzip')
        UMI_df.columns = ['UMI_sequence']
        UMI_summary_df = pd.read_csv(f"{input_path_1}/{sampleID}_UMI_summary.txt",  header=None, sep="\t", index_col=0)
        UMI_summary_df.columns = ['Count']
        
        print("###### (2) Separation UMI group (Unique / Duplicate)")
        uniqUMI_seq = UMI_summary_df.index[UMI_summary_df.Count == 1]
        dupUMI_seq = UMI_summary_df.index[UMI_summary_df.Count > 1]
        dup_twoUMI_seq = UMI_summary_df.index[UMI_summary_df.Count == 2]
        dup_threeUMI_seq = UMI_summary_df.index[UMI_summary_df.Count == 3]
        
        uniqUMI_idx = UMI_df[UMI_df.UMI_sequence.isin(uniqUMI_seq)].index
        dupUMI_idx = UMI_df[UMI_df.UMI_sequence.isin(dupUMI_seq)].index
        dup_twoUMI_idx = UMI_df[UMI_df.UMI_sequence.isin(dup_twoUMI_seq)].index
        dup_threeUMI_idx = UMI_df[UMI_df.UMI_sequence.isin(dup_threeUMI_seq)].index
        
        print("###### (3) Save survival read")
        uniqUMI_idx_file = f'{output_path}/{sampleID}_R{read_N}_uniqUMI_idx.csv'
        dupUMI_idx_file = f'{output_path}/{sampleID}_R{read_N}_dupUMI_idx.csv'
        dup_twoUMI_idx_file = f'{output_path}/{sampleID}_R{read_N}_dup2_UMI_idx.csv'
        dup_threeUMI_idx_file = f'{output_path}/{sampleID}_R{read_N}_dup3_UMI_idx.csv'

        np.savetxt(uniqUMI_idx_file, uniqUMI_idx, fmt='%d', delimiter=',')
        np.savetxt(dupUMI_idx_file, dupUMI_idx, fmt='%d', delimiter=',')
        np.savetxt(dup_twoUMI_idx_file, dup_twoUMI_idx, fmt='%d', delimiter=',')
        np.savetxt(dup_threeUMI_idx_file, dup_threeUMI_idx, fmt='%d', delimiter=',')

        print(f'> Read that have unique UMI (Count = 1):\t{len(uniqUMI_idx)}\t({round(len(uniqUMI_idx)/read_count*100, 5)} %)')
        print(f'> Read that have duplicate UMI (Count > 1):\t{len(dupUMI_idx)}\t({round(len(dupUMI_idx)/read_count*100, 5)} %)')
        print(f'> Read that have duplicate UMI (Count = 2):\t{len(dup_twoUMI_idx)}\t({round(len(dup_twoUMI_idx)/read_count*100, 5)} %)')
        print(f'> Read that have duplicate UMI (Count = 3):\t{len(dup_threeUMI_idx)}\t({round(len(dup_threeUMI_idx)/read_count*100, 5)} %)')
        print("--------------------- Save reads that have unique/duplicate UMI ------------------- ")
        print()
    
print("######################################################")
print(f"Summary")  
print("######################################################")
print()
print("### 1. Extract intersect reads in R1 & R2")
R1_Survive_idx_1 = pd.read_csv(f'{output_path}/{sampleID}_R1_Survive_idx_1_bySTdev.csv', index_col=0).index.tolist()
R2_Survive_idx_1 = pd.read_csv(f'{output_path}/{sampleID}_R2_Survive_idx_1_bySTdev.csv', index_col=0).index.tolist()
Survive_idx_1 = sorted(list(set(R1_Survive_idx_1) & set(R2_Survive_idx_1)))
print(f'> 1st Survival Reads (R1):\t{len(R1_Survive_idx_1)}\t({round(len(R1_Survive_idx_1)/read_count*100, 5)} %)')
print(f'> 1st Survival Reads (R2):\t{len(R2_Survive_idx_1)}\t({round(len(R2_Survive_idx_1)/read_count*100, 5)} %)')
print(f'> 1st Survival Reads (Overlap):\t{len(Survive_idx_1)}\t({round(len(Survive_idx_1)/read_count*100, 5)} %)')
Survive_idx_1_file = f'{output_path}/{sampleID}_Survive_idx_1_bySTdev.csv'
np.savetxt(Survive_idx_1_file, Survive_idx_1, fmt='%d', delimiter=',')
print("-------------- Save 1st filtering reads that are intersect in R1 & R2 ------------- ")
print()

R1_Survive_idx_2 = pd.read_csv(f'{output_path}/{sampleID}_R1_Survive_idx_2_byMut.csv', index_col=0).index.tolist()
R2_Survive_idx_2 = pd.read_csv(f'{output_path}/{sampleID}_R2_Survive_idx_2_byMut.csv', index_col=0).index.tolist()
Survive_idx_2 = sorted(list(set(R1_Survive_idx_2) & set(R2_Survive_idx_2)))
print(f'> 2nd Survival Reads (R1):\t{len(R1_Survive_idx_2)}\t({round(len(R1_Survive_idx_2)/read_count*100, 5)} %)')
print(f'> 2nd Survival Reads (R2):\t{len(R2_Survive_idx_2)}\t({round(len(R2_Survive_idx_2)/read_count*100, 5)} %)')
print(f'> 2nd Survival Reads (Overlap):\t{len(Survive_idx_2)}\t({round(len(Survive_idx_2)/read_count*100, 5)} %)')
Survive_idx_2_file = f'{output_path}/{sampleID}_Survive_idx_2_byMut.csv'
np.savetxt(Survive_idx_2_file, Survive_idx_2, fmt='%d', delimiter=',')
print("------------- Save 2nd filtering reads that are intersect in R1 & R2 -------------- ")
print()

print("### 2. Extract intersect reads in 1st and 2nd filtering steps")
Survive_idx_final = sorted(list(set(Survive_idx_1) & set(Survive_idx_2)))
print(f'> 1st Survival Reads:\t{len(Survive_idx_1)}\t({round(len(Survive_idx_1)/read_count*100, 5)} %)')
print(f'> 2nd Survival Reads:\t{len(Survive_idx_2)}\t({round(len(Survive_idx_2)/read_count*100, 5)} %)')
print(f'> Final Survival Reads:\t{len(Survive_idx_final)}\t({round(len(Survive_idx_final)/read_count*100, 5)} %)')
Survive_idx_final_file = f'{output_path}/{sampleID}_Survive_idx_final.csv'
np.savetxt(Survive_idx_final_file, Survive_idx_final, fmt='%d', delimiter=',')
print("--- Save final filtering reads that are intersect in 1st & 2nd filtering steps ---- ")
print()

print(f'-----------------------------------------------------------------------------------------------------')
print(f"End:\t{time.ctime()}"); end_time = time.time()
print(f'-----------------------------------------------------------------------------------------------------')
print()

duration_time = end_time - start_time
conversion_time = time_conversion(duration_time)
print(f"Time:\t{duration_time} ({conversion_time})")