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
parser = argparse.ArgumentParser(description='Filtering data by survived reads')

parser.add_argument('--sampleID')
parser.add_argument('--lenDF_file')
parser.add_argument('--read_file')
parser.add_argument('--filter_prefix')
parser.add_argument('--output_path', default=os.getcwd())
parser.add_argument('--polyA_1_len')
parser.add_argument('--Linker_len', default=np.nan)
parser.add_argument('--polyA_2_len', default=np.nan)

args = parser.parse_args()

sampleID = args.sampleID
lenDF_file = args.lenDF_file
read_file = args.read_file
filter_prefix = args.filter_prefix
output_path = args.output_path
polyA_1_len = int(args.polyA_1_len)
Linker_len = int(args.Linker_len) if args.Linker_len is not np.nan else np.nan
polyA_2_len = int(args.polyA_2_len) if args.polyA_2_len is not np.nan else np.nan

##################################
##### Define Custom Function #####
##################################
###################################################################################
def time_conversion(duration_time):
    hour, remainder = divmod(duration_time, 3600)
    minute, second = divmod(remainder, 60)
    
    return f'{int(hour)}h {int(minute)}m {int(second)}s'
###################################################################################
def rmse(p, a):  
    rmse = np.sqrt(np.mean(((p-a)**2)))
    return rmse

def rmse_val(predict_result,y):
    val = rmse(np.array(predict_result), np.array(y))
    return round(val, 5)

def mae(p, a): 
    mae = np.mean(np.abs(p-a))
    return mae

def mae_val(predict_result, y):
    val = mae(np.array(predict_result), np.array(y)) 
    return round(val, 5)
###################################################################################
sampleID = args.sampleID
lenDF_file = args.lenDF_file
read_file = args.read_file
filter_prefix = args.filter_prefix
output_path = args.output_path
polyA_1_len = int(args.polyA_1_len)
Linker_len = int(args.Linker_len) if args.Linker_len is not np.nan else np.nan
polyA_2_len = int(args.polyA_2_len) if args.polyA_2_len is not np.nan else np.nan


print(f'==================================================================================================================================================================')
print(f"Command :")
print(f"python 3-3_Filtering_Data_byReads.py --samplePrefix {sampleID} --filterPrefix {filter_prefix} --lenDF {lenDF_file} --readList {read_file} --output {output_path} --polyA_1 {polyA_1_len} --linker {Linker_len} --polyA_2 {polyA_2_len}")
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
print()
if polyA_2_len is np.nan:
    polyA_regionDic = {"Homopolymer": polyA_1_len}
else:
    polyA_regionDic = {"Homopolymer_1": polyA_1_len,
                       "Linker": Linker_len,
                       "Homopolymer_2": polyA_2_len}


### 1. Import data
print("### 1. Importing CSV (final polyA length / Survival reads)")
final_lenDF = pd.read_csv(lenDF_file, sep="\t", index_col=0)
survive_readDF = pd.read_csv(read_file, sep="\t", index_col=0, header=None)

final_lenDF_idx = final_lenDF.index
survive_readDF_idx = survive_readDF.index
print()

### 2. Read Count that filtered out in polyA length measured step
print("### 2. Read Count that filtered out in polyA length measured step")
overlap_idx = set(final_lenDF_idx) & set(survive_readDF_idx)
measuredStep_fail_readCount = len(survive_readDF_idx) - len(overlap_idx)
print(f'> Total count of survival read: {len(overlap_idx)}')
print(f'> Fail reads in polyA length measurement step: {measuredStep_fail_readCount} ({round(measuredStep_fail_readCount/len(survive_readDF_idx) * 100, 2)} %)')
print()

### 3. Extract polyA length by survival reads
print("### 3. Extract polyA length by survival reads")
final_lenDF_bySurvive = final_lenDF.loc[overlap_idx].sort_index()
final_lenDF_bySurvive.to_csv(f'{output_path}/{sampleID}_Final_polyA_lenDF.{filter_prefix}.csv', sep="\t")
print("----------------- Save poly(A) lengths filtered by survival reads ----------------- ")
print()

### 4. Calculalte RMSE & MAE
print("### 4. Calculalte RMSE & MAE")

i = 1
for region, ref_len in polyA_regionDic.items():
    print(f"###### ({i}) {region}"); i += 1
    
    rmse_previous = rmse_val(final_lenDF[region],[ref_len]*len(final_lenDF))
    mae_previous = mae_val(final_lenDF[region],[ref_len]*len(final_lenDF))
    
    rmse_final = rmse_val(final_lenDF_bySurvive[region],[ref_len]*len(final_lenDF_bySurvive))
    mae_final = mae_val(final_lenDF_bySurvive[region],[ref_len]*len(final_lenDF_bySurvive))
    print(f"> Type\tBefore\tAfter")
    print(f"> RMSE\t{rmse_previous}\t{rmse_final}")
    print(f"> MAE\t{mae_previous}\t{mae_final}")
    print()
    
print(f'-----------------------------------------------------------------------------------------------------')
print(f"End:\t{time.ctime()}"); end_time = time.time()
print(f'-----------------------------------------------------------------------------------------------------')
print()

duration_time = end_time - start_time
conversion_time = time_conversion(duration_time)
print(f"Time:\t{duration_time} ({conversion_time})")
