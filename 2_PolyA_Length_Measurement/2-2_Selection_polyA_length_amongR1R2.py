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
parser = argparse.ArgumentParser(description='Selection polyA length among R1 and R2')

parser.add_argument('--sampleID')
parser.add_argument('--input_path')
parser.add_argument('--output_path', default=os.getcwd())
parser.add_argument('--polyA_1_len')
parser.add_argument('--Linker_len', default=np.nan)
parser.add_argument('--polyA_2_len', default=np.nan)

args = parser.parse_args()

sampleID = args.sampleID
input_path = args.input_path
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
def Select_polyA_len_withLinker(R1_lenDF, R2_lenDF, polyA_lenDic):
        
    R1_lenDF['H1_diff'] = abs(R1_lenDF['Homopolymer_1'] - polyA_lenDic['Homopolymer_1'])
    R1_lenDF['L_diff'] = abs(R1_lenDF['Linker'] - polyA_lenDic['Linker'])
    R1_lenDF['H2_diff'] = abs(R1_lenDF['Homopolymer_2'] - polyA_lenDic['Homopolymer_2'])

    R2_lenDF['H1_diff'] = abs(R2_lenDF['Homopolymer_1'] - polyA_lenDic['Homopolymer_1'])
    R2_lenDF['L_diff'] = abs(R2_lenDF['Linker'] - polyA_lenDic['Linker'])
    R2_lenDF['H2_diff'] = abs(R2_lenDF['Homopolymer_2'] - polyA_lenDic['Homopolymer_2'])

    ## Homopolymer_1 ##
    H1_finalRead_none = R1_lenDF.loc[(R1_lenDF['Homopolymer_1'].isnull()) & (R2_lenDF['Homopolymer_1'].isnull())].index
    H1_finalRead_byR1_NaN = R1_lenDF.loc[R2_lenDF['Homopolymer_1'].isnull()].index
    H1_finalRead_byR2_NaN = R1_lenDF.loc[R1_lenDF['Homopolymer_1'].isnull()].index
    H1_finalRead_byR1_diff = R1_lenDF.loc[R1_lenDF['H1_diff'] <= R2_lenDF['H1_diff']].index
    H1_finalRead_byR2_diff = R1_lenDF.loc[R1_lenDF['H1_diff'] > R2_lenDF['H1_diff']].index

    H1_finalRead_byR1 = pd.Index(H1_finalRead_byR1_NaN.tolist() + H1_finalRead_byR1_diff.tolist())
    H1_finalRead_byR2 = pd.Index(H1_finalRead_byR2_NaN.tolist() + H1_finalRead_byR2_diff.tolist())

    ## Homopolymer_2 ##
    H2_finalRead_none = R1_lenDF.loc[(R1_lenDF['Homopolymer_2'].isnull()) & (R2_lenDF['Homopolymer_2'].isnull())].index
    H2_finalRead_byR1_NaN = R1_lenDF.loc[R2_lenDF['Homopolymer_2'].isnull()].index
    H2_finalRead_byR2_NaN = R1_lenDF.loc[R1_lenDF['Homopolymer_2'].isnull()].index
    H2_finalRead_byR1_diff = R1_lenDF.loc[R1_lenDF['H2_diff'] < R2_lenDF['H2_diff']].index
    H2_finalRead_byR2_diff = R1_lenDF.loc[R1_lenDF['H2_diff'] >= R2_lenDF['H2_diff']].index

    H2_finalRead_byR1 = pd.Index(H2_finalRead_byR1_NaN.tolist() + H2_finalRead_byR1_diff.tolist())
    H2_finalRead_byR2 = pd.Index(H2_finalRead_byR2_NaN.tolist() + H2_finalRead_byR2_diff.tolist())
    
    ## Linker ##
    L_finalRead_none = R1_lenDF.loc[(R1_lenDF['Linker'].isnull()) & (R2_lenDF['Linker'].isnull())].index
    L_finalRead_byR1_NaN = R1_lenDF.loc[R2_lenDF['Linker'].isnull()].index
    L_finalRead_byR2_NaN = R1_lenDF.loc[R1_lenDF['Linker'].isnull()].index
    L_finalRead_byR1_overlap = H1_finalRead_byR1.intersection(H2_finalRead_byR2)
    L_finalRead_byR1_H1 = H1_finalRead_byR1.difference(L_finalRead_byR1_overlap.union(L_finalRead_byR1_NaN).union(L_finalRead_byR2_NaN))
    L_finalRead_byR2_H2 = H2_finalRead_byR2.difference(L_finalRead_byR1_overlap.union(L_finalRead_byR2_NaN).union(L_finalRead_byR1_NaN))
    L_finalRead_byR1_remaining = R1_lenDF.index.difference(L_finalRead_byR1_overlap.union(L_finalRead_byR1_NaN).union(L_finalRead_byR2_NaN).union(L_finalRead_byR1_H1).union(L_finalRead_byR2_H2))

    L_finalRead_byR1 = pd.Index(L_finalRead_byR1_NaN.tolist() + L_finalRead_byR1_overlap.tolist() + L_finalRead_byR1_H1.tolist() + L_finalRead_byR1_remaining.tolist())
    L_finalRead_byR2 = pd.Index(L_finalRead_byR2_NaN.tolist() + L_finalRead_byR2_H2.tolist())
    
    ## Summary ##
    finalRead_none = H1_finalRead_none.union(H2_finalRead_none).union(L_finalRead_none)
    
    None_list, H1_lenDic, H2_lenDic, L_lenDic = [], {}, {}, {}
    
    with alive_bar(R1_lenDF.shape[0], force_tty=True) as bar:
        for idx in R1_lenDF.index:
            if idx in finalRead_none:
                None_list.append(idx)
            else:
                H1_len = R1_lenDF.loc[idx, 'Homopolymer_1'] if idx in H1_finalRead_byR1 else R2_lenDF.loc[idx, 'Homopolymer_1']
                H1_lenDic[idx] = [H1_len, 1 if idx in H1_finalRead_byR1 else 2]

                H2_len = R1_lenDF.loc[idx, 'Homopolymer_2'] if idx in H2_finalRead_byR1 else R2_lenDF.loc[idx, 'Homopolymer_2']
                H2_lenDic[idx] = [H2_len, 1 if idx in H2_finalRead_byR1 else 2]

                L_len = R1_lenDF.loc[idx, 'Linker'] if idx in L_finalRead_byR1 else R2_lenDF.loc[idx, 'Linker']
                L_lenDic[idx] = [L_len, 1 if idx in L_finalRead_byR1 else 2]
                
            bar()
            
    H1_lenDF = pd.DataFrame(H1_lenDic).transpose()
    H2_lenDF = pd.DataFrame(H2_lenDic).transpose()
    L_lenDF = pd.DataFrame(L_lenDic).transpose()

    H1_lenDF.columns = ['Length', 'Read_inf']
    H2_lenDF.columns = ['Length', 'Read_inf']
    L_lenDF.columns = ['Length', 'Read_inf']

    final_lenDF = H1_lenDF[['Length']].rename(columns={'Length': 'Homopolymer_1'})
    final_lenDF['Linker'] = L_lenDF['Length']
    final_lenDF['Homopolymer_2'] = H2_lenDF['Length']

    final_readInfDF = H1_lenDF[['Read_inf']].rename(columns={'Read_inf': 'Homopolymer_1'})
    final_readInfDF['Linker'] = L_lenDF['Read_inf']
    final_readInfDF['Homopolymer_2'] = H2_lenDF['Read_inf']
    
    return {'NA': None_list, 'Length':final_lenDF, 'Read_inf': final_readInfDF}
###################################################################################
def Select_polyA_len_withoutLinker(R1_lenDF, R2_lenDF, polyA_lenDic):
        
    R1_lenDF['H_diff'] = abs(R1_lenDF['Homopolymer'] - polyA_lenDic['Homopolymer'])
    R2_lenDF['H_diff'] = abs(R2_lenDF['Homopolymer'] - polyA_lenDic['Homopolymer'])

    ## Homopolymer ##
    H_finalRead_none = R1_lenDF.loc[(R1_lenDF['Homopolymer'].isnull()) & (R2_lenDF['Homopolymer'].isnull())].index
    H_finalRead_byR1_NaN = R1_lenDF.loc[R2_lenDF['Homopolymer'].isnull()].index
    H_finalRead_byR2_NaN = R1_lenDF.loc[R1_lenDF['Homopolymer'].isnull()].index
    H_finalRead_byR1_diff = R1_lenDF.loc[R1_lenDF['H_diff'] < R2_lenDF['H_diff']].index
    H_finalRead_byR2_diff = R1_lenDF.loc[R1_lenDF['H_diff'] >= R2_lenDF['H_diff']].index

    H_finalRead_byR1 = pd.Index(H_finalRead_byR1_NaN.tolist() + H_finalRead_byR1_diff.tolist())
    H_finalRead_byR2 = pd.Index(H_finalRead_byR2_NaN.tolist() + H_finalRead_byR2_diff.tolist())

    #### Summary ######
    finalRead_none = H_finalRead_none
    None_list, H_lenDic = [], {}
    
    with alive_bar(R1_lenDF.shape[0], force_tty=True) as bar:
        for idx in R1_lenDF.index:
            if idx in finalRead_none:
                None_list.append(idx)
            else:
                H_len = R1_lenDF.loc[idx, 'Homopolymer'] if idx in H_finalRead_byR1 else R2_lenDF.loc[idx, 'Homopolymer']
                H_lenDic[idx] = [H_len, 1 if idx in H_finalRead_byR1 else 2]

            bar()
            
    H_lenDF = pd.DataFrame(H_lenDic).transpose()
    H_lenDF.columns = ['Length', 'Read_inf']

    final_lenDF = H_lenDF[['Length']].rename(columns={'Length': 'Homopolymer'})
    final_readInfDF = H_lenDF[['Read_inf']].rename(columns={'Read_inf': 'Homopolymer'})
    
    return {'NA': None_list, 'Length':final_lenDF, 'Read_inf': final_readInfDF}
###################################################################################

print(f'==================================================================================================================================================================')
print(f"Command :")
print(f"python 2-2_Selection_polyA_length_amongR1R2.py --prefix {sampleID} --input {input_path} --output {output_path} --polyA_1 {polyA_1_len} --linker {Linker_len} --polyA_2 {polyA_2_len}")
print(f'==================================================================================================================================================================')
print()

try: 
    if not os.path.exists(f'{output_path}'): 
        os.makedirs(f'{output_path}')
except OSError: 
    print("Error: Failed to create the directory.")

if polyA_2_len is np.nan:
    withLinker, polyA_lenDic = False, {'Homopolymer': polyA_1_len}
    Select_polyA_len_func = Select_polyA_len_withoutLinker
else:
    withLinker, polyA_lenDic = True, {'Homopolymer_1': polyA_1_len, 'Linker': Linker_len, 'Homopolymer_2': polyA_2_len}
    Select_polyA_len_func = Select_polyA_len_withLinker
    
print(f'-----------------------------------------------------------------------------------------------------')
print(f"Start:\t{time.ctime()}"); start_time = time.time()
print(f'-----------------------------------------------------------------------------------------------------')
print()

### 1. Import data
print("### 1. Importing CSV (Measured Length among R1 and R2 data)")
R1_lenDF = pd.read_csv(f"{input_path}/{sampleID}_R1_polyA_lenDF.csv", sep="\t", index_col=0)
R2_lenDF = pd.read_csv(f"{input_path}/{sampleID}_R2_polyA_lenDF.csv", sep="\t", index_col=0)
print()

### 2. Selection polyA length
print("### 2. Selection PolyA Length")
Select_polyA_data = Select_polyA_len_func(R1_lenDF, R2_lenDF, polyA_lenDic)

NA_read_lst = Select_polyA_data['NA']
if len(NA_read_lst) > 0:
    idx_file = f'{output_path}/{sampleID}_selection_fail_idx.csv'
    np.savetxt(idx_file, list(NA_read_lst), fmt='%d', delimiter=',')
    print("--------------------- Save NA read list from selection step -----------------------")

final_lenDF = Select_polyA_data['Length']
final_lenDF.to_csv(f'{output_path}/{sampleID}_Final_polyA_lenDF.csv', sep="\t")
print("----------------------- Save selected poly(A) length table ------------------------")

final_readInfDF = Select_polyA_data['Read_inf']
final_readInfDF.to_csv(f'{output_path}/{sampleID}_Final_select_readInfDF.csv', sep="\t", encoding='utf-8')
print("---------------------- Save selected read information table -----------------------")
print()

print(f'-----------------------------------------------------------------------------------------------------')
print(f"End:\t{time.ctime()}"); end_time = time.time()
print(f'-----------------------------------------------------------------------------------------------------')
print()

duration_time = end_time - start_time
conversion_time = time_conversion(duration_time)
print(f"Time:\t{duration_time} ({conversion_time})")