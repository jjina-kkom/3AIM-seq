#!/usr/bin/env python
# coding: utf-8

# In[252]:


import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter1d


# In[254]:


def count_mutations_withoutL(seq, homopolymer_base):
        homopolymer_len = 0
        point = 0
        for base in seq:
            point += 1
            if base == homopolymer_base:
                homopolymer_len += 1
            else:
                if point < 10:
                    seq_window = seq[:10]
                    mutCount = sum(base != homopolymer_base for base in seq_window)
                    if mutCount > 1:
                        break
                    else:
                        homopolymer_len += 1
                        continue
                else:
                    seq_window = seq[point-10:point]
                    mutCount = sum(base != homopolymer_base for base in seq_window)
                    if mutCount > 1:
                        break
                    else:
                        homopolymer_len += 1
                        continue
        return homopolymer_len
    
def count_mutations_beforeL(seq, homopolymer_base):
    homopolymer_len = 0
    point = 0
    for base in seq:
        point += 1
        if base == homopolymer_base:
            homopolymer_len += 1
        else:
            if point < 10:
                seq_window = seq[:10]
                mutCount = sum(base != homopolymer_base for base in seq_window)
                if mutCount > 1:
                    break
                else:
                    homopolymer_len += 1
                    continue
            else:
                seq_window = seq[point-10:point]
                mutCount = sum(base != homopolymer_base for base in seq_window)
                if mutCount > 1:
                    break
                else:
                    homopolymer_len += 1
                    continue
    return homopolymer_len

def count_mutations_linker(seq, start, homopolymer_base):
    seq = seq[start:] 
    linker_len = 1
    for base in seq:
        if ''.join(seq[linker_len:linker_len+10]).count(homopolymer_base) >= 9:
            break
        else:
            linker_len += 1
            
    return linker_len

def count_mutations_afterL(seq, start, homopolymer_base):
    seq = seq[start:] 
    homopolymer_len = 0
    point = 0
    for base in seq:
        point += 1
        if base == homopolymer_base:
            homopolymer_len += 1
        else:
            if point < 10:
                seq_window = seq[:10]
                mutCount = sum(base != homopolymer_base for base in seq_window)
                if mutCount > 1:
                    break
                else:
                    homopolymer_len += 1
                    continue
            else:
                seq_window = seq[point-10:point]
                mutCount = sum(base != homopolymer_base for base in seq_window)
                if mutCount > 1:
                    break
                else:
                    homopolymer_len += 1
                    continue
    return homopolymer_len


# In[255]:


path='/BiO2/research/collab2/Catholic/JHNam/polyA_seq/240622'
info_df = pd.read_csv(f'{path}/info_summary_forpython.txt', sep="\t", index_col=0)


# In[256]:


sampleID_list = list(info_df.index)
sampleID_lenDic = {'STD_A40':[40], 'STD_A50':[50], 'STD_A60':[60], 
                   'Direct-lig_A50':[50], 'Lig-free_A50':[50],
                   'STD_A45L45':[45,10,45], 'STD_A50L50':[50,10,50], 'STD_A55L55':[55,10,55], 
                   'Direct-lig_A50L50':[50,10,50], 'Lig-free_A50L50':[50,10,50], 
                   'STD_A25L65':[25,10,65], 'STD_A30L70':[30,10,70], 'STD_A70L30':[70,10,30],'STD_A35L75':[35,10,75], 
                   'Direct-lig_A30L70':[30,10,70], 'Lig-free_A30L70':[30,10,70], 
                   'STD_A110':[110], 'STD_A120':[120], 'STD_A130':[130], 
                   'Direct-lig_A120':[120], 'Lig-free_A120':[120]}


# In[259]:


path = '/BiO2/research/collab2/Catholic/JHNam/polyA_seq/240622'
output_path = '/BiO2/research/collab2/Catholic/JHNam/polyA_seq/240622/6_Analysis/BaseCall'


# In[261]:


for sampleID in sampleID_list:
    lenList = sampleID_lenDic[sampleID]
    print('---------------------------------------------------------------------------')
    print(sampleID)
    if len(lenList) > 1:
        withL = True
    else:
        withL = False
    
    R1_KnownSeq = info_df.loc[sampleID]['R1_KnownSeq']
    R2_KnownSeq = info_df.loc[sampleID]['R2_KnownSeq']

    R1_homopolymer_base, R2_homopolymer_base = 'T', 'A'
    
    R1_split_Seq_df = pd.read_csv(f'{path}/1_preprocessing/{sampleID}/{sampleID}_Sequence_split.R1.txt.gz',
                               header=None, sep="\t", compression='gzip')
    R2_split_Seq_df = pd.read_csv(f'{path}/1_preprocessing/{sampleID}/{sampleID}_Sequence_split.R2.txt.gz',
                               header=None, sep="\t", compression='gzip')

    SurviveRead_byMut = pd.read_csv(f'{path}/3_filtering_reads/{sampleID}/whis1p5/{sampleID}_Survive_idx_2_byMut.csv',
                               index_col=0, header=None).index.to_list()

    R1_split_Seq_df_survive = R1_split_Seq_df.loc[SurviveRead_byMut]
    R1_split_Seq_df_survive_H_region = R1_split_Seq_df_survive.loc[:, len(R1_KnownSeq):]

    R2_split_Seq_df_survive = R2_split_Seq_df.loc[SurviveRead_byMut]
    R2_split_Seq_df_survive_H_region = R2_split_Seq_df_survive.loc[:, len(R2_KnownSeq):]
    
    try: 
        if not os.path.exists(f'{output_path}/{sampleID}'): 
            os.makedirs(f'{output_path}/{sampleID}')
    except OSError: 
        print("Error: Failed to create the directory.")
        

    R1_polyA_lenDic, R2_polyA_lenDic = {}, {}

    R1_seq_array = R1_split_Seq_df_survive_H_region.values
    R2_seq_array = R2_split_Seq_df_survive_H_region.values

    if withL:
        with alive_bar(len(R1_seq_array), force_tty=True) as bar:
            for idx, row in enumerate(R1_seq_array):
                R1_beforeL_len = count_mutations_beforeL(row, R1_homopolymer_base)
                R1_linker_len = count_mutations_linker(row, R1_beforeL_len, R1_homopolymer_base)
                R1_afterL_len = count_mutations_afterL(row, R1_beforeL_len+R1_linker_len, R1_homopolymer_base)
                R1_polyA_lenDic[idx] = [R1_beforeL_len, R1_linker_len, R1_afterL_len]
                bar()
        with alive_bar(len(R2_seq_array), force_tty=True) as bar:
            for idx, row in enumerate(R2_seq_array):
                R2_beforeL_len = count_mutations_beforeL(row, R2_homopolymer_base)
                R2_linker_len = count_mutations_linker(row, R2_beforeL_len, R2_homopolymer_base)
                R2_afterL_len = count_mutations_afterL(row, R2_beforeL_len+R2_linker_len, R2_homopolymer_base)
                R2_polyA_lenDic[idx] = [R2_beforeL_len, R2_linker_len, R2_afterL_len]
                bar()
                
        R1_polyA_lenDF = pd.DataFrame.from_dict(R1_polyA_lenDic).transpose().rename(columns={0:'Homopolymer_2', 1:'Linker', 2:'Homopolymer_1'})
        R2_polyA_lenDF = pd.DataFrame.from_dict(R2_polyA_lenDic).transpose().rename(columns={0:'Homopolymer_1', 1:'Linker', 2:'Homopolymer_2'})

    else:
        with alive_bar(len(R1_seq_array), force_tty=True) as bar:
            for idx, row in enumerate(R1_seq_array):
                R1_polyA_lenDic[idx] = [count_mutations_withoutL(row, R1_homopolymer_base)]
                bar()
        with alive_bar(len(R2_seq_array), force_tty=True) as bar:
            for idx, row in enumerate(R2_seq_array):
                R2_polyA_lenDic[idx] = [count_mutations_withoutL(row, R2_homopolymer_base)]
                bar()
                
        R1_polyA_lenDF = pd.DataFrame.from_dict(R1_polyA_lenDic).transpose().rename(columns={0:'Homopolymer'})
        R2_polyA_lenDF = pd.DataFrame.from_dict(R2_polyA_lenDic).transpose().rename(columns={0:'Homopolymer'})

    R1_polyA_lenDF.to_csv(f'{output_path}/{sampleID}/{sampleID}_R1_polyA_lenDF_byBaseCall.csv', sep="\t", encoding='utf-8')
    R2_polyA_lenDF.to_csv(f'{output_path}/{sampleID}/{sampleID}_R2_polyA_lenDF_byBaseCall.csv', sep="\t", encoding='utf-8')


# In[ ]:




