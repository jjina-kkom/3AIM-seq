#!/usr/bin/env python
# coding: utf-8

# In[38]:


import pandas as pd
import numpy as np
from scipy.ndimage import gaussian_filter1d
from multiprocessing import Pool, cpu_count
from alive_progress import alive_bar


# In[39]:


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

def Finding_Peak(positions, values, cutoff):
    continuousLoc_lst, collection_lst = [], []
    
    for pos, value in zip(positions, values):
        if value > 0:
            continuousLoc_lst.append(pos)
        else:
            if len(continuousLoc_lst) >= cutoff:
                collection_lst.append(continuousLoc_lst[-1])
            continuousLoc_lst = []

    if len(continuousLoc_lst) >= cutoff:
        collection_lst.append(continuousLoc_lst[-1])
    
    return collection_lst

def hist_maxPeak(data, readLen, binSize, cutoff):
    hist, bin_edges = np.histogram(data, range=(0, readLen), bins=binSize, density=True)
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    smoothed_hist = gaussian_filter1d(hist, sigma=1)
    
    summary_df = pd.DataFrame({'density': smoothed_hist, 'x': bin_centers})
    summary_df['diff'] = summary_df['density'].diff().fillna(0)
    
    pred_idx = Finding_Peak(summary_df.index, summary_df['diff'], cutoff)
    summary_df = summary_df.loc[pred_idx].sort_values(by='density', ascending=False)
    summary_dic = dict(zip(summary_df['x'], summary_df['density']))

    return summary_dic


# In[40]:


path='/BiO2/research/collab2/Catholic/JHNam/polyA_seq/240622'


# ### [1] Poly(A) without Linker

# In[ ]:


ref_len_list = [40,50,60,110,120,130]
spike_len_list = [40,50,60,110,120,130]

def process_combination(args):
    ref_len, spike_len, spike_prop, path, total_count, cutoff = args
    withoutL_mixDF = pd.read_csv(f'{path}/6_Analysis/STD_mixData/withoutL/STD_A{ref_len}_{spike_len}/STD_A{ref_len}_{spike_len}-{spike_prop}_withoutL_mixData.csv', 
                                 sep='\t', index_col=0)
    withoutL_mixDF = withoutL_mixDF.dropna(how='all')
    results = []

    for rep in range(1, 11):
        df = withoutL_mixDF[[str(rep)]].dropna()
        total_len = df.index.repeat(df[str(rep)].astype(int)).tolist()

        if len(total_len) == total_count:
            RMSE_mix = rmse_val(total_len, [ref_len] * len(total_len))
            MAE_mix = mae_val(total_len, [ref_len] * len(total_len))
            Peak_dic = hist_maxPeak(total_len, 300, 100, cutoff)
            maxPeak = max(Peak_dic, key=Peak_dic.get)
            maxPeak_density = Peak_dic[maxPeak]
            hist_maxPeak_dic = {k: v for k, v in Peak_dic.items() if v >= maxPeak_density * 0.01}
            results.append([rep, ref_len, spike_len, spike_prop / 100, RMSE_mix, MAE_mix, hist_maxPeak_dic])
        else:
            print('Wrong')
    
    return results

def main(max_workers=None):
    spike_len_list = [40, 50, 60, 110, 120, 130]
    ref_len_list = [40, 50, 60, 110, 120, 130]

    total_count = 140000
    cutoff = 2
    Final_data = []
    path = '/BiO2/research/collab2/Catholic/JHNam/polyA_seq/240622'

    if max_workers is None:
        max_workers = 50
        
    with Pool(processes=max_workers) as pool:
        args_list = [
            (ref_len, spike_len, spike_prop, path, total_count, cutoff)
            for ref_len in ref_len_list
            for spike_len in spike_len_list if spike_len != ref_len
            for spike_prop in range(0, 101)
        ]
        for results in pool.imap(process_combination, args_list):
            Final_data.extend(results)
            
    columns = ['Replicate', 'Ref', 'Spike', 'Spike%', 'RMSE', 'MAE', 'maxPeak']
    Final_data_df = pd.DataFrame(Final_data, columns=columns)
    
    Final_data_df.to_csv(f'{path}/6_Analysis/STD_mixData/withoutL_mixData_DF.csv', sep="\t", encoding='utf-8')

if __name__ == "__main__":
    main()


# ### [2-1] 1st Poly(A) with Linker

# In[41]:


ref_lenList = [[45,10,45], [50,10,50], [55,10,55], [25,10,65], [30,10,70], [70,10,30], [35,10,75]]
spike_lenList = [[45,10,45], [50,10,50], [55,10,55], [25,10,65], [30,10,70], [70,10,30], [35,10,75]]

def process_combination(args):
    ref_list, spike_list, spike_prop, path, total_count, cutoff = args
    ref_len, spike_len = ref_list[0], spike_list[0]
    
    withL_mixDF = pd.read_csv(f'{path}/6_Analysis/STD_mixData/withL_H1/STD_A{ref_list[0]}L{ref_list[-1]}_{spike_list[0]}L{spike_list[-1]}/STD_A{ref_list[0]}L{ref_list[-1]}_{spike_list[0]}L{spike_list[-1]}-{spike_prop}_withL_H1_mixData.csv', 
                                 sep='\t', index_col=0)
    withL_mixDF = withL_mixDF.dropna(how='all')
    results = []

    for rep in range(1, 11):
        df = withL_mixDF[[str(rep)]].dropna()
        total_len = df.index.repeat(df[str(rep)].astype(int)).tolist()

        if len(total_len) == total_count:
            RMSE_mix = rmse_val(total_len, [ref_len] * len(total_len))
            MAE_mix = mae_val(total_len, [ref_len] * len(total_len))
            Peak_dic = hist_maxPeak(total_len, 300, 100, cutoff)
            maxPeak = max(Peak_dic, key=Peak_dic.get)
            maxPeak_density = Peak_dic[maxPeak]
            hist_maxPeak_dic = {k: v for k, v in Peak_dic.items() if v >= maxPeak_density * 0.01}
            results.append([rep, ref_len, spike_len, spike_prop / 100, RMSE_mix, MAE_mix, hist_maxPeak_dic])
        else:
            print('Wrong')
    
    return results

def main(max_workers=None):
    spike_len_list = spike_lenList
    ref_len_list = ref_lenList
    total_count = 140000
    cutoff = 2
    Final_data = []
    path = '/BiO2/research/collab2/Catholic/JHNam/polyA_seq/240622'

    if max_workers is None:
        max_workers = 50
        
    with Pool(processes=max_workers) as pool:
        args_list = [
            (ref_list, spike_list, spike_prop, path, total_count, cutoff)
            for ref_list in ref_len_list
            for spike_list in spike_len_list if spike_list != ref_list
            for spike_prop in range(0, 101)
        ]
        for results in pool.imap(process_combination, args_list):
            Final_data.extend(results)
            
    columns = ['Replicate', 'Ref', 'Spike', 'Spike%', 'RMSE', 'MAE', 'maxPeak']
    Final_data_df = pd.DataFrame(Final_data, columns=columns)
    
    Final_data_df.to_csv(f'{path}/6_Analysis/STD_mixData/withL_H1_mixData_DF.csv', sep="\t", encoding='utf-8')

if __name__ == "__main__":
    main()


# ### [2-2] 2nd Poly(A) with Linker

# In[2]:


ref_lenList = [[45,10,45], [50,10,50], [55,10,55], [25,10,65], [30,10,70], [70,10,30], [35,10,75]]
spike_lenList = [[45,10,45], [50,10,50], [55,10,55], [25,10,65], [30,10,70], [70,10,30], [35,10,75]]

def process_combination(args):
    ref_list, spike_list, spike_prop, path, total_count, cutoff = args
    ref_len, spike_len = ref_list[-1], spike_list[-1]
    
    withL_mixDF = pd.read_csv(f'{path}/6_Analysis/STD_mixData/withL_H2/STD_A{ref_list[0]}L{ref_list[-1]}_{spike_list[0]}L{spike_list[-1]}/STD_A{ref_list[0]}L{ref_list[-1]}_{spike_list[0]}L{spike_list[-1]}-{spike_prop}_withL_H2_mixData.csv', 
                                 sep='\t', index_col=0)
    withL_mixDF = withL_mixDF.dropna(how='all')
    results = []

    for rep in range(1, 11):
        df = withL_mixDF[[str(rep)]].dropna()
        total_len = df.index.repeat(df[str(rep)].astype(int)).tolist()

        if len(total_len) == total_count:
            RMSE_mix = rmse_val(total_len, [ref_len] * len(total_len))
            MAE_mix = mae_val(total_len, [ref_len] * len(total_len))
            Peak_dic = hist_maxPeak(total_len, 300, 100, cutoff)
            maxPeak = max(Peak_dic, key=Peak_dic.get)
            maxPeak_density = Peak_dic[maxPeak]
            hist_maxPeak_dic = {k: v for k, v in Peak_dic.items() if v >= maxPeak_density * 0.01}
            results.append([rep, ref_len, spike_len, spike_prop / 100, RMSE_mix, MAE_mix, hist_maxPeak_dic])
        else:
            print('Wrong')
    
    return results

def main(max_workers=None):
    spike_len_list = spike_lenList
    ref_len_list = ref_lenList
    
    total_count = 140000
    cutoff = 2
    Final_data = []
    path = '/BiO2/research/collab2/Catholic/JHNam/polyA_seq/240622'

    if max_workers is None:
        max_workers = 50
        
    with Pool(processes=max_workers) as pool:
        args_list = [
            (ref_list, spike_list, spike_prop, path, total_count, cutoff)
            for ref_list in ref_len_list
            for spike_list in spike_len_list if spike_list != ref_list
            for spike_prop in range(0, 101)
        ]
        for results in pool.imap(process_combination, args_list):
            Final_data.extend(results)
            
    columns = ['Replicate', 'Ref', 'Spike', 'Spike%', 'RMSE', 'MAE', 'maxPeak']
    Final_data_df = pd.DataFrame(Final_data, columns=columns)
    
    Final_data_df.to_csv(f'{path}/6_Analysis/STD_mixData/withL_H2_mixData_DF.csv', sep="\t", encoding='utf-8')

if __name__ == "__main__":
    main()


# In[ ]:




