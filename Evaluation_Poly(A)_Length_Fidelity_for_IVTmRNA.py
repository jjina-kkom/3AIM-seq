#!/usr/bin/env python
# coding: utf-8

# In[227]:


import pandas as pd
import numpy as np
from scipy.ndimage import gaussian_filter1d
from alive_progress import alive_bar


# In[240]:


def mae(p, a): 
    mae = np.mean(np.abs(p-a))
    return mae

def mae_val(predict_result, y):
    val = mae(np.array(predict_result), np.array(y)) 
    return round(val, 2)

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
    summary_dic = dict(zip(np.ceil(summary_df['x']), summary_df['density']))

    return summary_dic

def errorRange(value, error, n):
    value_min = round(float(value) - float(error), n)
    value_max = round(float(value) + float(error), n)
    
    return (value_min, value_max)

def equation_execute(A,B,C, equation):
    return eval(equation)

def evalution_mRNA_polyA(ref_len_list, spike_len_list, spike_prop_list, pred_diff_m, pred_diff_M, equation):
    possible_values = {}

    if type(ref_len_list) == int:
        ref_len_list = [ref_len_list]
        
    for ref_len in ref_len_list:
        spike_possible_values = {}
        
        for spike_len in spike_len_list:
            Diff_len = int(spike_len - ref_len)
            possible_prop_list = []
            
            for spike_prop in spike_prop_list:
                y_diff = equation_execute(ref_len, Diff_len, spike_prop, equation)
                
                if pred_diff_m <= y_diff <= pred_diff_M:
                    possible_prop_list.append(spike_prop)

            if len(possible_prop_list) > 0:
                spike_possible_values[spike_len] = possible_prop_list
                
        if len(spike_possible_values) > 0:
            possible_values[ref_len] = spike_possible_values

    return possible_values



def Length_fromPeak_x(peak_x, coef, intercept):
    length = (peak_x - intercept)/coef
    
    return length


# In[229]:


path='/BiO2/research/collab2/Catholic/JHNam/polyA_seq/240622'
density_cutoff=0.25


# # 1. Calculated MAE & Peak Position

# ### [1] Poly(A) without Linker

# In[258]:


withoutL_ref_list = [50,120]
Final_data = []

for ref_len in withoutL_ref_list:
    cutoff, region = 2, 'Homopolymer'
    withoutL_mixDF = pd.read_csv(f'{path}/4_filtering_data/Lig-free_A{ref_len}/Lig-free_A{ref_len}_Final_polyA_lenDF.final_uniqUMI.csv', 
                                 sep='\t', index_col=0)
    results = [] 
    lenList = list(withoutL_mixDF[region])
    
    MAE_mix = mae_val(lenList, [ref_len] * len(lenList))
    Peak_dic = hist_maxPeak(lenList, 300, 100, cutoff)
    maxPeak = max(Peak_dic, key=Peak_dic.get)
    maxPeak_density = Peak_dic[maxPeak]
    hist_maxPeak_dic = {k: v for k, v in Peak_dic.items() if v >= maxPeak_density * density_cutoff}
    primaryPeak = list(hist_maxPeak_dic.keys())[0]
    primaryPeak_val = round(hist_maxPeak_dic[primaryPeak], 5)
    if len(hist_maxPeak_dic) > 1:
        secondaryPeak = list(hist_maxPeak_dic.keys())[1]
        secondaryPeak_val = round(hist_maxPeak_dic[secondaryPeak], 5)
    else:
        secondaryPeak, secondaryPeak_val = 0, 0
        
    Final_data.append([f'Lig-free_A{ref_len}', ref_len, MAE_mix, primaryPeak, primaryPeak_val, secondaryPeak, secondaryPeak_val])
    
columns = ['SampleID', 'Length', 'MAE', 'primaryPeak_x', 'primaryPeak_y', 'secondaryPeak_x', 'secondaryPeak_y']
withoutL_mRNA_sumamry_df = pd.DataFrame(Final_data, columns=columns)
withoutL_mRNA_sumamry_df.set_index('SampleID', inplace=True)

withoutL_mRNA_sumamry_df


# ### [2-1] 1st Poly(A) with Linker

# In[259]:


withL_H1_ref_list = [[30, 10, 70],[50, 10, 50]]
Final_data = []

for ref_list in withL_H1_ref_list:
    ref_len = ref_list[0]
    cutoff, region = 2, 'Homopolymer_1'
    withL_mixDF = pd.read_csv(f'{path}/4_filtering_data/Lig-free_A{ref_list[0]}L{ref_list[-1]}/Lig-free_A{ref_list[0]}L{ref_list[-1]}_Final_polyA_lenDF.final_uniqUMI.csv', 
                                 sep='\t', index_col=0)
    results = [] 
    lenList = list(withL_mixDF[region])
    
    MAE_mix = mae_val(lenList, [ref_len] * len(lenList))
    Peak_dic = hist_maxPeak(lenList, 300, 100, cutoff)
    maxPeak = max(Peak_dic, key=Peak_dic.get)
    maxPeak_density = Peak_dic[maxPeak]
    hist_maxPeak_dic = {k: v for k, v in Peak_dic.items() if v >= maxPeak_density * density_cutoff}
    primaryPeak = list(hist_maxPeak_dic.keys())[0]
    primaryPeak_val = round(hist_maxPeak_dic[primaryPeak], 5)
    if len(hist_maxPeak_dic) > 1:
        secondaryPeak = list(hist_maxPeak_dic.keys())[1]
        secondaryPeak_val = round(hist_maxPeak_dic[secondaryPeak], 5)
    else:
        secondaryPeak, secondaryPeak_val = 0, 0
        
    Final_data.append([f'Lig-free_A{ref_list[0]}L{ref_list[-1]}', ref_len, MAE_mix, primaryPeak, primaryPeak_val, secondaryPeak, secondaryPeak_val])
    
columns = ['SampleID', 'Length', 'MAE', 'primaryPeak_x', 'primaryPeak_y', 'secondaryPeak_x', 'secondaryPeak_y']
withL_H1_mRNA_sumamry_df = pd.DataFrame(Final_data, columns=columns)
withL_H1_mRNA_sumamry_df.set_index('SampleID', inplace=True)

withL_H1_mRNA_sumamry_df


# ### [2-2] 2nd Poly(A) with Linker

# In[260]:


withL_H2_ref_list = [[30, 10, 70],[50, 10, 50]]
Final_data = []

for ref_list in withL_H2_ref_list:
    ref_len = ref_list[-1]
    cutoff, region = 2, 'Homopolymer_2'
    withL_mixDF = pd.read_csv(f'{path}/4_filtering_data/Lig-free_A{ref_list[0]}L{ref_list[-1]}/Lig-free_A{ref_list[0]}L{ref_list[-1]}_Final_polyA_lenDF.final_uniqUMI.csv', 
                                 sep='\t', index_col=0)
    results = [] 
    lenList = list(withL_mixDF[region])
    
    MAE_mix = mae_val(lenList, [ref_len] * len(lenList))
    Peak_dic = hist_maxPeak(lenList, 300, 100, cutoff)
    maxPeak = max(Peak_dic, key=Peak_dic.get)
    maxPeak_density = Peak_dic[maxPeak]
    hist_maxPeak_dic = {k: v for k, v in Peak_dic.items() if v >= maxPeak_density * density_cutoff}
    primaryPeak = list(hist_maxPeak_dic.keys())[0]
    primaryPeak_val = round(hist_maxPeak_dic[primaryPeak], 5)
    if len(hist_maxPeak_dic) > 1:
        secondaryPeak = list(hist_maxPeak_dic.keys())[1]
        secondaryPeak_val = round(hist_maxPeak_dic[secondaryPeak], 5)
    else:
        secondaryPeak, secondaryPeak_val = 0, 0
        
    Final_data.append([f'Lig-free_A{ref_list[0]}L{ref_list[-1]}', ref_len, MAE_mix, primaryPeak, primaryPeak_val, secondaryPeak, secondaryPeak_val])

columns = ['SampleID', 'Length', 'MAE', 'primaryPeak_x', 'primaryPeak_y', 'secondaryPeak_x', 'secondaryPeak_y']
withL_H2_mRNA_sumamry_df = pd.DataFrame(Final_data, columns=columns)
withL_H2_mRNA_sumamry_df.set_index('SampleID', inplace=True)

withL_H2_mRNA_sumamry_df


# # 2. Loading Information of Polynomial Regression

# In[261]:


withoutL_polyReg_info_df = pd.read_csv(f'{path}/6_Analysis/STD_calibration/3_STD_mixData_PolyReg/withoutL_polyReg_info_df.csv',
                               sep='\t', index_col=0)
withL_H1_polyReg_info_df = pd.read_csv(f'{path}/6_Analysis/STD_calibration/3_STD_mixData_PolyReg/withL_H1_polyReg_info_df.csv',
                               sep='\t', index_col=0)
withL_H2_polyReg_info_df = pd.read_csv(f'{path}/6_Analysis/STD_calibration/3_STD_mixData_PolyReg/withL_H2_polyReg_info_df.csv',
                               sep='\t', index_col=0)


# In[262]:


withoutL_polyReg_MAE = withoutL_polyReg_info_df['polyReg'].loc['MAE']
withoutL_polyReg_equation = withoutL_polyReg_info_df['polyReg'].loc['Equation'] 

withL_H1_polyReg_MAE = withL_H1_polyReg_info_df['polyReg'].loc['MAE']
withL_H1_polyReg_equation = withL_H1_polyReg_info_df['polyReg'].loc['Equation'] 

withL_H2_polyReg_MAE = withL_H2_polyReg_info_df['polyReg'].loc['MAE']
withL_H2_polyReg_equation = withL_H2_polyReg_info_df['polyReg'].loc['Equation'] 


# In[263]:


result_dic = {}
for file_type in ['withoutL', 'withL_H1', 'withL_H2']:
    result_dic[file_type] = {}
    mRNA_list = {'withoutL':['Lig-free_A50', 'Lig-free_A120'], 
                'withL_H1':['Lig-free_A50L50', 'Lig-free_A30L70'], 
                'withL_H2':['Lig-free_A50L50', 'Lig-free_A30L70'], }
    mRNA_summaryDF = globals()[f'{file_type}_mRNA_sumamry_df']

    pred_DF = pd.read_csv(f'{path}/6_Analysis/STD_calibration/2_STD_Data_Equation/{file_type}_STD_predict_df.csv',sep = '\t', index_col=0)
    pred_equation_coef_DF = pd.read_csv(f'{path}/6_Analysis/STD_calibration/2_STD_Data_Equation/STD_predict_equation_coef_df.csv',sep = '\t', index_col=0)
    pred_equation_intercept_DF = pd.read_csv(f'{path}/6_Analysis/STD_calibration/2_STD_Data_Equation/STD_predict_equation_intercept_df.csv',sep = '\t', index_col=0)
    pred_error_DF = pd.read_csv(f'{path}/6_Analysis/STD_calibration/2_STD_Data_Equation/{file_type}_STD_predict_error_df.csv',sep = '\t', index_col=0)

    polyReg_MAE = globals()[f'{file_type}_polyReg_MAE']
    polyReg_equation = globals()[f'{file_type}_polyReg_equation']

    mixData = mRNA_summaryDF
    pred_df = pred_DF
    maxPeak_equation_a = pred_equation_coef_DF[[f'{file_type}']].loc['maxPeak_x'][0]
    maxPeak_equation_b = pred_equation_intercept_DF[[f'{file_type}']].loc['maxPeak_x'][0]

    pred_error_1 = pred_error_DF[['Max']]
    pred_error_diff_2 = float(polyReg_MAE)
    equation = polyReg_equation

    stat="MAE"
    maxLen = 150
    
    spike_prop_list = [i/100 for i in range(0,101,2)]

    for idx in mRNA_list[file_type]:
    ### Mix data
        row = mixData.loc[idx]

        designed_len = int(row['Length'])
        measured_value = round(row[stat], 1)

        measured_1stPeak_x, measured_1stPeak_y = row['primaryPeak_x'], row['primaryPeak_y']
        measured_2ndPeak_x, measured_2ndPeak_y = row['secondaryPeak_x'], row['secondaryPeak_y']    

    ### Predict data
        pred_row = pred_df.loc[pred_df['Length'] == designed_len]

        pred_1stPeak_x, pred_error_maxPeak_x = float(pred_row['maxPeak_x']), pred_error_1.loc['maxPeak_x'][0]
        pred_1stPeak_x_m, pred_1stPeak_x_M = errorRange(pred_1stPeak_x, pred_error_maxPeak_x, 5)

        pred_value, pred_error_value  = float(pred_row[stat]), pred_error_1.loc['MAE'][0]
        pred_value_m, pred_value_M = errorRange(pred_value, pred_error_value, 5)

    ### Other
        diff_value_m = (measured_value - pred_value_M - pred_error_diff_2)
        diff_value_M = (measured_value - pred_value_m + pred_error_diff_2)

        measured_1stLength = Length_fromPeak_x(measured_1stPeak_x, maxPeak_equation_a, maxPeak_equation_b)
        measured_2ndLength = Length_fromPeak_x(measured_2ndPeak_x, maxPeak_equation_a, maxPeak_equation_b)

    ### Check
        if (pred_value_m <= measured_value <= pred_value_M):
            step = "No Contamination (MAE)"
            dic = [measured_value, pred_value, round(pred_error_value,1)]

        else:
            if (pred_1stPeak_x_m <= measured_1stPeak_x <= pred_1stPeak_x_M):
                if (measured_2ndPeak_x > 0):
                    
                    measured_2ndLength_m = Length_fromPeak_x((measured_2ndPeak_x - pred_error_maxPeak_x), maxPeak_equation_a, maxPeak_equation_b)
                    measured_2ndLength_M = Length_fromPeak_x((measured_2ndPeak_x + pred_error_maxPeak_x), maxPeak_equation_a, maxPeak_equation_b)

                    spike_len_m, spike_len_M = int(round(measured_2ndLength_m,0)), int(round(measured_2ndLength_M, 0))
                    if spike_len_m!=spike_len_M: spike_len_list = [i for i in range(spike_len_m, spike_len_M+1)]
                    else: spike_len_list = [spike_len_m]

                    measured_Peak_ratio = int(measured_2ndPeak_y*100//measured_1stPeak_y)
                    spike_prop_m, spike_prop_M = measured_Peak_ratio-5, measured_Peak_ratio+5
                    spike_prop_list = [i/100 for i in range(spike_prop_m,spike_prop_M+1,1)]

                    dic = evalution_mRNA_polyA(designed_len, spike_len_list, spike_prop_list, diff_value_m, diff_value_M, equation)
                    step = "Step 1"

                else:
                    step = "Step 2"
                    dic = {}

            else:
                if (measured_2ndPeak_x > 0):
                    if (pred_1stPeak_x_m <= measured_2ndPeak_x <= pred_1stPeak_x_M):
                        measured_1stLength_m = Length_fromPeak_x((measured_1stPeak_x - pred_error_maxPeak_x), maxPeak_equation_a, maxPeak_equation_b)
                        measured_1stLength_M = Length_fromPeak_x((measured_1stPeak_x + pred_error_maxPeak_x), maxPeak_equation_a, maxPeak_equation_b)

                        spike_len_m, spike_len_M = int(round(measured_1stLength_m,0)), int(round(measured_1stLength_M,0))
                        if spike_len_m!=spike_len_M: spike_len_list = [i for i in range(spike_len_m, spike_len_M+1)]
                        else: spike_len_list = [spike_len_m]

                        measured_Peak_ratio = int(measured_2ndPeak_y*100//measured_1stPeak_y)
                        spike_prop_m, spike_prop_M = measured_Peak_ratio-5, measured_Peak_ratio+5
                        spike_prop_list = sorted([1.0 - i/100 for i in range(spike_prop_m,spike_prop_M+1,1)])

                        dic = evalution_mRNA_polyA(designed_len, spike_len_list, spike_prop_list, diff_value_m, diff_value_M, equation)
                        step = "Step 3"

                    else:
                        measured_1stLength_m = Length_fromPeak_x((measured_1stPeak_x - pred_error_maxPeak_x), maxPeak_equation_a, maxPeak_equation_b)
                        measured_1stLength_M = Length_fromPeak_x((measured_1stPeak_x + pred_error_maxPeak_x), maxPeak_equation_a, maxPeak_equation_b)

                        measured_2ndLength_m = Length_fromPeak_x((measured_2ndPeak_x - pred_error_maxPeak_x), maxPeak_equation_a, maxPeak_equation_b)
                        measured_2ndLength_M = Length_fromPeak_x((measured_2ndPeak_x + pred_error_maxPeak_x), maxPeak_equation_a, maxPeak_equation_b)

                        dominant_len_m, dominant_len_M = int(round(measured_1stLength_m,0)), int(round(measured_1stLength_M, 0))
                        next_len_m, next_len_M = int(round(measured_2ndLength_m,0)), int(round(measured_2ndLength_M, 0))
                        measured_Peak_ratio = int(measured_2ndPeak_y*100//measured_1stPeak_y)
                        prop_m, prop_M = measured_Peak_ratio-5, measured_Peak_ratio+5
                        dominant_prop_list = [i/100 for i in range(prop_m,prop_M+1,1)]
                        next_prop_list = [1-i/100 for i in range(prop_m,prop_M+1,1)]

                        if abs(designed_len - measured_1stLength_m) <= abs(designed_len - measured_2ndLength_m):
                            ref_len_m, ref_len_M = dominant_len_m, dominant_len_M
                            spike_len_m, spike_len_M = next_len_m, next_len_M
                            spike_prop_list = next_prop_list

                        else:
                            ref_len_m, ref_len_M = next_len_m, next_len_M
                            spike_len_m, spike_len_M = dominant_len_m, dominant_len_M
                            spike_prop_list = dominant_prop_list

                        if spike_len_m!=spike_len_M: spike_len_list = [i for i in range(spike_len_m, spike_len_M+1)]
                        else: spike_len_list = [spike_len_m]
                        if ref_len_m!=ref_len_M: ref_len_list = [i for i in range(ref_len_m, ref_len_M+1)]
                        else: ref_len_list = [ref_len_m]

                        dic = evalution_mRNA_polyA(ref_len_list, spike_len_list, spike_prop_list, diff_value_m, diff_value_M, equation)
                        step = "Step 4"


        if len(dic) < 1:
            comment = "Many similar length"
            dic = [comment, designed_len, measured_value, pred_value, round(pred_error_value,2)]

        result_dic[file_type][idx] = [step, designed_len, dic]


# In[264]:


result_dic

