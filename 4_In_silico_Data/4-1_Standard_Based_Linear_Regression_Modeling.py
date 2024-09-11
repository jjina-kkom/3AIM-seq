#!/usr/bin/env python
# coding: utf-8

# In[19]:


import pandas as pd
import numpy as np
from scipy.ndimage import gaussian_filter1d
from alive_progress import alive_bar
from scipy.stats import pearsonr
from sklearn.linear_model import LinearRegression
from sklearn.metrics import mean_squared_error, r2_score


# In[20]:


path='/BiO2/research/collab2/Catholic/JHNam/polyA_seq/240622'
withoutL_ref_list = [40,50,60,110,120,130]
withL_H1_ref_list = [[25, 10, 65], [30, 10, 70], [35, 10, 75], [45, 10, 45], [50, 10, 50], [55, 10, 55], [70, 10, 30]]
withL_H2_ref_list = [[70, 10, 30], [45, 10, 45], [50, 10, 50], [55, 10, 55], [25, 10, 65], [30, 10, 70], [35, 10, 75]]


# # 1. Calculated MAE & Peak Position

# In[21]:


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


# ### [1] Poly(A) without Linker

# In[22]:


Final_data = []

for ref_len in withoutL_ref_list:
    cutoff, region = 2, 'Homopolymer'
    withoutL_mixDF = pd.read_csv(f'{path}/4_filtering_data/STD_A{ref_len}/STD_A{ref_len}_Final_polyA_lenDF.final_uniqUMI.csv', 
                                 sep='\t', index_col=0)
    results = [] 
    lenList = list(withoutL_mixDF[region])
    
    MAE_mix = mae_val(lenList, [ref_len] * len(lenList))
    Peak_dic = hist_maxPeak(lenList, 300, 100, cutoff)
    maxPeak = max(Peak_dic, key=Peak_dic.get)
    maxPeak_density = Peak_dic[maxPeak]
    hist_maxPeak_dic = {k: v for k, v in Peak_dic.items() if v >= maxPeak_density * 0.15}
    primaryPeak = list(hist_maxPeak_dic.keys())[0]
    primaryPeak_val = round(hist_maxPeak_dic[primaryPeak], 5)
    if len(hist_maxPeak_dic) > 1:
        secondaryPeak = list(hist_maxPeak_dic.keys())[0]
        secondaryPeak_val = round(hist_maxPeak_dic[secondaryPeak], 5)
    else:
        secondaryPeak, secondaryPeak_val = 0, 0

    Final_data.append([f'STD_A{ref_len}', ref_len, MAE_mix, primaryPeak, primaryPeak_val])

columns = ['SampleID', 'Length', 'MAE', 'maxPeak_x', 'maxPeak_val']
withoutL_STD_summary_df = pd.DataFrame(Final_data, columns=columns)
withoutL_STD_summary_df.set_index('SampleID', inplace=True)

withoutL_STD_summary_df['logMAE'] = list(np.log10(withoutL_STD_summary_df['MAE']))
withoutL_STD_summary_df


# ### [2-1] 1st Poly(A) with Linker

# In[23]:


Final_data = []

for ref_list in withL_H1_ref_list:
    ref_len = ref_list[0]
    cutoff, region = 3, 'Homopolymer_1'
    withL_mixDF = pd.read_csv(f'{path}/4_filtering_data/STD_A{ref_list[0]}L{ref_list[-1]}/STD_A{ref_list[0]}L{ref_list[-1]}_Final_polyA_lenDF.final_uniqUMI.csv', 
                                 sep='\t', index_col=0)
    results = [] 
    lenList = list(withL_mixDF[region])
    
    MAE_mix = mae_val(lenList, [ref_len] * len(lenList))
    Peak_dic = hist_maxPeak(lenList, 300, 100, cutoff)
    maxPeak = max(Peak_dic, key=Peak_dic.get)
    maxPeak_density = Peak_dic[maxPeak]
    hist_maxPeak_dic = {k: v for k, v in Peak_dic.items() if v >= maxPeak_density * 0.15}
    primaryPeak = list(hist_maxPeak_dic.keys())[0]
    primaryPeak_val = round(hist_maxPeak_dic[primaryPeak], 5)
    if len(hist_maxPeak_dic) > 1:
        secondaryPeak = list(hist_maxPeak_dic.keys())[1]
        secondaryPeak_val = round(hist_maxPeak_dic[secondaryPeak], 5)
    else:
        secondaryPeak, secondaryPeak_val = 0, 0
        
    Final_data.append([f'STD_A{ref_list[0]}L{ref_list[-1]}', ref_len, MAE_mix, primaryPeak, primaryPeak_val])

columns = ['SampleID', 'Length', 'MAE', 'maxPeak_x', 'maxPeak_val']
withL_H1_STD_summary_df = pd.DataFrame(Final_data, columns=columns)
withL_H1_STD_summary_df.set_index('SampleID', inplace=True)

withL_H1_STD_summary_df['logMAE'] = list(np.log10(withL_H1_STD_summary_df['MAE']))
withL_H1_STD_summary_df


# ### [2-2] 2nd Poly(A) with Linker

# In[24]:


Final_data = []

for ref_list in withL_H2_ref_list:
    ref_len = ref_list[-1]
    cutoff, region = 3, 'Homopolymer_2'
    withL_mixDF = pd.read_csv(f'{path}/4_filtering_data/STD_A{ref_list[0]}L{ref_list[-1]}/STD_A{ref_list[0]}L{ref_list[-1]}_Final_polyA_lenDF.final_uniqUMI.csv', 
                                 sep='\t', index_col=0)
    results = [] 
    lenList = list(withL_mixDF[region])
    
    MAE_mix = mae_val(lenList, [ref_len] * len(lenList))
    Peak_dic = hist_maxPeak(lenList, 300, 100, cutoff)
    maxPeak = max(Peak_dic, key=Peak_dic.get)
    maxPeak_density = Peak_dic[maxPeak]
    hist_maxPeak_dic = {k: v for k, v in Peak_dic.items() if v >= maxPeak_density * 0.01}
    primaryPeak = list(hist_maxPeak_dic.keys())[0]
    primaryPeak_val = round(hist_maxPeak_dic[primaryPeak], 5)
    if len(hist_maxPeak_dic) > 1:
        secondaryPeak = list(hist_maxPeak_dic.keys())[1]
        secondaryPeak_val = round(hist_maxPeak_dic[secondaryPeak], 5)
    else:
        secondaryPeak, secondaryPeak_val = 0, 0
        
    Final_data.append([f'STD_A{ref_list[0]}L{ref_list[-1]}', ref_len, MAE_mix, primaryPeak, primaryPeak_val])

columns = ['SampleID', 'Length', 'MAE', 'maxPeak_x', 'maxPeak_val']
withL_H2_STD_summary_df = pd.DataFrame(Final_data, columns=columns)
withL_H2_STD_summary_df.set_index('SampleID', inplace=True)

withL_H2_STD_summary_df['logMAE'] = list(np.log10(withL_H2_STD_summary_df['MAE']))
withL_H2_STD_summary_df


# # 2. Linear Regression Modeling

# In[25]:


def linearReg(x, y, log_trans):
    
    if log_trans:
        x, y = np.log10(list(x)), np.log10(list(y))
    else:
        x, y = np.array(list(x)), np.array(list(y))
    
    model = LinearRegression()
    model.fit(x.reshape(-1, 1), y)
    predicted_y = model.predict(x.reshape(-1, 1))
    
    mse = mean_squared_error(y, predicted_y)
    rmse = np.sqrt(mse)
    mae = sum(abs(y - predicted_y))/len(predicted_y)
    r2 = r2_score(y, predicted_y)
    correlation, _ = pearsonr(x, y)
    intercept = model.intercept_
    coefficient = model.coef_[0]
    
    Equation=f'{coefficient} * x + {intercept}'

    return {'a':coefficient, 'b':intercept, 'equation': Equation,
           'correlation': correlation,
           'RMSE':rmse, 'MAE':mae, 'R^2':r2,
           'predicted_y': predicted_y, 'residuals': y - predicted_y}


# In[26]:


withoutL_STD_predict_df = withoutL_STD_summary_df[['Length']]
withL_H1_STD_predict_df = withL_H1_STD_summary_df[['Length']]
withL_H2_STD_predict_df = withL_H2_STD_summary_df[['Length']]


# ### [1] For Primary Peak

# In[27]:


x1, y1 = withoutL_STD_summary_df['Length'], withoutL_STD_summary_df['maxPeak_x']
x2, y2 = withL_H1_STD_summary_df['Length'], withL_H1_STD_summary_df['maxPeak_x']
x3, y3 = withL_H2_STD_summary_df['Length'], withL_H2_STD_summary_df['maxPeak_x']

withoutL_maxPeakX = linearReg(x1, y1, False)
withL_H1_maxPeakX = linearReg(x2, y2, False)
withL_H2_maxPeakX = linearReg(x3, y3, False)

withoutL_maxPeakX_equation, withoutL_maxPeakX_a, withoutL_maxPeakX_b, withoutL_maxPeakX_r2 = withoutL_maxPeakX['equation'], withoutL_maxPeakX['a'], withoutL_maxPeakX['b'], withoutL_maxPeakX['R^2']
withL_H1_maxPeakX_equation, withL_H1_maxPeakX_a, withL_H1_maxPeakX_b, withL_H1_maxPeakX_r2 = withL_H1_maxPeakX['equation'], withL_H1_maxPeakX['a'], withL_H1_maxPeakX['b'], withL_H1_maxPeakX['R^2']
withL_H2_maxPeakX_equation, withL_H2_maxPeakX_a, withL_H2_maxPeakX_b, withL_H2_maxPeakX_r2 = withL_H2_maxPeakX['equation'], withL_H2_maxPeakX['a'], withL_H2_maxPeakX['b'], withL_H2_maxPeakX['R^2']

withoutL_maxPeakX_rmse, withoutL_maxPeakX_mae = withoutL_maxPeakX['RMSE'], withoutL_maxPeakX['MAE']
withL_H1_maxPeakX_rmse, withL_H1_maxPeakX_mae = withL_H1_maxPeakX['RMSE'], withL_H1_maxPeakX['MAE']
withL_H2_maxPeakX_rmse, withL_H2_maxPeakX_mae = withL_H2_maxPeakX['RMSE'], withL_H2_maxPeakX['MAE']

print(f"Equation: {withoutL_maxPeakX['equation']}")
print(f"Equation: {withL_H1_maxPeakX['equation']}")
print(f"Equation: {withL_H2_maxPeakX['equation']}")
print()
print(f"RMSE: {withoutL_maxPeakX['RMSE']}")
print(f"RMSE: {withL_H1_maxPeakX['RMSE']}")
print(f"RMSE: {withL_H2_maxPeakX['RMSE']}")


# In[28]:


####################################
### Poly(A) Without Linker
####################################
df = withoutL_STD_summary_df
pred_y_list, equation = [], withoutL_maxPeakX_equation
residuals_list = withoutL_maxPeakX['residuals']

for idx in df.index:
    row = df.loc[idx]
    x = row['Length']
    pred_y = eval(equation)
    pred_y_list.append(np.ceil(pred_y))
    
withoutL_STD_predict_df.loc[:,'pred_x_residual'] = residuals_list
withoutL_STD_predict_df.loc[:,'maxPeak_x'] = pred_y_list

####################################
### 1st Poly(A) With Linker
####################################
df = withL_H1_STD_summary_df
pred_y_list, equation = [], withL_H1_maxPeakX_equation
residuals_list = withL_H1_maxPeakX['residuals']

for idx in df.index:
    row = df.loc[idx]
    x = row['Length']
    pred_y = eval(equation)
    pred_y_list.append(np.ceil(pred_y))
    
withL_H1_STD_predict_df.loc[:,'pred_x_residual'] = residuals_list
withL_H1_STD_predict_df.loc[:,'maxPeak_x'] = pred_y_list

####################################
### 2nd Poly(A) With Linker
####################################
df = withL_H2_STD_summary_df
pred_y_list, equation = [], withL_H2_maxPeakX_equation
residuals_list = withL_H2_maxPeakX['residuals']

for idx in df.index:
    row = df.loc[idx]
    x = row['Length']
    pred_y = eval(equation)
    pred_y_list.append(np.ceil(pred_y))
    
withL_H2_STD_predict_df.loc[:,'pred_x_residual'] = residuals_list
withL_H2_STD_predict_df.loc[:,'maxPeak_x'] = pred_y_list


# ### [2] For MAE

# In[29]:


x1, y1 = withoutL_STD_summary_df['Length'], withoutL_STD_summary_df['MAE']
x2, y2 = withL_H1_STD_summary_df['Length'], withL_H1_STD_summary_df['MAE']
x3, y3 = withL_H2_STD_summary_df['Length'], withL_H2_STD_summary_df['MAE']

withoutL_MAE = linearReg(x1, y1, True)
withL_H1_MAE = linearReg(x2, y2, True)
withL_H2_MAE = linearReg(x3, y3, True)

withoutL_MAE_equation, withoutL_MAE_a, withoutL_MAE_b, withoutL_MAE_r2 = withoutL_MAE['equation'], withoutL_MAE['a'], withoutL_MAE['b'], withoutL_MAE['R^2']
withL_H1_MAE_equation, withL_H1_MAE_a, withL_H1_MAE_b, withL_H1_MAE_r2 = withL_H1_MAE['equation'], withL_H1_MAE['a'], withL_H1_MAE['b'], withL_H1_MAE['R^2']
withL_H2_MAE_equation, withL_H2_MAE_a, withL_H2_MAE_b, withL_H2_MAE_r2 = withL_H2_MAE['equation'], withL_H2_MAE['a'], withL_H2_MAE['b'], withL_H2_MAE['R^2']

withoutL_MAE_rmse, withoutL_MAE_mae = withoutL_MAE['RMSE'], withoutL_MAE['MAE']
withL_H1_MAE_rmse, withL_H1_MAE_mae = withL_H1_MAE['RMSE'], withL_H1_MAE['MAE']
withL_H2_MAE_rmse, withL_H2_MAE_mae = withL_H2_MAE['RMSE'], withL_H2_MAE['MAE']

print(f"correlation: {withoutL_MAE['correlation']}")
print(f"correlation: {withL_H1_MAE['correlation']}")
print(f"correlation: {withL_H2_MAE['correlation']}")
print()
print(f"RMSE: {withoutL_MAE['RMSE']}")
print(f"RMSE: {withL_H1_MAE['RMSE']}")
print(f"RMSE: {withL_H2_MAE['RMSE']}")


# In[30]:


####################################
### Poly(A) Without Linker
####################################
df = withoutL_STD_predict_df
pred_y_list, log_pred_y_list = [], []
for idx in df.index:
    row = df.loc[idx]
    x =  row['Length']
    
    log_pred_y_add = withoutL_MAE_a * np.log10(x) + withoutL_MAE_b
    log_pred_y_list.append(log_pred_y_add)
    pred_y_add = 10**(log_pred_y_add)
    pred_y_list.append(round(pred_y_add, 2))
    
withoutL_STD_predict_df.loc[:,'MAE'] = pred_y_list
withoutL_STD_predict_df.loc[:,'logMAE'] = log_pred_y_list

####################################
### 1st Poly(A) With Linker
####################################
df = withL_H1_STD_predict_df
pred_y_list, log_pred_y_list = [], []
for idx in df.index:
    row = df.loc[idx]
    x =  row['Length']

    log_pred_y_add = withL_H1_MAE_a * np.log10(x) + withL_H1_MAE_b
    log_pred_y_list.append(log_pred_y_add)
    pred_y_add = 10**(log_pred_y_add)
    pred_y_list.append(round(pred_y_add, 2))
    
withL_H1_STD_predict_df.loc[:,'MAE'] = pred_y_list
withL_H1_STD_predict_df.loc[:,'logMAE'] = log_pred_y_list

####################################
### 2nd Poly(A) With Linker
####################################
df = withL_H2_STD_predict_df
pred_y_list, log_pred_y_list = [], []
for idx in df.index:
    row = df.loc[idx]
    x =  row['Length']
    
    log_pred_y_add = withL_H2_MAE_a * np.log10(x) + withL_H2_MAE_b
    log_pred_y_list.append(log_pred_y_add)
    pred_y_add = 10**(log_pred_y_add)
    pred_y_list.append(round(pred_y_add, 2))
    
withL_H2_STD_predict_df.loc[:,'MAE'] = pred_y_list
withL_H2_STD_predict_df.loc[:,'logMAE'] = log_pred_y_list


# ### [3] Summary

# In[31]:


columns = ['Length', 'maxPeak_x', 'MAE', 'logMAE']
withoutL_STD_summary_df = withoutL_STD_summary_df.loc[:, columns]
withL_H1_STD_summary_df = withL_H1_STD_summary_df.loc[:, columns]
withL_H2_STD_summary_df = withL_H2_STD_summary_df.loc[:, columns]

columns = ['Length', 'maxPeak_x', 'MAE', 'logMAE']
withoutL_STD_predict_df = withoutL_STD_predict_df.loc[:, columns]
withL_H1_STD_predict_df = withL_H1_STD_predict_df.loc[:, columns]
withL_H2_STD_predict_df = withL_H2_STD_predict_df.loc[:, columns]


# In[32]:


equation_a_dict = {'maxPeak_x':[withoutL_maxPeakX_a, withL_H1_maxPeakX_a, withL_H2_maxPeakX_a], 
                   'MAE':[withoutL_MAE_a, withL_H1_MAE_a, withL_H2_MAE_a]}

equation_b_dict = {'maxPeak_x':[withoutL_maxPeakX_b, withL_H1_maxPeakX_b, withL_H2_maxPeakX_b], 
                   'MAE':[withoutL_MAE_b, withL_H1_MAE_b, withL_H2_MAE_b]}

equation_r2_dict = {'maxPeak_x':[withoutL_maxPeakX_r2, withL_H1_maxPeakX_r2, withL_H2_maxPeakX_r2], 
                   'MAE':[withoutL_MAE_r2, withL_H1_MAE_r2, withL_H2_MAE_r2]}

STD_predict_equation_coef_df = pd.DataFrame.from_dict(equation_a_dict, orient='index').rename(columns={0:'withoutL', 1:'withL_H1', 2:'withL_H2'})
STD_predict_equation_intercept_df = pd.DataFrame.from_dict(equation_b_dict, orient='index').rename(columns={0:'withoutL', 1:'withL_H1', 2:'withL_H2'})
STD_predict_equation_r2_df = pd.DataFrame.from_dict(equation_r2_dict, orient='index').rename(columns={0:'withoutL', 1:'withL_H1', 2:'withL_H2'})


# ### [4] Set Error Range for Predict Model

# In[33]:


withoutL_STD_diff_df = abs(withoutL_STD_predict_df-withoutL_STD_summary_df)
withL_H1_STD_diff_df = abs(withL_H1_STD_predict_df-withL_H1_STD_summary_df)
withL_H2_STD_diff_df = abs(withL_H2_STD_predict_df-withL_H2_STD_summary_df)


# In[34]:


withoutL_dict = {'maxPeak_x':[max(withoutL_STD_diff_df['maxPeak_x'])], 
                 'MAE':[max(withoutL_STD_diff_df['MAE'])],
                 'logMAE':[max(withoutL_STD_diff_df['logMAE'])]}

witL_H1_dict = {'maxPeak_x':[max(withL_H1_STD_diff_df['maxPeak_x'])], 
                 'MAE':[max(withL_H1_STD_diff_df['MAE'])],
                 'logMAE':[max(withL_H1_STD_diff_df['logMAE'])]}

witL_H2_dict = {'maxPeak_x':[max(withL_H2_STD_diff_df['maxPeak_x'])], 
                 'MAE':[max(withL_H2_STD_diff_df['MAE'])],
                 'logMAE':[max(withL_H2_STD_diff_df['logMAE'])]}

withoutL_STD_predict_error_df = pd.DataFrame.from_dict(withoutL_dict, orient='index').rename(columns={0:'Max'})
withL_H1_STD_predict_error_df = pd.DataFrame.from_dict(witL_H1_dict, orient='index').rename(columns={0:'Max'})
withL_H2_STD_predict_error_df = pd.DataFrame.from_dict(witL_H2_dict, orient='index').rename(columns={0:'Max'})


# ### [5] Save Data

# In[ ]:


output=f'{path}/6_Analysis/STD_calibration/2_STD_Data_Equation'

withoutL_STD_summary_df.to_csv(f'{output}/withoutL_STD_summary_df.csv', sep='\t')
withL_H1_STD_summary_df.to_csv(f'{output}/withL_H1_STD_summary_df.csv', sep='\t')
withL_H2_STD_summary_df.to_csv(f'{output}/withL_H2_STD_summary_df.csv', sep='\t')

withoutL_STD_predict_df.to_csv(f'{output}/withoutL_STD_predict_df.csv', sep='\t')
withL_H1_STD_predict_df.to_csv(f'{output}/withL_H1_STD_predict_df.csv', sep='\t')
withL_H2_STD_predict_df.to_csv(f'{output}/withL_H2_STD_predict_df.csv', sep='\t')

withoutL_STD_predict_error_df.to_csv(f'{output}/withoutL_STD_predict_error_df.csv', sep='\t')
withL_H1_STD_predict_error_df.to_csv(f'{output}/withL_H1_STD_predict_error_df.csv', sep='\t')
withL_H2_STD_predict_error_df.to_csv(f'{output}/withL_H2_STD_predict_error_df.csv', sep='\t')

STD_predict_equation_coef_df.to_csv(f'{output}/STD_predict_equation_coef_df.csv', sep='\t')
STD_predict_equation_intercept_df.to_csv(f'{output}/STD_predict_equation_intercept_df.csv', sep='\t')
STD_predict_equation_r2_df.to_csv(f'{output}/STD_predict_equation_r2_df.csv', sep='\t')

