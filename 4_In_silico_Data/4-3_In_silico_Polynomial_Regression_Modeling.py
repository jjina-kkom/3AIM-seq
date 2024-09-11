#!/usr/bin/env python
# coding: utf-8

# In[45]:


import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import PolynomialFeatures
from sklearn.linear_model import LinearRegression
from sklearn.pipeline import make_pipeline
from sklearn.metrics import mean_squared_error
import matplotlib.pyplot as plt
from sklearn.metrics import r2_score
from alive_progress import alive_bar


# In[46]:


path='/BiO2/research/collab2/Catholic/JHNam/polyA_seq/240622'
output=f'{path}/6_Analysis/STD_calibration/3_STD_mixData_PolyReg'


# ### [1] Poly(A) without Linker

# In[47]:


file_type='withoutL'
mixDF = pd.read_csv(f'{path}/6_Analysis/STD_mixData/{file_type}_mixData_DF.csv', sep = '\t', index_col=0)

pred_DF = pd.read_csv(f'{path}/6_Analysis/STD_calibration/2_STD_Data_Equation/{file_type}_STD_predict_df.csv',sep = '\t', index_col=0)
pred_equation_coef_DF = pd.read_csv(f'{path}/6_Analysis/STD_calibration/2_STD_Data_Equation/STD_predict_equation_coef_df.csv',sep = '\t', index_col=0)
pred_equation_intercept_DF = pd.read_csv(f'{path}/6_Analysis/STD_calibration/2_STD_Data_Equation/STD_predict_equation_intercept_df.csv',sep = '\t', index_col=0)
pred_error_DF = pd.read_csv(f'{path}/6_Analysis/STD_calibration/2_STD_Data_Equation/{file_type}_STD_predict_error_df.csv',sep = '\t', index_col=0)


# In[119]:


df = mixDF
pred_df = pred_DF

primaryPeak_x_list, secondaryPeak_x_list = [], []

diff_len_list = []
diff_mae_list = []

with alive_bar(df.shape[0], force_tty=True) as bar:
    for idx in df.index:
    
    ################
    ### Max Peak ###
    ################
        maxPeak_dic = eval(df.loc[idx]['maxPeak'])
    
        primaryPeak_x = list(maxPeak_dic.keys())[0]

        if len(maxPeak_dic) > 1: secondaryPeak_x = list(maxPeak_dic.keys())[1]
        else: secondaryPeak_x = np.nan

        primaryPeak_x_list.append(np.ceil(primaryPeak_x))
        secondaryPeak_x_list.append(np.ceil(secondaryPeak_x))
        
    ############
    ### Diff ###
    ############
        row = df.loc[idx]
        
        ref_len, spike_len, mae = int(row['Ref']), int(row['Spike']), round(row['MAE'], 2)
        diff_len = spike_len - ref_len
        diff_len_list.append(diff_len)
        
        pred_mae = float(pred_df.loc[pred_df['Length']==ref_len]['MAE'])
        
        diff_mae = mae - pred_mae
        diff_mae_list.append(diff_mae)
        
        bar()


# In[120]:


input_mixData_DF = mixDF.loc[:, ['Ref', 'Spike', 'Spike%', 'MAE']]
input_mixData_DF.loc[:, 'Diff_Length'] = diff_len_list
input_mixData_DF.loc[:, 'Diff_MAE'] = diff_mae_list
input_mixData_DF.loc[:, 'primaryPeak_x'] = primaryPeak_x_list
input_mixData_DF.loc[:, 'secondaryPeak_x'] = secondaryPeak_x_list


# In[122]:


input_df = input_mixData_DF
X = input_df.loc[:, ['Ref', 'Diff_Length', 'Spike%']].rename(columns={'Ref':'A', 'Diff_Length':'B', 'Spike%':'C', })
y = input_df.loc[:, ['Diff_MAE']].rename(columns={'Diff_MAE':'y'})
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.1, random_state=1209)

# Polynomial regression model
degree = 9
polyreg_model = make_pipeline(PolynomialFeatures(degree), LinearRegression())
polyreg_model.fit(X_train, y_train)

y_pred = polyreg_model.predict(X_test)
mse = mean_squared_error(y_test, y_pred)
rmse = np.sqrt(mse)
y_pred_list = [i[0] for i in y_pred]
mae = sum(abs(np.array(y_test['y']) - np.array(y_pred_list)))/len(y_pred_list)
print(f"RMSE: {rmse}")
print(f"MAE: {mae}")

feature_names = PolynomialFeatures(degree=degree, include_bias=False).fit(X_train).get_feature_names_out(['A','B','C'])

coefficients = polyreg_model.named_steps['linearregression'].coef_
intercept = polyreg_model.named_steps['linearregression'].intercept_

polynomial_equation = f"{intercept[0]}" 
for coef, name in zip(coefficients[0][1:], feature_names):
    polynomial_equation += f"+({coef})*{name}" 

for i in range(degree,1,-1):
    polynomial_equation = polynomial_equation.replace(' ', '*')
    polynomial_equation = polynomial_equation.replace(f'A^{i}',f'(A^{i})').replace(f'B^{i}',f'(B^{i})').replace(f'C^{i}',f'(C^{i})')
    
polynomial_equation = polynomial_equation.replace('^','**')


# In[123]:


polyReg_RMSE, polyReg_MAE = rmse, mae
polyReg_equation = polynomial_equation

polyReg_list = [polyReg_RMSE, polyReg_MAE, degree, polyReg_equation]
polyReg_df = pd.DataFrame(polyReg_list).rename(index={0:'RMSE', 1:'MAE', 2:'Degree', 3:'Equation'}).rename(columns={0:'polyReg'})
polyReg_df.to_csv(f'{output}/{file_type}_polyReg_info_df.csv', sep='\t')
polyReg_df


# ### [2-1] 1st Poly(A) with Linker

# In[124]:


file_type='withL_H1'
mixDF = pd.read_csv(f'{path}/6_Analysis/STD_mixData/{file_type}_mixData_DF.csv', sep = '\t', index_col=0)

pred_DF = pd.read_csv(f'{path}/6_Analysis/STD_calibration/2_STD_Data_Equation/{file_type}_STD_predict_df.csv',sep = '\t', index_col=0)
pred_equation_coef_DF = pd.read_csv(f'{path}/6_Analysis/STD_calibration/2_STD_Data_Equation/STD_predict_equation_coef_df.csv',sep = '\t', index_col=0)
pred_equation_intercept_DF = pd.read_csv(f'{path}/6_Analysis/STD_calibration/2_STD_Data_Equation/STD_predict_equation_intercept_df.csv',sep = '\t', index_col=0)
pred_error_DF = pd.read_csv(f'{path}/6_Analysis/STD_calibration/2_STD_Data_Equation/{file_type}_STD_predict_error_df.csv',sep = '\t', index_col=0)


# In[125]:


df = mixDF
pred_df = pred_DF

primaryPeak_x_list, secondaryPeak_x_list = [], []

diff_len_list = []
diff_mae_list = []

with alive_bar(df.shape[0], force_tty=True) as bar:
    for idx in df.index:
    
    ################
    ### Max Peak ###
    ################
        maxPeak_dic = eval(df.loc[idx]['maxPeak'])
    
        primaryPeak_x = list(maxPeak_dic.keys())[0]

        if len(maxPeak_dic) > 1: secondaryPeak_x = list(maxPeak_dic.keys())[1]
        else: secondaryPeak_x = np.nan

        primaryPeak_x_list.append(np.ceil(primaryPeak_x))
        secondaryPeak_x_list.append(np.ceil(secondaryPeak_x))
        
    ############
    ### Diff ###
    ############
        row = df.loc[idx]
        
        ref_len, spike_len, mae = int(row['Ref']), int(row['Spike']), round(row['MAE'], 2)
        diff_len = spike_len - ref_len
        diff_len_list.append(diff_len)
        
        pred_mae = float(pred_df.loc[pred_df['Length']==ref_len]['MAE'])
        
        diff_mae = mae - pred_mae
        diff_mae_list.append(diff_mae)
        
        bar()


# In[126]:


input_mixData_DF = mixDF.loc[:, ['Ref', 'Spike', 'Spike%', 'MAE']]
input_mixData_DF.loc[:, 'Diff_Length'] = diff_len_list
input_mixData_DF.loc[:, 'Diff_MAE'] = diff_mae_list
input_mixData_DF.loc[:, 'primaryPeak_x'] = primaryPeak_x_list
input_mixData_DF.loc[:, 'secondaryPeak_x'] = secondaryPeak_x_list


# In[128]:


input_df = input_mixData_DF
X = input_df.loc[:, ['Ref', 'Diff_Length', 'Spike%']].rename(columns={'Ref':'A', 'Diff_Length':'B', 'Spike%':'C', })
y = input_df.loc[:, ['Diff_MAE']].rename(columns={'Diff_MAE':'y'})
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.1, random_state=1209)

# Polynomial regression model
degree = 9
polyreg_model = make_pipeline(PolynomialFeatures(degree), LinearRegression())
polyreg_model.fit(X_train, y_train)

y_pred = polyreg_model.predict(X_test)
mse = mean_squared_error(y_test, y_pred)
rmse = np.sqrt(mse)
y_pred_list = [i[0] for i in y_pred]
mae = sum(abs(np.array(y_test['y']) - np.array(y_pred_list)))/len(y_pred_list)
print(f"RMSE: {rmse}")
print(f"MAE: {mae}")

feature_names = PolynomialFeatures(degree=degree, include_bias=False).fit(X_train).get_feature_names_out(['A','B','C'])

coefficients = polyreg_model.named_steps['linearregression'].coef_
intercept = polyreg_model.named_steps['linearregression'].intercept_

polynomial_equation = f"{intercept[0]}" 
for coef, name in zip(coefficients[0][1:], feature_names):
    polynomial_equation += f"+({coef})*{name}" 

for i in range(degree,1,-1):
    polynomial_equation = polynomial_equation.replace(' ', '*')
    polynomial_equation = polynomial_equation.replace(f'A^{i}',f'(A^{i})').replace(f'B^{i}',f'(B^{i})').replace(f'C^{i}',f'(C^{i})')
    
polynomial_equation = polynomial_equation.replace('^','**')


# In[129]:


polyReg_RMSE, polyReg_MAE = rmse, mae
polyReg_equation = polynomial_equation

polyReg_list = [polyReg_RMSE, polyReg_MAE, degree, polyReg_equation]
polyReg_df = pd.DataFrame(polyReg_list).rename(index={0:'RMSE', 1:'MAE', 2:'Degree', 3:'Equation'}).rename(columns={0:'polyReg'})
polyReg_df.to_csv(f'{output}/{file_type}_polyReg_info_df.csv', sep='\t')
polyReg_df


# ### [2-2] 2nd Poly(A) with Linker

# In[130]:


file_type='withL_H2'
mixDF = pd.read_csv(f'{path}/6_Analysis/STD_mixData/{file_type}_mixData_DF.csv', sep = '\t', index_col=0)

pred_DF = pd.read_csv(f'{path}/6_Analysis/STD_calibration/2_STD_Data_Equation/{file_type}_STD_predict_df.csv',sep = '\t', index_col=0)
pred_equation_coef_DF = pd.read_csv(f'{path}/6_Analysis/STD_calibration/2_STD_Data_Equation/STD_predict_equation_coef_df.csv',sep = '\t', index_col=0)
pred_equation_intercept_DF = pd.read_csv(f'{path}/6_Analysis/STD_calibration/2_STD_Data_Equation/STD_predict_equation_intercept_df.csv',sep = '\t', index_col=0)
pred_error_DF = pd.read_csv(f'{path}/6_Analysis/STD_calibration/2_STD_Data_Equation/{file_type}_STD_predict_error_df.csv',sep = '\t', index_col=0)


# In[131]:


df = mixDF
pred_df = pred_DF

primaryPeak_x_list, secondaryPeak_x_list = [], []

diff_len_list = []
diff_mae_list = []

with alive_bar(df.shape[0], force_tty=True) as bar:
    for idx in df.index:
    
    ################
    ### Max Peak ###
    ################
        maxPeak_dic = eval(df.loc[idx]['maxPeak'])
    
        primaryPeak_x = list(maxPeak_dic.keys())[0]

        if len(maxPeak_dic) > 1: secondaryPeak_x = list(maxPeak_dic.keys())[1]
        else: secondaryPeak_x = np.nan

        primaryPeak_x_list.append(np.ceil(primaryPeak_x))
        secondaryPeak_x_list.append(np.ceil(secondaryPeak_x))
        
    ############
    ### Diff ###
    ############
        row = df.loc[idx]
        
        ref_len, spike_len, mae = int(row['Ref']), int(row['Spike']), round(row['MAE'], 2)
        diff_len = spike_len - ref_len
        diff_len_list.append(diff_len)
        
        pred_mae = float(pred_df.loc[pred_df['Length']==ref_len]['MAE'])
        
        diff_mae = mae - pred_mae
        diff_mae_list.append(diff_mae)
        
        bar()


# In[132]:


input_mixData_DF = mixDF.loc[:, ['Ref', 'Spike', 'Spike%', 'MAE']]
input_mixData_DF.loc[:, 'Diff_Length'] = diff_len_list
input_mixData_DF.loc[:, 'Diff_MAE'] = diff_mae_list
input_mixData_DF.loc[:, 'primaryPeak_x'] = primaryPeak_x_list
input_mixData_DF.loc[:, 'secondaryPeak_x'] = secondaryPeak_x_list


# In[134]:


input_df = input_mixData_DF
X = input_df.loc[:, ['Ref', 'Diff_Length', 'Spike%']].rename(columns={'Ref':'A', 'Diff_Length':'B', 'Spike%':'C', })
y = input_df.loc[:, ['Diff_MAE']].rename(columns={'Diff_MAE':'y'})
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.1, random_state=1209)

# Polynomial regression model
degree = 10
polyreg_model = make_pipeline(PolynomialFeatures(degree), LinearRegression())
polyreg_model.fit(X_train, y_train)

y_pred = polyreg_model.predict(X_test)
mse = mean_squared_error(y_test, y_pred)
rmse = np.sqrt(mse)
y_pred_list = [i[0] for i in y_pred]
mae = sum(abs(np.array(y_test['y']) - np.array(y_pred_list)))/len(y_pred_list)
print(f"RMSE: {rmse}")
print(f"MAE: {mae}")

feature_names = PolynomialFeatures(degree=degree, include_bias=False).fit(X_train).get_feature_names_out(['A','B','C'])

coefficients = polyreg_model.named_steps['linearregression'].coef_
intercept = polyreg_model.named_steps['linearregression'].intercept_

polynomial_equation = f"{intercept[0]}" 
for coef, name in zip(coefficients[0][1:], feature_names):
    polynomial_equation += f"+({coef})*{name}" 

for i in range(degree,1,-1):
    polynomial_equation = polynomial_equation.replace(' ', '*')
    polynomial_equation = polynomial_equation.replace(f'A^{i}',f'(A^{i})').replace(f'B^{i}',f'(B^{i})').replace(f'C^{i}',f'(C^{i})')
    
polynomial_equation = polynomial_equation.replace('^','**')


# In[135]:


polyReg_RMSE, polyReg_MAE = rmse, mae
polyReg_equation = polynomial_equation

polyReg_list = [polyReg_RMSE, polyReg_MAE, degree, polyReg_equation]
polyReg_df = pd.DataFrame(polyReg_list).rename(index={0:'RMSE', 1:'MAE', 2:'Degree', 3:'Equation'}).rename(columns={0:'polyReg'})
polyReg_df.to_csv(f'{output}/{file_type}_polyReg_info_df.csv', sep='\t')
polyReg_df

