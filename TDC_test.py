# -*- coding: utf-8 -*-
"""
Created on Fri Dec  9 07:54:11 2022

@author: JohnDoe
"""

import pandas as pd
import numpy as np
import time
from hierarchical_filtering import fs
import coracle
from tdc import tdc
directory = "/Users/JohnDoe/Desktop/data/datasets/"

#%%

######### CBASS dataset




### CBASS dataset
cbass_meta = pd.read_csv(directory + "CBASS_gradient_metafile.txt", sep= "\s+")
cbass_ASV = pd.read_csv(directory + "ASV_table_QCfiltered.txt", sep= "\s+")
cbass_ED50 = pd.read_csv(directory + "ED50_per_colony.txt", sep= "\s+", index_col=(0))
cbass_tax = cbass_ASV.iloc[:, -7:-1] #get taxonomic information
cbass_ASV = cbass_ASV.iloc[:, :-8]

#cbass_ASV = cbass_ASV.transpose()
cbass_ASV = cbass_ASV.loc[:, cbass_ASV.columns.str.contains('30')]
cbass_ASV.columns = cbass_ASV.columns.to_series().astype(str).str.replace(r'_30.+','')


#change meta
cbass_meta.set_index(["Sample"], inplace = True)
cbass_meta = cbass_meta.loc[cbass_meta.index.str.contains('30'), :]
cbass_meta.index = cbass_meta.index.to_series().astype(str).str.replace(r'_30.+','')
cbass_meta.drop(["Temperature", "Colony"], axis = 1,  inplace = True)


#merge
cbass_ASV = cbass_ASV.transpose()
cbass = cbass_ED50.merge(cbass_meta, left_index=True, right_index=True)
cbass = cbass.merge(cbass_ASV, left_index=True, right_index=True)

#split by species
acropora = cbass.loc[cbass.Species=='Acropora']
pocillopora = cbass.loc[cbass.Species=='Pocillopora']
porites = cbass.loc[cbass.Species=='Porites']
stylophora = cbass.loc[cbass.Species=='Stylophora']

acropora.drop(["Site", "Species"], axis = 1,  inplace = True)
pocillopora.drop(["Site", "Species"], axis = 1,  inplace = True)
porites.drop(["Site", "Species"], axis = 1,  inplace = True)
stylophora.drop(["Site", "Species"], axis = 1,  inplace = True)

acropora_ED50 = pd.DataFrame(acropora["ED50"])
acropora_ASV = acropora.loc[:, (acropora != 0).any(axis=0)]
acropora_ASV = acropora.iloc[:, 1:].transpose()


pocillopora_ED50 = pd.DataFrame(pocillopora["ED50"])
pocillopora_ASV = pocillopora.loc[:, (pocillopora != 0).any(axis=0)]
pocillopora_ASV = pocillopora.iloc[:, 1:].transpose()

porites_ED50 = pd.DataFrame(porites["ED50"])
porites_ASV = porites.loc[:, (porites != 0).any(axis=0)]
porites_ASV = porites.iloc[:, 1:].transpose()

stylophora_ED50 = pd.DataFrame(stylophora["ED50"])
stylophora_ASV = stylophora.loc[:, (stylophora != 0).any(axis=0)]
stylophora_ASV = stylophora.iloc[:, 1:].transpose()


#%%
#prepare tax:
values = cbass_tax.index.to_series()
for column in cbass_tax.columns:
    cbass_tax[column] = cbass_tax[column].fillna(values)


### test on cbass data with small threshold
#%%
start = time.time()
result_acropora = tdc(acropora_ASV, acropora_ED50, cbass_tax, threshold=0.2)
print(time.time()- start)
#%%
start = time.time()
result_pocillopora = td_coracle(pocillopora_ASV, pocillopora_ED50, cbass_tax, threshold=0.2)
print(time.time()- start)
#%%
start = time.time()
result_porites = td_coracle(porites_ASV, porites_ED50, cbass_tax, threshold=0.2)
print(time.time()- start)
#%%
start = time.time()
result_stylophora = td_coracle(stylophora_ASV, stylophora_ED50, cbass_tax, threshold=0.2)
print(time.time()- start)


#%%
result_acropora.to_csv("C:/Users/JohnDoe/Desktop/data/datasets/test_result/result_tdc_acropora.csv", float_format="%.4g")
result_pocillopora.to_csv("C:/Users/JohnDoe/Desktop/data/datasets/test_result/result_tdc_pocillopora.csv", float_format="%.4g")
result_porites.to_csv("C:/Users/JohnDoe/Desktop/data/datasets/test_result/result_tdc_porites.csv", float_format="%.4g")
result_stylophora.to_csv("C:/Users/JohnDoe/Desktop/data/datasets/test_result/result_tdc_stylophora.csv", float_format="%.4g")

#%%

from td_coracle import td_coracle_wo_uc

### test on cbass data with small threshold
#%%
start = time.time()
result_acropora_wo_uc = td_coracle_wo_uc(acropora_ASV, acropora_ED50, cbass_tax, threshold=0.2)
print(time.time()- start)
#%%
start = time.time()
result_pocillopora_wo_uc = td_coracle_wo_uc(pocillopora_ASV, pocillopora_ED50, cbass_tax, threshold=0.2)
print(time.time()- start)
#%%
start = time.time()
result_porites_wo_uc = td_coracle_wo_uc(porites_ASV, porites_ED50, cbass_tax, threshold=0.2)
print(time.time()- start)
#%%
start = time.time()
result_stylophora_wo_uc = td_coracle_wo_uc(stylophora_ASV, stylophora_ED50, cbass_tax, threshold=0.2)
print(time.time()- start)


#%%
result_acropora_wo_uc.to_csv("C:/Users/JohnDoe/Desktop/data/datasets/test_result/result_tdc_acropora_wo_uc.csv", float_format="%.4g")
result_pocillopora_wo_uc.to_csv("C:/Users/JohnDoe/Desktop/data/datasets/test_result/result_tdc_pocillopora_wo_uc.csv", float_format="%.4g")
result_porites_wo_uc.to_csv("C:/Users/JohnDoe/Desktop/data/datasets/test_result/result_tdc_porites_wo_uc.csv", float_format="%.4g")
result_stylophora_wo_uc.to_csv("C:/Users/JohnDoe/Desktop/data/datasets/test_result/result_tdc_stylophora_wo_uc.csv", float_format="%.4g")