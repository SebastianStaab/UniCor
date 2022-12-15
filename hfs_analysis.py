# -*- coding: utf-8 -*-
"""
Created on Thu Sep 29 11:22:59 2022

@author: JohnDoe
"""
import pandas as pd
import numpy as np
import time
from hierarchical_filtering import fs

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

#%%
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

acropora_filtered_ASV = fs(acropora_ASV, acropora_ED50, cbass_tax)
pocillopora_filtered_ASV = fs(pocillopora_ASV, pocillopora_ED50, cbass_tax)
porites_filtered_ASV = fs(porites_ASV, porites_ED50, cbass_tax)
stylophora_filtered_ASV = fs(stylophora_ASV, stylophora_ED50, cbass_tax)

acropora_ASV = acropora_ASV.merge(cbass_tax, left_index=True, right_index=True)
pocillopora_ASV = pocillopora_ASV.merge(cbass_tax, left_index=True, right_index=True)
porites_ASV = porites_ASV.merge(cbass_tax, left_index=True, right_index=True)
stylophora_ASV = stylophora_ASV.merge(cbass_tax, left_index=True, right_index=True)

#%%

#Acropora
#change ASV
acropora_ASV_genus = acropora_ASV.groupby(["Genus"]).sum().transpose()
acropora_ASV_family = acropora_ASV.groupby(["Family"]).sum().transpose()
acropora_ASV_order = acropora_ASV.groupby(["Order"]).sum().transpose()
acropora_ASV_class = acropora_ASV.groupby(["Class"]).sum().transpose()
acropora_ASV_phylum = acropora_ASV.groupby(["Phylum"]).sum().transpose()

acropora_ASV_genus_f = acropora_filtered_ASV.groupby(["Genus"]).sum().transpose()
acropora_ASV_family_f = acropora_filtered_ASV.groupby(["Family"]).sum().transpose()
acropora_ASV_order_f = acropora_filtered_ASV.groupby(["Order"]).sum().transpose()
acropora_ASV_class_f = acropora_filtered_ASV.groupby(["Class"]).sum().transpose()
acropora_ASV_phylum_f = acropora_filtered_ASV.groupby(["Phylum"]).sum().transpose()

#merge
acropora_genus = acropora_ED50.merge(acropora_ASV_genus, left_index=True, right_index=True)
acropora_family = acropora_ED50.merge(acropora_ASV_family, left_index=True, right_index=True)
acropora_order = acropora_ED50.merge(acropora_ASV_order, left_index=True, right_index=True)
acropora_class = acropora_ED50.merge(acropora_ASV_class, left_index=True, right_index=True)
acropora_phylum = acropora_ED50.merge(acropora_ASV_phylum, left_index=True, right_index=True)

acropora_genus_f = acropora_ED50.merge(acropora_ASV_genus_f, left_index=True, right_index=True)
acropora_family_f = acropora_ED50.merge(acropora_ASV_family_f, left_index=True, right_index=True)
acropora_order_f = acropora_ED50.merge(acropora_ASV_order_f, left_index=True, right_index=True)
acropora_class_f = acropora_ED50.merge(acropora_ASV_class_f, left_index=True, right_index=True)
acropora_phylum_f = acropora_ED50.merge(acropora_ASV_phylum_f, left_index=True, right_index=True)

#%%
import coracle

x = acropora_genus.iloc[:,3:]
y = acropora_genus[['ED50']]
start = time.time()
result_g_acropora = coracle.coracle(x, y, alpha_l1 = 10**(-2.9), alpha_clr = 10**(-0.4))
print(time.time()- start)

x_f = acropora_genus_f.iloc[:,3:]
y_f = acropora_genus_f[['ED50']]
start = time.time()
result_g_acropora_f = coracle.coracle(x_f, y_f, alpha_l1 = 10**(-2.9), alpha_clr = 10**(-0.4))
print(time.time()- start)

x = acropora_family.iloc[:,3:]
y = acropora_family[['ED50']]
start = time.time()
result_f_acropora = coracle.coracle(x, y, alpha_l1 = 10**(-2.9), alpha_clr = 10**(-0.4))
print(time.time()- start)

x_f = acropora_family_f.iloc[:,3:]
y_f = acropora_family_f[['ED50']]
start = time.time()
result_f_acropora_f = coracle.coracle(x_f, y_f, alpha_l1 = 10**(-2.9), alpha_clr = 10**(-0.4))
print(time.time()- start)

x = acropora_order.iloc[:,3:]
y = acropora_order[['ED50']]
start = time.time()
result_o_acropora = coracle.coracle(x, y, alpha_l1 = 10**(-2.9), alpha_clr = 10**(-0.4))
print(time.time()- start)

x_f = acropora_order_f.iloc[:,3:]
y_f = acropora_order_f[['ED50']]
start = time.time()
result_o_acropora_f = coracle.coracle(x_f, y_f, alpha_l1 = 10**(-2.9), alpha_clr = 10**(-0.4))
print(time.time()- start)

x = acropora_class.iloc[:,3:]
y = acropora_class[['ED50']]
start = time.time()
result_c_acropora = coracle.coracle(x, y, alpha_l1 = 10**(-2.9), alpha_clr = 10**(-0.4))
print(time.time()- start)

x_f = acropora_class_f.iloc[:,3:]
y_f = acropora_class_f[['ED50']]
start = time.time()
result_c_acropora_f = coracle.coracle(x_f, y_f, alpha_l1 = 10**(-2.9), alpha_clr = 10**(-0.4))
print(time.time()- start)

x = acropora_phylum.iloc[:,3:]
y = acropora_phylum[['ED50']]
start = time.time()
result_p_acropora = coracle.coracle(x, y, alpha_l1 = 10**(-2.9), alpha_clr = 10**(-0.4))
print(time.time()- start)

x_f = acropora_phylum_f.iloc[:,3:]
y_f = acropora_phylum_f[['ED50']]
start = time.time()
result_p_acropora_f = coracle.coracle(x_f, y_f, alpha_l1 = 10**(-2.9), alpha_clr = 10**(-0.4))
print(time.time()- start)

#%%
result_g_acropora.to_csv("C:/Users/JohnDoe/Desktop/data/datasets/test_result/result_g_acropora.csv", float_format="%.4g")
result_g_acropora_f.to_csv("C:/Users/JohnDoe/Desktop/data/datasets/test_result/result_g_acropora_f.csv", float_format="%.4g")
result_f_acropora.to_csv("C:/Users/JohnDoe/Desktop/data/datasets/test_result/result_f_acropora.csv", float_format="%.4g")
result_f_acropora_f.to_csv("C:/Users/JohnDoe/Desktop/data/datasets/test_result/result_f_acropora_f.csv", float_format="%.4g")
result_o_acropora.to_csv("C:/Users/JohnDoe/Desktop/data/datasets/test_result/result_o_acropora.csv", float_format="%.4g")
result_o_acropora_f.to_csv("C:/Users/JohnDoe/Desktop/data/datasets/test_result/result_o_acropora_f.csv", float_format="%.4g")
result_c_acropora.to_csv("C:/Users/JohnDoe/Desktop/data/datasets/test_result/result_c_acropora.csv", float_format="%.4g")
result_c_acropora_f.to_csv("C:/Users/JohnDoe/Desktop/data/datasets/test_result/result_c_acropora_f.csv", float_format="%.4g")
result_p_acropora.to_csv("C:/Users/JohnDoe/Desktop/data/datasets/test_result/result_p_acropora.csv", float_format="%.4g")
result_p_acropora_f.to_csv("C:/Users/JohnDoe/Desktop/data/datasets/test_result/result_p_acropora_f.csv", float_format="%.4g")


#%%

#pocillopora
#change ASV
pocillopora_ASV_genus = pocillopora_ASV.groupby(["Genus"]).sum().transpose()
pocillopora_ASV_family = pocillopora_ASV.groupby(["Family"]).sum().transpose()
pocillopora_ASV_order = pocillopora_ASV.groupby(["Order"]).sum().transpose()
pocillopora_ASV_class = pocillopora_ASV.groupby(["Class"]).sum().transpose()
pocillopora_ASV_phylum = pocillopora_ASV.groupby(["Phylum"]).sum().transpose()

pocillopora_ASV_genus_f = pocillopora_filtered_ASV.groupby(["Genus"]).sum().transpose()
pocillopora_ASV_family_f = pocillopora_filtered_ASV.groupby(["Family"]).sum().transpose()
pocillopora_ASV_order_f = pocillopora_filtered_ASV.groupby(["Order"]).sum().transpose()
pocillopora_ASV_class_f = pocillopora_filtered_ASV.groupby(["Class"]).sum().transpose()
pocillopora_ASV_phylum_f = pocillopora_filtered_ASV.groupby(["Phylum"]).sum().transpose()

#merge
pocillopora_genus = pocillopora_ED50.merge(pocillopora_ASV_genus, left_index=True, right_index=True)
pocillopora_family = pocillopora_ED50.merge(pocillopora_ASV_family, left_index=True, right_index=True)
pocillopora_order = pocillopora_ED50.merge(pocillopora_ASV_order, left_index=True, right_index=True)
pocillopora_class = pocillopora_ED50.merge(pocillopora_ASV_class, left_index=True, right_index=True)
pocillopora_phylum = pocillopora_ED50.merge(pocillopora_ASV_phylum, left_index=True, right_index=True)

pocillopora_genus_f = pocillopora_ED50.merge(pocillopora_ASV_genus_f, left_index=True, right_index=True)
pocillopora_family_f = pocillopora_ED50.merge(pocillopora_ASV_family_f, left_index=True, right_index=True)
pocillopora_order_f = pocillopora_ED50.merge(pocillopora_ASV_order_f, left_index=True, right_index=True)
pocillopora_class_f = pocillopora_ED50.merge(pocillopora_ASV_class_f, left_index=True, right_index=True)
pocillopora_phylum_f = pocillopora_ED50.merge(pocillopora_ASV_phylum_f, left_index=True, right_index=True)

#%%
import coracle

x = pocillopora_genus.iloc[:,3:]
y = pocillopora_genus[['ED50']]
start = time.time()
result_g_pocillopora = coracle.coracle(x, y, alpha_l1 = 10**(-2.9), alpha_clr = 10**(-0.4))
print(time.time()- start)

x_f = pocillopora_genus_f.iloc[:,3:]
y_f = pocillopora_genus_f[['ED50']]
start = time.time()
result_g_pocillopora_f = coracle.coracle(x_f, y_f, alpha_l1 = 10**(-2.9), alpha_clr = 10**(-0.4))
print(time.time()- start)

x = pocillopora_family.iloc[:,3:]
y = pocillopora_family[['ED50']]
start = time.time()
result_f_pocillopora = coracle.coracle(x, y, alpha_l1 = 10**(-2.9), alpha_clr = 10**(-0.4))
print(time.time()- start)

x_f = pocillopora_family_f.iloc[:,3:]
y_f = pocillopora_family_f[['ED50']]
start = time.time()
result_f_pocillopora_f = coracle.coracle(x_f, y_f, alpha_l1 = 10**(-2.9), alpha_clr = 10**(-0.4))
print(time.time()- start)

x = pocillopora_order.iloc[:,3:]
y = pocillopora_order[['ED50']]
start = time.time()
result_o_pocillopora = coracle.coracle(x, y, alpha_l1 = 10**(-2.9), alpha_clr = 10**(-0.4))
print(time.time()- start)

x_f = pocillopora_order_f.iloc[:,3:]
y_f = pocillopora_order_f[['ED50']]
start = time.time()
result_o_pocillopora_f = coracle.coracle(x_f, y_f, alpha_l1 = 10**(-2.9), alpha_clr = 10**(-0.4))
print(time.time()- start)

x = pocillopora_class.iloc[:,3:]
y = pocillopora_class[['ED50']]
start = time.time()
result_c_pocillopora = coracle.coracle(x, y, alpha_l1 = 10**(-2.9), alpha_clr = 10**(-0.4))
print(time.time()- start)

x_f = pocillopora_class_f.iloc[:,3:]
y_f = pocillopora_class_f[['ED50']]
start = time.time()
result_c_pocillopora_f = coracle.coracle(x_f, y_f, alpha_l1 = 10**(-2.9), alpha_clr = 10**(-0.4))
print(time.time()- start)

x = pocillopora_phylum.iloc[:,3:]
y = pocillopora_phylum[['ED50']]
start = time.time()
result_p_pocillopora = coracle.coracle(x, y, alpha_l1 = 10**(-2.9), alpha_clr = 10**(-0.4))
print(time.time()- start)

x_f = pocillopora_phylum_f.iloc[:,3:]
y_f = pocillopora_phylum_f[['ED50']]
start = time.time()
result_p_pocillopora_f = coracle.coracle(x_f, y_f, alpha_l1 = 10**(-2.9), alpha_clr = 10**(-0.4))
print(time.time()- start)

#%%
result_g_pocillopora.to_csv("C:/Users/JohnDoe/Desktop/data/datasets/test_result/result_g_pocillopora.csv", float_format="%.4g")
result_g_pocillopora_f.to_csv("C:/Users/JohnDoe/Desktop/data/datasets/test_result/result_g_pocillopora_f.csv", float_format="%.4g")
result_f_pocillopora.to_csv("C:/Users/JohnDoe/Desktop/data/datasets/test_result/result_f_pocillopora.csv", float_format="%.4g")
result_f_pocillopora_f.to_csv("C:/Users/JohnDoe/Desktop/data/datasets/test_result/result_f_pocillopora_f.csv", float_format="%.4g")
result_o_pocillopora.to_csv("C:/Users/JohnDoe/Desktop/data/datasets/test_result/result_o_pocillopora.csv", float_format="%.4g")
result_o_pocillopora_f.to_csv("C:/Users/JohnDoe/Desktop/data/datasets/test_result/result_o_pocillopora_f.csv", float_format="%.4g")
result_c_pocillopora.to_csv("C:/Users/JohnDoe/Desktop/data/datasets/test_result/result_c_pocillopora.csv", float_format="%.4g")
result_c_pocillopora_f.to_csv("C:/Users/JohnDoe/Desktop/data/datasets/test_result/result_c_pocillopora_f.csv", float_format="%.4g")
result_p_pocillopora.to_csv("C:/Users/JohnDoe/Desktop/data/datasets/test_result/result_p_pocillopora.csv", float_format="%.4g")
result_p_pocillopora_f.to_csv("C:/Users/JohnDoe/Desktop/data/datasets/test_result/result_p_pocillopora_f.csv", float_format="%.4g")

#%%

#porites
#change ASV
porites_ASV_genus = porites_ASV.groupby(["Genus"]).sum().transpose()
porites_ASV_family = porites_ASV.groupby(["Family"]).sum().transpose()
porites_ASV_order = porites_ASV.groupby(["Order"]).sum().transpose()
porites_ASV_class = porites_ASV.groupby(["Class"]).sum().transpose()
porites_ASV_phylum = porites_ASV.groupby(["Phylum"]).sum().transpose()

porites_ASV_genus_f = porites_filtered_ASV.groupby(["Genus"]).sum().transpose()
porites_ASV_family_f = porites_filtered_ASV.groupby(["Family"]).sum().transpose()
porites_ASV_order_f = porites_filtered_ASV.groupby(["Order"]).sum().transpose()
porites_ASV_class_f = porites_filtered_ASV.groupby(["Class"]).sum().transpose()
porites_ASV_phylum_f = porites_filtered_ASV.groupby(["Phylum"]).sum().transpose()

#merge
porites_genus = porites_ED50.merge(porites_ASV_genus, left_index=True, right_index=True)
porites_family = porites_ED50.merge(porites_ASV_family, left_index=True, right_index=True)
porites_order = porites_ED50.merge(porites_ASV_order, left_index=True, right_index=True)
porites_class = porites_ED50.merge(porites_ASV_class, left_index=True, right_index=True)
porites_phylum = porites_ED50.merge(porites_ASV_phylum, left_index=True, right_index=True)

porites_genus_f = porites_ED50.merge(porites_ASV_genus_f, left_index=True, right_index=True)
porites_family_f = porites_ED50.merge(porites_ASV_family_f, left_index=True, right_index=True)
porites_order_f = porites_ED50.merge(porites_ASV_order_f, left_index=True, right_index=True)
porites_class_f = porites_ED50.merge(porites_ASV_class_f, left_index=True, right_index=True)
porites_phylum_f = porites_ED50.merge(porites_ASV_phylum_f, left_index=True, right_index=True)

#%%
import coracle

x = porites_genus.iloc[:,3:]
y = porites_genus[['ED50']]
start = time.time()
result_g_porites = coracle.coracle(x, y, alpha_l1 = 10**(-2.9), alpha_clr = 10**(-0.4))
print(time.time()- start)

x_f = porites_genus_f.iloc[:,3:]
y_f = porites_genus_f[['ED50']]
start = time.time()
result_g_porites_f = coracle.coracle(x_f, y_f, alpha_l1 = 10**(-2.9), alpha_clr = 10**(-0.4))
print(time.time()- start)

x = porites_family.iloc[:,3:]
y = porites_family[['ED50']]
start = time.time()
result_f_porites = coracle.coracle(x, y, alpha_l1 = 10**(-2.9), alpha_clr = 10**(-0.4))
print(time.time()- start)

x_f = porites_family_f.iloc[:,3:]
y_f = porites_family_f[['ED50']]
start = time.time()
result_f_porites_f = coracle.coracle(x_f, y_f, alpha_l1 = 10**(-2.9), alpha_clr = 10**(-0.4))
print(time.time()- start)

x = porites_order.iloc[:,3:]
y = porites_order[['ED50']]
start = time.time()
result_o_porites = coracle.coracle(x, y, alpha_l1 = 10**(-2.9), alpha_clr = 10**(-0.4))
print(time.time()- start)

x_f = porites_order_f.iloc[:,3:]
y_f = porites_order_f[['ED50']]
start = time.time()
result_o_porites_f = coracle.coracle(x_f, y_f, alpha_l1 = 10**(-2.9), alpha_clr = 10**(-0.4))
print(time.time()- start)

x = porites_class.iloc[:,3:]
y = porites_class[['ED50']]
start = time.time()
result_c_porites = coracle.coracle(x, y, alpha_l1 = 10**(-2.9), alpha_clr = 10**(-0.4))
print(time.time()- start)

x_f = porites_class_f.iloc[:,3:]
y_f = porites_class_f[['ED50']]
start = time.time()
result_c_porites_f = coracle.coracle(x_f, y_f, alpha_l1 = 10**(-2.9), alpha_clr = 10**(-0.4))
print(time.time()- start)

x = porites_phylum.iloc[:,3:]
y = porites_phylum[['ED50']]
start = time.time()
result_p_porites = coracle.coracle(x, y, alpha_l1 = 10**(-2.9), alpha_clr = 10**(-0.4))
print(time.time()- start)

x_f = porites_phylum_f.iloc[:,3:]
y_f = porites_phylum_f[['ED50']]
start = time.time()
result_p_porites_f = coracle.coracle(x_f, y_f, alpha_l1 = 10**(-2.9), alpha_clr = 10**(-0.4))
print(time.time()- start)

#%%
result_g_porites.to_csv("C:/Users/JohnDoe/Desktop/data/datasets/test_result/result_g_porites.csv", float_format="%.4g")
result_g_porites_f.to_csv("C:/Users/JohnDoe/Desktop/data/datasets/test_result/result_g_porites_f.csv", float_format="%.4g")
result_f_porites.to_csv("C:/Users/JohnDoe/Desktop/data/datasets/test_result/result_f_porites.csv", float_format="%.4g")
result_f_porites_f.to_csv("C:/Users/JohnDoe/Desktop/data/datasets/test_result/result_f_porites_f.csv", float_format="%.4g")
result_o_porites.to_csv("C:/Users/JohnDoe/Desktop/data/datasets/test_result/result_o_porites.csv", float_format="%.4g")
result_o_porites_f.to_csv("C:/Users/JohnDoe/Desktop/data/datasets/test_result/result_o_porites_f.csv", float_format="%.4g")
result_c_porites.to_csv("C:/Users/JohnDoe/Desktop/data/datasets/test_result/result_c_porites.csv", float_format="%.4g")
result_c_porites_f.to_csv("C:/Users/JohnDoe/Desktop/data/datasets/test_result/result_c_porites_f.csv", float_format="%.4g")
result_p_porites.to_csv("C:/Users/JohnDoe/Desktop/data/datasets/test_result/result_p_porites.csv", float_format="%.4g")
result_p_porites_f.to_csv("C:/Users/JohnDoe/Desktop/data/datasets/test_result/result_p_porites_f.csv", float_format="%.4g")

#%%



#stylophora
#change ASV
stylophora_ASV_genus = stylophora_ASV.groupby(["Genus"]).sum().transpose()
stylophora_ASV_family = stylophora_ASV.groupby(["Family"]).sum().transpose()
stylophora_ASV_order = stylophora_ASV.groupby(["Order"]).sum().transpose()
stylophora_ASV_class = stylophora_ASV.groupby(["Class"]).sum().transpose()
stylophora_ASV_phylum = stylophora_ASV.groupby(["Phylum"]).sum().transpose()

stylophora_ASV_genus_f = stylophora_filtered_ASV.groupby(["Genus"]).sum().transpose()
stylophora_ASV_family_f = stylophora_filtered_ASV.groupby(["Family"]).sum().transpose()
stylophora_ASV_order_f = stylophora_filtered_ASV.groupby(["Order"]).sum().transpose()
stylophora_ASV_class_f = stylophora_filtered_ASV.groupby(["Class"]).sum().transpose()
stylophora_ASV_phylum_f = stylophora_filtered_ASV.groupby(["Phylum"]).sum().transpose()

#merge
stylophora_genus = stylophora_ED50.merge(stylophora_ASV_genus, left_index=True, right_index=True)
stylophora_family = stylophora_ED50.merge(stylophora_ASV_family, left_index=True, right_index=True)
stylophora_order = stylophora_ED50.merge(stylophora_ASV_order, left_index=True, right_index=True)
stylophora_class = stylophora_ED50.merge(stylophora_ASV_class, left_index=True, right_index=True)
stylophora_phylum = stylophora_ED50.merge(stylophora_ASV_phylum, left_index=True, right_index=True)

stylophora_genus_f = stylophora_ED50.merge(stylophora_ASV_genus_f, left_index=True, right_index=True)
stylophora_family_f = stylophora_ED50.merge(stylophora_ASV_family_f, left_index=True, right_index=True)
stylophora_order_f = stylophora_ED50.merge(stylophora_ASV_order_f, left_index=True, right_index=True)
stylophora_class_f = stylophora_ED50.merge(stylophora_ASV_class_f, left_index=True, right_index=True)
stylophora_phylum_f = stylophora_ED50.merge(stylophora_ASV_phylum_f, left_index=True, right_index=True)

#%%
import coracle

x = stylophora_genus.iloc[:,3:]
y = stylophora_genus[['ED50']]
start = time.time()
result_g_stylophora = coracle.coracle(x, y, alpha_l1 = 10**(-2.9), alpha_clr = 10**(-0.4))
print(time.time()- start)

x_f = stylophora_genus_f.iloc[:,3:]
y_f = stylophora_genus_f[['ED50']]
start = time.time()
result_g_stylophora_f = coracle.coracle(x_f, y_f, alpha_l1 = 10**(-2.9), alpha_clr = 10**(-0.4))
print(time.time()- start)

x = stylophora_family.iloc[:,3:]
y = stylophora_family[['ED50']]
start = time.time()
result_f_stylophora = coracle.coracle(x, y, alpha_l1 = 10**(-2.9), alpha_clr = 10**(-0.4))
print(time.time()- start)

x_f = stylophora_family_f.iloc[:,3:]
y_f = stylophora_family_f[['ED50']]
start = time.time()
result_f_stylophora_f = coracle.coracle(x_f, y_f, alpha_l1 = 10**(-2.9), alpha_clr = 10**(-0.4))
print(time.time()- start)

x = stylophora_order.iloc[:,3:]
y = stylophora_order[['ED50']]
start = time.time()
result_o_stylophora = coracle.coracle(x, y, alpha_l1 = 10**(-2.9), alpha_clr = 10**(-0.4))
print(time.time()- start)

x_f = stylophora_order_f.iloc[:,3:]
y_f = stylophora_order_f[['ED50']]
start = time.time()
result_o_stylophora_f = coracle.coracle(x_f, y_f, alpha_l1 = 10**(-2.9), alpha_clr = 10**(-0.4))
print(time.time()- start)

x = stylophora_class.iloc[:,3:]
y = stylophora_class[['ED50']]
start = time.time()
result_c_stylophora = coracle.coracle(x, y, alpha_l1 = 10**(-2.9), alpha_clr = 10**(-0.4))
print(time.time()- start)

x_f = stylophora_class_f.iloc[:,3:]
y_f = stylophora_class_f[['ED50']]
start = time.time()
result_c_stylophora_f = coracle.coracle(x_f, y_f, alpha_l1 = 10**(-2.9), alpha_clr = 10**(-0.4))
print(time.time()- start)

x = stylophora_phylum.iloc[:,3:]
y = stylophora_phylum[['ED50']]
start = time.time()
result_p_stylophora = coracle.coracle(x, y, alpha_l1 = 10**(-2.9), alpha_clr = 10**(-0.4))
print(time.time()- start)

x_f = stylophora_phylum_f.iloc[:,3:]
y_f = stylophora_phylum_f[['ED50']]
start = time.time()
result_p_stylophora_f = coracle.coracle(x_f, y_f, alpha_l1 = 10**(-2.9), alpha_clr = 10**(-0.4))
print(time.time()- start)

#%%
result_g_stylophora.to_csv("C:/Users/JohnDoe/Desktop/data/datasets/test_result/result_g_stylophora.csv", float_format="%.4g")
result_g_stylophora_f.to_csv("C:/Users/JohnDoe/Desktop/data/datasets/test_result/result_g_stylophora_f.csv", float_format="%.4g")
result_f_stylophora.to_csv("C:/Users/JohnDoe/Desktop/data/datasets/test_result/result_f_stylophora.csv", float_format="%.4g")
result_f_stylophora_f.to_csv("C:/Users/JohnDoe/Desktop/data/datasets/test_result/result_f_stylophora_f.csv", float_format="%.4g")
result_o_stylophora.to_csv("C:/Users/JohnDoe/Desktop/data/datasets/test_result/result_o_stylophora.csv", float_format="%.4g")
result_o_stylophora_f.to_csv("C:/Users/JohnDoe/Desktop/data/datasets/test_result/result_o_stylophora_f.csv", float_format="%.4g")
result_c_stylophora.to_csv("C:/Users/JohnDoe/Desktop/data/datasets/test_result/result_c_stylophora.csv", float_format="%.4g")
result_c_stylophora_f.to_csv("C:/Users/JohnDoe/Desktop/data/datasets/test_result/result_c_stylophora_f.csv", float_format="%.4g")
result_p_stylophora.to_csv("C:/Users/JohnDoe/Desktop/data/datasets/test_result/result_p_stylophora.csv", float_format="%.4g")
result_p_stylophora_f.to_csv("C:/Users/JohnDoe/Desktop/data/datasets/test_result/result_p_stylophora_f.csv", float_format="%.4g")













"""


#%%%
### CBASS84 dataset
#cbass84_meta = pd.read_csv("cbass84_metadata.txt", sep= "[ \t \t]+", error_bad_lines=False) #, usecols=(["Sample", "Site", "Reef", "Species"]))
cbass84_ASV = pd.read_csv(directory + "cbass84_ASVs.txt", sep= "\s+")
cbass84_ED50 = pd.read_csv(directory + "cbass84_ED50s.txt", sep= "\s+")
cbass84_tax = cbass84_ASV.iloc[:, -5:] #get taxonomic information
cbass84_ASV = cbass84_ASV.iloc[:, :-5] #get rest
cbass84_ASV = cbass84_ASV.transpose()
cbass84_ASV = cbass84_ASV.loc[:, (cbass84_ASV != 0).any(axis=0)]
cbass84_ASV = cbass84_ASV.transpose()
cbass84_ASV = cbass84_ASV.merge(cbass84_tax, left_index=True, right_index=True)
cbass84_tax = cbass84_ASV.iloc[:, -5:] #get taxonomic information

#change ASV
cbass84_ASV_genus = cbass84_ASV.groupby(["Genus"]).sum()
cbass84_ASV_genus = cbass84_ASV_genus.transpose()
cbass84_ASV_family = cbass84_ASV.groupby(["Family"]).sum()
cbass84_ASV_family = cbass84_ASV_family.transpose()
cbass84_ASV_order = cbass84_ASV.groupby(["Order"]).sum()
cbass84_ASV_order = cbass84_ASV_order.transpose()
cbass84_ASV_class = cbass84_ASV.groupby(["Class"]).sum()
cbass84_ASV_class = cbass84_ASV_class.transpose()


#change ED50
cbass84_ED50["Site"] = None
cbass84_ED50["Site"][cbass84_ED50["Sample"].str.contains("AF")] = "Al Fahal (AF)"
cbass84_ED50["Site"][cbass84_ED50["Sample"].str.contains("ExT")] = "Tahala (ExT)"
cbass84_ED50["Site"][cbass84_ED50["Sample"].str.contains("PrT")] = "Tahala (PrT)"
cbass84_ED50["Site"][cbass84_ED50["Sample"].str.contains("ICN")] = "Interuniversity Institute for Marine Science (IUI)"
cbass84_ED50["Species"] = "Stylophora pistillata"
cbass84_ED50.set_index(["Sample"], inplace = True)
#specific for correlation analysis:
cbass84_ED50 = pd.DataFrame(cbass84_ED50["ED50"])




ASV = cbass84_ASV.iloc[:, :-5]
cbass84_filtered_ASV = fs(ASV, cbass84_ED50, cbass84_tax)

#%%



#change ASV
cbass84_ASV_genus = cbass84_ASV.groupby(["Genus"]).sum()
cbass84_ASV_genus = cbass84_ASV_genus.transpose()
cbass84_ASV_family = cbass84_ASV.groupby(["Family"]).sum()
cbass84_ASV_family = cbass84_ASV_family.transpose()
cbass84_ASV_order = cbass84_ASV.groupby(["Order"]).sum()
cbass84_ASV_order = cbass84_ASV_order.transpose()
cbass84_ASV_class = cbass84_ASV.groupby(["Class"]).sum()
cbass84_ASV_class = cbass84_ASV_class.transpose()
cbass84_ASV_phylum = cbass84_ASV.groupby(["Phylum"]).sum()
cbass84_ASV_phylum = cbass84_ASV_phylum.transpose()


cbass84_ASV_genus_f = cbass84_filtered_ASV.groupby(["Genus"]).sum()
cbass84_ASV_genus_f = cbass84_ASV_genus_f.transpose()
cbass84_ASV_family_f = cbass84_filtered_ASV.groupby(["Family"]).sum()
cbass84_ASV_family_f = cbass84_ASV_family_f.transpose()
cbass84_ASV_order_f = cbass84_filtered_ASV.groupby(["Order"]).sum()
cbass84_ASV_order_f = cbass84_ASV_order_f.transpose()
cbass84_ASV_class_f = cbass84_filtered_ASV.groupby(["Class"]).sum()
cbass84_ASV_class_f = cbass84_ASV_class_f.transpose()
cbass84_ASV_phylum_f = cbass84_filtered_ASV.groupby(["Phylum"]).sum()
cbass84_ASV_phylum_f = cbass84_ASV_phylum_f.transpose()



#%%

#merge
cbass84_genus = cbass84_ED50.merge(cbass84_ASV_genus, left_index=True, right_index=True)
cbass84_family = cbass84_ED50.merge(cbass84_ASV_family, left_index=True, right_index=True)
cbass84_order = cbass84_ED50.merge(cbass84_ASV_order, left_index=True, right_index=True)
cbass84_class = cbass84_ED50.merge(cbass84_ASV_class, left_index=True, right_index=True)
cbass84_phylum = cbass84_ED50.merge(cbass84_ASV_phylum, left_index=True, right_index=True)


cbass84_genus_f = cbass84_ED50.merge(cbass84_ASV_genus_f, left_index=True, right_index=True)
cbass84_family_f = cbass84_ED50.merge(cbass84_ASV_family_f, left_index=True, right_index=True)
cbass84_order_f = cbass84_ED50.merge(cbass84_ASV_order_f, left_index=True, right_index=True)
cbass84_class_f = cbass84_ED50.merge(cbass84_ASV_class_f, left_index=True, right_index=True)
cbass84_phylum_f = cbass84_ED50.merge(cbass84_ASV_phylum_f, left_index=True, right_index=True)


#%%

import coracle

x = cbass84_genus.iloc[:,3:]
y = cbass84_genus[['ED50']]

start = time.time()
result_g84 = coracle.coracle(x, y, alpha_l1 = 10**(-2.9), alpha_clr = 10**(-0.4))
print(time.time()- start)

x_f = cbass84_genus_f.iloc[:,3:]
y_f = cbass84_genus_f[['ED50']]

start = time.time()
result_g84_f = coracle.coracle(x_f, y_f, alpha_l1 = 10**(-2.9), alpha_clr = 10**(-0.4))
print(time.time()- start)

#%%

x = cbass84_family.iloc[:,3:]
y = cbass84_family[['ED50']]

start = time.time()
result_f84 = coracle.coracle(x, y, alpha_l1 = 10**(-2.9), alpha_clr = 10**(-0.4))
print(time.time()- start)

x_f = cbass84_family_f.iloc[:,3:]
y_f = cbass84_family_f[['ED50']]

start = time.time()
result_f84_f = coracle.coracle(x_f, y_f, alpha_l1 = 10**(-2.9), alpha_clr = 10**(-0.4))
print(time.time()- start)

#%%

x = cbass84_order.iloc[:,3:]
y = cbass84_order[['ED50']]

start = time.time()
result_o84 = coracle.coracle(x, y, alpha_l1 = 10**(-2.9), alpha_clr = 10**(-0.4))
print(time.time()- start)

x_f = cbass84_order_f.iloc[:,3:]
y_f = cbass84_order_f[['ED50']]

start = time.time()
result_o84_f = coracle.coracle(x_f, y_f, alpha_l1 = 10**(-2.9), alpha_clr = 10**(-0.4))
print(time.time()- start)

#%%

x = cbass84_class.iloc[:,3:]
y = cbass84_class[['ED50']]

start = time.time()
result_c84 = coracle.coracle(x, y, alpha_l1 = 10**(-2.9), alpha_clr = 10**(-0.4))
print(time.time()- start)

x_f = cbass84_class_f.iloc[:,3:]
y_f = cbass84_class_f[['ED50']]

start = time.time()
result_c84_f = coracle.coracle(x_f, y_f, alpha_l1 = 10**(-2.9), alpha_clr = 10**(-0.4))
print(time.time()- start)


#%%

x = cbass84_phylum.iloc[:,3:]
y = cbass84_phylum[['ED50']]

start = time.time()
result_p84 = coracle.coracle(x, y, alpha_l1 = 10**(-2.9), alpha_clr = 10**(-0.4))
print(time.time()- start)

x_f = cbass84_phylum_f.iloc[:,3:]
y_f = cbass84_phylum_f[['ED50']]

start = time.time()
result_p84_f = coracle.coracle(x_f, y_f, alpha_l1 = 10**(-2.9), alpha_clr = 10**(-0.4))
print(time.time()- start)

#%%

result_g84.to_csv("C:/Users/JohnDoe/Desktop/data/datasets/result_g84.csv", float_format="%.2g")
result_g84_f.to_csv("C:/Users/JohnDoe/Desktop/data/datasets/result_g84_f.csv", float_format="%.2g")
result_f84.to_csv("C:/Users/JohnDoe/Desktop/data/datasets/result_f84.csv", float_format="%.2g")
result_f84_f.to_csv("C:/Users/JohnDoe/Desktop/data/datasets/result_f84_f.csv", float_format="%.2g")
result_o84.to_csv("C:/Users/JohnDoe/Desktop/data/datasets/result_o84.csv", float_format="%.2g")
result_o84_f.to_csv("C:/Users/JohnDoe/Desktop/data/datasets/result_o84_f.csv", float_format="%.2g")
result_c84.to_csv("C:/Users/JohnDoe/Desktop/data/datasets/result_c84.csv", float_format="%.2g")
result_c84_f.to_csv("C:/Users/JohnDoe/Desktop/data/datasets/result_c84_f.csv", float_format="%.2g")
result_p84.to_csv("C:/Users/JohnDoe/Desktop/data/datasets/result_p84.csv", float_format="%.2g")
result_p84_f.to_csv("C:/Users/JohnDoe/Desktop/data/datasets/result_p84_f.csv", float_format="%.2g")


#%%

directory = "/Users/JohnDoe/Desktop/data/datasets/"


#Roseo Dataset
Roseo_ASV = pd.read_csv(directory + "Roseo_ASVs.txt", sep= "\s+")
Roseo_ED50 = pd.read_csv(directory + "Roseo_ED50s.txt", sep= "\s+", index_col=(0))
Roseo_tax = Roseo_ASV.iloc[:, -6:-1] #get taxonomic information (ignore species as to many NaNs)
Roseo_ASV = Roseo_ASV.iloc[:, :-6] #get rest
######filter out zero columns within the algorithm and INCLUDE taxonomy filtering!!!
Roseo_ASV = Roseo_ASV.transpose()
Roseo_ASV = Roseo_ASV.loc[:, (Roseo_ASV != 0).any(axis=0)]
Roseo_ASV = Roseo_ASV.transpose()
Roseo_ASV = Roseo_ASV.merge(Roseo_tax, left_index=True, right_index=True)
Roseo_tax = Roseo_ASV.iloc[:, -5:] #get taxonomic information
Roseo_ASV = Roseo_ASV.iloc[:, :-5] #get rest
#######
Roseo_ASV_f = fs(Roseo_ASV, Roseo_ED50, Roseo_tax)
Roseo_ASV = Roseo_ASV.merge(Roseo_tax, left_index=True, right_index=True)

#%%



#change ASV
Roseo_ASV_genus = Roseo_ASV.groupby(["Genus"]).sum()
Roseo_ASV_genus = Roseo_ASV_genus.transpose()
Roseo_ASV_family = Roseo_ASV.groupby(["Family"]).sum()
Roseo_ASV_family = Roseo_ASV_family.transpose()
Roseo_ASV_order = Roseo_ASV.groupby(["Order"]).sum()
Roseo_ASV_order = Roseo_ASV_order.transpose()
Roseo_ASV_class = Roseo_ASV.groupby(["Class"]).sum()
Roseo_ASV_class = Roseo_ASV_class.transpose()
Roseo_ASV_phylum = Roseo_ASV.groupby(["Phylum"]).sum()
Roseo_ASV_phylum = Roseo_ASV_phylum.transpose()

#change ASV with hfs
Roseo_ASV_genus_f = Roseo_ASV_f.groupby(["Genus"]).sum()
Roseo_ASV_genus_f = Roseo_ASV_genus_f.transpose()
Roseo_ASV_family_f = Roseo_ASV_f.groupby(["Family"]).sum()
Roseo_ASV_family_f = Roseo_ASV_family_f.transpose()
Roseo_ASV_order_f = Roseo_ASV_f.groupby(["Order"]).sum()
Roseo_ASV_order_f = Roseo_ASV_order_f.transpose()
Roseo_ASV_class_f = Roseo_ASV_f.groupby(["Class"]).sum()
Roseo_ASV_class_f = Roseo_ASV_class_f.transpose()
Roseo_ASV_phylum_f = Roseo_ASV_f.groupby(["Phylum"]).sum()
Roseo_ASV_phylum_f = Roseo_ASV_phylum_f.transpose()

#change ED50
Roseo_ED50 = pd.read_csv(directory + "Roseo_ED50s.txt", sep= "\s+")
Roseo_ED50["Site"] = "Al Fahal"
Roseo_ED50["Species"] = None
Roseo_ED50["Species"][Roseo_ED50["Sample"].str.contains("humilis")] = "Acropora humilis"
Roseo_ED50["Species"][Roseo_ED50["Sample"].str.contains("cytherea")] = "Acropora cytherea"
Roseo_ED50["Species"][Roseo_ED50["Sample"].str.contains("hemprichii")] = "Acropora hemprichii"
Roseo_ED50.set_index(["Sample"], inplace = True)

#merge
roseo_genus = Roseo_ED50.merge(Roseo_ASV_genus, left_index=True, right_index=True)
roseo_family = Roseo_ED50.merge(Roseo_ASV_family, left_index=True, right_index=True)
roseo_order = Roseo_ED50.merge(Roseo_ASV_order, left_index=True, right_index=True)
roseo_class = Roseo_ED50.merge(Roseo_ASV_class, left_index=True, right_index=True)
roseo_phylum = Roseo_ED50.merge(Roseo_ASV_phylum, left_index=True, right_index=True)
#merge with hfs
roseo_genus_f = Roseo_ED50.merge(Roseo_ASV_genus_f, left_index=True, right_index=True)
roseo_family_f = Roseo_ED50.merge(Roseo_ASV_family_f, left_index=True, right_index=True)
roseo_order_f = Roseo_ED50.merge(Roseo_ASV_order_f, left_index=True, right_index=True)
roseo_class_f = Roseo_ED50.merge(Roseo_ASV_class_f, left_index=True, right_index=True)
roseo_phylum_f = Roseo_ED50.merge(Roseo_ASV_phylum_f, left_index=True, right_index=True)


#%%
import coracle


x = roseo_genus.iloc[:,3:]
y = roseo_genus[['ED50']]

start = time.time()
result_gr = coracle.coracle(x, y, alpha_l1 = 10**(-2.9), alpha_clr = 10**(-0.4))
print(time.time()- start)

x_f = roseo_genus_f.iloc[:,3:]
y_f = roseo_genus_f[['ED50']]

start = time.time()
result_gr_f = coracle.coracle(x_f, y_f, alpha_l1 = 10**(-2.9), alpha_clr = 10**(-0.4))
print(time.time()- start)

#%%

x = roseo_family.iloc[:,3:]
y = roseo_family[['ED50']]

start = time.time()
result_fr = coracle.coracle(x, y, alpha_l1 = 10**(-2.9), alpha_clr = 10**(-0.4))
print(time.time()- start)

x_f = roseo_family_f.iloc[:,3:]
y_f = roseo_family_f[['ED50']]

start = time.time()
result_fr_f = coracle.coracle(x_f, y_f, alpha_l1 = 10**(-2.9), alpha_clr = 10**(-0.4))
print(time.time()- start)

#%%

x = roseo_order.iloc[:,3:]
y = roseo_order[['ED50']]

start = time.time()
result_or = coracle.coracle(x, y, alpha_l1 = 10**(-2.9), alpha_clr = 10**(-0.4))
print(time.time()- start)

x_f = roseo_order_f.iloc[:,3:]
y_f = roseo_order_f[['ED50']]

start = time.time()
result_or_f = coracle.coracle(x_f, y_f, alpha_l1 = 10**(-2.9), alpha_clr = 10**(-0.4))
print(time.time()- start)

#%%

x = roseo_class.iloc[:,3:]
y = roseo_class[['ED50']]

start = time.time()
result_cr = coracle.coracle(x, y, alpha_l1 = 10**(-2.9), alpha_clr = 10**(-0.4))
print(time.time()- start)

x_f = roseo_class_f.iloc[:,3:]
y_f = roseo_class_f[['ED50']]

start = time.time()
result_cr_f = coracle.coracle(x_f, y_f, alpha_l1 = 10**(-2.9), alpha_clr = 10**(-0.4))
print(time.time()- start)


#%%

x = roseo_phylum.iloc[:,3:]
y = roseo_phylum[['ED50']]

start = time.time()
result_pr = coracle.coracle(x, y, alpha_l1 = 10**(-2.9), alpha_clr = 10**(-0.4))
print(time.time()- start)

x_f = roseo_phylum_f.iloc[:,3:]
y_f = roseo_phylum_f[['ED50']]

start = time.time()
result_pr_f = coracle.coracle(x_f, y_f, alpha_l1 = 10**(-2.9), alpha_clr = 10**(-0.4))
print(time.time()- start)

#%%

result_gr.to_csv("C:/Users/JohnDoe/Desktop/data/datasets/result_gr.csv", float_format="%.2g")
result_gr_f.to_csv("C:/Users/JohnDoe/Desktop/data/datasets/result_gr_f.csv", float_format="%.2g")
result_fr.to_csv("C:/Users/JohnDoe/Desktop/data/datasets/result_fr.csv", float_format="%.2g")
result_fr_f.to_csv("C:/Users/JohnDoe/Desktop/data/datasets/result_fr_f.csv", float_format="%.2g")
result_or.to_csv("C:/Users/JohnDoe/Desktop/data/datasets/result_or.csv", float_format="%.2g")
result_or_f.to_csv("C:/Users/JohnDoe/Desktop/data/datasets/result_or_f.csv", float_format="%.2g")
result_cr.to_csv("C:/Users/JohnDoe/Desktop/data/datasets/result_cr.csv", float_format="%.2g")
result_cr_f.to_csv("C:/Users/JohnDoe/Desktop/data/datasets/result_cr_f.csv", float_format="%.2g")
result_pr.to_csv("C:/Users/JohnDoe/Desktop/data/datasets/result_pr.csv", float_format="%.2g")
result_pr_f.to_csv("C:/Users/JohnDoe/Desktop/data/datasets/result_pr_f.csv", float_format="%.2g")

"""