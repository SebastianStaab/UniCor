# -*- coding: utf-8 -*-
"""
Created on Wed Mar  1 15:43:34 2023

@author: JohnDoe2Go
"""

import pandas as pd
from uni_cor import uniCor

directory = "/Users/JohnDoe2Go/Desktop/PhD/10_data/datasets/"

#%%
### CBASS84 dataset
#cbass84_meta = pd.read_csv("cbass84_metadata.txt", sep= "[ \t \t]+", error_bad_lines=False) #, usecols=(["Sample", "Site", "Reef", "Species"]))
cbass84_ASV = pd.read_csv(directory + "cbass84_ASVs.txt", sep= "\s+")

#for index in cbass84_ASV.index.values:
#    cbass84_ASV.loc[index] = cbass84_ASV.loc[index].fillna(str(index))

cbass84_ED50 = pd.read_csv(directory + "cbass84_ED50s.txt", sep= "\s+")
cbass84_tax = cbass84_ASV.iloc[:, -5:] #get taxonomic information
cbass84_ASV = cbass84_ASV.iloc[:, :-5] #get rest
#cbass84_ASV = cbass84_ASV.transpose()
#cbass84_ASV = cbass84_ASV.loc[:, (cbass84_ASV != 0).any(axis=0)]
#cbass84_ASV = cbass84_ASV.transpose()
cbass84_ASV = cbass84_ASV.merge(cbass84_tax, left_index=True, right_index=True)
cbass84_tax = cbass84_ASV.iloc[:, -5:] #get taxonomic information





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

cbass84_ASV = cbass84_ASV.iloc[:, :-5].transpose()

result = uniCor(cbass84_ASV, cbass84_ED50, cbass84_tax)





#%%
import unittest


#unit testing
class TestWrongInput(unittest.TestCase):
    def test_wrong_input_asv(self): #tests if TypeError is thrown if asv file is not a pandas dataframe
        self.assertRaises(TypeError, uniCor, 1, cbass84_ED50, cbass84_tax)
    def test_wrong_input_target(self): #tests if TypeError is thrown if target file is not a pandas dataframe, or series
        self.assertRaises(TypeError, uniCor, cbass84_ASV, 1, cbass84_tax)
    def test_wrong_input_tax(self): #tests if TypeError is thrown if tax file is not a pandas dataframe
        self.assertRaises(TypeError, uniCor, cbass84_ASV, cbass84_ED50, 1)
    def test_wrong_input_threshold(self): #tests if TypeError is thrown for non-float/integer input for threshold
        self.assertRaises(TypeError, uniCor, cbass84_ASV, cbass84_ED50, cbass84_tax, "1")
    def test_wrong_threshol1(self): #tests ValueError is thrown if threshold not between 0 and 1
        self.assertRaises(ValueError, uniCor, cbass84_ASV, cbass84_ED50, cbass84_tax, 1.5)
    def test_wrong_threshold2(self): #tests ValueError is thrown if threshold not between 0 and 1
        self.assertRaises(ValueError, uniCor, cbass84_ASV, cbass84_ED50, cbass84_tax, -1)
    def test_wrong_dimensions1(self): #tests ValueError is thrown if asv is not two dimensional
        self.assertRaises(ValueError, uniCor, cbass84_ED50, cbass84_ED50, cbass84_tax)
    def test_wrong_dimensions2(self): #tests ValueError is thrown if tax is not two dimesnional
        self.assertRaises(ValueError, uniCor, cbass84_ASV, cbass84_ED50, cbass84_ED50)
    def test_wrong_dimensions3(self): #tests ValueError is thrown if asv and tax don't match in their dimensions (number of features)
        self.assertRaises(ValueError, uniCor, cbass84_ASV, cbass84_ED50, cbass84_ASV)
    def test_wrong_dimensions4(self): #tests ValueError is thrown if asv and target don't match in their dimensions (number of samples)
        self.assertRaises(ValueError, uniCor, cbass84_ASV, cbass84_tax, cbass84_tax)
    
"""
#check values    
if aSV.ndim != 2 or tax.ndim != 2:  
    raise ValueError("ValueError exception thrown. Expected ASV and tax to have two dimensions")
if aSV.shape[1] != tax.shape[0]:
    raise ValueError("ValueError exception thrown. Expected ASV and tax to have the same number of samples")
if aSV.shape[0] != target.shape[0]:
    raise ValueError("ValueError exception thrown. Expected ASV and target to have the same number of samples")
"""

#test runner
if __name__ == "__main__":
    unittest.main()
    

#%%
cbass84 = cbass84_ASV.transpose().merge(cbass84_tax, left_index=True, right_index=True)

cbass84_ASV_genus = cbass84.groupby(["Genus"]).sum()
cbass84_ASV_genus = cbass84_ASV_genus.transpose()
cbass84_ASV_genus = cbass84_ED50.merge(cbass84_ASV_genus, left_index=True, right_index=True)
cbass84_ASV_family = cbass84.groupby(["Family"]).sum()
cbass84_ASV_family = cbass84_ASV_family.transpose()
cbass84_ASV_family = cbass84_ED50.merge(cbass84_ASV_family, left_index=True, right_index=True)
cbass84_ASV_order = cbass84.groupby(["Order"]).sum()
cbass84_ASV_order = cbass84_ASV_order.transpose()
cbass84_ASV_order = cbass84_ED50.merge(cbass84_ASV_order, left_index=True, right_index=True)
cbass84_ASV_class = cbass84.groupby(["Class"]).sum()
cbass84_ASV_class = cbass84_ASV_class.transpose()
cbass84_ASV_class = cbass84_ED50.merge(cbass84_ASV_class, left_index=True, right_index=True)
cbass84_ASV_phylum = cbass84.groupby(["Phylum"]).sum()
cbass84_ASV_phylum = cbass84_ASV_phylum.transpose()
cbass84_ASV_phylum = cbass84_ED50.merge(cbass84_ASV_phylum, left_index=True, right_index=True)

#%%

uc_cbass84_ASV_genus = result.groupby(["Genus"]).sum()
uc_cbass84_ASV_genus = uc_cbass84_ASV_genus.transpose()
uc_cbass84_ASV_genus = cbass84_ED50.merge(uc_cbass84_ASV_genus, left_index=True, right_index=True)
uc_cbass84_ASV_family = result.groupby(["Family"]).sum()
uc_cbass84_ASV_family = uc_cbass84_ASV_family.transpose()
uc_cbass84_ASV_family = cbass84_ED50.merge(uc_cbass84_ASV_family, left_index=True, right_index=True)
uc_cbass84_ASV_order = result.groupby(["Order"]).sum()
uc_cbass84_ASV_order = uc_cbass84_ASV_order.transpose()
uc_cbass84_ASV_order = cbass84_ED50.merge(uc_cbass84_ASV_order, left_index=True, right_index=True)
uc_cbass84_ASV_class = result.groupby(["Class"]).sum()
uc_cbass84_ASV_class = uc_cbass84_ASV_class.transpose()
uc_cbass84_ASV_class = cbass84_ED50.merge(uc_cbass84_ASV_class, left_index=True, right_index=True)
uc_cbass84_ASV_phylum = result.groupby(["Phylum"]).sum()
uc_cbass84_ASV_phylum = uc_cbass84_ASV_phylum.transpose()
uc_cbass84_ASV_phylum = cbass84_ED50.merge(uc_cbass84_ASV_phylum, left_index=True, right_index=True)


#%%
import coracle

### without unicor
#same for normal and uc version since it is the lowest level
x = cbass84_ASV_genus.iloc[:,1:]
y = cbass84_ASV_genus[['ED50']]
result_g84 = coracle.coracle(x, y, alpha_l1 = 10**(-2.9), alpha_clr = 10**(-0.4))
x = cbass84_ASV_family.iloc[:,1:]
y = cbass84_ASV_family[['ED50']]
result_f84 = coracle.coracle(x, y, alpha_l1 = 10**(-2.9), alpha_clr = 10**(-0.4))
x = cbass84_ASV_order.iloc[:,1:]
y = cbass84_ASV_order[['ED50']]
result_o84 = coracle.coracle(x, y, alpha_l1 = 10**(-2.9), alpha_clr = 10**(-0.4))
x = cbass84_ASV_class.iloc[:,1:]
y = cbass84_ASV_class[['ED50']]
result_c84 = coracle.coracle(x, y, alpha_l1 = 10**(-2.9), alpha_clr = 10**(-0.4))
x = cbass84_ASV_phylum.iloc[:,1:]
y = cbass84_ASV_phylum[['ED50']]
result_p84 = coracle.coracle(x, y, alpha_l1 = 10**(-2.9), alpha_clr = 10**(-0.4))

### with unicor
x = uc_cbass84_ASV_genus.iloc[:,1:]
y = uc_cbass84_ASV_genus[['ED50']]
result_uc_g84 = coracle.coracle(x, y, alpha_l1 = 10**(-2.9), alpha_clr = 10**(-0.4))
x = uc_cbass84_ASV_family.iloc[:,1:]
y = uc_cbass84_ASV_family[['ED50']]
result_uc_f84 = coracle.coracle(x, y, alpha_l1 = 10**(-2.9), alpha_clr = 10**(-0.4))
x = uc_cbass84_ASV_order.iloc[:,1:]
y = uc_cbass84_ASV_order[['ED50']]
result_uc_o84 = coracle.coracle(x, y, alpha_l1 = 10**(-2.9), alpha_clr = 10**(-0.4))
x = uc_cbass84_ASV_class.iloc[:,1:]
y = uc_cbass84_ASV_class[['ED50']]
result_uc_c84 = coracle.coracle(x, y, alpha_l1 = 10**(-2.9), alpha_clr = 10**(-0.4))
x = uc_cbass84_ASV_phylum.iloc[:,1:]
y = uc_cbass84_ASV_phylum[['ED50']]
result_uc_p84 = coracle.coracle(x, y, alpha_l1 = 10**(-2.9), alpha_clr = 10**(-0.4))


#%%

result_g84.to_csv("C:/Users/JohnDoe2Go/Desktop/PhD/10_data/datasets/unicor_result/g84.csv", float_format="%.4g")
result_f84.to_csv("C:/Users/JohnDoe2Go/Desktop/PhD/10_data/datasets/unicor_result/f84.csv", float_format="%.4g")
result_o84.to_csv("C:/Users/JohnDoe2Go/Desktop/PhD/10_data/datasets/unicor_result/o84.csv", float_format="%.4g")
result_c84.to_csv("C:/Users/JohnDoe2Go/Desktop/PhD/10_data/datasets/unicor_result/c84.csv", float_format="%.4g")
result_p84.to_csv("C:/Users/JohnDoe2Go/Desktop/PhD/10_data/datasets/unicor_result/p84.csv", float_format="%.4g")
result_uc_f84.to_csv("C:/Users/JohnDoe2Go/Desktop/PhD/10_data/datasets/unicor_result/uc_f84.csv", float_format="%.4g")
result_uc_o84.to_csv("C:/Users/JohnDoe2Go/Desktop/PhD/10_data/datasets/unicor_result/uc_o84.csv", float_format="%.4g")
result_uc_c84.to_csv("C:/Users/JohnDoe2Go/Desktop/PhD/10_data/datasets/unicor_result/uc_c84.csv", float_format="%.4g")
result_uc_p84.to_csv("C:/Users/JohnDoe2Go/Desktop/PhD/10_data/datasets/unicor_result/uc_p84.csv", float_format="%.4g")