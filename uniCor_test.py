# -*- coding: utf-8 -*-
"""
Created on Wed Mar  1 15:43:34 2023

@author: JohnDoe2Go
"""

import pandas as pd
from uni_cor import uniCor

directory = "C:/Users/JohnDoe2Go/Downloads/" #use your path to the variables

### CBASS84 dataset
ASV = pd.read_csv(directory + "cbass.csv", index_col=0) #read in cbass ASVs
tax = pd.read_csv(directory + "cbass_tax.csv", index_col=0) #read in cbass taxonomic information
### 2. Split Target Variable
y = ASV["ED50"].to_frame()
### 3. Combine ASV with tax
x = ASV.iloc[:,3:]


result = uniCor(x, y, tax)
#%%
import unittest


#unit testing
class TestWrongInput(unittest.TestCase):
    def test_wrong_input_asv(self): #tests if TypeError is thrown if asv file is not a pandas dataframe
        self.assertRaises(TypeError, uniCor, 1, y, tax)
    def test_wrong_input_target(self): #tests if TypeError is thrown if target file is not a pandas dataframe, or series
        self.assertRaises(TypeError, uniCor, x, 1, tax)
    def test_wrong_input_tax(self): #tests if TypeError is thrown if tax file is not a pandas dataframe
        self.assertRaises(TypeError, uniCor, x, y, 1)
    def test_wrong_input_threshold(self): #tests if TypeError is thrown for non-float/integer input for threshold
        self.assertRaises(TypeError, uniCor, x, y, tax, "1")
    def test_wrong_threshold1(self): #tests ValueError is thrown if threshold not between 0 and 1
        self.assertRaises(ValueError, uniCor, x, y, tax, 1.5)
    def test_wrong_threshold2(self): #tests ValueError is tax if threshold not between 0 and 1
        self.assertRaises(ValueError, uniCor, x, y, tax, -1)
    def test_wrong_dimensions1(self): #tests ValueError is thrown if asv is not two dimensional
        self.assertRaises(ValueError, uniCor, y, y, tax)
    def test_wrong_dimensions2(self): #tests ValueError is thrown if tax is not two dimesnional
        self.assertRaises(ValueError, uniCor, x, y, y)
    def test_wrong_dimensions3(self): #tests ValueError is thrown if asv and tax don't match in their dimensions (number of features)
        self.assertRaises(ValueError, uniCor, x, y, x)
    def test_wrong_dimensions4(self): #tests ValueError is thrown if asv and target don't match in their dimensions (number of samples)
        self.assertRaises(ValueError, uniCor, x, tax, tax)

#test runner
if __name__ == "__main__":
    unittest.main()
    
