# -*- coding: utf-8 -*-
"""
Created on Wed Mar  1 15:43:34 2023

@author: JohnDoe2Go
"""


#%%
import pandas as pd
#from uni_cor import uniCorP, unicor_metric
from unicor import unicorp, unicor_metric
directory = "C:/Users/JohnDoe2Go/Downloads/" #use your path to the variables

### CBASS84 dataset
ASV = pd.read_csv(directory + "cbass.csv", index_col=0) #read in cbass ASVs
tax = pd.read_csv(directory + "cbass_tax.csv", index_col=0) #read in cbass taxonomic information
### 2. Split Target Variable
y = ASV["ED50"].to_frame()
### 3. Combine ASV with tax
x = ASV.iloc[:,3:]

#%%
import unittest


#unit testing
class TestUnicorMetricFunctionality(unittest.TestCase):
    def test_wrong_input_features2(self): #tests if TypeError is thrown if features file is not a pandas dataframe
        self.assertRaises(TypeError, unicor_metric, 1, y)
    def test_wrong_input_target2(self): #tests if TypeError is thrown if target file is not a pandas dataframe, or series
        self.assertRaises(TypeError, unicor_metric, x, 1)
    def test_wrong_dimensions5(self): #tests ValueError is thrown if features and target don't match in their dimensions (number of samples)
        self.assertRaises(ValueError, unicor_metric, y, tax)
    def test_output_type(self):
        result = unicor_metric(x, y)
        self.assertIsInstance(result, dict)
    def test_contains_all_features(self):
        result = unicor_metric(x, y)
        self.assertTrue(set(result.keys()).issubset(set(x.columns)))
    def test_metric_values_reasonable(self):
        result = unicor_metric(x, y)
        for val in result.values():
            self.assertIsInstance(val, float)
    def test_wrong_method_unicor_metric(self):  # invalid correlation method
        self.assertRaises(ValueError, unicor_metric, x, y, method='kendall')
    def test_wrong_transformation_unicor_metric(self):  # invalid transformation
        self.assertRaises(ValueError, unicor_metric, x, y, transformation='lognorm')



class TestUnicorPFunctionality(unittest.TestCase):
    def test_output_type(self):
        result = unicorp(x, y, tax, threshold=0.01)
        self.assertIsInstance(result, pd.DataFrame)
    def test_includes_hierarchy_columns(self):
        result = unicorp(x, y, tax, threshold=0.01)
        for col in tax.columns:
            self.assertIn(col, result.columns)
    def test_propagation_changes_hierarchy(self):
        result_thresh = unicorp(x, y, tax, threshold=0.2)
        result_topx = unicorp(x, y, tax, top_k=20)
        self.assertFalse(result_thresh.equals(result_topx))  # top_k and threshold should give different results
    def test_wrong_input_features(self): #tests if TypeError is thrown if features file is not a pandas dataframe
        self.assertRaises(TypeError, unicorp, 1, y, tax, top_k=1)
    def test_wrong_input_target(self): #tests if TypeError is thrown if target file is not a pandas dataframe, or series
        self.assertRaises(TypeError, unicorp, x, 1, tax, top_k=1)
    def test_wrong_input_hierarchy(self): #tests if TypeError is thrown if tax file is not a pandas dataframe
        self.assertRaises(TypeError, unicorp, x, y, 1, top_k=1)
    def test_wrong_input_threshold(self): #tests if TypeError is thrown for non-float/integer input for threshold
        self.assertRaises(TypeError, unicorp, x, y, tax, threshold="1")
    def test_wrong_threshold1(self): #tests ValueError is thrown if threshold not between 0 and 1
        self.assertRaises(ValueError, unicorp, x, y, tax, threshold=1.5)
    def test_wrong_threshold2(self): #tests ValueError is thrown if threshold not between 0 and 1
        self.assertRaises(ValueError, unicorp, x, y, tax, threshold=-1)
    def test_wrong_dimensions1(self): #tests ValueError is thrown if features is not two dimensional
        self.assertRaises(ValueError, unicorp, y, y, tax, top_k=1)
    def test_wrong_dimensions2(self): #tests ValueError is thrown if hierarchy is not two dimesnional
        self.assertRaises(ValueError, unicorp, x, y, y, top_k=1)
    def test_wrong_dimensions3(self): #tests ValueError is thrown if features and hierarchy don't match in their dimensions (number of features)
        self.assertRaises(ValueError, unicorp, x, y, x, top_k=1)
    def test_wrong_dimensions4(self): #tests ValueError is thrown if features and target don't match in their dimensions (number of samples)
        self.assertRaises(ValueError, unicorp, x, tax, tax, top_k=1)
    def test_wrong_dimensions6(self): #tests ValueError is thrown if features is not two dimensional
        self.assertRaises(ValueError, unicorp, y, y, tax, top_k=1)
    def test_both_topx_and_threshold_given(self):  # should raise error if both are provided
        self.assertRaises(ValueError, unicorp, x, y, tax, threshold=0.5, top_k=5)
    def test_neither_topx_nor_threshold_given(self):  # should raise error if neither is provided
        self.assertRaises(ValueError, unicorp, x, y, tax)
    def test_wrong_method_unicorp(self):  # invalid correlation method
        self.assertRaises(ValueError, unicorp, x, y, tax, threshold=0.1, method='kendall')
    def test_wrong_transformation_unicorp(self):  # invalid transformation
        self.assertRaises(ValueError, unicorp, x, y, tax, threshold=0.1, transformation='lognorm')
    def test_wrong_topx_type(self):  # top_k must be integer
        self.assertRaises(TypeError, unicorp, x, y, tax, top_k='five')
    def test_wrong_topx_value(self):  # top_k must be >= 1
        self.assertRaises(ValueError, unicorp, x, y, tax, top_k=0)


#test runner
if __name__ == "__main__":
    unittest.main()
    
