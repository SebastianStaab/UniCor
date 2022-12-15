# -*- coding: utf-8 -*-
"""
Created on Wed Dec 14 11:35:53 2022

@author: JohnDoe
"""
from hierarchical_filtering import fs
import coracle
import pandas as pd


def tdc(ASV, target, tax, threshold = 0.2, uc = True, uc_threshold = 0.15):
    """
    Algorithm to analyse hierarchical ASV datasets and identify bacteria/ASV that are associated with a continuous target variable
    First it uses the UniCor algorithm to propagate the most unique and valuable features to the next higher level (preserve information). Then a top-down skimming approach uses Coracle to select the [threshold] top partition consecutively at each level. The lowest level is once again anylized with Coracle and given back as the final result.

    Parameters
    ----------
    ASV : pd.dataframe
        ASV dataset
    target_var : pd.dataframe
        Continuous target variable
    tax : pd.dataframe
        Taxonomic hierarchy
    threshold : bool, optional
        Propagation threshold. Has to be between 0 and 1. Gives the proportion of the best performing groups that are propagated to the next lower level. Default is 0.2. In general you want this value to be as small as possible for runtime efficiency reasons.
    uc_threshold : bool, optional
        UniCor metric threshold. Can only be between 0 and 1. Optimal values depend on the dataset. The default is 0.15.

    Raises
    ------
    TypeError
        check types of input
    ValueError
        check allowed ranges for the thresholds, dimensions for the dataframes

    Returns
    -------
    resultO : pd.dataframe
        Final Coracle result on the lowest level
 """       

    
    
    ### check input
    ###########################################################################
    #check type
    if not isinstance(ASV, (pd.DataFrame)):
        raise TypeError("TypeError exception thrown. Expected pandas dataframe for ASV")
    if not isinstance(target, (pd.DataFrame)):
        raise TypeError("TypeError exception thrown. Expected pandas dataframe for target")
    if not isinstance(tax, (pd.DataFrame, pd.Series)):
        raise TypeError("TypeError exception thrown. Expected pandas dataframe for tax")
    if not isinstance(threshold, (float)):
        raise TypeError("TypeError exception thrown. Expected float for threshold")
    if not isinstance(uc, (bool)):
        raise TypeError("TypeError exception thrown. Expected bool for uc")
    if not isinstance(uc_threshold, (float)):
        raise TypeError("TypeError exception thrown. Expected float uc_threshold")
    #check values    
    if threshold > 1 or threshold < 0:
        raise ValueError("ValueError exception thrown. threshold is expected to be a float between 0 and 1")
    if uc_threshold > 1 or uc_threshold < 0:
        raise ValueError("ValueError exception thrown. uc_threshold is expected to be a float between 0 and 1")
    #check dimensions
    if ASV.ndim != 2 or tax.ndim != 2:  
        raise ValueError("ValueError exception thrown. Expected ASV and tax to have two dimensions")
    if ASV.shape[0] != tax.shape[0]:
        raise ValueError("ValueError exception thrown. Expected ASV and tax to have the same number of samples")
    if ASV.shape[1] != target.shape[0]:
        raise ValueError("ValueError exception thrown. Expected ASV and target to have the same number of samples")
    ###########################################################################
    
    
    ### clean zero columns/rows (ASVs that are not present in any samples)
    ###########################################################################
    ASV = ASV.transpose()
    ASV = ASV.loc[:, (ASV != 0).any(axis=0)]
    ASV = ASV.transpose()
    ###########################################################################
    
    
    ### check order/hierarchy of taxonomy
    ###########################################################################
    tax_levels = list(tax.columns) #get tax level names
    ASV_tax = ASV.merge(tax, how="left", left_index=True, right_index=True) #merge to compare size
    size = {} #dict to check sizes
    last_size = 0 #helper variable to save the last size
    count = 0 #helper variable
    
    for i in tax_levels: #for every taxonomic level
        print(i)
        size[i] = len(ASV_tax[i].unique()) #get number of unique entries in that level
        print(size[i])
        
        if size[i] < last_size: #if number of unique entries in current level is smaller than in the last level
            count += 1 #increase count
        last_size = size[i] #set current size to future last size
    
    #check if hierarchical order is fulfilled
    if count == len(tax_levels)-1: #check if hierarchical order is ascending from left to right(-1 because first comparison is with count=0)
        tax_levels.reverse() #if that is the case, reverse the order
        print("taxonomic order has been reversed") #print a note
    elif count > 0 and count < len(tax_levels)-1: #if order is not unambiguous: (-1 because first comparison is with count=0)
        print("Warning: unclear hierarchy as number of unique entries is neither clearly ascending nor descending! The order will be unchanged and used as if descending from left to right") #print a warning
    ###########################################################################
        
        
    ### UniCor to filter ASV
    ###########################################################################
    if uc == True: #if UniCor is wanted
        filtered_ASV = fs(ASV, target, tax, threshold=uc_threshold)
    else: #if not
        filtered_ASV = ASV.merge(tax, how="left", left_index=True, right_index=True)
    ###########################################################################
    
    ### Top-Down Coracle propagation
    ###########################################################################
    for i in range(len(tax_levels)): #for every taxonomic level
        selected = [] #list of selected groups to propagate
        highest_level = filtered_ASV.groupby([tax_levels[i]]).sum().transpose() #start at the highest level
        highest_level = target.merge(highest_level, left_index=True, right_index=True) #merge with target variable
        x = highest_level.iloc[:,1:] 
        y = highest_level.iloc[:,0].to_frame()
        
        resultO = coracle.coracle(x, y, alpha_l1 = 10**(-2.9), alpha_clr = 10**(-0.4)) #run coracle
        result = resultO.iloc[3:,0] #select only the groups and their respective score
        for j in range(int(result.size * threshold)+1): #for every group within the threshold partition (always round up to prevent zero features)
            if result[j] > 0: #if score is not null
                #print(result.index[j])
                selected.append(result.index[j]) #add to selected-list
            else: #else print a warning that no more features will be added
                print("Warning: as the score reached zero, threshold partition will not be fulfilled but will include the top", ((100*j)/(result.size * threshold)), "percent")
                break
        
        
        filtered_ASV = filtered_ASV[filtered_ASV[tax_levels[i]].isin(selected)] #propagate them to the next level
        
        
        
    ###########################################################################
    
    return resultO #return final coracle result