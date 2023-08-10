# -*- coding: utf-8 -*-
"""
Created on Tue Aug 16 08:05:35 2022

@author: JohnDoe
"""

import pandas as pd
import numpy as np





def bottom_up_propagation(tax, upper, lower, aSV_acc, tv):
    """
    create dictionary of upper taxonomic levels that holds the correlation matrix of all included lower taxonomic level taxes

    Parameters
    ----------
    tax : pd.dataframe
        taxonomy
    upper : string
        upper taxonomic levels
    lower : string
        lower taxonomic levels
    aSV_acc : pd.dataframe
        accumulated ASV
    tv : pd.Index
        name of target variable

    Returns
    -------
    final : dictionary of pd.dataframes
        dictionary of grouped features according to taxonomic hierarchy

    """   
    
    final = {} #final 
    gr_upper = tax[upper].unique() #unique entries of the upper tax
    for i in gr_upper: #for every entry
        new_list = tax.loc[tax[upper] == i][lower].unique() #insert all included unique lower taxa
        new_list = np.append(tv, new_list) #append target variable... for general purpose: CHANGE!!!
        new_list = [x for x in new_list if not pd.isnull(x)] #delete NaNs, they won't be propagated to the next level (as there might be problems with NaN subgroups from other groups) but will appear with their lowest defined taxonomic level for the first time!
        
        Corr = aSV_acc[new_list].corr() #get correlation
        final[i] = Corr
    return final


def uc_metric(corr_dict, target):
    """
    function that computes the UniCor metric given a correlation dictionary

    Parameters
    ----------
    corr_dict : pd.dataframe
        dictionary that gives correlations between features and correlation of feature to target variable
    target : pd.Index
        name of target variable

    Returns
    -------
    metrics : dict
        the UniCor metric for each feature

    """
    metrics = {}
    strains = {}
    
    # go through all entries in the correlation dictionary to access the correlation matrix
    for i in corr_dict:
        corr_matr = corr_dict[i]
        for column in corr_matr:
            if column == target:
                continue
            else:
                strain = corr_matr[column]
                
                n = len(strain)
                fc_corr = strain.iloc[0] #feature-target correlation
                ff_corr = (strain.iloc[1:].sum()-1)/(n-1.99) #feature-feature correlation (-own corr and -ED50 corr + 0.01 to prevent zero division)
                unicor = (0.5*abs(fc_corr)) - (0.5*ff_corr)
                
                
                if column in strains:
                    print("found it")
                    print(metrics[column])
                    print(unicor)
                    if  metrics[column] > unicor:
                        unicor = metrics[column]            
                    print(unicor)
                metrics[column], strains[column] = unicor, strain
    return metrics, strains



#%%
def uniCor(aSV, target, tax, threshold = 0.15):
    """
    UniCor algortihm: Takes in a hierarchical continous dataset and propagates significant features (UniCor metric) to higher taxonomic levels in order to preserve crucial information while significantly reducing the feature set in a biologically meaningfull way

    Parameters
    ----------
    aSV : pd.dataframe
        ASV dataset
    target : pd.dataframe
        continuous target variable
    tax : pd.dataframe
        taxonomic hierarchy
    threshold : bool, optional
        UniCor metric threshold. Can only be between 0 and 1. Optimal values depend on the dataset. The default is 0.15.

    Raises
    ------
    ValueError
        check allowed ranges for the thresholds, dimensions for the dataframes
    TypeError
        check types of input

    Returns
    -------
    new_ASV : pd.dataframe
        returns ASV with changed taxonomic levels, still every level can be accessed and used

    """
    ### check input
    ###########################################################################
    #check type
    if not isinstance(aSV, (pd.DataFrame)):
        raise TypeError("TypeError exception thrown. Expected pandas dataframe for ASV")
    if not isinstance(target, (pd.DataFrame, pd.Series)):
        raise TypeError("TypeError exception thrown. Expected pandas dataframe for target")
    if not isinstance(tax, (pd.DataFrame)):
        raise TypeError("TypeError exception thrown. Expected pandas dataframe for tax")
    if not isinstance(threshold, (float, int)):
        raise TypeError("TypeError exception thrown. Expected float for threshold")
    #check values    
    if threshold > 1 or threshold < 0:
        raise ValueError("ValueError exception thrown. threshold is expected to be a float between 0 and 1")
    #check dimensions
    if aSV.ndim != 2 or tax.ndim != 2:  
        raise ValueError("ValueError exception thrown. Expected ASV and tax to have two dimensions")
    if aSV.shape[1] != tax.shape[0]:
        raise ValueError("ValueError exception thrown. Expected ASV and tax to have the same number of samples")
    if aSV.shape[0] != target.shape[0]:
        raise ValueError("ValueError exception thrown. Expected ASV and target to have the same number of samples")
    ###########################################################################
    
    
    ### clean zero columns/rows (ASVs that are not present in any samples)
    ###########################################################################
    aSV = aSV.loc[:, (aSV != 0).any(axis=0)]
    aSV = aSV.transpose()
    ###########################################################################
        
        
    ### check order/hierarchy of taxonomy
    ###########################################################################
    tax_levels = list(tax.columns) #get tax level names
    aSV_tax = aSV.merge(tax, how="left", left_index=True, right_index=True) #merge to compare size
    size = {} #dict to check sizes
    last_size = 0 #helper variable to save the last size
    count = 0 #helper variable
    
    for i in tax_levels: #for every taxonomic level
        print(i)
        size[i] = len(aSV_tax[i].unique()) #get number of unique entries in that level
        print(size[i])
        
        if size[i] > last_size: #if number of unique entries in current level is smaller than in the last level
            count += 1 #increase count
        last_size = size[i] #set current size to future last size
    
    #check if hierarchical order is fulfilled
    if count == len(tax_levels): #check if hierarchical order is ascending from left to right
        tax_levels.reverse() #if that is the case, reverse the order
        print("taxonomic order has been reversed to start with the lowest level from left") #print a note
    elif count > 0 and count < len(tax_levels): #if order is not unambiguous:
        print("Warning: unclear hierarchy as number of unique entries is neither clearly ascending nor descending! The order will be unchanged and used as if number of unique features is descending from left to right") #print a warning
    ###########################################################################
    
    ### check if hierarchical order is strict
    ###########################################################################
    strict = True # helper variable
    for m in range(len(tax_levels)-1): #for every tax level
        for n in tax[tax_levels[m]].unique(): #for every entry in that tax level
            proxy = tax[tax_levels[m+1]].loc[tax[tax_levels[m]] == n].unique()#create list of all groups that include this entry
            
            if len(proxy) > 1: #if more the one parent group
                print(n, proxy) #print out the potential problem child
                strict = False #hierarchy is not strict anymore
    
    if strict == False: #if hierarchy is not strict print warning
        print("hierarchy is not strict")
    ###########################################################################
    
    
    ### bottom up UniCor propagation
    ###########################################################################
    # create dictionarys
    aSVdict = {} #for full ASV tables per tax level
    tax_groups = {} #taxonomic groups
    hg = {} #new hierarchy groups
    metrics = {} #UniCor hierarchy groups
    strains = {}
    new_ASV = aSV_tax.copy() #to keep old and new ASV table distinct
    tv = target.columns
    
    # propagation level by level
    for j in range(len(tax_levels)-1): #for every taxonomic level
        i = tax_levels[j] #access specific tax level
        aSV_acc = new_ASV.groupby(i).sum() #accumulate
        aSV_acc = aSV_acc.transpose()
        aSV_acc = target.merge(aSV_acc, left_index=True, right_index=True) #merge with target variable
        aSVdict[i] = aSV_acc #save ASV in dictionary - might be interesting to see later
        taxonomy = new_ASV[tax_levels] #get back the taxonomy
        tax_groups[i] = taxonomy[i].unique() #get the tax groups of this level
        
        
        # create hierarchy with selected features
        if j <= len(tax_levels)-2: #-1 due to index starting at zero, -1 as we don't need the last computation
            hg[i] = bottom_up_propagation(taxonomy, tax_levels[j+1], tax_levels[j], aSV_acc, tv) #get hierarchy groups
            #if nan -> treat as same group and just use next higher tax (no selection/propagation to next level)
            metrics[i], strains[i] = uc_metric(hg[i], tv) #get metrics for hierarchy groups
        
        # propagate applicable group to next level
        for q in metrics[i]: #for specific strain metric in metrics
            strain_metric = metrics[i][q] #assign strain metric to variable
            if strain_metric > threshold: #check if bigger then defined threshold
                new_ASV.loc[new_ASV[i] == q,tax_levels[j+1]] = q #add relevant strain to next higher level
    ###########################################################################

    return new_ASV #return the adapted asv incuding the taxonomic information



