# -*- coding: utf-8 -*-
"""
Created on Tue Aug 16 08:05:35 2022

@author: JohnDoe
"""

import pandas as pd
import numpy as np





def helper_bottom_up_propagation(hir, upper, lower, f_acc, target_name):
    """
    create dictionary of upper hierarchical levels that holds the correlation 
    matrix of all included lower hierarchical level features

    Parameters
    ----------
    hir : pd.dataframe
        hierarchy
    upper : string
        upper hierarchical levels
    lower : string
        lower hierarchical levels
    f_acc : pd.dataframe
        accumulated feature file
    target_name : pd.Index
        name of target variable

    Returns
    -------
    final : dictionary of pd.dataframes
        dictionary of grouped features according to hierarchy

    """   
    
    final = {} #final 
    gr_upper = hir[upper].unique() #unique entries of the upper hierarchy level
    
    for i in gr_upper: #for every entry
        new_list = hir.loc[hir[upper] == i][lower].unique() #insert all included unique lower hierarchy features
        new_list = np.append(target_name, new_list) #append target variable... for general purpose: change
        new_list = [x for x in new_list if not pd.isnull(x)] #delete NaNs, they won't be propagated to the next level (as there might be problems with NaN subgroups from other groups) but will appear with their lowest defined hierarchical level for the first time!
        
        Corr = f_acc[new_list].corr() #get correlation
        final[i] = Corr #save correlation
        
    return final


def helper_hierarchical_unicor_metric(corr_dict, target_name):
    """
    function that computes the UniCor metric given a correlation dictionary

    Parameters
    ----------
    corr_dict : pd.dataframe
        dictionary that gives correlations between features and correlation of 
        feature to target variable
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
            if column == target_name:
                continue
            else:
                strain = corr_matr[column]
                
                n = len(strain)
                fc_corr = strain.iloc[0] #feature-target correlation
                ff_corr = (strain.iloc[1:].sum()-1)/(n-1.99) #feature-feature correlation (-own corr and -ED50 corr + 0.01 to prevent zero division)
                unicor = (0.5*abs(fc_corr)) - (0.5*ff_corr)
                
                
                if column in strains:
                    #print("found it")                    
                    #print(metrics[column])
                    #print(unicor)
                    if  metrics[column] > unicor:
                        unicor = metrics[column]            
                    #print(unicor)
                metrics[column], strains[column] = unicor, strain
    return metrics, strains





def unicor_metric(features, target):
    """
    function that computes the unicor metric for a simple featureset without a hierarchy

    Parameters
    ----------
    features : pd.Dataframe
        continuous feature set, index is sample number, columns are the features.
    target : pd.Dataframe
        continuous target variable, one column with the data, index is the sample number.

    Returns
    -------
    metric : dictionary
        features and respective unicor metrics

    """
    ### check input
    ###########################################################################
    #check type
    if not isinstance(target, (pd.DataFrame, pd.Series)):
        raise TypeError("TypeError exception thrown. Expected pandas dataframe for target")
    if not isinstance(features, (pd.DataFrame)):
        raise TypeError("TypeError exception thrown. Expected pandas dataframe for features")
    #check dimensions
    if features.ndim != 2:  
        raise ValueError("ValueError exception thrown. Expected features to have two dimensions")
    if features.shape[0] != target.shape[0]:
        raise ValueError("ValueError exception thrown. Expected features and target to have the same number of samples")
    ###########################################################################
    
    
    ### clean zero columns/rows (features that are not present in any samples)
    ###########################################################################
    features = features.loc[:, (features != 0).any(axis=0)]
    ###########################################################################
    
    ### compute unicor metric for all features
    ###########################################################################
    table = target.merge(features, left_index=True, right_index=True) #merge features and target variable
    cor_matrix = table.corr() #get correlation matrix
    metric = {} #dictionary to save results to
    
    
    for feature in cor_matrix: #for each feature
        if feature in target.columns: #skip target varialbe
            continue
        else: 
            n = len(cor_matrix[feature])
            fc_corr = cor_matrix[feature].iloc[0] #feature-target correlation
            ff_corr = (cor_matrix[feature].iloc[1:].sum()-1)/(n-1.99) #feature-feature correlation (-own corr and -target variable corr + 0.01 to prevent zero division)
            unicor = (0.5*abs(fc_corr)) - (0.5*ff_corr) #unicor metric computation
            metric[feature] = unicor #save in dictionary
    ###########################################################################
    return metric




def uniCor(features, target, hir, threshold = 0.15):
    """
    UniCor algortihm: Takes in a hierarchical continous dataset and propagates 
    significant features (UniCor metric) to higher hierarchical levels in order 
    to preserve crucial information while significantly reducing the feature set 
    in a biologically meaningfull way

    Parameters
    ----------
    features : pd.dataframe
        features dataset
    target : pd.dataframe
        continuous target variable
    hir : pd.dataframe
        table that contains the  hierarchy
    threshold : bool, optional
        UniCor metric threshold. Can only be between 0 and 1. 
        Optimal values depend on the dataset. The default is 0.15.

    Raises
    ------
    ValueError
        check allowed ranges for the thresholds, dimensions for the dataframes
    TypeError
        check types of input

    Returns
    -------
    new_features : pd.dataframe
        returns features with changed hierarchical levels, still every level 
        can be accessed and used

    """
    ### check input
    ###########################################################################
    #check type
    if not isinstance(features, (pd.DataFrame)):
        raise TypeError("TypeError exception thrown. Expected pandas dataframe for features")
    if not isinstance(target, (pd.DataFrame, pd.Series)):
        raise TypeError("TypeError exception thrown. Expected pandas dataframe for target")
    if not isinstance(hir, (pd.DataFrame)):
        raise TypeError("TypeError exception thrown. Expected pandas dataframe for hir")
    if not isinstance(threshold, (float, int)):
        raise TypeError("TypeError exception thrown. Expected float for threshold")
    #check values    
    if threshold > 1 or threshold < 0:
        raise ValueError("ValueError exception thrown. threshold is expected to be a float between 0 and 1")
    #check dimensions
    if features.ndim != 2 or hir.ndim != 2:  
        raise ValueError("ValueError exception thrown. Expected features and hir to have two dimensions")
    if features.shape[1] != hir.shape[0]:
        raise ValueError("ValueError exception thrown. Expected features and hir to have the same number of features")
    if features.shape[0] != target.shape[0]:
        raise ValueError("ValueError exception thrown. Expected features and target to have the same number of samples")
    ###########################################################################
    
    
    ### clean zero columns/rows (features that are not present in any samples)
    ###########################################################################
    features = features.loc[:, (features != 0).any(axis=0)]
    features = features.transpose()
    ###########################################################################
        
        
    ### check order of hierarchy
    ###########################################################################
    hir_levels = list(hir.columns) #get hierarchy level names
    features_hir = features.merge(hir, how="left", left_index=True, right_index=True) #merge to compare size
    size = {} #dict to check sizes
    last_size = 0 #helper variable to save the last size
    count = 0 #helper variable
    
    for i in hir_levels: #for every hierarchical level
        print(i)
        size[i] = len(features_hir[i].unique()) #get number of unique entries in that level
        print(size[i])
        
        if size[i] > last_size: #if number of unique entries in current level is smaller than in the last level
            count += 1 #increase count
        last_size = size[i] #set current size to future last size
    
    #check if hierarchical order is fulfilled
    if count == len(hir_levels): #check if hierarchical order is ascending from left to right
        hir_levels.reverse() #if that is the case, reverse the order
        print("hierarchical order has been reversed to start with the lowest level from left") #print a note
    elif count > 0 and count < len(hir_levels): #if order is not unambiguous:
        print("Warning: unclear hierarchy as number of unique entries is neither clearly ascending nor descending! The order will be unchanged and used as if number of unique features is descending from left to right") #print a warning
    ###########################################################################
    
    ### check if hierarchical order is strict
    ###########################################################################
    strict = True # helper variable
    for m in range(len(hir_levels)-1): #for every hierarchy level
        for n in hir[hir_levels[m]].unique(): #for every entry in that hir level
            proxy = hir[hir_levels[m+1]].loc[hir[hir_levels[m]] == n].unique()#create list of all groups that include this entry
            
            if len(proxy) > 1: #if more the one parent group
                print(n, proxy) #print out the potential problem child
                strict = False #hierarchy is not strict anymore
    
    if strict == False: #if hierarchy is not strict print warning
        print("hierarchy is not strict")
    ###########################################################################
    
    
    ### bottom up UniCor propagation
    ###########################################################################
    # create dictionarys
    f_dict = {} #for full feature tables per hierarchy level
    hir_groups = {} #hierarchical groups
    hg = {} #new hierarchy groups
    metrics = {} #UniCor hierarchy groups
    strains = {}
    new_features = features_hir.copy() #to keep old and new features table distinct
    tv = target.columns
    
    # propagation level by level
    for j in range(len(hir_levels)-1): #for every hierarchical level
        i = hir_levels[j] #access specific hierarchy level
        f_acc = new_features.groupby(i).sum() #accumulate
        f_acc = f_acc.transpose()
        f_acc = target.merge(f_acc, left_index=True, right_index=True) #merge with target variable
        f_dict[i] = f_acc #save features in dictionary - might be interesting to see later
        hierarchy = new_features[hir_levels] #get back the hierarchy
        hir_groups[i] = hierarchy[i].unique() #get the hierarchy groups of this level
        
        
        # create hierarchy with selected features
        if j <= len(hir_levels)-2: #-1 due to index starting at zero, -1 as we don't need the last computation
            hg[i] = helper_bottom_up_propagation(hierarchy, hir_levels[j+1], hir_levels[j], f_acc, tv) #get hierarchy groups
            #if nan -> treat as same group and just use next higher hierarchy (no selection/propagation to next level)
            metrics[i], strains[i] = helper_hierarchical_unicor_metric(hg[i], tv) #get metrics for hierarchy groups
        
        # propagate applicable group to next level
        for q in metrics[i]: #for specific strain metric in metrics
            strain_metric = metrics[i][q] #assign strain metric to variable
            if strain_metric > threshold: #check if bigger then defined threshold
                new_features.loc[new_features[i] == q,hir_levels[j+1]] = q #add relevant strain to next higher level
    ###########################################################################

    return new_features #return the adapted features incuding the hierarchical information



