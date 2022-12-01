# -*- coding: utf-8 -*-
"""
Created on Tue Aug 16 08:05:35 2022

@author: JohnDoe
"""

import pandas as pd
import numpy as np





def hfsbu(tax, upper, lower, corr_matrix):
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
    corr_matrix : pd.dataframe
        correlation matrices

    Returns
    -------
    final : dictionary of pd.dataframes
        dictionary of grouped features according to taxonomic hierarchy

    """
    
    final = {} #final 
    gr_upper = tax[upper].unique() #unique entries of the upper tax
    for i in gr_upper: #for every entry
        new_list = tax.loc[tax[upper] == i][lower].unique() #insert all included unique lower taxa
        new_list = np.append(["ED50"], new_list) #append target variable... for general purpose: CHANGE!!!
        new_list = [x for x in new_list if not pd.isnull(x)] #delete NaNs, they won't be propagated to the next level (as there might be problems with NaN subgroups from other groups) but will appear with their lowest defined taxonomic level for the first time!
        
        
        
        new_corr = corr_matrix.loc[corr_matrix.columns.isin(new_list)]
        new_corr = new_corr[new_list]
        final[i] = new_corr
    return final


def metric(corr_dict):
    """
    function that computes the UniCor metric given a correlation dictionary

    Parameters
    ----------
    corr_dict : pd.dataframe
        dictionary that gives correlations between features and correlation of feature to target variable

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
            if column == "ED50":
                continue
            else:
                strain = corr_matr[column]
                
                n = len(strain)
                fc_corr = strain.iloc[0]
                #print(fc_corr)
                ff_corr = (strain.iloc[1:].sum()-1)/(n-1.99) #-own corr and -ED50 corr + 0.01 to prevent zero division
                #print(ff_corr)
                comp = (0.5*abs(fc_corr)) - (0.5*ff_corr)
                #print(strain)
                #print(metric)
                #print("----")
                metrics[column], strains[column] = comp, strain
    return metrics



#%%
def fs(ASV, target_var, taxonomy, threshold = 0.15):
    """
    UniCor algortihm: Takes in a hierarchical continous dataset and propagates significant features (UniCor metric) to higher taxonomic levels in order to preserve crucial information while significantly reducing the feature set in a biologically meaningfull way

    Parameters
    ----------
    ASV : pd.dataframe
        ASV dataset
    target_var : pd.dataframe
        continuous target variable
    taxonomy : pd.dataframe
        taxonomic hierarchy
    threshold : bool, optional
        UniCor metric threshold. Can only be between 0 and 1. Optimal values depend on the dataset. The default is 0.15.

    Raises
    ------
    ValueError
        Taxonomy is expected to have a clear hierarchy

    Returns
    -------
    new_ASV : pd.dataframe
        returns ASV with changed taxonomic levels, still every level can be accessed and used

    """
    
    
    #clean zero columns/rows
    ASV = ASV.transpose()
    ASV = ASV.loc[:, (ASV != 0).any(axis=0)]
    ASV = ASV.transpose()
    #get ASV including tax
    tax_size = taxonomy.columns
    ASV = ASV.merge(taxonomy, how="left", left_index=True, right_index=True)
    # get tax levels
    tax_levels = list(taxonomy.columns)
    #print(tax_levels)
    
    #check sizes
    size = {}
    for i in tax_levels:
        size[i] = len(ASV[i].unique())
        print(i)
        print(size[i])
    
    
    
    ###order & check taxonomic levels
    #order them correctly if wrong (e.g. unique entries in phylum vs class level )
    if tax_levels[0] > tax_levels[1]: 
        tax_levels.reverse()
    """
    unique_entries = 1000000    # absurdely high baselevel
    #check if there is a clear hierarchy
    for i in tax_levels: #for every level
        
        if len(ASV[i].unique()) > unique_entries: #if number of unique values bigger than next higher level
            raise ValueError("ValueError exception thrown. Taxonomy is expected to have a clear hierarchy (lower number of unique strains in higher taxonomic brackets)")
        unique_entries = len(ASV[i].unique()) #set current number of unique entries as new max number for the next level
    """    

    # create empty dictionarys
    ASVdict = {} #for full ASV tables per tax level
    Corr = {} #for full correlation matrix per tax level
    tax_groups = {} #create tax groups
    hg = {} #new hierarchy groups
    metrics = {} #UniCor hierarchy groups
    new_ASV = ASV.copy()
    
    for j in range(len(tax_levels)-1): #for every taxonomic level
        i = tax_levels[j] #access specific tax level
        #print(i)
        #print("-------------")
        x = new_ASV.groupby(i).sum() #accumulate
        x = x.transpose()
        x = target_var.merge(x, left_index=True, right_index=True) #merge with targetvar
        ASVdict[i] = x #save in dictionary
        Corr[i] = x.corr() #save correlations in dictionary
        taxonomy = new_ASV[tax_size]
        tax_groups[i] = taxonomy[i].unique()
        
        
        #create hierarchy with selected features
        if j <= len(tax_levels)-2: #-1 due to index starting at zero, -1 as we don't need the last computation
            
            hg[i] = hfsbu(taxonomy, tax_levels[j+1], tax_levels[j], Corr[i]) #get hierarchy groups
            #if nan -> treat as same group and just use next higher tax (no selection/propagation to next level)
            metrics[i] = metric(hg[i]) #get metrics for hierarchy groups
            #return hg
        
        
        for q in metrics[i]: #for specific strain metric in metrics
            strain_metric = metrics[i][q] #assign strain metric to variable
            if strain_metric > threshold: #check if bigger then defined threshold
                print(q)
                #print(strain_metric)
                #new_ASV[new_ASV[tax_levels[o] == q]] = 
                new_ASV.loc[new_ASV[i] == q,tax_levels[j+1]] = q #add relevant strain to next higher level
     

    return new_ASV

