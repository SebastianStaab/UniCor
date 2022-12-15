# -*- coding: utf-8 -*-
"""
Created on Tue Aug 16 08:05:35 2022

@author: JohnDoe
"""

import pandas as pd
import numpy as np
import time




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
def fs(ASV, target, tax, threshold = 0.15):
    """
    UniCor algortihm: Takes in a hierarchical continous dataset and propagates significant features (UniCor metric) to higher taxonomic levels in order to preserve crucial information while significantly reducing the feature set in a biologically meaningfull way

    Parameters
    ----------
    ASV : pd.dataframe
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
    if not isinstance(ASV, (pd.DataFrame)):
        raise TypeError("TypeError exception thrown. Expected pandas dataframe for ASV")
    if not isinstance(target, (pd.DataFrame)):
        raise TypeError("TypeError exception thrown. Expected pandas dataframe for target")
    if not isinstance(tax, (pd.DataFrame, pd.Series)):
        raise TypeError("TypeError exception thrown. Expected pandas dataframe for tax")
    if not isinstance(threshold, (float)):
        raise TypeError("TypeError exception thrown. Expected float for threshold")
    #check values    
    if threshold > 1 or threshold < 0:
        raise ValueError("ValueError exception thrown. threshold is expected to be a float between 0 and 1")
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
    
    
    ### bottom up UniCor propagation
    ###########################################################################
    # create dictionarys
    ASVdict = {} #for full ASV tables per tax level
    Corr = {} #for full correlation matrix per tax level
    tax_groups = {} #taxonomic groups
    hg = {} #new hierarchy groups
    metrics = {} #UniCor hierarchy groups
    new_ASV = ASV_tax.copy() #to keep old and new ASV table distinct
    
    # propagation level by level
    for j in range(len(tax_levels)-1): #for every taxonomic level
        i = tax_levels[j] #access specific tax level
        x = new_ASV.groupby(i).sum() #accumulate
        x = x.transpose()
        x = target.merge(x, left_index=True, right_index=True) #merge with target variable
        ASVdict[i] = x #save ASV in dictionary
        Corr[i] = x.corr() #save correlations in dictionary
        taxonomy = new_ASV[tax_levels] #get back the taxonomy
        tax_groups[i] = taxonomy[i].unique() #get the tax groups of this level
        
        
        # create hierarchy with selected features
        if j <= len(tax_levels)-2: #-1 due to index starting at zero, -1 as we don't need the last computation
            hg[i] = hfsbu(taxonomy, tax_levels[j+1], tax_levels[j], Corr[i]) #get hierarchy groups
            #if nan -> treat as same group and just use next higher tax (no selection/propagation to next level)
            metrics[i] = metric(hg[i]) #get metrics for hierarchy groups
        
        # propagate applicable group to next level
        for q in metrics[i]: #for specific strain metric in metrics
            strain_metric = metrics[i][q] #assign strain metric to variable
            if strain_metric > threshold: #check if bigger then defined threshold
                new_ASV.loc[new_ASV[i] == q,tax_levels[j+1]] = q #add relevant strain to next higher level
    ###########################################################################
    

    return new_ASV #return the adapted asv incuding the taxonomic information




#%%%
### CBASS84 dataset
directory = "/Users/JohnDoe/Desktop/data/datasets/"
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

#%%


ASV = cbass84_ASV.iloc[:, :-5]
cbass84_filtered_ASV = fs(ASV, cbass84_ED50, cbass84_tax)



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

