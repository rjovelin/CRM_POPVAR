# -*- coding: utf-8 -*-
"""
Created on Wed Jun 10 19:12:34 2015

@author: Richard
"""

from scipy.stats import rankdata
import numpy as np


# In the Benjamini Hochberg method, the P-values are first sorted and ranked.
# The smallest value gets rank 1, the second rank 2, and the largest gets rank N.
# Then, each P-value is multiplied by N and divided by its assigned rank to give
# the adjusted P-values. In order to restrict the false discovery rate to
# (say) 0.05, all the genes with adjusted P-values less than 0.05 are selected.



def Benjamini_Hochberg_correction(p_vals):
    '''
    (list) -> dict
    Take the list of tuples (p-value, associated identifier) and return a dict
    with identifier as key and adjusted p_values
    
    
    '''
    # create a list with the p_values
    pvalues = [i[0] for i in p_vals]
    # create a list of identifiers in the same order as their associated p-values in pvalue list
    identifiers = [i[1] for i in p_vals]
    
    # create a dictionnary with index of p_value in list as key and identifier as value
    identifier_position = {}
    for i in range(len(identifiers)):
        identifier_position[i] = identifiers[i]
    
    # create a dictionnary with {pvalue: [list of index]} in list pairs
    pvalue_position = {}
    for i in range(len(pvalues)):
        if pvalues[i] in pvalue_position:
            pvalue_position[pvalues[i]].append(i)
        else:
            pvalue_position[pvalues[i]] = [i]
    
    # sort P-values
    pvalues.sort()    
    
    # create a an array with p_values
    Pval_arr = np.array(pvalues)
    
    # create  a array with corrected p_values 
    corrected_pvalues = Pval_arr * len(Pval_arr) / rankdata(Pval_arr)
    
    # assign to 1 p-values > 1
    corrected_pvalues[corrected_pvalues > 1] = 1

    # create a dict {identifier: adjusted-pvalue}
    ID_pval = {}
    
    # loop over pvalues
    for i in range(len(pvalues)):
        # get the index of pvalue in non-sorted list
        for j in pvalue_position[pvalues[i]]:
            # find the identifier
            # populate dict ID_pval with identifier and adjusted p-value 
            ID_pval[identifier_position[j]] = corrected_pvalues[i]
            
    return ID_pval
   




    
    
