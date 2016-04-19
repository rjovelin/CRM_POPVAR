# -*- coding: utf-8 -*-
"""
Created on Wed Aug  5 09:53:48 2015

@author: RJovelin
"""

import math
import numpy as np


# use this function to compute theta in a sliding window in a single sequence 
def sequence_sliding_window(chromo, start, end, orientation, upstream, window_size, step, SNP_sites, threshold):
    '''
    (str, int, int, str, bool, int, int, dict, int) -> dict
    Take a the chromo, the start and end position of a given sequence feature
    on chromo (0-based indices), the feature orientation, a boolean indicating if the focal
    region is upstream of the feature, the window size and the incremental step, 
    a dictionnary with all site position (0-based indices) with coverage in the genome
    : list of allele counts pairs and the threshold of missing site to accept per window
    Preconditions: start and end are 0-based indices, and the dictionary
    of SNP/no-SNP for each site has been cleaned (ie. minimum sample size applied)
    '''
    
    if upstream == True:
        if orientation == '+':
            orientation = '-'
        elif orientation == '-':
            orientation = '+'
        
    # create a dict with {window position: theta}     
    theta_windows = {}
    
    # initiate key
    j = -1
    
    # make a list of indices for each site in the sequence
    if orientation  == '+':
        positions = [i for i in range(start, end - window_size + step, step)]
    elif orientation == '-':
        positions = [i for i in range(end -1, start + window_size - step + 1, -step)]
    
    # loop over the positions
    for i in positions:
        # get the indices in the window
        if orientation == '+':
            w = [k for k in range(i, i + window_size, 1)]
        elif orientation == '-':
            w = [k for k in range(i, i - window_size, -1)]
        # get the position of the window
        j += 1
        # compute the number of polymorphic sites
        polym = 0
        # make a list to store the sample size of polymorphic site
        sample_size = []
        # count the number of sites
        N = 0
        # loop over indices in window
        for k in w:
            # check if site in dict
            if k in SNP_sites[chromo]:
                # adjust site counter
                N += 1
                # get reference, aternative alleles and their counts
                ref, alt, ref_count, alt_count = SNP_sites[chromo][k][0], SNP_sites[chromo][k][1], SNP_sites[chromo][k][2], SNP_sites[chromo][k][3]
                # check if site polymorphism or not
                # need to check counts instead of alleles, because a monoallelic
                # site in KSR+PX may be polymorphic in the PB strains
                # and so ref and alt will be different
                if ref_count != 0 and alt_count != 0:
                    # site is polymorphic
                    # double check
                    assert ref != alt, 'counts of both allele different than 0 but ref and alt are the same'
                    # adjust polymorophic site counter
                    polym += 1
                    # add sample size of polymorphic site to list
                    sample_size.append(ref_count + alt_count)
        # check that a maximum of threshold sites are missing in window
        if N >= window_size - threshold:
            # compute theta
            theta = polym / N
            if theta != 0:
                # apply correction for varying sample size
                theta = theta / math.log(np.mean(sample_size) - 1)
            # populate dict with window position : theta pairs
            theta_windows[j] = [theta]
        else:
            # populate dict with window position : empty list
            # (to keep track of undefined windows)
            theta_windows[j] = []
            
    return theta_windows
        
            
        
             
    