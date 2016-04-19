# -*- coding: utf-8 -*-
"""
Created on Tue Aug 11 16:04:01 2015

@author: RJovelin
"""

import numpy as np
import math



# use this function to compute theta per site for a sequence of interest
# does not consider synonymous or replacement sites
def compute_theta_non_coding(chromo_sites, chromo, start, end, threshold):
    '''
    (dict, str, int, int, int) -> float
    Take the dictionary with allele counts in the genome, a chromo where the 
    sequence of interest is located, the coordinates of the sequence, and the
    maximum number of accptable missing sites and return theta per site
    adjusted for varying sample size among sites
    Precondition: all positions in the dict and start and end are 0-based    
    '''
    
    # get the positions of the sequence of interest
    positions = [i for i in range(start, end)]
    
    # count the number of polymorphic sites
    diff = 0
    
    # create a list to store the sample size of the polymorphic sites    
    sample_size = []
    
    # count the number of sites
    N = 0
    
    # loop over indices
    for i in positions:
        # check if site in chromo_sites
        if i in chromo_sites[chromo]:
            # adjust site counter
            N += 1
            # get reference, alternative allele and their counts
            ref, alt = chromo_sites[chromo][i][0], chromo_sites[chromo][i][1]
            ref_count, alt_count = chromo_sites[chromo][i][2], chromo_sites[chromo][i][3]
            # check if site is polymorphic
            if ref_count != 0 and alt_count != 0:
                # site is polymorphich
                # double check
                assert ref != alt, 'counts of both allele different than 0 but ref and alt are the same'
                # adjust polym site counter
                diff += 1
                # add sample size of polymorphich site to list
                sample_size.append(ref_count + alt_count)
    # check that a maximum of threshold site are missing
    if N >= len(positions) - threshold:
        # compute theta per site
        theta = diff / N
        if theta != 0:
            # apply a correction for varying sample size
            theta = theta / math.log(np.mean(sample_size) -1)
        
        return theta
        
        
