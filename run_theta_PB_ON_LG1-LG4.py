# -*- coding: utf-8 -*-
"""
Created on Fri Sep 11 21:10:55 2015

@author: Richard
"""



from accessories import *
from piRNAs import *
from miRNA_target import *
from repeats_TEs import *
from sliding_windows import *
from sites_with_coverage import *
from divergence import *
import numpy as np
from scipy import stats
import math


# convert genome fasta to dict
genome = convert_fasta('noamb_356_v1_4.txt')
print('genome converted to fasta dict')

# compute theta per 50 Kb window in all strains
# get the allele counts for all sites with coverage, keep all sites 
chromo_sites = get_all_strains_snps('../PB_ON_SNP_files/', 0)
print('got allele counts at all sites')

# create a list of starting positions of each window on chromo
LG1_positions = [i for i in range(0, len(genome['linkage_group_1']), 50000)]
LG2_positions = [i for i in range(0, len(genome['linkage_group_2']), 50000)]
LG4_positions = [i for i in range(0, len(genome['linkage_group_4']), 50000)]

# open file for writing
newfile = open('LG1_theta_PB+ON_50Kb.txt', 'w')
# loop over start position in LG1:
for i in LG1_positions:
    # compute theta per window
    # use a large threshold > window length to accept any number of missing site  
    theta = compute_theta_non_coding(chromo_sites, 'linkage_group_1', i, i + 50000, 100000)
    newfile.write(str(i) + '\t' + str(theta) + '\n')
# close file
newfile.close()
print('theta per windows saved to file for LG1')

# open file for writing
newfile = open('LG2_theta_PB+ON_50Kb.txt', 'w')
# loop over start positions in LG2
for i in LG2_positions:
    # compute theta per window
    # use a large threshold > window length to accept any number of missing site 
    theta = compute_theta_non_coding(chromo_sites, 'linkage_group_2', i, i + 50000, 100000)
    newfile.write(str(i) + '\t' + str(theta) + '\n')
# close file
newfile.close()
print('theta per windows saved to file for LG2')

# open file for writing
newfile = open('LG4_theta_PB+ON_50Kb.txt', 'w')
# loop over start positions in LG4:
for i in LG4_positions:
    # compute theta per window
    # use a large threshold > window length to accept any number of missing site  
    theta = compute_theta_non_coding(chromo_sites, 'linkage_group_4', i, i + 50000, 100000)
    newfile.write(str(i) + '\t' + str(theta) + '\n')
# close file
newfile.close()
print('theta per windows saved to file for LG4')




