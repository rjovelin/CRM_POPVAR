# -*- coding: utf-8 -*-
"""
Created on Thu Aug  6 22:00:54 2015

@author: Richard
"""

from sites_with_coverage import *
from sliding_windows import *
from accessories import *
import numpy as np
import math
from scipy import stats


# convert genome fasta file to dict
genome = convert_fasta('../CREM_CLA_protein_divergence/noamb_356_v1_4.txt')

print('fasta converted to dict')

# create a dictionary with all sites with coverage in the genome
# consider only site with a minimum sample size of 10
chromo_sites = get_non_coding_snps('../SNP_files/', 10)

print('got sites with coverage')
print(len(chromo_sites))

# get the mirnas coordinates
mirnas_coord = {}
# openfile for reading
infile = open('crm_miRBase21_premina_coordinates.txt', 'r')
infile.readline()
for line in infile:
    line = line.rstrip()
    if line != '':
        line = line.split()
        mir = line[0]
        chromo = line[1]
        start = int(line[2]) -1
        end = int(line[3])
        orientation = line[4]
        mirnas_coord[mir] = [chromo, start, end, orientation]
infile.close()

print('got miRNAs coordinates')

# create a dict for upstream downstream and mirnas regions
# {window_position : [list of thetas across all regions at that position]}
upstream_pos, mirna_pos, downstream_pos = {}, {}, {}

# loop over chromo
# compute theta for each mirna
# populate the mirna dict with window positions : list of theta
for mir in mirnas_coord:
    # get coordinates
    chromo = mirnas_coord[mir][0]
    start = mirnas_coord[mir][1]
    end = mirnas_coord[mir][2]
    orientation = mirnas_coord[mir][3]
    # compute the sliding window of theta for each pirna
    theta_windows = sequence_sliding_window(chromo, start, end, orientation, False, 10, 3, chromo_sites, 2)
    # add theta at each position in the mirna dict
    for j in theta_windows:
        # check if j is key in dict
        if j in mirna_pos:
            # add theta if theta is defined
            if len(theta_windows[j]) != 0:
                mirna_pos[j].append(theta_windows[j][0])
        else:
            # initiate list value and add theta if theta is defined
            mirna_pos[j] = []
            if len(theta_windows[j]) != 0:
                mirna_pos[j].append(theta_windows[j][0])
                    
print('computed sliding windows for pirnas')

# compute sliding window theta for each downstream region
for mir in mirnas_coord:
    # loop over each mirna 
    # get pirna coordinates
    chromo = mirnas_coord[mir][0]
    start = mirnas_coord[mir][1]
    end = mirnas_coord[mir][2]
    orientation = mirnas_coord[mir][3]
    # compute the sliding window of theta for each 500 bp of downstream region
    # compute downstream coordinates
    # check orientation
    if orientation == '+':
        down_start = end
        if down_start + 500 < len(genome[chromo]):
            down_end = down_start + 500
        else:
            down_end = len(genome[chromo])
    elif orientation == '-':
        down_end = start
        if start - 500 < 0:
            down_start = 0
        else:
            down_start = start - 500
    # compute the sliding windows of theta for each downstream seq
    theta_windows = sequence_sliding_window(chromo, down_start, down_end, orientation, False, 10, 3, chromo_sites, 2)
    # add theta at each position in the downstream dict
    for j in theta_windows:
        # check if j is key in dict
        if j in downstream_pos:
            # add theta if theta is defined
            if len(theta_windows[j]) != 0:
                downstream_pos[j].append(theta_windows[j][0])
        else:
            # initiate list value and add theta if theta is defined
            downstream_pos[j] = []
            if len(theta_windows[j]) != 0:
                downstream_pos[j].append(theta_windows[j][0])

print('computed sliding windows for downstream')


# loop over mir
# compute sliding window theta for each upstream region
for mir in mirnas_coord:
    # get pirna coordinates
    chromo = mirnas_coord[mir][0]
    start = mirnas_coord[mir][1]
    end = mirnas_coord[mir][2]
    orientation = mirnas_coord[mir][3]
    # compute the upstream coordinates
    # check orientation
    if orientation == '+':
        up_end = start
        if start - 500 < 0:
            up_start = 0
        else:
            up_start = start - 500
    elif orientation == '-':
        up_start = end
        if up_start + 500 < len(genome[chromo]):
            up_end = up_start + 500
        else:
            up_end = len(genome[chromo])
    # compute the sliding windows of theta for each upstream seq
    theta_windows = sequence_sliding_window(chromo, up_start, up_end, orientation, True, 10, 3, chromo_sites, 2)
    # add theta at each position in the upstream dict
    for j in theta_windows:
        # check if j in key in dict
        if j in upstream_pos:
            # add theta if theta is defined
            if len(theta_windows[j]) != 0:
                upstream_pos[j].append(theta_windows[j][0])
        else:
            # initiate list value and add theta if theta is defined
            upstream_pos[j] = []
            if len(theta_windows[j]) != 0:
                upstream_pos[j].append(theta_windows[j][0])
    
print('computed sliding windows for upstream')                           
    
# create lists of lists to store the mean theta at each position with 95% CI
# [sample size, mean, stderror, std, low_CI, high_CI]
    
# lower index in upstream_pos corresponds to end the upstream region
# need to invert positions
# create a list of keys in upstream
up_pos = [i for i in upstream_pos]
# sort keys
up_pos.sort()
# reverse sort keys
up_pos.reverse()
# create list of lists
up_theta = []
    
    
# loop over reverse soreted keys in up_pos
for i in up_pos:
    # compute the mean theta at that position
    mean_theta = np.mean(upstream_pos[i])
    # compute standard error
    stderror = np.std(upstream_pos[i]) / math.sqrt(len(upstream_pos[i]))
    # compute margin error (critical value = 1.96 for a 95% CI)
    margin = 1.96 * stderror
    lCI = mean_theta - margin
    hCI = mean_theta + margin
    up_theta.append([len(upstream_pos[i]), mean_theta, stderror, np.std(upstream_pos[i]), lCI, hCI])
        
# create list of keys in mirnas
mir_pos = [i for i in mirna_pos]
# sort list
mir_pos.sort()
# create list of lists
mir_theta = []
# loop over the sorted keys in pi_pos
for i in mir_pos:
    # compute the mean theta at that position
    mean_theta = np.mean(mirna_pos[i])
    # compute the standard error
    stderror = np.std(mirna_pos[i]) / math.sqrt(len(mirna_pos[i]))
    # compute the margin error (critical value = 1.96 for a 95% CI)
    margin = 1.96 * stderror
    lCI = mean_theta - margin
    hCI = mean_theta + margin
    mir_theta.append([len(mirna_pos[i]), mean_theta, stderror, np.std(mirna_pos[i]), lCI, hCI])
        
        
# create a list of keys in downstream
down_pos = [i for i in downstream_pos]
# sort list
down_pos.sort()
# create list of lists
down_theta = []
# loop over sorted keys in down_pos
for i in down_pos:
    # compute the mean theta at that position
    mean_theta = np.mean(downstream_pos[i])
    # compyte the standard error
    stderror = np.std(downstream_pos[i]) / math.sqrt(len(downstream_pos[i]))
    # compute the margin error (critical value = 1.96 for a 95% CI)
    margin = 1.96 * stderror
    lCI = mean_theta - margin
    hCI = mean_theta + margin
    down_theta.append([len(downstream_pos[i]), mean_theta, stderror, np.std(downstream_pos[i]), lCI, hCI])
        
        
# create lambda function
Gstr = lambda x: str(x)
    
# open file for writing
newfile = open('miRNAs_sliding_windows_theta.txt', 'w')
# write header
newfile.write('\t'.join(['window_number', 'window_position_region', 'region', 'N', 'mean_theta', 'SEM', 'STD', 'low_95%_CI', 'high_95%_CI']) + '\n')
    
# set up window counter
j = 0
    
# loop over values in up_theta
for i in range(len(up_theta)):
    # update window counter        
    j += 1
    # write to file
    newfile.write(str(j) + '\t' + str(i) + '\t' + 'upstream' + '\t')
    newfile.write('\t'.join(list(map(Gstr, up_theta[i]))) + '\n')
# loop over values in pi_theta
for i in range(len(mir_theta)):
    # update window counter
    j += 1
    # write to file
    newfile.write(str(j) + '\t' + str(i) + '\t' + 'miRNA' + '\t')
    newfile.write('\t'.join(list(map(Gstr, mir_theta[i]))) + '\n')
# loop over values in down_theta
for i in range(len(down_theta)):
    # update window counter
    j += 1
    # write to file
    newfile.write(str(j) + '\t' + str(i) + '\t' + 'downstream' + '\t')
    newfile.write('\t'.join(list(map(Gstr, down_theta[i]))) + '\n')
        
        
# close file after writing
newfile.close()
            
    
        
        
        
        
        
        
        
        
        
        
    
    