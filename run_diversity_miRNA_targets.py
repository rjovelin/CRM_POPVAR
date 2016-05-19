# -*- coding: utf-8 -*-
"""
Created on Sat Aug 15 23:09:46 2015

@author: Richard
"""

# use this script to plot a graph of nucleotide diversity at miRNA target sites and other types of sites


# use Agg backend on server without X server
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from manipulate_sequences import *
from sliding_windows import *
from sites_with_coverage import *
from miRNA_target import *
from divergence import *
from genomic_coordinates import *
import numpy as np
from scipy import stats
import math


# convert genome fasta to dict
genome = convert_fasta('../Genome_Files/noamb_356_v1_4.txt')
print('converted genome fasta to dict')

# get the allele counts for all sites with coverage, exclude sites with sample size < 10
chromo_sites = get_non_coding_snps('../SNP_files/', 10) 
print('got allele counts at all sites')

# get mirna target coordinates 
# {gene: [[chromo, start, end, orientation, seed, N_mirnas, site_type, conservation, utr]]}
target_coord = get_miRNA_target_coord('Cremanei_miRNA_sites.txt')
print('got miRNA target coordinates')

# get the set of valid transcripts
valid_transcripts = get_valid_transcripts('../Genome_Files/unique_transcripts.txt')
print('made list of valid transcripts')

# remove genes that are not not valid (ie to keep a single transcript per gene)
to_remove = [gene for gene in target_coord if gene not in valid_transcripts]
for gene in to_remove:
    del target_coord[gene]
print('deleted {0} non-valid genes'.format(len(to_remove)))
print('filtered non-valid transcripts') 

# compute theta at synonymous sites
theta_syn = compute_theta_diversity('../Genome_Files/CDS_SNP_DIVERG.txt', valid_transcripts, 'SYN', 10)
# make a list of theta values
SYN_theta = []
for gene in theta_syn:
    SYN_theta.append(theta_syn[gene])
print('sites\tmean\tmin\tmax\t')
print('SYN', np.mean(SYN_theta), min(SYN_theta), max(SYN_theta))
print('done computing theta at synonymous sites')

# compute theta at nonsynonymous sites
theta_rep = compute_theta_diversity('../Genome_Files/CDS_SNP_DIVERG.txt', valid_transcripts, 'REP', 10)
# make a list of theta values
REP_theta = []
for gene in theta_rep:
    REP_theta.append(theta_rep[gene])
print('sites\tmean\tmin\tmax\t')
print('REP', np.mean(REP_theta), min(REP_theta), max(REP_theta))
print('done computing theta at replacement sites')


# compute diversity for all mirna targets
# create a list theta at target sites
targets_theta =  [] 
# loop over gene in target cood
for gene in target_coord:
    # loop over targets in gene
    for i in range(len(target_coord[gene])):
        # get chromo
        chromo = target_coord[gene][i][0]
        # get start, end position
        start = target_coord[gene][i][1]
        end = target_coord[gene][i][2]
        # get conservation
        conservation = target_coord[gene][i][7]
        # get utr type
        utr = target_coord[gene][i][8]
        # check that chromo is in chromo_sites
        if chromo in chromo_sites:
            # compute theta, accepting 0 missing sites
            theta = compute_theta_non_coding(chromo_sites, chromo, start, end, 0)
            # check that theta is defined
            if theta != 'NA':
                # add to list 
                targets_theta.append(theta)

# compare diversity between targetsm REP and SYN
# make a list of theta lists
alldata = [REP_theta, SYN_theta, targets_theta]
# make a liust of site categories
site_types = ['Rep', 'Syn', 'targets']

# write header in file
summary_file.write('comparison of mean theta between protein coding genes and miRNA target sites\n')
summary_file.write('-' * 77 + '\n')
summary_file.write('\n')
summary_file.write('sites' + '\t' + 'N' + '\t' + 'mean_theta' + '\t' + 'SEM' + '\n')

for i in range(len(site_theta)):
    # compute mean and SEM
    summary_file.write('\t'.join([site_names[i], str(len(site_theta[i])), str(np.mean(site_theta[i])), str(np.std(site_theta[i]) / math.sqrt(len(site_theta[i])))]) + '\n')

summary_file.write('\n')

summary_file.write('sites' + '\t' + 'wilcoxon' + '\t' + 'P-val' + '\n')

for i in range(len(site_theta)-1):
    for j in range(i+1, len(site_theta)):
        # compare mean differences
        wilcoxon, p = stats.ranksums(site_theta[i], site_theta[j])
        # write results to file
        summary_file.write('\t'.join([site_names[i] + '_vs_' + site_names[j], str(wilcoxon), str(p)]) + '\n')
        
summary_file.write('\n')

# open file to store theta values
newfile = open('theta_values_all_target_sites.txt', 'w')
for i in range(len(site_theta)):
    for theta in site_theta[i]:
        newfile.write(site_names[i] + '\t' + str(theta) + '\n')
# close file
newfile.close()        

        
# make a list of theta for different site categories
theta_cons = [crm_targets_theta, crmcla_targets_theta, crmclacel_targets_theta]
cons_names = ['crm-only', 'crm-cla', 'crm-cla-cel']

summary_file.write('targets sites with varying conservation level\n')
summary_file.write('-' *46 + '\n')
summary_file.write('\n')
summary_file.write('sites' + '\t' + 'N' + '\t' + 'mean_theta' + '\t' + 'SEM' + '\n')

for i in range(len(theta_cons)):
    # compute mean and sem
   summary_file.write('\t'.join([cons_names[i], str(len(theta_cons[i])), str(np.mean(theta_cons[i])), str(np.std(theta_cons[i]) / math.sqrt(len(theta_cons[i])))]) + '\n')

summary_file.write('\n')

summary_file.write('sites' + '\t' + 'wilcoxon' + '\t' + 'P-val' + '\n')

for i in range(len(theta_cons)-1):
    for j in range(i+1, len(theta_cons)):
        wilcoxon, p = stats.ranksums(theta_cons[i], theta_cons[j])
        summary_file.write('\t'.join([cons_names[i] + '_vs_' + cons_names[j], str(wilcoxon), str(p)]) + '\n')
        
summary_file.write('\n')

# open file to store theta values
newfile = open('theta_target_sites_conservation.txt', 'w')
for i in range(len(theta_cons)):
    for theta in theta_cons[i]:
        newfile.write(cons_names[i] + '\t' + str(theta) + '\n')
# close file
newfile.close()        

# compare theta for target sites from different downstream sequences
summary_file.write('targets sites with annotated UTRs\n')
summary_file.write('-' * 34 + '\n')
summary_file.write('\n')
summary_file.write('sites' + '\t' + 'N' + '\t' + 'mean_theta' + '\t' + 'SEM' + '\n')

summary_file.write('\t'.join(['downtream', str(len(nonUTR_targets_theta)), str(np.mean(nonUTR_targets_theta)), str(np.std(nonUTR_targets_theta) / math.sqrt(len(nonUTR_targets_theta)))]) + '\n')
summary_file.write('\t'.join(['UTR', str(len(UTR_targets_theta)), str(np.mean(UTR_targets_theta)), str(np.std(UTR_targets_theta) / math.sqrt(len(UTR_targets_theta)))]) + '\n')

summary_file.write('\n')
summary_file.write('sites' + '\t' + 'wilcoxon' + '\t' + 'P-val' + '\n')
wilcoxon, p = stats.ranksums(nonUTR_targets_theta, UTR_targets_theta)
summary_file.write('downtream_vs_UTR' + '\t' + str(wilcoxon) + '\t' + str(p) + '\n')

# open file to store the theta value os downand UTR sequences
newfile = open('theta_target_sites_UTR.txt', 'w')
for theta in nonUTR_targets_theta:
    newfile.write('downstream' + '\t' + str(theta) + '\n')
for theta in UTR_targets_theta:
    newfile.write('UTR' + '\t' + str(theta) + '\n')
newfile.close()    

print('theta values written to files')

# compare theta and flanking sites

# create a dict for upstream and downstream flanking sites
# {window_position : [list of thetas across all regions at that position]}
upstream_pos, downstream_pos = {}, {}

# compute theta for each non-overlapping windows upstream of target sites

# loop over targets
# compute theta theta for each upstream region
for gene in target_coord:
    # loop over target on gene
    for i in range(len(target_coord[gene])):
        # get target coordinates
        chromo = target_coord[gene][i][0]
        start = target_coord[gene][i][1]
        end = target_coord[gene][i][2]
        orientation = target_coord[gene][i][3]
        # check that chromo in chromo-sites
        if chromo in chromo_sites:
            # determine the length of the target sites
            target_size = end - start
            upstream_size = target_size * 3
            # compute upstream coordinates
            # take 3 windows on each side of the target
            if orientation == '+':
                up_end = start
                if start - upstream_size < 0:
                    up_start = 0
                else:
                    up_start = start - upstream_size
            elif orientation == '-':
                up_start = end
                if up_start + upstream_size < len(genome[chromo]):
                    up_end = up_start + upstream_size
                else:
                    up_end = up_start + upstream_size
            # compute theta for each non-verlapping window
            theta_windows = sequence_sliding_window(chromo, up_start, up_end, orientation, True, target_size, target_size, chromo_sites, 0)
            # add theta at each position in the upstream dict
            for j in theta_windows:
                # check if j in dict
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
            
# loop over targets
for gene in target_coord:
    # loop over target on gene
    for i in range(len(target_coord[gene])):
        # get target coordinates
        chromo = target_coord[gene][i][0]
        start = target_coord[gene][i][1]
        end = target_coord[gene][i][2]
        orientation = target_coord[gene][i][3]
        # check that chromo in chromo-sites
        if chromo in chromo_sites:
            # determine the length of the target sites
            target_size = end - start
            downstream_size = target_size * 3
            # compute downstream coordinates
            # take 3 windows on each side of the target
            if orientation == '+':
                down_start = end
                if down_start + downstream_size < len(genome[chromo]):
                    down_end = down_start + downstream_size
                else:
                    down_end = len(genome[chromo])
            elif orientation == '-':
                down_end = start
                if start - downstream_size < 0:
                    down_start = 0
                else:
                    down_start = start - downstream_size
            theta_windows = sequence_sliding_window(chromo, down_start, down_end, orientation, False, target_size, target_size, chromo_sites, 0)
            # add theta at each position in the downstream dict
            for j in theta_windows:
                # check if j is key
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

summary_file.write('\n')

summary_file.write('Mean theta at all target sites and flanking sites\n')
summary_file.write('-' * 50 + '\n')
summary_file.write('sites' + '\t' +  'N' + '\t' + 'mean_theta' + '\t' + 'SEM' + '\n')

# create lists of lists to store the mean theta at each position 
# [sample size, mean, stderror]
    
# lower index in upstream_pos corresponds to end the upstream region
# need to invert positions
# create a list of keys in upstream
up_pos = [i for i in upstream_pos]
# sort keys
up_pos.sort()
# reverse sort keys
up_pos.reverse()
    
# loop over reverse soreted keys in up_pos
for i in up_pos:
    # compute the mean theta at that position
    mean_theta = np.mean(upstream_pos[i])
    # compute standard error
    stderror = np.std(upstream_pos[i]) / math.sqrt(len(upstream_pos[i]))
    summary_file.write('\t'.join(['upstream_' + str(i), str(len(upstream_pos[i])), str(mean_theta), str(stderror)]) + '\n')
    
# write mean and SEM for targets
summary_file.write('\t'.join(['targets', str(len(targets_theta)), str(np.mean(targets_theta)), str(np.std(targets_theta) / math.sqrt(len(targets_theta)))]) + '\n')   
    
# create a list of keys in downstream
down_pos = [i for i in downstream_pos]
# sort list
down_pos.sort()

# loop over sorted keys in down_pos
for i in down_pos:
    # compute the mean theta at that position
    mean_theta = np.mean(downstream_pos[i])
    # compyte the standard error
    stderror = np.std(downstream_pos[i]) / math.sqrt(len(downstream_pos[i]))
    summary_file.write('\t'.join(['downstream_' + str(i), str(len(downstream_pos[i])), str(mean_theta), str(stderror)]) + '\n')

# save theta value for targets and flaning sites to a single file
newfile = open('theta_target_flanking_sites.txt', 'w')
for i in up_pos:
    for theta in upstream_pos[i]:
        newfile.write('upstream_' + str(i) + '\t' + str(theta) + '\n')
for theta in targets_theta:
    newfile.write('targets' + '\t' + str(theta) + '\n')
for i in down_pos:
    for theta in downstream_pos[i]:
        newfile.write('downstream_' + str(i) + '\t' + str(theta) + '\n')
newfile.close()

print('saved theta values to file')

# test mean differences between theta at target sites and flanking sites

summary_file.write('\n')
summary_file.write('Mean differences between all targets and flanking sites\n')
summary_file.write('-' * 56 + '\n')
summary_file.write('sites' + '\t' + 'wilcoxon' + 'P-val' + '\n')
for i in up_pos:
    wilcoxon, p = stats.ranksums(upstream_pos[i], targets_theta)
    summary_file.write('\t'.join(['upstream_' + str(i) + '_vs_targets', str(wilcoxon), str(p)]) + '\n')
for i in down_pos:
    wilcoxon, p = stats.ranksums(downstream_pos[i], targets_theta)
    summary_file.write('\t'.join(['downstream_' + str(i) + '_vs_targets', str(wilcoxon), str(p)]) + '\n')

print('tested mean differences between targets and flanking sites')

summary_file.close()

