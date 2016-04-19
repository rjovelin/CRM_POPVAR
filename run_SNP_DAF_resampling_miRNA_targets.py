# -*- coding: utf-8 -*-
"""
Created on Thu Aug 27 12:21:18 2015

@author: RJovelin
"""


import math
import numpy as np
from premature_stops import *
from scipy import stats
from accessories import *
from miRNA_target import *
from piRNAs import *
from sites_with_coverage import *
from Cel_UTR import *
from cel_UTR_length import *
from randomize_SNPs import *
from get_coding_sequences import *


Gstr = lambda x : str(x)

# get the allele counts for all sites with coverage, exclude sites with sample size < 10
chromo_sites = get_non_coding_snps('../SNP_files/', 10)
print('got allele counts in genome')

# get conservation level of all remanei targets (excluding targets of non-valid transcripts)
conservation_scores = get_conservation_score('Cremanei_miRNA_sites.txt', '../CREM_CLA_protein_divergence/unique_transcripts.txt')
print('got conservation level of miRNA target sites')

# get DAF all targets
DAF_all_targets = get_DAF_miRNA_targets(chromo_sites, '../CREM_CLA_protein_divergence/noamb_356_v1_4.txt', 'Cremanei_Clatens_miRNA_sites.txt', 'Crm_Cla_UTR_sequences/', conservation_scores, 'all', '../CREM_CLA_protein_divergence/unique_transcripts.txt')
# get DAF remanei-specific targets
DAF_crm_targets = get_DAF_miRNA_targets(chromo_sites, '../CREM_CLA_protein_divergence/noamb_356_v1_4.txt', 'Cremanei_Clatens_miRNA_sites.txt', 'Crm_Cla_UTR_sequences/', conservation_scores, 'crm', '../CREM_CLA_protein_divergence/unique_transcripts.txt')
# get DAF remanei-latens conserved targets
DAF_crmcla_targets = get_DAF_miRNA_targets(chromo_sites, '../CREM_CLA_protein_divergence/noamb_356_v1_4.txt', 'Cremanei_Clatens_miRNA_sites.txt', 'Crm_Cla_UTR_sequences/', conservation_scores, 'crm-cla', '../CREM_CLA_protein_divergence/unique_transcripts.txt')
# get DAF remanei-latens-elegans targets
DAF_crmclacel_targets = get_DAF_miRNA_targets(chromo_sites, '../CREM_CLA_protein_divergence/noamb_356_v1_4.txt', 'Cremanei_Clatens_miRNA_sites.txt', 'Crm_Cla_UTR_sequences/', conservation_scores, 'crm-cla-cel', '../CREM_CLA_protein_divergence/unique_transcripts.txt')
print('got DAF for miRNA targets')

# get DAF of sites in UTR
DAF_UTR = get_DAF_UTR_non_target(chromo_sites, 'Cremanei_Clatens_miRNA_sites.txt', 'Crm_Cla_UTR_sequences/', '../CREM_CLA_protein_divergence/noamb_356_v1_4.txt')
print('got DAF of UTR sites')

# get the positions of miRNA target sites on each chromo
target_pos = get_all_miRNA_target_sites('Cremanei_miRNA_sites.txt')
print('got positions of miRNA target sites on each chromo')

# remove UTR positions corresponding to miRNA target sites
for chromo in target_pos:
    if chromo in DAF_UTR:
        to_remove = [i for i in target_pos[chromo]]
        for i in to_remove:
            if i in DAF_UTR[chromo]:
                del DAF_UTR[chromo][i]
# remove chromo without positions
to_remove = []
for chromo in DAF_UTR:
    if len(DAF_UTR[chromo]) == 0:
        to_remove.append(chromo)
if len(to_remove) != 0:
    for chromo in to_remove:
        del DAF_UTR[chromo]
print('got DAF in UTR exluding miRNA targets')                

     
# get the total number of snps in non-target UTRs
non_target_snps = 0
for chromo in DAF_UTR:
    non_target_snps += len(DAF_UTR[chromo])
print('SNPs at non-target UTR sites: ', non_target_snps)

# resample SNPs 
UTR_resampled_DAF = SNP_DAF_randomization(DAF_UTR, 10000, 1000)

# get the proportions of SNPs in each MAF bin
UTR_DAF_proportions = get_DAF_distribution_from_replicates(UTR_resampled_DAF)

# get the SNP proportions for the observed SNPs in each DAF bin
empirical_all_DAF = SNP_proportions_DAF_bin(DAF_all_targets)
empirical_crm_DAF = SNP_proportions_DAF_bin(DAF_crm_targets)
empirical_crmcla_DAF = SNP_proportions_DAF_bin(DAF_crmcla_targets)
empirical_crmclacel_DAF = SNP_proportions_DAF_bin(DAF_crmclacel_targets)

print(len(empirical_all_DAF))
print(len(empirical_crm_DAF))
print(len(empirical_crmcla_DAF))
print(len(empirical_crmclacel_DAF))

# creat a list of keys, being the DAF lower bound in the dicts with DAF proportions
daf_limit = [i for i in empirical_all_DAF]
# sort list
daf_limit.sort()

# open file to save results of Z-test
newfile = open('summary_DAF_resampled_SNPs_miRNA_targets.txt', 'w')

newfile.write('Z-test of SNP proportions in DAF bins\n')
newfile.write('random resampling of 10000 SNPs among {0} non-target UTR SNPs with 1000 replicates\n'.format(non_target_snps)) 

# make a list of empirical MAF 
empirical = [empirical_all_DAF, empirical_crm_DAF, empirical_crmcla_DAF, empirical_crmclacel_DAF]
# make a list of target names
names = ['all_targets', 'crm-specific-targets', 'crm-cla-conserved-targets', 'crm-cla-cel-conserved-targets']


# loop over names
for j in range(len(names)):
    newfile.write('\n')
    newfile.write('DAF for ' + names[j] + '\n')
    newfile.write('\t'.join(['DAF', 'mean_sample', 'margin', 'stdev_sample', 'l95', 'h95', 'observed', 'z-score', 'alpha_0.05', 'alpha_0.01', 'alpha_0.001']) + '\n')
    # loop over sorted keys
    for i in daf_limit:
        # compute mean
        average = np.mean(UTR_DAF_proportions[i])
        # get standard deviation
        stdev = np.std(UTR_DAF_proportions[i])
        # compute 95% CI
        stderror = stdev / math.sqrt(len(UTR_DAF_proportions[i]))
        # compute the margin error (critical value = 1.96 for 95% CI)
        margin = 1.96 * stderror
        lCI = average - margin
        hCI = average + margin
        observed = empirical[j][i]
        z_score = (observed - average) / stdev
        # critical values for 1-sample 2-tailed z-test: 0.05: +- 1.96, 0.01: +- 2.58, 0.001: +-3.27
        # Ho : obsvered == mean, H1: observed != mean
        if z_score < -1.96 or z_score > 1.96:
            P5 = '*'
        elif -1.96 <= z_score <= 1.96:
            P5 = 'NS'
        if z_score < -2.58 or z_score > 2.58:
            P1 = '*'
        elif -2.58 <= z_score <= 2.58:
            P1 = 'NS'
        if z_score < -3.27 or z_score > 3.27:
            P01 = '*'
        elif -3.27 <= z_score <= 3.27:
            P01 = 'NS'
        newfile.write('\t'.join([str(i), str(average), str(margin), str(stdev), str(lCI), str(hCI), str(observed), str(z_score), str(P5), str(P1), str(P01)]) + '\n')
    
# close file after writing
newfile.close()
    
