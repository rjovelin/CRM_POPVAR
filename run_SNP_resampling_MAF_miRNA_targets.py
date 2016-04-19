# -*- coding: utf-8 -*-
"""
Created on Mon Aug 24 22:53:20 2015

@author: Richard
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


# get the coordinates of all target sites
all_target_coord = get_miRNA_target_loci('Cremanei_miRNA_sites.txt', '../CREM_CLA_protein_divergence/unique_transcripts.txt', 'all')
# get the coordinates of remanei-specific targets
crm_target_coord = get_miRNA_target_loci('Cremanei_miRNA_sites.txt', '../CREM_CLA_protein_divergence/unique_transcripts.txt', 'crm')
# get the coordinates of the remanei-latens conserved targets
crmcla_target_coord = get_miRNA_target_loci('Cremanei_miRNA_sites.txt', '../CREM_CLA_protein_divergence/unique_transcripts.txt', 'crm-cla')
# get the coordinates of the remanei-latens-elegans targets
crmclacel_target_coord = get_miRNA_target_loci('Cremanei_miRNA_sites.txt', '../CREM_CLA_protein_divergence/unique_transcripts.txt', 'crm-cla-cel')
print('got miRNA target site coordinates')


# get the allele counts for all targets
all_target_sites = get_feature_sites(chromo_sites, all_target_coord)
# get allele counts for the remanei-specific targets
crm_target_sites = get_feature_sites(chromo_sites, crm_target_coord)
# get allele counts for the remanei-latens conserved targets
crmcla_target_sites = get_feature_sites(chromo_sites, crmcla_target_coord)
# get allele counts for the remanei-latens-elegans conserved sites
crmclacel_target_sites = get_feature_sites(chromo_sites, crmclacel_target_coord)
print('got allele counts for target sites')


# compute MAF for all targets
MAF_all_targets = MAF_non_coding(all_target_sites)
# compute MAF for remanei-specific targets
MAF_crm_targets = MAF_non_coding(crm_target_sites)
# compute MAF for remanei-latens conserved targets
MAF_crmcla_targets = MAF_non_coding(crmcla_target_sites)
# compute MAF for remanei-latens-elegans conserved targets
MAF_crmclacel_targets = MAF_non_coding(crmclacel_target_sites)
print('MAF for miRNA targets done')


# get SNPs in non-miRNA target UTRs
UTR_snps = get_UTR_SNPs(chromo_sites, '../CREM_CLA_protein_divergence/356_10172014.gff3', 'c_elegans.PRJNA13758.WS248.annotations.gff3', '../CREM_CLA_protein_divergence/noamb_356_v1_4.txt', 99, '../CREM_CLA_protein_divergence/unique_transcripts.txt', 'Cremanei_miRNA_sites.txt')
print('got allele counts for non-target UTR sites')

# get the total number of snps in non-target UTRs
non_target_snps = 0
for chromo in UTR_snps:
    non_target_snps += len(UTR_snps[chromo])
print('SNPs at non-target UTR sites: ', non_target_snps)


# resample SNPs and compute MAF
UTR_resampled_MAF = SNP_MAF_randomization(UTR_snps, 10000, 1000)

# get the proportions of SNPs in each MAF bin
UTR_MAF_proportions = get_MAF_distribution_from_replicates(UTR_resampled_MAF)

# get the SNP proportions for the observed SNPs
empirical_all_MAF = SNP_proportions_MAF_bin(MAF_all_targets)
empirical_crm_MAF = SNP_proportions_MAF_bin(MAF_crm_targets)
empirical_crmcla_MAF = SNP_proportions_MAF_bin(MAF_crmcla_targets)
empirical_crmclacel_MAF = SNP_proportions_MAF_bin(MAF_crmclacel_targets)

print(len(empirical_all_MAF))
print(len(empirical_crm_MAF))
print(len(empirical_crmcla_MAF))
print(len(empirical_crmclacel_MAF))

# creat a list of keys, being the MAF lower bound in the dicts with MAF proportions
maf_limit = [i for i in empirical_all_MAF]
# sort list
maf_limit.sort()

# open file to save results of Z-test
newfile = open('summary_MAF_resampled_SNPs_miRNA_targets.txt', 'w')

newfile.write('Z-test of SNP proportions in MAF bins\n')
newfile.write('random resampling of 10000 SNPs among {0} non-target UTR SNPs with 1000 replicates\n'.format(non_target_snps)) 


# make a list of empirical MAF 
empirical = [empirical_all_MAF, empirical_crm_MAF, empirical_crmcla_MAF, empirical_crmclacel_MAF]
# make a list of target names
names = ['all_targets', 'crm-specific-targets', 'crm-cla-conserved-targets', 'crm-cla-cel-conserved-targets']


# loop over names
for j in range(len(names)):
    newfile.write('\n')
    newfile.write('MAF for ' + names[j] + '\n')
    newfile.write('\t'.join(['MAF', 'mean_sample', 'margin', 'stdev_sample', 'l95', 'h95', 'observed', 'z-score', 'alpha_0.05', 'alpha_0.01', 'alpha_0.001']) + '\n')
    # loop over sorted keys
    for i in maf_limit:
        # compute mean
        average = np.mean(UTR_MAF_proportions[i])
        # get standard deviation
        stdev = np.std(UTR_MAF_proportions[i])
        # compute 95% CI
        stderror = stdev / math.sqrt(len(UTR_MAF_proportions[i]))
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
    
