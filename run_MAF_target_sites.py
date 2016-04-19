# -*- coding: utf-8 -*-
"""
Created on Mon Aug 24 14:24:58 2015

@author: RJovelin
"""


import math
import numpy as np
from premature_stops import *
from scipy import stats
from accessories import *
from miRNA_target import *
from sites_with_coverage import *
from Cel_UTR import *
from cel_UTR_length import *
from randomize_SNPs import *

Gstr = lambda x : str(x)

# make a list with MAF of replacement SNPs, excluding sites with sample size < 10
MAF_REP = MAF_SNP('../CREM_CLA_protein_divergence/CDS_SNP_DIVERG.txt', '../CREM_CLA_protein_divergence/unique_transcripts.txt', 
                  '../Premature_stops/transcripts_indels_CDS.txt', 'REP', 10)
print('MAF for replacement sites done')

# make a list with MAF of synonymous sites, excluding sites with sample size < 10
MAF_SYN = MAF_SNP('../CREM_CLA_protein_divergence/CDS_SNP_DIVERG.txt', '../CREM_CLA_protein_divergence/unique_transcripts.txt', 
                  '../Premature_stops/transcripts_indels_CDS.txt', 'SYN', 10)
print('MAF for synonymous sites done')

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


# express frequencies in %
for i in range(len(MAF_REP)):
    MAF_REP[i] = MAF_REP[i] * 100
for i in range(len(MAF_SYN)):
    MAF_SYN[i] = MAF_SYN[i] * 100
for i in range(len(MAF_all_targets)):
    MAF_all_targets[i] = MAF_all_targets[i] * 100
for i in range(len(MAF_crm_targets)):
    MAF_crm_targets[i] = MAF_crm_targets[i] * 100
for i in range(len(MAF_crmcla_targets)):
    MAF_crmcla_targets[i] = MAF_crmcla_targets[i] * 100
for i in range(len(MAF_crmclacel_targets)):
    MAF_crmclacel_targets[i] = MAF_crmclacel_targets[i] * 100
    
print('conversion to % frequencies done')

# make histograms
MAF_REP_hist = np.histogram(MAF_REP, range(0, 51, 10))
MAF_SYN_hist = np.histogram(MAF_SYN, range(0, 51, 10))
MAF_all_targets_hist = np.histogram(MAF_all_targets, range(0, 51, 10))
MAF_crm_targets_hist = np.histogram(MAF_crm_targets, range(0, 51, 10))
MAF_crmcla_targets_hist = np.histogram(MAF_crmcla_targets, range(0, 51, 10))
MAF_crmclacel_targets_hist = np.histogram(MAF_crmclacel_targets, range(0, 51, 10))

print('histograms done')

# compare MAF distributions
diff_MAF_REP_SYN = stats.chi2_contingency([MAF_REP_hist[0], MAF_SYN_hist[0]])
diff_MAF_REP_all = stats.chi2_contingency([MAF_REP_hist[0], MAF_all_targets_hist[0]])
diff_MAF_REP_crm = stats.chi2_contingency([MAF_REP_hist[0], MAF_crm_targets_hist[0]])
diff_MAF_REP_crmcla = stats.chi2_contingency([MAF_REP_hist[0], MAF_crmcla_targets_hist[0]])
diff_MAF_REP_crmclacel = stats.chi2_contingency([MAF_REP_hist[0], MAF_crmclacel_targets_hist[0]])
diff_MAF_SYN_all = stats.chi2_contingency([MAF_SYN_hist[0], MAF_all_targets_hist[0]])
diff_MAF_SYN_crm = stats.chi2_contingency([MAF_SYN_hist[0], MAF_crm_targets_hist[0]])
diff_MAF_SYN_crmcla = stats.chi2_contingency([MAF_SYN_hist[0], MAF_crmcla_targets_hist[0]])
diff_MAF_SYN_crmclacel = stats.chi2_contingency([MAF_SYN_hist[0], MAF_crmclacel_targets_hist[0]])
diff_MAF_all_crm = stats.chi2_contingency([MAF_all_targets_hist[0], MAF_crm_targets_hist[0]])
diff_MAF_all_crmcla = stats.chi2_contingency([MAF_all_targets_hist[0], MAF_crm_targets_hist[0]])
diff_MAF_all_crmclacel = stats.chi2_contingency([MAF_all_targets_hist[0], MAF_crmclacel_targets_hist[0]])
diff_MAF_crm_crmcla = stats.chi2_contingency([MAF_crm_targets_hist[0], MAF_crmcla_targets_hist[0]])
diff_MAF_crm_crmclacel = stats.chi2_contingency([MAF_crm_targets_hist[0], MAF_crmclacel_targets_hist[0]])
diff_MAF_crmcla_crmclacel = stats.chi2_contingency([MAF_crmcla_targets_hist[0], MAF_crmclacel_targets_hist[0]])

print('chi2 MAF comparisons done')

# make a list of test results
test_results = [diff_MAF_REP_SYN, diff_MAF_REP_all, diff_MAF_REP_crm,
                diff_MAF_REP_crmcla, diff_MAF_REP_crmclacel, diff_MAF_SYN_all,
                diff_MAF_SYN_crm, diff_MAF_SYN_crmcla, diff_MAF_SYN_crmclacel,
                diff_MAF_all_crm, diff_MAF_all_crmcla, diff_MAF_all_crmclacel,
                diff_MAF_crm_crmcla, diff_MAF_crm_crmclacel, diff_MAF_crmcla_crmclacel]

# make a list of strings for pairwise comparisons
pairwise_comp = ['REP_SYN', 'REP_all', 'REP_crm', 'REP_crmcla', 'REP_crmclacel',
                 'SYN_all', 'SYN_crm', 'SYN_crmcla', 'SYN_crmclacel', 'all_crm',
                 'all_crmcla', 'all_crmclacel', 'crm_crmcla', 'crm_crmclacel',
                 'crmcla_crmclacel']

# open summary file
summary_file = open('summary_MAF_miRNA_targets.txt', 'w')

# write results to summary file
summary_file.write('comparison of the MAF distributions\n')
summary_file.write('-' * 36 + '\n')
summary_file.write('distributions' + '\t' + 'chi2' + '\t' + 'p-val' + '\t' + 'dof' + '\n')

# loop over list of string comparisons
# write comp to file and write results of chi2 to file
for i in range(len(pairwise_comp)):
    summary_file.write(pairwise_comp[i] + '\t')
    content = '\t'.join([str(test_results[i][j]) for j in range(3)])
    summary_file.write(content + '\n')

print('chi2 results of MAF differences written to file')

# write tables to file from histograms
# open file to write the MAF frequencies of the different SNPs
newfile = open('MAF_SNPs_miRNA_targets.txt', 'w')
newfile.write('\t'.join(['MAF', 'REP', 'SYN', 'all', 'crm', 'crm-cla', 'crm-cla-cel']) + '\n')

print(len(MAF_REP_hist[0]))
print(len(MAF_REP_hist[1]))

for i in range(len(MAF_REP_hist[0])):
    newfile.write(str(MAF_REP_hist[1][i]) + ':' + str(MAF_REP_hist[1][i] + 10) + '\t' +
    str(MAF_REP_hist[0][i]) + '\t' + str(MAF_SYN_hist[0][i]) + '\t' + str(MAF_all_targets_hist[0][i]) + '\t' +
    str(MAF_crm_targets_hist[0][i]) + '\t' + str(MAF_crmcla_targets_hist[0][i]) + '\t' + str(MAF_crmclacel_targets_hist[0][i]) + '\n')
newfile.close()

print('MAF distribution written to file')


# compute pairwise difference between average MAF

# make a list of MAF lists
MAF_list = [MAF_REP, MAF_SYN, MAF_all_targets, MAF_crm_targets, MAF_crmcla_targets, MAF_crmclacel_targets]
# make a list of MAF names
MAF_names = ['REP', 'SYN', 'all', 'crm', 'crm-cla', 'crm-cla-cel']

summary_file.write('\n')
summary_file.write('comparison of mean MAF differences\n')
summary_file.write('-' * 35 + '\n')
summary_file.write('distributions' + '\t' + 'wilcoxon' + '\t' + 'p-val' + '\n')

# loop over MAF_list, comparre average differences
for i in range(0, len(MAF_list) - 1):
    for j in range(i+1, len(MAF_list)):
        # test mean differences    
        wilcoxon, p_val = stats.ranksums(MAF_list[i], MAF_list[j])
        summary_file.write(MAF_names[i] + '_' + MAF_names[j] + '\t' + str(wilcoxon) + '\t' + str(p_val) + '\n')

summary_file.write('\n')
summary_file.write('Mean MAF at sites\n')
summary_file.write('-' * 18 + '\n')
summary_file.write('sites' + '\t' + 'mean_MAF' + '\t' + 'SEM' + '\n')
for i in range(len(MAF_list)):
    summary_file.write('\t'.join([MAF_names[i], str(np.mean(MAF_list[i])), str(np.std(MAF_list[i]) / math.sqrt(len(MAF_list[i])))]) + '\n') 

# close summary file
summary_file.close()















