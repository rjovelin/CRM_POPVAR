# -*- coding: utf-8 -*-
"""
Created on Tue Aug 25 10:45:56 2015

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


# get the allele counts for all sites with coverage, exclude sites with sample size < 10
chromo_sites = get_non_coding_snps('../SNP_files/', 10)
print('got allele counts in genome')

# make a list with DAF of replacement SNPs, excluding sites with sample size , 10
DAF_REP = DAF_SNP('../CREM_CLA_protein_divergence/CDS_SNP_DIVERG.txt', '../CREM_CLA_protein_divergence/unique_transcripts.txt', '../Premature_stops/transcripts_indels_CDS.txt', 'REP', 10)
print('got DAF for replacement SNPs')
# make a list with DAF of synonymous sites, excluding sites with sample size < 10
DAF_SYN = DAF_SNP('../CREM_CLA_protein_divergence/CDS_SNP_DIVERG.txt', '../CREM_CLA_protein_divergence/unique_transcripts.txt', '../Premature_stops/transcripts_indels_CDS.txt', 'SYN', 10)
print('got DAF for synonymous SNPs')


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


print('number of snps in:')
print('all_targets: ', len(DAF_all_targets))
print('crm_targets: ', len(DAF_crm_targets))
print('crmcla_targets: ', len(DAF_crmcla_targets))
print('crmclacel_targets: ', len(DAF_crmclacel_targets))


# express frequencies in %
for i in range(len(DAF_REP)):
    DAF_REP[i] = DAF_REP[i] * 100
for i in range(len(DAF_SYN)):
    DAF_SYN[i] = DAF_SYN[i] * 100
for i in range(len(DAF_all_targets)):
    DAF_all_targets[i] = DAF_all_targets[i] * 100
for i in range(len(DAF_crm_targets)):
    DAF_crm_targets[i] = DAF_crm_targets[i] * 100
for i in range(len(DAF_crmcla_targets)):
    DAF_crmcla_targets[i] = DAF_crmcla_targets[i] * 100
for i in range(len(DAF_crmclacel_targets)):
    DAF_crmclacel_targets[i] = DAF_crmclacel_targets[i] * 100
print('done converting frequencies to %')    

# make histograms
DAF_REP_hist = np.histogram(DAF_REP, range(0, 101, 10))
DAF_SYN_hist = np.histogram(DAF_SYN, range(0, 101, 10))
DAF_all_hist = np.histogram(DAF_all_targets, range(0, 101, 10))
DAF_crm_hist = np.histogram(DAF_crm_targets, range(0, 101, 10))
DAF_crmcla_hist = np.histogram(DAF_crmcla_targets, range(0, 101, 10))
DAF_crmclacel_hist = np.histogram(DAF_crmclacel_targets, range(0, 101, 10))
print('done making histograms')

# compare DAF distributions
diff_DAF_REP_SYN = stats.chi2_contingency([DAF_REP_hist[0], DAF_SYN_hist[0]])
diff_DAF_REP_all = stats.chi2_contingency([DAF_REP_hist[0], DAF_all_hist[0]])
diff_DAF_REP_crm = stats.chi2_contingency([DAF_REP_hist[0], DAF_crm_hist[0]])
diff_DAF_REP_crmcla = stats.chi2_contingency([DAF_REP_hist[0], DAF_crmcla_hist[0]])
diff_DAF_REP_crmclacel = stats.chi2_contingency([DAF_REP_hist[0], DAF_crmclacel_hist[0]])
diff_DAF_SYN_all = stats.chi2_contingency([DAF_SYN_hist[0], DAF_all_hist[0]])
diff_DAF_SYN_crm = stats.chi2_contingency([DAF_SYN_hist[0], DAF_crm_hist[0]])
diff_DAF_SYN_crmcla = stats.chi2_contingency([DAF_SYN_hist[0], DAF_crmcla_hist[0]])
diff_DAF_SYN_crmclacel = stats.chi2_contingency([DAF_SYN_hist[0], DAF_crmclacel_hist[0]])
diff_DAF_all_crm = stats.chi2_contingency([DAF_all_hist[0], DAF_crm_hist[0]])
diff_DAF_all_crmcla = stats.chi2_contingency([DAF_all_hist[0], DAF_crmcla_hist[0]])
diff_DAF_all_crmclacel = stats.chi2_contingency([DAF_all_hist[0], DAF_crmclacel_hist[0]])
diff_DAF_crm_crmcla = stats.chi2_contingency([DAF_crm_hist[0], DAF_crmcla_hist[0]])
diff_DAF_crm_crmclacel = stats.chi2_contingency([DAF_crm_hist[0], DAF_crmclacel_hist[0]])
diff_DAF_crmcla_crmclacel = stats.chi2_contingency([DAF_crmcla_hist[0], DAF_crmclacel_hist[0]])
print('done comparing DAF distributions with chi2')


# make a list of test results
test_results = [diff_DAF_REP_SYN, diff_DAF_REP_all, diff_DAF_REP_crm,
                diff_DAF_REP_crmcla, diff_DAF_REP_crmclacel, diff_DAF_SYN_all,
                diff_DAF_SYN_crm, diff_DAF_SYN_crmcla, diff_DAF_SYN_crmclacel,
                diff_DAF_all_crm, diff_DAF_all_crmcla, diff_DAF_all_crmclacel,
                diff_DAF_crm_crmcla, diff_DAF_crm_crmclacel, diff_DAF_crmcla_crmclacel]

# make a list of strings for pairwise comparisons
pairwise_comp = ['REP_SYN', 'REP_all', 'REP_crm', 'REP_crmcla', 'REP_crmclacel',
                 'SYN_all', 'SYN_crm', 'SYN_crmcla', 'SYN_crmclacel', 'all_crm',
                 'all_crmcla', 'all_crmclacel', 'crm_crmcla', 'crm_crmclacel',
                 'crmcla_crmclacel']

# open summary file
summary_file = open('summary_DAF_miRNA_targets.txt', 'w')

# write results to summary file
summary_file.write('comparison of the DAF distributions\n')
summary_file.write('-' * 36 + '\n')
summary_file.write('distributions' + '\t' + 'chi2' + '\t' + 'p-val' + '\t' + 'dof' + '\n')

# loop over list of string comparisons
# write comp to file and write results of chi2 to file
for i in range(len(pairwise_comp)):
    summary_file.write(pairwise_comp[i] + '\t')
    content = '\t'.join([str(test_results[i][j]) for j in range(3)])
    summary_file.write(content + '\n')

print('chi2 results of DAF differences written to file')

# write tables to file from histograms
# open file to write the MAF frequencies of the different SNPs
newfile = open('DAF_SNPs_miRNA_targets.txt', 'w')
newfile.write('\t'.join(['DAF', 'REP', 'SYN', 'all', 'crm', 'crm-cla', 'crm-cla-cel']) + '\n')

print(len(DAF_REP_hist[0]))
print(len(DAF_REP_hist[1]))

for i in range(len(DAF_REP_hist[0])):
    newfile.write(str(DAF_REP_hist[1][i]) + ':' + str(DAF_REP_hist[1][i] + 10) + '\t' +
    str(DAF_REP_hist[0][i]) + '\t' + str(DAF_SYN_hist[0][i]) + '\t' + str(DAF_all_hist[0][i]) + '\t' +
    str(DAF_crm_hist[0][i]) + '\t' + str(DAF_crmcla_hist[0][i]) + '\t' + str(DAF_crmclacel_hist[0][i]) + '\n')
newfile.close()

print('DAF distribution written to file')


# compute pairwise difference between average DAF
# make a list of DAF lists
DAF_list = [DAF_REP, DAF_SYN, DAF_all_targets, DAF_crm_targets, DAF_crmcla_targets, DAF_crmclacel_targets]
# make a list of MAF names
DAF_names = ['REP', 'SYN', 'all', 'crm', 'crm-cla', 'crm-cla-cel']

summary_file.write('\n')
summary_file.write('comparison of mean DAF differences\n')
summary_file.write('-' * 35 + '\n')
summary_file.write('distributions' + '\t' + 'wilcoxon' + '\t' + 'p-val' + '\n')

# loop over DAF_list, compare average differences
for i in range(0, len(DAF_list) - 1):
    for j in range(i+1, len(DAF_list)):
        # test mean differences    
        wilcoxon, p_val = stats.ranksums(DAF_list[i], DAF_list[j])
        summary_file.write(DAF_names[i] + '_' + DAF_names[j] + '\t' + str(wilcoxon) + '\t' + str(p_val) + '\n')

summary_file.write('\n')
summary_file.write('Mean DAF at sites\n')
summary_file.write('-' * 18 + '\n')
summary_file.write('sites' + '\t' + 'mean_DAF' + '\t' + 'SEM' + '\n')
for i in range(len(DAF_list)):
    summary_file.write('\t'.join([DAF_names[i], str(np.mean(DAF_list[i])), str(np.std(DAF_list[i]) / math.sqrt(len(DAF_list[i])))]) + '\n') 

# close summary file
summary_file.close()

