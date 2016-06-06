# -*- coding: utf-8 -*-
"""
Created on Sun Jun  5 23:31:07 2016

@author: Richard
"""





# use this script to plot the Pn/Ps ratio and Dn/Ds ratio for chemo and non chemo genes
# and the distribution of these 2 ratios per gene


# use Agg backend on server without X server
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib import rc
rc('mathtext', default='regular')
# import modules
import numpy as np
from scipy import stats
import math
import os
import sys
import random
# import custom modules
from chemoreceptors import *
from manipulate_sequences import *
from mk_test import *
from premature_stops import *


# create a set of valid transcripts
transcripts = get_valid_transcripts('../Genome_Files/unique_transcripts.txt')
print('got list of valid transcripts')

# get the set of chemoreceptors from the iprscan outputfile
chemo = get_chemoreceptors('../Genome_Files//PX356_protein_seq.tsv') 
print('parsed chemo genes')

# create a set of valid chemoreceptors
GPCRs = set(gene for gene in chemo if gene in transcripts)
print('GPCRS', len(GPCRs))
print('made set of valid chemo genes')

# create a set of non-chemo gene
NonGPCRs = set(gene for gene in transcripts if gene not in chemo)
print('NonGPCRs', len(NonGPCRs))
print('made set of valid non-chemo genes')








# use this funtion to find the 5' most upstream PTC SNPs
def position_first_PTC(snp_file, GeneList, indel_transcripts, CDS_fasta):




# make a list of position of 5' most upstream PTC SNP relative to CDS length
upstream_PTC = position_first_PTC('../CREM_CLA_protein_divergence/CDS_SNP_DIVERG.txt', '../CREM_CLA_protein_divergence/unique_transcripts.txt', 'transcripts_indels_CDS.txt', '../CREM_CLA_protein_divergence/noamb_PX356_all_CDS.fasta')
# express relative position in % of the CDS length
for i in range(len(upstream_PTC)):
    upstream_PTC[i] = upstream_PTC[i] * 100
# make a histogram
upstream_hist = np.histogram(upstream_PTC, range(0, 101, 10))
# open file to write the CDS truncations
newfile = open('CDS_truncation.txt', 'w')
newfile.write('PTC_relative_position' + '\t' + 'PTC_counts' + '\n')
for i in range(len(upstream_hist[0])):
    newfile.write(str(upstream_hist[1][i]) + ':' + str(upstream_hist[1][i] + 10) + '\t' + str(upstream_hist[0][i]) + '\n')
newfile.close()



# test uniformity of the distribution of the 5' most upstream PTC
test_uniform_upstream = stats.chisquare(upstream_hist[0])
# write results of test to summary file
summary_file.write('\n')
summary_file.write('test of uniformity of the distribution of 5\' most upstream PTC\n') 
summary_file.write('-' * 63 + '\n')
summary_file.write('chi2' + '\t' + 'p-val' + '\n')
summary_file.write(str(test_uniform_upstream[0]) + '\t' + str(test_uniform_upstream[1]) + '\n')





# compare expression level between PTC and non-PTC genes
# use all PTC genes, including genes with indels

# get the set of genes with PTC
PTC_genes = count_genes_with_indels_premature_stops('../CREM_CLA_protein_divergence/CDS_SNP_DIVERG.txt', '../CREM_CLA_protein_divergence/unique_transcripts.txt')
# get the set of valid transcripts
transcripts = get_valid_transcripts('../CREM_CLA_protein_divergence/unique_transcripts.txt')
# create a set of non-PTC genes
non_PTC_genes = set(gene for gene in transcripts if gene not in PTC_genes)
# create a dict of {gene: expression leve;}
crm_expression = expression_developmental_stages('../CREM_CLA_protein_divergence/WBPaper00041190.cre.mr.csv', 'c_remanei.PRJNA53967.WS248.geneIDs.txt')
# create lists with expression level for PTC and non-PTC genes
PTC_expression = [crm_expression[gene] for gene in PTC_genes if gene in crm_expression]
non_PTC_expression = [crm_expression[gene] for gene in non_PTC_genes if gene in crm_expression]
# test differences in gene expression between PTC and non-PTC genes
wilcoxon, p_val = stats.ranksums(PTC_expression, non_PTC_expression)
# write mean expression, SEM and differences to file

summary_file.write('\n')
summary_file.write('Expression level (includes PTC genes with indels):\n')
summary_file.write('-' * 51 + '\n')
summary_file.write('\t'.join(['genes', 'N', 'mean_expression', 'SEM']) + '\n')
# write the sample size, mean expression, SEM
summary_file.write('\t'.join(['PTC_genes', str(len(PTC_expression)), str(np.mean(PTC_expression)), str(np.std(PTC_expression) / math.sqrt(len(PTC_expression)))]) + '\n')
summary_file.write('\t'.join(['non_PTC_genes', str(len(non_PTC_expression)), str(np.mean(non_PTC_expression)), str(np.std(non_PTC_expression) / math.sqrt(len(non_PTC_expression)))]) + '\n')
# write results of statistical test
summary_file.write('Wilcoxon_sum_rank_test:' + '\t' + str(wilcoxon) + '\t' + str(p_val) + '\n')





# use this funtion to get the relative position of all PTC SNPs, including genes with alleles
def RelativePositionPTC(snp_file, GeneList, CDS_fasta):




