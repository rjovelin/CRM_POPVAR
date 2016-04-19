# -*- coding: utf-8 -*-
"""
Created on Sun Jun 14 22:48:10 2015

@author: Richard
"""

import math
import numpy as np
from premature_stops import *
from scipy import stats
from accessories import *
from crm_expression import *
from caeno_chromosomes import *
from miRNA_target import *


# write summary to file
summary_file = open('summary_PTC_SNPs.txt', 'w')


summary_file.write('genes with PTC SNPs, including genes with indels:\n')
summary_file.write('-' * 50 + '\n')

# count the number of genes PTC SNPs, including genes with indels
all_PTC_genes = count_genes_with_indels_premature_stops('../CREM_CLA_protein_divergence/CDS_SNP_DIVERG.txt', '../CREM_CLA_protein_divergence/unique_transcripts.txt')
summary_file.write('total number of genes with PTC:' + '\t' + str(len(all_PTC_genes)) + '\n')

# count the total number of PTC SNPs
total_SNPs = count_total_PTC_SNP('../CREM_CLA_protein_divergence/CDS_SNP_DIVERG.txt', '../CREM_CLA_protein_divergence/unique_transcripts.txt')
summary_file.write('total number of PTC SNPs:' + '\t' + str(total_SNPs) + '\n')

# count the number of genes with PTC SNPs excluding genes with indels
PTC_genes = genes_with_premature_stops('../CREM_CLA_protein_divergence/CDS_SNP_DIVERG.txt', '../CREM_CLA_protein_divergence/unique_transcripts.txt', 'transcripts_indels_CDS.txt')

# create lists to count the number of premature stop codons and PTC SNPs
PTC_codons = []
PTC_SNPs = []
for gene in PTC_genes:
    PTC_SNPs.append(PTC_genes[gene][0])
    PTC_codons.append(PTC_genes[gene][1])

summary_file.write('\n')
summary_file.write('genes with SNPs-induced PTCs, exlcuding genes with indels:\n')
summary_file.write('-' * 59 + '\n')
summary_file.write('number of genes with PTCs:' + '\t' + str(len(PTC_genes)) + '\n')
summary_file.write('number of SNPs causing PTCs:' + '\t' + str(sum(PTC_SNPs)) + '\n')
summary_file.write('number of premature stop codons:' + '\t' + str(sum(PTC_codons)) + '\n')

# get the number of PTC codons per gene for genes with PTC
mean_codons = np.mean(PTC_codons)
# get the mean number of PTC SNPs per gene for genes with PTC
mean_SNPs = np.mean(PTC_SNPs)
# compute standard error of the mean PTC codons
SEM_codons = np.std(PTC_codons) / math.sqrt(len(PTC_codons))
# compute standard error of the mean PTC SNPs
SEM_SNPs = np.std(PTC_SNPs) / math.sqrt(len(PTC_SNPs))

summary_file.write('mean +- SEM (PTC codons / gene for genes with PTC):' + '\t' + str(mean_codons) + '\t' + '+-' + '\t' + str(SEM_codons) + '\n')
summary_file.write('mean += SEM (PTC SNPs / gene for genes with PTC):' + '\t' + str(mean_SNPs) + '\t' + '+-' + '\t' + str(SEM_SNPs) + '\n')

# make a list with the MAF of PTC SNPs, excluding genes with indels, and excluding site with sample size < 10
MAF_PTC = MAF_SNP('../CREM_CLA_protein_divergence/CDS_SNP_DIVERG.txt', '../CREM_CLA_protein_divergence/unique_transcripts.txt', 'transcripts_indels_CDS.txt', 'PTC', 10)

# make a list with MAF of replacement SNPs, excluding sites with sample size < 10
MAF_REP = MAF_SNP('../CREM_CLA_protein_divergence/CDS_SNP_DIVERG.txt', '../CREM_CLA_protein_divergence/unique_transcripts.txt', 'transcripts_indels_CDS.txt', 'REP', 10)

# make a list with MAF of synonymous sites, excluding sites with sample size , 10
MAF_SYN = MAF_SNP('../CREM_CLA_protein_divergence/CDS_SNP_DIVERG.txt', '../CREM_CLA_protein_divergence/unique_transcripts.txt', 'transcripts_indels_CDS.txt', 'SYN', 10)


# express frequencies in %
for i in range(len(MAF_PTC)):
    MAF_PTC[i] = MAF_PTC[i] * 100
for i in range(len(MAF_REP)):
    MAF_REP[i] = MAF_REP[i] * 100
for i in range(len(MAF_SYN)):
    MAF_SYN[i] = MAF_SYN[i] * 100
    

# make histograms
MAF_PTC_hist = np.histogram(MAF_PTC, range(0, 51, 10))
MAF_REP_hist = np.histogram(MAF_REP, range(0, 51, 10))
MAF_SYN_hist = np.histogram(MAF_SYN, range(0, 51, 10))

# compare MAF distributions
diff_MAF_PTC_REP = stats.chi2_contingency([MAF_PTC_hist[0], MAF_REP_hist[0]])
diff_MAF_PTC_SYN = stats.chi2_contingency([MAF_PTC_hist[0], MAF_SYN_hist[0]])
diff_MAF_REP_SYN = stats.chi2_contingency([MAF_REP_hist[0], MAF_SYN_hist[0]])

# write results to summary file
summary_file.write('\n')
summary_file.write('comparison of the MAF distributions\n')
summary_file.write('-' * 36 + '\n')
summary_file.write('distributions' + '\t' + 'chi2' + '\t' + 'p-val' + '\t' + 'dof' + '\n')
summary_file.write('PTC_vs_REP' + '\t')
content = '\t'.join([str(diff_MAF_PTC_REP[i]) for i in range(3)])
summary_file.write(content + '\n')

summary_file.write('PTC_vs_SYN' +'\t')
content = '\t'.join([str(diff_MAF_PTC_SYN[i]) for i in range(3)])
summary_file.write(content + '\n')

summary_file.write('REP_vs_SYN' + '\t')
content = '\t'.join([str(diff_MAF_REP_SYN[i]) for i in range(3)])
summary_file.write(content + '\n')


# make a list with the DAF of PTC SNPs, excluding genes with indels, excluding sites with sample size , 10
DAF_PTC = DAF_SNP('../CREM_CLA_protein_divergence/CDS_SNP_DIVERG.txt', '../CREM_CLA_protein_divergence/unique_transcripts.txt', 'transcripts_indels_CDS.txt', 'PTC', 10)

# make a list with DAF of replacement SNPs, excluding sites with sample size , 10
DAF_REP = DAF_SNP('../CREM_CLA_protein_divergence/CDS_SNP_DIVERG.txt', '../CREM_CLA_protein_divergence/unique_transcripts.txt', 'transcripts_indels_CDS.txt', 'REP', 10)

# make a list with DAF of synonymous sites, excluding sites with sample size < 10
DAF_SYN = DAF_SNP('../CREM_CLA_protein_divergence/CDS_SNP_DIVERG.txt', '../CREM_CLA_protein_divergence/unique_transcripts.txt', 'transcripts_indels_CDS.txt', 'SYN', 10)


# express frequencies in %
for i in range(len(DAF_PTC)):
    DAF_PTC[i] = DAF_PTC[i] * 100
for i in range(len(DAF_REP)):
    DAF_REP[i] = DAF_REP[i] * 100
for i in range(len(DAF_SYN)):
    DAF_SYN[i] = DAF_SYN[i] * 100
    
# make histograms
DAF_PTC_hist = np.histogram(DAF_PTC, range(0, 101, 10))
DAF_REP_hist = np.histogram(DAF_REP, range(0, 101, 10))
DAF_SYN_hist = np.histogram(DAF_SYN, range(0, 101, 10))

# compare DAF distributions
diff_DAF_PTC_REP = stats.chi2_contingency([DAF_PTC_hist[0], DAF_REP_hist[0]])
diff_DAF_PTC_SYN = stats.chi2_contingency([DAF_PTC_hist[0], DAF_SYN_hist[0]])
diff_DAF_REP_SYN = stats.chi2_contingency([DAF_REP_hist[0], DAF_SYN_hist[0]])

# write results to summary file
summary_file.write('\n')
summary_file.write('comparison of the DAF distributions\n')
summary_file.write('-' * 36 + '\n')
summary_file.write('distributions' + '\t' + 'chi2' + '\t' + 'p-val' + '\t' + 'dof' + '\n')
summary_file.write('PTC_vs_REP' + '\t')
content = '\t'.join([str(diff_DAF_PTC_REP[i]) for i in range(3)])
summary_file.write(content + '\n')

summary_file.write('PTC_vs_SYN' +'\t')
content = '\t'.join([str(diff_DAF_PTC_SYN[i]) for i in range(3)])
summary_file.write(content + '\n')

summary_file.write('REP_vs_SYN' + '\t')
content = '\t'.join([str(diff_DAF_REP_SYN[i]) for i in range(3)])
summary_file.write(content + '\n')


# write tables to file from histograms
# open file to write the MAF frequencies of the different SNPs
newfile = open('MAF_SNPs_CDS.txt', 'w')
newfile.write('MAF' + '\t' + 'PTC' + '\t' + 'REP' + '\t' + 'SYN' + '\n')
for i in range(len(MAF_PTC_hist[0])):
    newfile.write(str(MAF_PTC_hist[1][i]) + ':' + str(MAF_PTC_hist[1][i] + 10) + '\t' + str(MAF_PTC_hist[0][i]) + '\t' + str(MAF_REP_hist[0][i]) + '\t' + str(MAF_SYN_hist[0][i]) + '\n')
newfile.close()

newfile = open('DAF_SNPs_CDS.txt', 'w')
newfile.write('DAF' +'\t' + 'PTC' + '\t' + 'REP' + '\t' + 'SYN' + '\n')
for i in range(len(DAF_PTC_hist[0])):
    newfile.write(str(DAF_PTC_hist[1][i]) + ':' + str(DAF_PTC_hist[1][i] + 10) + '\t' + str(DAF_PTC_hist[0][i]) + '\t' + str(DAF_REP_hist[0][i]) + '\t' + str(DAF_SYN_hist[0][i]) + '\n')
newfile.close()


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

# get the mean allelic frequencies of the 5' most upstream PTC SNP in each decile of CDS length
PTC_frequency = freq_first_PTC('../CREM_CLA_protein_divergence/CDS_SNP_DIVERG.txt', '../CREM_CLA_protein_divergence/unique_transcripts.txt', 'transcripts_indels_CDS.txt', '../CREM_CLA_protein_divergence/noamb_PX356_all_CDS.fasta')

# write table to file with mean and SEM in each decile length
newfile = open('PTC_average_frequencies.txt', 'w')
newfile.write('\t'.join(['CDS_length_decile', 'average_allele_frequency', 'SEM_freq']) + '\n')
# compute mean and SEM
j = 0
for i in range(len(PTC_frequency)):
    mean_freq = np.mean(PTC_frequency[i])
    SEM_freq = np.std(PTC_frequency[i]) / math.sqrt(len(PTC_frequency[i]))
    newfile.write(str(j) + ':' + str(j + 10) + '\t' + str(mean_freq) + '\t' + str(SEM_freq) + '\n')
    j += 10
newfile.close()


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


# compare distribution of PTC genes osn the X and autosomes

# generate lists of PTC X-linked and autosomal genes
PTC_X, PTC_autosomal = X_autosomal_genes(PTC_genes, '../CREM_CLA_protein_divergence/356_10172014.gff3', '../CREM_CLA_protein_divergence/scaffold_chromosome.csv')
non_PTC_X, non_PTC_autosomal = X_autosomal_genes(non_PTC_genes, '../CREM_CLA_protein_divergence/356_10172014.gff3', '../CREM_CLA_protein_divergence/scaffold_chromosome.csv')     
# get the p_value of the Chi2 test
chi2, p_val_chi2, dof, expected = stats.chi2_contingency([[len(PTC_X), len(PTC_autosomal)], [len(non_PTC_X), len(non_PTC_autosomal)]])
# get the p_value of fisher's exact test
p_val_fisher = stats.fisher_exact([[len(PTC_X), len(PTC_autosomal)], [len(non_PTC_X), len(non_PTC_autosomal)]])[1]

# write results to file
summary_file.write('\n')
summary_file.write('Chromosomal bias for PTC genes:\n')
summary_file.write('-' * 32 + '\n')
summary_file.write('\t' + '\t' + 'PTC' + '\t' + 'non_PTC' + '\n')
summary_file.write('X_linked' + '\t' + str(len(PTC_X)) + '\t' + str(len(non_PTC_X)) + '\n')
summary_file.write('autosomal' + '\t' + str(len(PTC_autosomal)) + '\t' + str(len(non_PTC_autosomal)) + '\n')
# write results of statistical test
# ask which test to report
a = [PTC_X, PTC_autosomal,non_PTC_X, non_PTC_autosomal]
use_fisher = False
for i in range(len(a)):
    if len(a[i]) <=10:
        # use fisher test
        use_fisher = True
        break
if use_fisher  == True:
    summary_file.write('fisher_exact_test:' + '\t' + str(p_val_fisher) + '\n')
else:
    summary_file.write('Chi2_contingency_test:' + '\t' + str(chi2) + '\t' + str(p_val_chi2) + '\t' + str(dof) + '\n')


# compare sex-bias expression between PTC and non-PTC genes

# create a dict with expression in males and females {gene : (female, male, p-val)} 
sex_bias_expression = expression_male_female('../CREM_CLA_protein_divergence/GSE41367_remanei_gns.expr.txt')
# create lists of female and male biased expressed genes
PTC_male, PTC_female, non_PTC_male, non_PTC_female = [], [], [], []
# loop over genes in sex_bias_expression
for gene in sex_bias_expression:
    # record only sex-biased expressed genes
    if sex_bias_expression[gene][-1] < 0.05:
        # gene is sex-biased at alpha = 5%
        # determine if gene is male bias or female bias
        if sex_bias_expression[gene][0] > sex_bias_expression[gene][1]:
            # gene is female biased
            if gene in PTC_genes:
                PTC_female.append(gene)
            elif gene in non_PTC_genes:
                non_PTC_female.append(gene)
        elif sex_bias_expression[gene][0] < sex_bias_expression[gene][1]:
            # gene is male biased
            if gene in PTC_genes:
                PTC_male.append(gene)
            elif gene in non_PTC_genes:
                non_PTC_male.append(gene)
# write results to file
summary_file.write('\n')
summary_file.write('Sex-expression bias for PTC genes:\n')
summary_file.write('-' * 35 + '\n')
summary_file.write('\t' + '\t' + 'PTC' + '\t' + 'non_PTC' + '\n')
summary_file.write('male_biased' + '\t' + str(len(PTC_male)) + '\t' + str(len(non_PTC_male)) + '\n')
summary_file.write('female_biased' + '\t' + str(len(PTC_female)) + '\t' + str(len(non_PTC_female)) + '\n')
# get the p_value of the Chi2 test
chi2, p_val_chi2, dof, expected = stats.chi2_contingency([[len(PTC_male), len(PTC_female)], [len(non_PTC_male), len(non_PTC_female)]])
# get the p_value of fisher's exact test
p_val_fisher = stats.fisher_exact([[len(PTC_male), len(PTC_female)], [len(non_PTC_male), len(non_PTC_female)]])[1]
# write results of statistical test
# ask which test to report
a = [PTC_male, PTC_female, non_PTC_male, non_PTC_female]
use_fisher = False
for i in range(len(a)):
    if len(a[i]) <=10:
        # use fisher test
        use_fisher = True
        break
if use_fisher  == True:
    summary_file.write('fisher_exact_test:' + '\t' + str(p_val_fisher) + '\n')
else:
    summary_file.write('Chi2_contingency_test:' + '\t' + str(chi2) + '\t' + str(p_val_chi2) + '\t' + str(dof) + '\n')


# compare nucleotide divergence between PTC and non-PTC genes

# create a dict to store nucleotide divergene for each gene
divergence = {}
infile = open('../CREM_CLA_protein_divergence/CRM_CLA_prot_diverg_filtered.txt', 'r')
infile.readline()
for line in infile:
    if line.rstrip() != '':
        line = line.rstrip().split()
        gene = line[0]
        dN = float(line[4])
        dS = float(line[5])
        if line[6] != 'NA':
            omega = float(line[6])
        else:
            omega = line[6]
        divergence[gene] = [dN, dS, omega]
infile.close()

# set up counter variable
i = 0

# write header
summary_file.write('\n')
summary_file.write('Nucleotide divergence in PTC and non-PTC genes:\n')
summary_file.write('-' * 48 + '\n')
summary_file.write('\t'.join(['\t', 'N_PTC', 'PTC_mean', 'PTC_SEM', 'N_non_PTC', 'non_PTC_mean', 'non_PTC_SEM', 'wilcoxon', 'p_val']) + '\n')

# loop over each gene, make lists of divergence for PTC and non_PTC genes for each nucleotide variable
while i != 3:
    # create list of divergence ofr PTC and non_PTC genes
    PTC_diverg, non_PTC_diverg = [], []
    # loop over each gene
    for gene in divergence:
        # check that divergence is defined (!= 'NA'), and < 5
        if divergence[gene][i] != 'NA' and divergence[gene][i] < 5:
            # ask if gene in PTC or non_PTC genes
            if gene in PTC_genes:
                PTC_diverg.append(divergence[gene][i])
            elif gene in non_PTC_genes:
                non_PTC_diverg.append(divergence[gene][i])
    # done populating lists, write resukts to file
    if i == 0:
        summary_file.write('dN' + '\t')
    elif i == 1:
        summary_file.write('dS' + '\t')
    elif i == 2:
        summary_file.write('dN/dS' + '\t')
    # test mean differences between PTC and non-PTC divergence   
    wilcoxon, p_val = stats.ranksums(PTC_diverg, non_PTC_diverg)
    # write results to file 
    summary_file.write('\t'.join([str(len(PTC_diverg)), str(np.mean(PTC_diverg)), str(np.std(PTC_diverg) / math.sqrt(len(PTC_diverg))), str(len(non_PTC_diverg)), 
                                  str(np.mean(non_PTC_diverg)), str(np.std(non_PTC_diverg) / math.sqrt(len(non_PTC_diverg))), str(wilcoxon), str(p_val)]) + '\n')
    # update counter
    i += 1


# compare MAF and DAF distribution between X-linked and autosomal genes

# get the lists of MAF frequencies for PTC SNPs for X-linked  and autosomal genes, excluding sites with sample size < 10
MAF_X, MAF_auto = MAF_PTC_chromo('../CREM_CLA_protein_divergence/CDS_SNP_DIVERG.txt', '../CREM_CLA_protein_divergence/unique_transcripts.txt', 'transcripts_indels_CDS.txt', '../CREM_CLA_protein_divergence/356_10172014.gff3', '../CREM_CLA_protein_divergence/scaffold_chromosome.csv', 10)
# get the lists of DAF frequencies for PTC SNPs for X-linked and autosomal genes, excluding sites with sample size < 10
DAF_X, DAF_auto = DAF_PTC_chromo('../CREM_CLA_protein_divergence/CDS_SNP_DIVERG.txt', '../CREM_CLA_protein_divergence/unique_transcripts.txt', 'transcripts_indels_CDS.txt', '../CREM_CLA_protein_divergence/356_10172014.gff3', '../CREM_CLA_protein_divergence/scaffold_chromosome.csv', 10)


# express frequencies in %
for i in range(len(MAF_X)):
    MAF_X[i] = MAF_X[i] * 100
for i in range(len(MAF_auto)):
    MAF_auto[i] = MAF_auto[i] * 100
for i in range(len(DAF_X)):
    DAF_X[i] = DAF_X[i] * 100
for i in range(len(DAF_auto)):
    DAF_auto[i] = DAF_auto[i] * 100
    
# make histograms
MAF_X_hist = np.histogram(MAF_X, range(0, 51, 10))
MAF_auto_hist = np.histogram(MAF_auto, range(0, 51, 10))
DAF_X_hist = np.histogram(DAF_X, range(0, 101, 10))
DAF_auto_hist = np.histogram(DAF_auto, range(0, 101, 10))

# compare DAF distributions
diff_MAF_X_auto = stats.chi2_contingency([MAF_X_hist[0], MAF_auto_hist[0]])
diff_DAF_X_auto = stats.chi2_contingency([DAF_X_hist[0], DAF_auto_hist[0]])

# write results to summary file
summary_file.write('\n')
summary_file.write('Freq distributions between X-linked autosomal PTC genes:\n')
summary_file.write('-' * 57 + '\n')
summary_file.write('distributions' + '\t' + 'chi2' + '\t' + 'p-val' + '\t' + 'dof' + '\n')
summary_file.write('MAF_PTC_X_vs_auto' + '\t')
content = '\t'.join([str(diff_MAF_X_auto[i]) for i in range(3)])
summary_file.write(content + '\n')

summary_file.write('DAF_PTC_X_vs_auto' + '\t')
content = '\t'.join([str(diff_DAF_X_auto[i]) for i in range(3)])
summary_file.write(content + '\n')


# write tables to file from histograms
# open file to write the MAF frequencies of the PTC SNPs
newfile = open('MAF_PTC_SNPs_X_auto.txt', 'w')
newfile.write('MAF' + '\t' + 'X_linked' + '\t' + 'autosomal' + '\n')
for i in range(len(MAF_X_hist[0])):
    newfile.write(str(MAF_X_hist[1][i]) + ':' + str(MAF_X_hist[1][i] + 10) + '\t' + str(MAF_X_hist[0][i]) + '\t' + str(MAF_auto_hist[0][i]) + '\n')
newfile.close()
    
newfile = open('DAF_PTC_SNPs_X_auto.txt', 'w')
newfile.write('DAF' + '\t' + 'X_linked' + '\t' + 'autosomal' + '\n')
for i in range(len(DAF_X_hist[0])):
    newfile.write(str(DAF_X_hist[1][i]) + ':' + str(DAF_X_hist[1][i] + 10) + '\t' + str(DAF_X_hist[0][i]) + '\t' + str(DAF_auto_hist[0][i]) + '\n')
newfile.close()

# close summary file
summary_file.close()

