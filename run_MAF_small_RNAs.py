# -*- coding: utf-8 -*-
"""
Created on Fri Aug  7 16:49:15 2015

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

# get the coordinates of the miRNA loci
mirnas_coord = get_mirna_loci('../miRNA_Target_sites/crm_miRBase21_premina_coordinates.txt')

# get the allele counts for miRNA sites
mirna_sites = get_feature_sites(chromo_sites, mirnas_coord)

# compute MAF for mirna sites (sites with sample size < 10 are already excluded)
MAF_mirna = MAF_non_coding(mirna_sites)

print('MAF for miRNA sites done')


# get the coordinates of the piRNA loci
pirnas_coord = get_pirna_loci('PX356_piRNA_coord.txt')

# get the allele counts for piRNA sites
pirna_sites = get_feature_sites(chromo_sites, pirnas_coord)

# compute MAF for piRNA sites, (sites with sample size < 10 are already excluded)
MAF_pirna = MAF_non_coding(pirna_sites)

print('MAF for piRNA sites done')


# get all the miRNA positions in genome
mirna_pos = get_small_rna_sites(mirnas_coord)

print('got miRNA positions')

# get all the piRNAs positions in the genome
pirna_pos = get_small_rna_sites(pirnas_coord)

print('got piRNA positions')


# compute threshold based on the distribution of elegans UTR length
UTR_length = celegans_three_prime_UTR_length('../Genome_Files/c_elegans.PRJNA13758.WS248.annotations.gff3')
threshold = get_percentile(UTR_length, 99)
# get UTR coord {TS1 : [chromo, start, end, orientation]}
three_prime = get_three_prime_UTR_positions('../Genome_Files/356_10172014.gff3', '../Genome_Files/noamb_356_v1_4.txt', threshold)
# get all the predicted UTR positions in the genome
UTR_pos = get_UTR_sites(three_prime)
print('got UTR positions')


# get all the gene positions in the genome
gene_pos = get_gene_sites('../CREM_CLA_protein_divergence/356_10172014.gff3')

print('got gene positions')

# get the coodinates of the intergenic sites
# this modifies the dict of allele counts for all sites

intergenic_sites = get_intergenic_sites(chromo_sites, gene_pos, pirna_pos, mirna_pos, UTR_pos, True)

print('got intergenic positions')

# get the MAF for intergenic sites, (sites with sample size < 10 are already excluded)
MAF_intergenic = MAF_non_coding(intergenic_sites)

print('MAF for intergenic sites done')


# express frequencies in %
for i in range(len(MAF_REP)):
    MAF_REP[i] = MAF_REP[i] * 100
for i in range(len(MAF_SYN)):
    MAF_SYN[i] = MAF_SYN[i] * 100
for i in range(len(MAF_mirna)):
    MAF_mirna[i] = MAF_mirna[i] * 100
for i in range(len(MAF_pirna)):
    MAF_pirna[i] = MAF_pirna[i] * 100
for i in range(len(MAF_intergenic)):
    MAF_intergenic[i] = MAF_intergenic[i] * 100
    
print('conversion to % frequencies done')

# make histograms
MAF_REP_hist = np.histogram(MAF_REP, range(0, 51, 10))
MAF_SYN_hist = np.histogram(MAF_SYN, range(0, 51, 10))
MAF_mirna_hist = np.histogram(MAF_mirna, range(0, 51, 10))
MAF_pirna_hist = np.histogram(MAF_pirna, range(0, 51, 10))
MAF_intergenic_hist = np.histogram(MAF_intergenic, range(0, 51, 10))

print('histograms done')

# compare MAF distributions
diff_MAF_REP_SYN = stats.chi2_contingency([MAF_REP_hist[0], MAF_SYN_hist[0]])
diff_MAF_REP_mirna = stats.chi2_contingency([MAF_REP_hist[0], MAF_mirna_hist[0]])
diff_MAF_REP_pirna = stats.chi2_contingency([MAF_REP_hist[0], MAF_pirna_hist[0]])
diff_MAF_REP_intergenic = stats.chi2_contingency([MAF_REP_hist[0], MAF_intergenic_hist[0]])
diff_MAF_SYN_mirna = stats.chi2_contingency([MAF_SYN_hist[0], MAF_mirna_hist[0]])
diff_MAF_SYN_pirna = stats.chi2_contingency([MAF_SYN_hist[0], MAF_pirna_hist[0]])
diff_MAF_SYN_intergenic = stats.chi2_contingency([MAF_SYN_hist[0], MAF_intergenic_hist[0]])
diff_MAF_mirna_pirna = stats.chi2_contingency([MAF_mirna_hist[0], MAF_pirna_hist[0]])
diff_MAF_mirna_intergenic = stats.chi2_contingency([MAF_mirna_hist[0], MAF_intergenic_hist[0]])
diff_MAF_pirna_intergenic = stats.chi2_contingency([MAF_pirna_hist[0], MAF_intergenic_hist[0]])

print('chi2 MAF comparisons done')

# make a list of test results
test_results = [diff_MAF_REP_SYN, diff_MAF_REP_mirna, diff_MAF_REP_pirna, diff_MAF_REP_intergenic,
                diff_MAF_SYN_mirna, diff_MAF_SYN_pirna, diff_MAF_SYN_intergenic,
                diff_MAF_mirna_pirna, diff_MAF_mirna_intergenic, diff_MAF_pirna_intergenic]

# make a list of strings for pairwise comparisons
pairwise_comp = ['REP_vs_SYN', 'REP_vs_mirna', 'REP_vs_pirna', 'REP_vs_intergenic', 'SYN_vs_mirna', 'SYN_vs_pirna', 'SYN_vs_intergenic',
'mirna_vs_pirna', 'mirna_vs_intergenic', 'pirna_vs_intergenic']

# open summary file
summary_file = open('summary_MAF_small_RNAs.txt', 'w')

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
newfile = open('MAF_SNPs_small_RNAs.txt', 'w')
newfile.write('MAF' + '\t' + 'REP' + '\t' + 'SYN' + '\t' + 'intergenic' + '\t' + 'miRNAs' + '\t' + 'piRNAs' + '\n')

print(len(MAF_REP_hist[0]))
print(len(MAF_REP_hist[1]))

for i in range(len(MAF_REP_hist[0])):
    newfile.write(str(MAF_REP_hist[1][i]) + ':' + str(MAF_REP_hist[1][i] + 10) + '\t' +
    str(MAF_REP_hist[0][i]) + '\t' + str(MAF_SYN_hist[0][i]) + '\t' + str(MAF_intergenic_hist[0][i]) + '\t' +
    str(MAF_mirna_hist[0][i]) + '\t' + str(MAF_pirna_hist[0][i]) + '\n')
newfile.close()

print('MAF distribution written to file')

# compute pairwise difference between average MAF

# make a list of MAF lists
MAF_list = [MAF_REP, MAF_SYN, MAF_intergenic, MAF_mirna, MAF_pirna]
# make a list of MAF names
MAF_names = ['REP', 'SYN', 'intergenic', 'mirna', 'pirna']


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

