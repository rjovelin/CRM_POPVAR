# -*- coding: utf-8 -*-
"""
Created on Tue May 24 14:19:12 2016

@author: RJovelin
"""

# use this script to plot the MAF distribution of miRNAs and other sites

# use Agg backend on server without X server
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import rc
rc('mathtext', default='regular')
# import modules
import numpy as np
from scipy import stats
import math
import sys
# import custom modules
from manipulate_sequences import *
from miRNA_target import *
from repeats_TEs import *
from sliding_windows import *
from sites_with_coverage import *
from divergence import *
from premature_stops import *



# make a list with MAF of replacement SNPs, excluding sites with sample size < 10
MAF_REP = MAF_SNP('../Genome_Files/CDS_SNP_DIVERG.txt', '../Genome_Files/unique_transcripts.txt', 
                  '../Genome_Files/transcripts_indels_CDS.txt', 'REP', 10)
print('REP', len(MAF_REP))
print('MAF for replacement sites done')

# make a list with MAF of synonymous sites, excluding sites with sample size < 10
MAF_SYN = MAF_SNP('../Genome_Files/CDS_SNP_DIVERG.txt', '../Genome_Files/unique_transcripts.txt', 
                  '../Genome_Files/transcripts_indels_CDS.txt', 'SYN', 10)
print('SNP', len(MAF_SYN))
print('MAF for synonymous sites done')

# create figure
fig = plt.figure(1, figsize = (4, 2))
# add a plot to figure (1 row, 1 column, 1 plot)
ax = fig.add_subplot(1, 1, 1)  

# plot the repeat of gene density per window
ax.hist([MAF_REP,MAF_SYN])

#ax.hist([MAF_REP,MAF_SYN], list(map(lambda x : x / 100, range(0, 100, 10))))





fig.savefig('testfile.pdf', bbox_inches = 'tight')

## get the allele counts for all sites with coverage, exclude sites with sample size < 10
#chromo_sites = get_non_coding_snps('../SNP_files/', 10)
#
#print('got allele counts in genome')
#
## get the coordinates of the miRNA loci
#mirnas_coord = get_mirna_loci('../miRNA_Target_sites/crm_miRBase21_premina_coordinates.txt')
#
## get the allele counts for miRNA sites
#mirna_sites = get_feature_sites(chromo_sites, mirnas_coord)
#
## compute MAF for mirna sites (sites with sample size < 10 are already excluded)
#MAF_mirna = MAF_non_coding(mirna_sites)
#
#print('MAF for miRNA sites done')
#
#
## get the coordinates of the piRNA loci
#pirnas_coord = get_pirna_loci('PX356_piRNA_coord.txt')
#
## get the allele counts for piRNA sites
#pirna_sites = get_feature_sites(chromo_sites, pirnas_coord)
#
## compute MAF for piRNA sites, (sites with sample size < 10 are already excluded)
#MAF_pirna = MAF_non_coding(pirna_sites)
#
#print('MAF for piRNA sites done')
#
## get all the miRNA positions in genome
#mirna_pos = get_small_rna_sites('../miRNA_Target_sites/crm_miRBase21_premina_coordinates.txt', 'miRNA')
#
#print('got miRNA positions')
#
## get all the piRNAs positions in the genome
#pirna_pos = get_small_rna_sites('PX356_piRNA_coord.txt', 'piRNA')
#
#print('got piRNA positions')
#
## get all the predicted UTR positions in the genome
#UTR_pos = get_UTR_sites('../CREM_CLA_protein_divergence/356_10172014.gff3', '../miRNA_Target_sites/c_elegans.PRJNA13758.WS248.annotations.gff3', '../CREM_CLA_protein_divergence/noamb_356_v1_4.txt', 99)
#
#print('got UTR positions')
#
## get all the gene positions in the genome
#gene_pos = get_gene_sites('../CREM_CLA_protein_divergence/356_10172014.gff3')
#
#print('got gene positions')
#
## get the coodinates of the intergenic sites
## this modifies the dict of allele counts for all sites
#
#intergenic_sites = get_intergenic_sites(chromo_sites, gene_pos, pirna_pos, mirna_pos, UTR_pos, True)
#
#print('got intergenic positions')
#
## get the MAF for intergenic sites, (sites with sample size < 10 are already excluded)
#MAF_intergenic = MAF_non_coding(intergenic_sites)
#
#print('MAF for intergenic sites done')
#
#
## express frequencies in %
#for i in range(len(MAF_REP)):
#    MAF_REP[i] = MAF_REP[i] * 100
#for i in range(len(MAF_SYN)):
#    MAF_SYN[i] = MAF_SYN[i] * 100
#for i in range(len(MAF_mirna)):
#    MAF_mirna[i] = MAF_mirna[i] * 100
#for i in range(len(MAF_pirna)):
#    MAF_pirna[i] = MAF_pirna[i] * 100
#for i in range(len(MAF_intergenic)):
#    MAF_intergenic[i] = MAF_intergenic[i] * 100
#    
#print('conversion to % frequencies done')
#
## make histograms
#MAF_REP_hist = np.histogram(MAF_REP, range(0, 51, 10))
#MAF_SYN_hist = np.histogram(MAF_SYN, range(0, 51, 10))
#MAF_mirna_hist = np.histogram(MAF_mirna, range(0, 51, 10))
#MAF_pirna_hist = np.histogram(MAF_pirna, range(0, 51, 10))
#MAF_intergenic_hist = np.histogram(MAF_intergenic, range(0, 51, 10))
#
#print('histograms done')
#
## compare MAF distributions
#diff_MAF_REP_SYN = stats.chi2_contingency([MAF_REP_hist[0], MAF_SYN_hist[0]])
#diff_MAF_REP_mirna = stats.chi2_contingency([MAF_REP_hist[0], MAF_mirna_hist[0]])
#diff_MAF_REP_pirna = stats.chi2_contingency([MAF_REP_hist[0], MAF_pirna_hist[0]])
#diff_MAF_REP_intergenic = stats.chi2_contingency([MAF_REP_hist[0], MAF_intergenic_hist[0]])
#diff_MAF_SYN_mirna = stats.chi2_contingency([MAF_SYN_hist[0], MAF_mirna_hist[0]])
#diff_MAF_SYN_pirna = stats.chi2_contingency([MAF_SYN_hist[0], MAF_pirna_hist[0]])
#diff_MAF_SYN_intergenic = stats.chi2_contingency([MAF_SYN_hist[0], MAF_intergenic_hist[0]])
#diff_MAF_mirna_pirna = stats.chi2_contingency([MAF_mirna_hist[0], MAF_pirna_hist[0]])
#diff_MAF_mirna_intergenic = stats.chi2_contingency([MAF_mirna_hist[0], MAF_intergenic_hist[0]])
#diff_MAF_pirna_intergenic = stats.chi2_contingency([MAF_pirna_hist[0], MAF_intergenic_hist[0]])
#
#print('chi2 MAF comparisons done')
#
## make a list of test results
#test_results = [diff_MAF_REP_SYN, diff_MAF_REP_mirna, diff_MAF_REP_pirna, diff_MAF_REP_intergenic,
#                diff_MAF_SYN_mirna, diff_MAF_SYN_pirna, diff_MAF_SYN_intergenic,
#                diff_MAF_mirna_pirna, diff_MAF_mirna_intergenic, diff_MAF_pirna_intergenic]
#
## make a list of strings for pairwise comparisons
#pairwise_comp = ['REP_vs_SYN', 'REP_vs_mirna', 'REP_vs_pirna', 'REP_vs_intergenic', 'SYN_vs_mirna', 'SYN_vs_pirna', 'SYN_vs_intergenic',
#'mirna_vs_pirna', 'mirna_vs_intergenic', 'pirna_vs_intergenic']
#
## open summary file
#summary_file = open('summary_MAF_small_RNAs.txt', 'w')
#
## write results to summary file
#summary_file.write('comparison of the MAF distributions\n')
#summary_file.write('-' * 36 + '\n')
#summary_file.write('distributions' + '\t' + 'chi2' + '\t' + 'p-val' + '\t' + 'dof' + '\n')
#
## loop over list of string comparisons
## write comp to file and write results of chi2 to file
#for i in range(len(pairwise_comp)):
#    summary_file.write(pairwise_comp[i] + '\t')
#    content = '\t'.join([str(test_results[i][j]) for j in range(3)])
#    summary_file.write(content + '\n')
#
#print('chi2 results of MAF differences written to file')
#
## write tables to file from histograms
## open file to write the MAF frequencies of the different SNPs
#newfile = open('MAF_SNPs_small_RNAs.txt', 'w')
#newfile.write('MAF' + '\t' + 'REP' + '\t' + 'SYN' + '\t' + 'intergenic' + '\t' + 'miRNAs' + '\t' + 'piRNAs' + '\n')
#
#print(len(MAF_REP_hist[0]))
#print(len(MAF_REP_hist[1]))
#
#for i in range(len(MAF_REP_hist[0])):
#    newfile.write(str(MAF_REP_hist[1][i]) + ':' + str(MAF_REP_hist[1][i] + 10) + '\t' +
#    str(MAF_REP_hist[0][i]) + '\t' + str(MAF_SYN_hist[0][i]) + '\t' + str(MAF_intergenic_hist[0][i]) + '\t' +
#    str(MAF_mirna_hist[0][i]) + '\t' + str(MAF_pirna_hist[0][i]) + '\n')
#newfile.close()
#
#print('MAF distribution written to file')
#
## compute pairwise difference between average MAF
#
## make a list of MAF lists
#MAF_list = [MAF_REP, MAF_SYN, MAF_intergenic, MAF_mirna, MAF_pirna]
## make a list of MAF names
#MAF_names = ['REP', 'SYN', 'intergenic', 'mirna', 'pirna']
#
#
#summary_file.write('\n')
#summary_file.write('comparison of mean MAF differences\n')
#summary_file.write('-' * 35 + '\n')
#summary_file.write('distributions' + '\t' + 'wilcoxon' + '\t' + 'p-val' + '\n')
#
## loop over MAF_list, comparre average differences
#for i in range(0, len(MAF_list) - 1):
#    for j in range(i+1, len(MAF_list)):
#        # test mean differences    
#        wilcoxon, p_val = stats.ranksums(MAF_list[i], MAF_list[j])
#        summary_file.write(MAF_names[i] + '_' + MAF_names[j] + '\t' + str(wilcoxon) + '\t' + str(p_val) + '\n')
#
#
#summary_file.write('\n')
#summary_file.write('Mean MAF at sites\n')
#summary_file.write('-' * 18 + '\n')
#summary_file.write('sites' + '\t' + 'mean_MAF' + '\t' + 'SEM' + '\n')
#for i in range(len(MAF_list)):
#    summary_file.write('\t'.join([MAF_names[i], str(np.mean(MAF_list[i])), str(np.std(MAF_list[i]) / math.sqrt(len(MAF_list[i])))]) + '\n') 
#
## close summary file
#summary_file.close()
#
#
#
#
#
#
#
#
#
#
#############################
#
#
#
#
## -*- coding: utf-8 -*-
#"""
#Created on Tue Aug 18 22:50:48 2015
#
#@author: Richard
#"""
#
#
#
#import math
#import numpy as np
#from premature_stops import *
#from scipy import stats
#from accessories import *
#from miRNA_target import *
#from piRNAs import *
#from sites_with_coverage import *
#from Cel_UTR import *
#from cel_UTR_length import *
#from randomize_SNPs import *
#from get_coding_sequences import *
#
#
#Gstr = lambda x : str(x)
#
#
## get the allele counts for all sites with coverage, exclude sites with sample size < 10
#chromo_sites = get_non_coding_snps('../SNP_files/', 10)
#
#print('got allele counts in genome')
#
#
## get CDS_coord {TS1: [chromo, sense, [(s1, end1), (s2, end2)]]}
#CDS_coord = get_CDS_positions('../CREM_CLA_protein_divergence/356_10172014.gff3')
## create new dict in the form {chromo: {set of positions}}
#CDS_pos = {}
#for gene in CDS_coord:
#    # get chromo
#    chromo = CDS_coord[gene][0]
#    # loop over cds coordinates 
#    for i in range(len(CDS_coord[gene][2])):
#        # get start and end positions in a list
#        # convert to 0-based index
#        start  = CDS_coord[gene][2][i][0] - 1
#        end = CDS_coord[gene][2][i][1]
#        # check if chromo in CDS_pos
#        if chromo not in CDS_pos:
#            CDS_pos[chromo] = set()
#        for j in range(start, end):
#            CDS_pos[chromo].add(j)
#        
#print('got CDS coord')
#
## get the coordinates of the miRNA loci
#mirnas_coord = get_mirna_loci('crm_miRBase21_premina_coordinates.txt')
#
## get the allele counts for miRNA sites
#mirna_sites = get_feature_sites(chromo_sites, mirnas_coord)
#
## compute MAF for mirna sites (sites with sample size < 10 are already excluded)
#MAF_mirna = MAF_non_coding(mirna_sites)
#
#print('MAF for miRNA sites done')
#
## get SNPs flanking miRNAs within 500 bp of the miRNAs
#mirna_flanking_snps = get_small_rna_flanking_SNPs(chromo_sites, 'crm_miRBase21_premina_coordinates.txt',
#                                                  'miRNA', 500)
#
## remove positions falling in coding sequences
#for chromo in CDS_pos:
#    # check if chromo in flanking sites
#    if chromo in mirna_flanking_snps:
#        # loop over CDS positions
#        for i in CDS_pos[chromo]:
#            # check if site in flanking
#            if i in mirna_flanking_snps[chromo]:
#                # remove site
#                del mirna_flanking_snps[chromo][i]
#                
#print('got SNPs flanking miRNAs')
#
#mirna_snps = 0
#for chromo in mirna_flanking_snps:
#    mirna_snps += len(mirna_flanking_snps[chromo])
#print('SNPs within 500 bp of miRNAs: ', mirna_snps)
#
## resample SNPs and compute MAF
#mirna_resampled_MAF = SNP_MAF_randomization(mirna_flanking_snps, 5000, 1000)
#
## get the proportions of SNPs in each MAF bin
#mirna_MAF_proportions = get_MAF_distribution_from_replicates(mirna_resampled_MAF)
#
## get the SNP proportions for the observed SNPs
#empirical_mirna_MAF = SNP_proportions_MAF_bin(MAF_mirna)
#
## creat a list of keys, being the MAF lower bound in the dicts with MAF proportions
#maf_limit = [i for i in empirical_mirna_MAF]
## sort list
#maf_limit.sort()
#
## open file to save results of Z-test
#newfile = open('summary_MAF_resampled_SNPs_miRNAs.txt', 'w')
#
#newfile.write('Z-test of SNP proportions in MAF bins\n')
#newfile.write('random resampling of 5000 SNPs within 500 bp of miRNAs, exluding miRNA and CDS sites with 1000 replicates\n')
#
#newfile.write('\t'.join(['MAF', 'mean_sample', 'margin', 'stdev_sample', 'l95', 'h95', 'observed', 'z-score', 'alpha_0.05', 'alpha_0.01', 'alpha_0.001']) + '\n')
#
## loop over sorted keys
#for i in maf_limit:
#    # compute mean
#    average = np.mean(mirna_MAF_proportions[i])
#    # get standard deviation
#    stdev = np.std(mirna_MAF_proportions[i])
#    # compute 95% CI
#    stderror = stdev / math.sqrt(len(mirna_MAF_proportions[i]))
#    # compute the margin error (critical value = 1.96 for 95% CI)
#    margin = 1.96 * stderror
#    lCI = average - margin
#    hCI = average + margin
#    observed = empirical_mirna_MAF[i]
#    z_score = (observed - average) / stdev
#    # critical values for 1-sample 2-tailed z-test: 0.05: +- 1.96, 0.01: +- 2.58, 0.001: +-3.27
#    # Ho : obsvered == mean, H1: observed != mean
#    if z_score < -1.96 or z_score > 1.96:
#        P5 = '*'
#    elif -1.96 <= z_score <= 1.96:
#        P5 = 'NS'
#    if z_score < -2.58 or z_score > 2.58:
#        P1 = '*'
#    elif -2.58 <= z_score <= 2.58:
#        P1 = 'NS'
#    if z_score < -3.27 or z_score > 3.27:
#        P01 = '*'
#    elif -3.27 <= z_score <= 3.27:
#        P01 = 'NS'
#    newfile.write('\t'.join([str(i), str(average), str(margin), str(stdev), str(lCI), str(hCI), str(observed), str(z_score), str(P5), str(P1), str(P01)]) + '\n')
#    
## close file after writing
#newfile.close()
#    
#
#
#
