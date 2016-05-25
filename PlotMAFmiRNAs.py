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
import json
# import custom modules
from manipulate_sequences import *
from miRNA_target import *
from genomic_coordinates import *
#from repeats_TEs import *
#from sliding_windows import *
from sites_with_coverage import *
from divergence import *
from premature_stops import *
from randomize_SNPs import *
from get_coding_sequences import *


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

# express frequencies in %
for i in range(len(MAF_REP)):
    MAF_REP[i] = MAF_REP[i] * 100
for i in range(len(MAF_SYN)):
    MAF_SYN[i] = MAF_SYN[i] * 100
MAF_REP_hist = np.histogram(MAF_REP, range(0, 51, 10))
MAF_SYN_hist = np.histogram(MAF_SYN, range(0, 51, 10))

# get the allele counts for all sites with coverage, exclude sites with sample size < 10
chromo_sites = get_non_coding_snps('../SNP_files/', 10)
print('got allele counts in genome')

# get the coordinates of the miRNA loci
mirnas_coord = get_mirna_loci('CRM_miRNAsCoordinatesFinal.txt')
print('got miRNA coordinates')

# get the allele counts for miRNA sites
mirna_sites = get_feature_sites(chromo_sites, mirnas_coord)
print('got allele counts for miRNAs')

# compute MAF for mirna sites (sites with sample size < 10 are already excluded)
MAF_mirna = MAF_non_coding(mirna_sites)
print('miRNAs', len(MAF_mirna))
print('MAF for miRNA sites done')

# get all the miRNA positions in genome
mirna_pos = get_small_rna_sites(mirnas_coord)
print('got miRNA positions')

# get UTR coordinates from json file
infile = open('CremUTRCoordsNo.json')
three_prime = json.load(infile)
infile.close()

## compute threshold based on the distribution of elegans UTR length
#UTR_length = celegans_three_prime_UTR_length('../Genome_Files/c_elegans.PRJNA13758.WS248.annotations.gff3')
#threshold = get_percentile(UTR_length, 99)
## get UTR coord {TS1 : [chromo, start, end, orientation]}
#three_prime = get_three_prime_UTR_positions('../Genome_Files/356_10172014.gff3', '../Genome_Files/noamb_356_v1_4.txt', threshold)
# get all the predicted UTR positions in the genome
UTR_pos = get_UTR_sites(three_prime)
print('got UTR positions')

# get all the gene positions in the genome
gene_pos = get_gene_sites('../Genome_Files/356_10172014.gff3')
print('got gene positions')

# get SNPs flanking miRNAs within 500 bp of the miRNAs
mirna_flanking_snps = get_small_rna_flanking_SNPs(chromo_sites, mirnas_coord, 500)
print('got allele counts in miRNA flanking regions')
# remove positions falling in coding sequences
# get CDS_coord {TS1: [chromo, sense, [(s1, end1), (s2, end2)]]}
CDS_coord = get_CDS_positions('../Genome_Files/356_10172014.gff3')
# create new dict in the form {chromo: {set of positions}}
CDS_pos = {}
for gene in CDS_coord:
    # get chromo
    chromo = CDS_coord[gene][0]
    # loop over cds coordinates 
    for i in range(len(CDS_coord[gene][2])):
        # get start and end positions in a list
        # convert to 0-based index
        start  = CDS_coord[gene][2][i][0] - 1
        end = CDS_coord[gene][2][i][1]
        # check if chromo in CDS_pos
        if chromo not in CDS_pos:
            CDS_pos[chromo] = set()
        for j in range(start, end):
            CDS_pos[chromo].add(j)
print('got CDS indices')
# remove positions falling in coding sequences
for chromo in CDS_pos:
    # check if chromo in flanking sites
    if chromo in mirna_flanking_snps:
        # loop over CDS positions
        for i in CDS_pos[chromo]:
            # check if site in flanking
            if i in mirna_flanking_snps[chromo]:
                # remove site
                del mirna_flanking_snps[chromo][i]
print('removed coding positions from miRNA flanking sites')

# count the number of snps in the vicinty of mirnas
mirna_snps = 0
for chromo in mirna_flanking_snps:
    mirna_snps += len(mirna_flanking_snps[chromo])
print('SNPs within 500 bp of miRNAs: ', mirna_snps)

# resample SNPs and compute MAF {replicate_number : [list of MAF values]}
# sampled 5000 SNPs 1000 times among the number of SNPs near miRNAs
mirna_resampled_MAF = SNP_MAF_randomization(mirna_flanking_snps, 5000, 1000)

###################### CONTNUE HERE

## obsolete code
## get the proportions of SNPs in each MAF bin 
#mirna_MAF_proportions = get_MAF_distribution_from_replicates(mirna_resampled_MAF)
## get the SNP proportions for the observed SNPs
#empirical_mirna_MAF = SNP_proportions_MAF_bin(MAF_mirna)
## creat a list of keys, being the MAF lower bound in the dicts with MAF proportions
#maf_limit = [i for i in empirical_mirna_MAF]
## sort list
#maf_limit.sort()
#
#
#average = np.mean(mirna_MAF_proportions[i])
#stdev = np.std(mirna_MAF_proportions[i])
#observed = empirical_mirna_MAF[i]
#z_score = (observed - average) / stdev
## critical values for 1-sample 2-tailed z-test: 0.05: +- 1.96, 0.01: +- 2.58, 0.001: +-3.27
## Ho : obsvered == mean, H1: observed != mean
#if z_score < -1.96 or z_score > 1.96:
#     P5 = '*'
#elif -1.96 <= z_score <= 1.96:
#     P5 = 'NS'
#if z_score < -2.58 or z_score > 2.58:
#     P1 = '*'
#elif -2.58 <= z_score <= 2.58:
#     P1 = 'NS'
#if z_score < -3.27 or z_score > 3.27:
#    P01 = '*'
#elif -3.27 <= z_score <= 3.27:
#    P01 = 'NS'


#####################

# get the coodinates of the intergenic sites
# this modifies the dict of allele counts for all sites
intergenic_sites = get_intergenic_sites(chromo_sites, gene_pos, pirna_pos, mirna_pos, UTR_pos, True)
print('got intergenic positions')

# get the MAF for intergenic sites, (sites with sample size < 10 are already excluded)
MAF_intergenic = MAF_non_coding(intergenic_sites)
print('MAF for intergenic sites done')

# express MAF frequencies in %
for i in range(len(MAF_REP)):
    MAF_REP[i] = MAF_REP[i] * 100
for i in range(len(MAF_SYN)):
    MAF_SYN[i] = MAF_SYN[i] * 100
for i in range(len(MAF_mirna)):
    MAF_mirna[i] = MAF_mirna[i] * 100
for i in range(len(MAF_intergenic)):
    MAF_intergenic[i] = MAF_intergenic[i] * 100
for i in mirna_resampled_MAF:
    for j in range(len(mirna_resampled_MAF[i])):
        mirna_resampled_MAF[i][j] = mirna_resampled_MAF[i][j] * 100
print('conversion to % frequencies done')


# make a list of maximum frequencies
maxfreq = []
for i in mirna_resampled_MAF:
    maxfreq.append(max(mirna_resampled_MAF))
for i in [MAF_REP, MAF_SYN, MAF_mirna, MAF_intergenic]:
    maxfreq.append(max(i))
print(max(maxfreq))
assert max(maxfreq) <= 50, 'MAF should be lower than 50%'








#
## make histograms
#MAF_REP_hist = np.histogram(MAF_REP, range(0, 51, 10))
#MAF_SYN_hist = np.histogram(MAF_SYN, range(0, 51, 10))
#MAF_mirna_hist = np.histogram(MAF_mirna, range(0, 51, 10))
#MAF_intergenic_hist = np.histogram(MAF_intergenic, range(0, 51, 10))
## create a dict to store counts in MAF bins {replicate_number : [MAF counts]}
#MAF_resampling = {}
#for i in mirna_resampled_MAF:
#    MAF_resampling[i] = np.histogram(mirna_resampled_MAF[i], range(0, 51, ))
#
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
################################################### SAVE FIG
#
#
## create figure
#fig = plt.figure(1, figsize = (4, 2))
## add a plot to figure (1 row, 1 column, 1 plot)
#ax = fig.add_subplot(1, 1, 1)  
#
#
##colorscheme = ['#810f7c', '#8856a7', '#8c96c6', '#b3cde3', '#edf8fb']
#
#
#colorscheme = ['#810f7c', '#8856a7']
#
#
#
#repfreq = []
#for i in MAF_REP_hist[0]:
#    repfreq.append(i / sum(MAF_REP_hist[0]))
#synfreq = []
#for i in MAF_SYN_hist[0]:
#    synfreq.append(i / sum(MAF_SYN_hist[0]))
#
## plot the repeat of gene density per window
##ax.hist([MAF_REP,MAF_SYN], color = colorscheme)
#
#width = 0.2
#
#print('REP', MAF_REP_hist[0], MAF_REP_hist[1])
#
#ax.bar([0, 0.4, 0.8, 1.2, 1.6], repfreq, width, yerr = [0.1, 0.2, 0.05, 0, 0.3], color = '#810f7c',
#       error_kw=dict(elinewidth=1, ecolor='black', markeredgewidth = 1))
#ax.bar([0.2, 0.6, 1, 1.4, 1.8], synfreq, width, color = '#8856a7')
#
#
#
#
#ax.set_ylabel('Proportion of SNPs', size = 10, ha = 'center', fontname = 'Arial')
#
#
#
# 
### determine tick position on x axis
##xpos =  [j for j in range(0, len(positions) + 50, 50)]
### convert interval windows numbers to genomic positions
##xtext = list(map(lambda x : (x * 50000) / 1000000, xpos))
##xtext = list(map(lambda x : str(x), xtext))
### set up tick positions and labels
##plt.xticks(xpos, xtext, fontsize = 10, fontname = 'Arial')
##plt.yticks(fontsize = 0)
#
## set x axis label
#ax.set_xlabel('Minor Allele Frequency', size = 10, ha = 'center', fontname = 'Arial')
#
#
## do not show lines around figure, keep bottow line  
#ax.spines["top"].set_visible(False)    
#ax.spines["bottom"].set_visible(True)    
#ax.spines["right"].set_visible(False)    
#ax.spines["left"].set_visible(False)      
## offset the spines
#for spine in ax.spines.values():
#  spine.set_position(('outward', 5))
#  
## add a light grey horizontal grid to the plot, semi-transparent, 
#ax.yaxis.grid(True, linestyle='--', which='major', color='lightgrey', alpha=0.5, linewidth = 0.5)  
## hide these grids behind plot objects
#ax.set_axisbelow(True)
#
#
## do not show ticks on 1st graph
#ax.tick_params(
#    axis='x',       # changes apply to the x-axis and y-axis (other option : x, y)
#    which='both',      # both major and minor ticks are affected
#    bottom='on',      # ticks along the bottom edge are off
#    top='off',         # ticks along the top edge are off
#    right = 'off',
#    left = 'off',          
#    labelbottom='on', # labels along the bottom edge are off 
#    colors = 'black',
#    labelsize = 10,
#    direction = 'out') # ticks are outside the frame when bottom = 'on
##
##
## do not show ticks
#ax.tick_params(
#    axis='y',       # changes apply to the x-axis and y-axis (other option : x, y)
#    which='both',      # both major and minor ticks are affected
#    bottom='off',      # ticks along the bottom edge are off
#    top='off',         # ticks along the top edge are off
#    right = 'off',
#    left = 'off',          
#    labelbottom='off', # labels along the bottom edge are off 
#    colors = 'black',
#    labelsize = 10,
#    direction = 'out') # ticks are outside the frame when bottom = 'on
#
#
#for label in ax.get_yticklabels():
#    label.set_fontname('Arial')
#
### add lines
##lns = graph1+graph2
### get labels
##if density == 'genes':
##    labs = ['Genes', 'Diversity']
##elif density == 'repeats':
##    labs = ['Repeats', 'Diversity']
### plot legend
##ax2.legend(lns, labs, loc=2, fontsize = 8, frameon = False)
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
#
#fig.savefig('testfile.pdf', bbox_inches = 'tight')
