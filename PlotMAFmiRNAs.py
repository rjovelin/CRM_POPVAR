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
import matplotlib.patches as mpatches
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
print('resampled SNPs near miRNAs')

# get target coordinates from json file
infile = open('AllmiRNATargetsCoords.json')
target_coord = json.load(infile)
infile.close()
print('target coords', len(target_coord))
print('got miRNA target site coordinates')
# get the allele counts for all targets
target_sites = get_feature_sites(chromo_sites, target_coord)
print('got allele counts for target sites')
# compute MAF for all targets
MAF_targets = MAF_non_coding(target_sites)
print('MAF for miRNA targets done')

# load UTR coordinates from json file
infile = open('CremUTRCoordsNo.json')
UTR_coord = json.load(infile)
infile.close()
print('UTR coords', len(UTR_coord))
print('got UTR coordinates from json file')
# get set of valid transcripts
valid_transcripts = get_valid_transcripts('../Genome_Files/unique_transcripts.txt')
# remove UTR of non-valid transcripts from UTR coords
to_remove = [gene for gene in UTR_coord if gene not in valid_transcripts]
if len(to_remove) != 0:
    for gene in to_remove:
        del UTR_coord[gene]
print('removed {0} genes from UTR coords'.format(len(to_remove)))
print('UTR coords', len(UTR_coord))
# get SNPs in non-miRNA target UTRs
UTR_snps = get_UTR_SNPs(chromo_sites, UTR_coord, target_coord)
print('got allele counts for non-target UTR sites')

# get the total number of snps in non-target UTRs
non_target_snps = 0
for chromo in UTR_snps:
    non_target_snps += len(UTR_snps[chromo])
print('SNPs at non-target UTR sites: ', non_target_snps)

# resample SNPs and compute MAF {replicate_number : [list of MAF values]}
# sampled 10000 SNPs 1000 times among the number of SNPs in UTRs excluding target sites
UTR_resampled_MAF = SNP_MAF_randomization(UTR_snps, 10000, 1000)
print('reseampled SNPs in UTRs')

# express MAF frequencies in %
for i in range(len(MAF_REP)):
    MAF_REP[i] = MAF_REP[i] * 100
for i in range(len(MAF_SYN)):
    MAF_SYN[i] = MAF_SYN[i] * 100
for i in range(len(MAF_mirna)):
    MAF_mirna[i] = MAF_mirna[i] * 100
for i in range(len(MAF_targets)):
    MAF_targets[i] = MAF_targets[i] * 100    
for i in mirna_resampled_MAF:
    for j in range(len(mirna_resampled_MAF[i])):
        mirna_resampled_MAF[i][j] = mirna_resampled_MAF[i][j] * 100
for i in UTR_resampled_MAF:
    for j in range(len(UTR_resampled_MAF[i])):
        UTR_resampled_MAF[i][j] = UTR_resampled_MAF[i][j] * 100
print('conversion to % frequencies done')


# make histograms to get the counts in each MAF bins
MAF_REP_hist = np.histogram(MAF_REP, range(0, 51, 10))
MAF_SYN_hist = np.histogram(MAF_SYN, range(0, 51, 10))
MAF_mirna_hist = np.histogram(MAF_mirna, range(0, 51, 10))
MAF_targets_hist = np.histogram(MAF_targets, range(0, 51, 10))
# create dicts to store counts in MAF bins {replicate_number : [[MAF counts], [MAF freq]]}
MAF_mirna_resampling = {}
for i in mirna_resampled_MAF:
    MAF_mirna_resampling[i] = np.histogram(mirna_resampled_MAF[i], range(0, 51, 10))
MAF_UTR_resampling = {}
for i in UTR_resampled_MAF:
    MAF_UTR_resampling[i] = np.histogram(UTR_resampled_MAF[i], range(0, 51, 10))
print('got SNP counts in each MAF bin')

# get the proportions of SNPs in each MAF bins
REP_freq = [i / sum(MAF_REP_hist[0]) for i in MAF_REP_hist[0]]
SYN_freq = [i / sum(MAF_SYN_hist[0]) for i in MAF_SYN_hist[0]]
mirna_freq = [i / sum(MAF_mirna_hist[0]) for i in MAF_mirna_hist[0]]
targets_freq = [i / sum(MAF_targets_hist[0]) for i in MAF_targets_hist[0]]
# check length of lists of SNP proportions
for i in [REP_freq, SYN_freq, mirna_freq, targets_freq]:
    assert len(i) == 5, 'there should be 5 values in list of SNP proportions'
mirna_resampling, targets_resampling = {}, {}
# get the snps proprtions for each replicate
for i in MAF_mirna_resampling:
    mirna_resampling[i] = [j / sum(MAF_mirna_resampling[i][0]) for j in MAF_mirna_resampling[i][0]]
for i in MAF_UTR_resampling:
    targets_resampling[i] = [j / sum(MAF_UTR_resampling[i][0]) for j in MAF_UTR_resampling[i][0]]
# get the mean and SEM SNP proportions for resampling
# create dicts with MAF bin index as value and list of SNP proportion as value    
mirna_prop, targets_prop = {}, {}
for i in mirna_resampling:
    # check length of list of SNP proportions
    assert len(mirna_resampling[i]) == 5, 'Number of SNP proportions is different than 5'
    for j in range(len(mirna_resampling[i])):
        if j not in mirna_prop:
            mirna_prop[j] = []
        mirna_prop[j].append(mirna_resampling[i][j])
for i in targets_resampling:
    # check length of list of SNP proportions
    assert len(targets_resampling[i]) == 5, 'Number of SNP proportions is different than 5'
    for j in range(len(targets_resampling[i])):
        if j not in targets_prop:
            targets_prop[j] = []
        targets_prop[j].append(targets_resampling[i][j])
assert len(targets_prop) == len(mirna_prop), 'should have same size'
assert len(targets_prop) == 5, 'should have size 5'
# take mean and SEM for each MAF bin
nums = [i for i in mirna_prop]
nums.sort()
mirna_samplefreq, mirna_sem, target_samplefreq, target_sem = [], [], [], []
for i in nums:
    mirna_samplefreq.append(np.mean(mirna_prop[i]))
    target_samplefreq.append(np.mean(targets_prop[i]))
    mirna_sem.append(np.std(mirna_prop[i]) / math.sqrt(len(mirna_prop[i])))
    target_sem.append(np.std(targets_prop[i]) / math.sqrt(len(targets_prop[i])))
print('got SNP proportions in each MAF bin')


# compare MAF miRNA to other sites
# save results in summary file
newfile = open('TestMAFDistribution.txt', 'w')
newfile.write('Comparison of MAF distributions\n')
newfile.write('=' * 32 + '\n')
newfile.write('\t'.join(['miRNA_sites', 'sites', 'Chi2', 'P']) + '\n')

sites = ['REP', 'SYN', 'Targets']
maf = [MAF_REP_hist[0], MAF_SYN_hist[0], MAF_targets_hist[0]]
for i in range(len(maf)):
    diff = stats.chi2_contingency(MAF_mirna_hist[0], maf[i])    
    newfile.write('\t'.join(['miRNA', sites[i], str(diff[0]), str(diff[1])]) + '\n')

sites = ['REP', 'SYN']
maf = [MAF_REP_hist[0], MAF_SYN_hist[0]]
for i in range(len(maf)):
    diff = stats.chi2_contingency(MAF_targets_hist[0], maf[i])
    newfile.write('\t'.join(['targets', sites[i], str(diff[0]), str(diff[1])]) + '\n')

newfile.write('\n' * 3)
newfile.write('Comparison of SNP proportions to random resampling\n')
newfile.write('=' * 31 + '\n')
newfile.write('\t'.join(['sites', 'MAF', 'observed', 'mean_sample', 'stdev_sample', 'z-score', 'P']) + '\n')

# write function to check significance of zscore
def ZScoreTest(zscore):
    # critical values for 1-sample 2-tailed z-test: 0.05: +- 1.96, 0.01: +- 2.58, 0.001: +-3.27
    # Ho : obsvered == mean, H1: observed != mean
    if zscore < -3.27 or zscore > 3.27:
        P = 0.001
    elif -3.27 <= zscore <= 3.27:
        P = 'NS'
    # if not significant, ask if significant at different alpha
    if P == 'Ns':
        if zscore < -2.58 or zscore > 2.58:
            P = 0.01
        elif -2.58 <= zscore <= 2.58:
            P = 'NS'
    if P == 'NS':
        if zscore < -1.96 or zscore > 1.96:
            P = 0.05
        elif -1.96 <= zscore <= 1.96:
            P = 'NS'
    return P

j = 0
# perform z-test for each MAF bin between mirnas and resampled SNPs near miRNAs
for i in nums:
    observed = mirna_freq[i]
    average = mirna_samplefreq[i]
    stdev = np.std(mirna_prop[i])
    zscore = (observed - average) / stdev
    P = ZScoreTest(zscore)        
    newfile.write('\t'.join(['miRNAs', ':'.join([str(j), str(j+10)]), str(observed), str(average), str(zscore), str(P)]) + '\n')    
    j += 10        

j = 0
# perform z-test for each MAF bin between targets and resampled SNPs in UTR
for i in nums:
    observed = targets_freq[i]
    average = target_samplefreq[i]
    stdev = np.std(targets_prop[i])
    zscore = (observed - average) / stdev
    P = ZScoreTest(zscore)        
    newfile.write('\t'.join(['targets', ':'.join([str(j), str(j+10)]), str(observed), str(average), str(zscore), str(P)]) + '\n')    
    j += 10        

# close file after writing
newfile.close()
print('results of MAF differences written to file')


# create figure
fig = plt.figure(1, figsize = (4, 2))
# add a plot to figure (1 row, 1 column, 1 plot)
ax = fig.add_subplot(1, 1, 1)  

# set width of bar
width = 0.2


# plot SNP proportions SYN
graph1 = ax.bar([0, 1.2, 2.4, 3.6, 4.8], SYN_freq, width, color = '#810f7c', edgecolor = 'black', linewidth = 1)
# plot SNP proportions REP
graph2 = ax.bar([0.2, 1.4, 2.6, 3.8, 5], REP_freq, width, color = '#8856a7', edgecolor = 'black', linewidth = 1)
# plot SNP proportions miRNAs
graph3 = ax.bar([0.4, 1.6, 2.8, 4, 5.2], mirna_freq, width, color = '#8c96c6', edgecolor = 'black', linewidth = 1)
# plot SNP proportions for resampled SNPs near mirnas
graph4 = ax.bar([0.6, 1.8, 3, 4.2, 5.4], mirna_samplefreq, width, yerr = mirna_sem, color = '#9ebcda', 
                edgecolor = 'black', linewidth = 1,
                error_kw=dict(elinewidth=1, ecolor='black', markeredgewidth = 1))
# plot SNP proportions for resampled SNPs near miRNAs
graph5 = ax.bar([0.8, 2, 3.2, 4.4, 5.6], targets_freq, width, color = '#bfd3e6', edgecolor = 'black', linewidth = 1)
# plot SNP proportions for resampled SNPs in UTRs
graph6 = ax.bar([1, 2.2, 3.4, 4.6, 5.8], target_samplefreq, width, yerr = target_sem, color = '#edf8fb',
                edgecolor = 'black', linewidth = 1,
                error_kw=dict(elinewidth=1, ecolor='black', markeredgewidth = 1))



#
#width = 0.2
#
#print('REP', MAF_REP_hist[0], MAF_REP_hist[1])
#
#ax.bar([0, 0.4, 0.8, 1.2, 1.6], repfreq, width, yerr = [0.1, 0.2, 0.05, 0, 0.3], color = '#810f7c',
#       error_kw=dict(elinewidth=1, ecolor='black', markeredgewidth = 1))
#ax.bar([0.2, 0.6, 1, 1.4, 1.8], synfreq, width, color = '#8856a7')

ax.set_ylabel('Proportion of SNPs', size = 10, ha = 'center', fontname = 'Arial')

# determine tick position on x axis
xpos =  [0, 1.2, 2.4, 3.6, 4.8, 6]
xtext = [0, 10, 20, 30, 40, 50]
xtext = list(map(lambda x : str(x), xtext))
# set up tick positions and labels
plt.xticks(xpos, xtext, fontsize = 10, fontname = 'Arial')
plt.yticks(fontsize = 0)

# set x axis label
ax.set_xlabel('Minor Allele Frequency (%)', size = 10, ha = 'center', fontname = 'Arial')

# do not show lines around figure, keep bottow line  
ax.spines["top"].set_visible(False)    
ax.spines["bottom"].set_visible(True)    
ax.spines["right"].set_visible(False)    
ax.spines["left"].set_visible(False)      
# offset the spines
for spine in ax.spines.values():
  spine.set_position(('outward', 5))
  
# add a light grey horizontal grid to the plot, semi-transparent, 
ax.yaxis.grid(True, linestyle='--', which='major', color='lightgrey', alpha=0.5, linewidth = 0.5)  
# hide these grids behind plot objects
ax.set_axisbelow(True)

# do not show ticks on 1st graph
ax.tick_params(
    axis='x',       # changes apply to the x-axis and y-axis (other option : x, y)
    which='both',      # both major and minor ticks are affected
    bottom='on',      # ticks along the bottom edge are off
    top='off',         # ticks along the top edge are off
    right = 'off',
    left = 'off',          
    labelbottom='on', # labels along the bottom edge are off 
    colors = 'black',
    labelsize = 10,
    direction = 'out') # ticks are outside the frame when bottom = 'on

# do not show ticks
ax.tick_params(
    axis='y',       # changes apply to the x-axis and y-axis (other option : x, y)
    which='both',      # both major and minor ticks are affected
    bottom='off',      # ticks along the bottom edge are off
    top='off',         # ticks along the top edge are off
    right = 'off',
    left = 'off',          
    labelbottom='off', # labels along the bottom edge are off 
    colors = 'black',
    labelsize = 10,
    direction = 'out') # ticks are outside the frame when bottom = 'on

for label in ax.get_yticklabels():
    label.set_fontname('Arial')

# create legend
syn = mpatches.Patch(color = '#810f7c' , edgecolor = 'black', fontsize = 8, label= 'Syn')
rep = mpatches.Patch(color = '#8856a7', edgecolor = 'black', fontsize = 8, label = 'Rep')
mirna = mpatches.Patch(color = '#8c96c6', edgecolor = 'black', fontsize = 8, label = 'miRNA')
nearmirna = mpatches.Patch(color = '#9ebcda', edgecolor = 'black', fontsize = 8, label = 'near miRNA')
targets = mpatches.Patch(color = '#bfd3e6', edgecolor = 'black', fontsize = 8, label = 'target')
utr = mpatches.Patch(color = '#edf8fb', edgecolor = 'black', fontsize = 8, label = 'UTR')
plt.legend(handles=[syn, rep, mirna, nearmirna, targets, utr], loc = 1, fontsize = 8, frameon = False)

# add margin on the x-axis
plt.margins(0.05)

fig.savefig('MAFmiRNADistribution.pdf', bbox_inches = 'tight')
