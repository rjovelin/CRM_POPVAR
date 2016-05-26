# -*- coding: utf-8 -*-
"""
Created on Thu May 26 07:21:34 2016

@author: Richard
"""



# use this script to plot the cumulative distribution function of the MAF distribution for miRNAs and other sites

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

# sort MAF values
MAF_REP.sort()
MAF_SYN.sort()
MAF_mirna.sort()
MAF_targets.sort()
print('values are sorted')

# create figure
fig = plt.figure(1, figsize = (2.5, 1.5))
# add a plot to figure (1 row, 1 column, 1 plot)
ax = fig.add_subplot(1, 1, 1)  

# plot MAF synonymous sites
graph1 = ax.step(MAF_SYN, np.linspace(0, 1, len(MAF_SYN), endpoint=False), linewidth = 1.2, color = '#33a02c', alpha = 0.5)
# plot MAF replacement sites
graph2 = ax.step(MAF_REP, np.linspace(0, 1, len(MAF_REP), endpoint=False), linewidth = 1.2, color = '#b2df8a', aplha = 0.5)
# plot MAF miRNAs
graph3 = ax.step(MAF_mirna, np.linspace(0, 1, len(MAF_mirna), endpoint=False), linewidth = 1.2, color = '#1f78b4', alpha = 0.5)
# plot MAF targets
graph4 = ax.step(MAF_targets, np.linspace(0, 1, len(MAF_targets), endpoint=False), linewidth = 1.2, color = '#a6cee3', alpha = 0.5)
print('plotted CDF')

# add label for the Y axis
ax.set_ylabel('Proportion of SNPs', size = 10, ha = 'center', fontname = 'Arial')
# set x axis label
ax.set_xlabel('Minor Allele Frequency', size = 10, ha = 'center', fontname = 'Arial')

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

# add lines
lns = graph1+graph2+graph3+graph4
# get labels
labs = ['Synonymous', 'Replacement', 'miRNAs', 'targets']
# plot legend
ax.legend(lns, labs, loc=2, fontsize = 8, frameon = False)

fig.savefig('CDFMAFmiRNAs.pdf', bbox_inches = 'tight')



