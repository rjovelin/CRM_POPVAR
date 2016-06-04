# -*- coding: utf-8 -*-
"""
Created on Sat Jun  4 17:54:41 2016

@author: Richard
"""


# use this script to compare divergence dS between chemoreceptor genes of different families 

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
# import custom modules
from chemoreceptors import *
from manipulate_sequences import *
from genomic_coordinates import *
from multiple_testing import *


# use this function to get the chemoreceptor genes in each chemoreceptor family
Families = chemo_families('../Genome_Files/PX356_protein_seq.tsv')
print('assigned GPCRs to gene families')

# parse protein divergence
ProtDiverg = {}
infile = open('../CREM_CLA_protein_divergence/CRM_CLA_ProtDiverg_FILTERED.txt')
infile.readline()
for line in infile:
    line = line.rstrip()
    if line != '':
        line = line.split('\t')
        gene = line[0]
        dN = float(line[4])
        dS = float(line[5])
        if dS == 0:
            omega = 'NA'
        else:
            omega = float(line[6])
        ProtDiverg[gene] = [dN, dS, omega]
print('parse divergence table')

# create a set of valid transcripts
transcripts = get_valid_transcripts('../Genome_Files/unique_transcripts.txt')
print('got list of valid transcripts')

# remove non valid genes
for family in Families:
    to_remove = []
    for gene in Families[family]:
        if gene not in transcripts:
            to_remove.append(gene)
    for gene in to_remove:
        Families[family].remove(gene)
print('removed non genes in chemo families')
    
# create dicts with divergence with {family name : [list of divergence]}
dS = {}

missing = {}

# populate dicts 
for family in Families:
    if family not in dS:
        # initialize list
        dS[family] = []
    # add dS values for all chemo genes in that family
    for gene in Families[family]:
        if gene in ProtDiverg:
            dS[family].append(ProtDiverg[gene][1])
        else:
            if family not in missing:
                missing[family] = []
            missing[family].append(gene)

print('missing', len(missing))
for family in missing:
    print(family, len(missing[family]))

for family in dS:
    print('dS', family, np.mean(dS[family]), max(dS[family]))

  
# create a list of [mean, SEM, family] for each family
MeanFam = []
for family in Families:
    MeanFam.append([np.mean(dS[family]), np.std(dS[family]) / math.sqrt(len(dS[family])), family])
# sort according to mean from highest to lowest mean
MeanFam.sort()
MeanFam.reverse()
    
# create parallel lists of means, SEN amd family names, sorted accoring to mean values
Means, SEM, FamNames = [], [], []
for i in range(len(MeanFam)):
    Means.append(MeanFam[i][0])
    SEM.append(MeanFam[i][1])
    FamNames.append(MeanFam[i][2])
print('created mean and SEM lists sorted according to means')

# create a list of list of DNN values, according to the family in FamNames
a = []
for family in FamNames:
    a.append(dN[family])


# compare divergence between chemo family wilcoxon sum rank tests
# create a dict of {family-family : P value} for each pairwise comparison
Pval = {}
for i in range(0, len(a) -1):
    for j in range(i+1, len(a)):
        P = stats.ranksums(a[i], a[j])[1]
        comparison = FamNames[i] + '_' + FamNames[j]
        Pval[comparison] = P
# make a list of [P, comparison]
Pvalues = [[val, key] for key, val in Pval.items()]
print('computed P values for pairwise differences')

# apply Benjamini_Hochberg correction, get a dict of comparison: adjusted Pval
CorrectedPval = Benjamini_Hochberg_correction(Pvalues)
print('corrected P values with Benjamini_Hochberg')

for comparison in CorrectedPval:
    # consider significance with FDR = 0.05
    if CorrectedPval[comparison] < 0.05:
        print(comparison, CorrectedPval[comparison])


# create figure
fig = plt.figure(1, figsize = (5, 2))
# add a plot to figure (1 row, 1 column, 1 plot)
ax = fig.add_subplot(1, 1, 1)  

# set width of bar
width = 0.2
# set colors
colorscheme = ['#9ecae1']

# plot nucleotide divergence
ax.bar([i / 10 for i in range(0, 38, 2)], Means, width, yerr = SEM, color = colorscheme, 
                edgecolor = 'black', linewidth = 1,
                error_kw=dict(elinewidth=1, ecolor='black', markeredgewidth = 1))

ax.set_ylabel('Nucleotide divergence (dN)', size = 10, ha = 'center', fontname = 'Arial')

# set y limits
plt.ylim([0, 0.20])

# set x axis label
ax.set_xlabel('Chemoreceptor families', size = 10, ha = 'center', fontname = 'Arial')

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
    labelsize = 8,
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
    labelsize = 8,
    direction = 'out') # ticks are outside the frame when bottom = 'on

# determine tick position on x axis
xpos =  [(i / 10) + 0.1 for i in range(0, 38, 2)]
xtext = FamNames
# set up tick positions and labels
plt.xticks(xpos, xtext, rotation = 30, ha = 'right', fontsize = 8, fontname = 'Arial')

# change font of y ticks
for label in ax.get_yticklabels():
    label.set_fontname('Arial')

# add margin on the x-axis
plt.margins(0.05)


fig.savefig('DivergencedSChemoFamilies.pdf', bbox_inches = 'tight')
