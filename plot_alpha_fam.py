# -*- coding: utf-8 -*-
"""
Created on Tue Jun  7 19:23:10 2016

@author: RJovelin
"""


# use this script to plot alpha Smith-EyreWalker 2002 for each family 



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
# import custom modules
from chemoreceptors import *
from manipulate_sequences import *
from genomic_coordinates import *
from multiple_testing import *
from mk_test import *


# get the chemoreceptor genes in each chemoreceptor family
Families = chemo_families('../Genome_Files/PX356_protein_seq.tsv')
print('assigned GPCRs to gene families')

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


# count divergence and polymorphisms in Ontario strains
# eliminate singleton polymorphisms
# dict in form {gene : [PN, PS, DN, DS]}        
MK = count_polym_diverg('../Genome_Files/CDS_SNP_DIVERG.txt', 'KSR_PX', 'count', 0, 1)

# remove genes that are not valid transcripts
to_remove = [gene for gene in MK if gene not in transcripts]
for gene in to_remove:
    del MK[gene]
print('removed {0} genes from MK analysis'.format(len(to_remove)))

# perform MK test using Fisher exact test
mk = MK_test(MK, 'fisher')

# make a list of genes with significant MK test
significant = [gene for gene in mk if mk[gene][-1] < 0.05]
# determine genes with significant MK that are under positive selection
negative, positive = NaturalSelection(MK, significant)
print('sorted significant genes based on selection pattern')

# count the genes under positive selection for each family
PositiveFam = {}
for family in Families:
    for gene in Families[family]:
        if gene in positive:
            if family in PositiveFam:
                PositiveFam[family] += 1
            else:
                PositiveFam[family] = 1
print('N genes positively selected without correction')
for family in PositiveFam:
    print(family, PositiveFam[family])

# apply a Bejamini-Hochberg correction for multiple testing {gene: corrected_Pval}
mk_p = {}
for gene in mk:
    mk_p[gene] = mk[gene][-1]
pvals = [(p, gene) for gene, p in mk_p.items()]
corrected = Benjamini_Hochberg_correction(pvals)

FDR = 0.1
signifcorr = [gene for gene in corrected if corrected[gene] < FDR]
# determine genes with significant MK that are under positive or negative selection
negative, positive = NaturalSelection(MK, signifcorr)

# count the genes under positive selection for each family
PositiveFam = {}
for family in Families:
    for gene in Families[family]:
        if gene in positive:
            if family in PositiveFam:
                PositiveFam[family] += 1
            else:
                PositiveFam[family] = 1
print('N genes positively selected after after BH with 10% FDR')
for family in PositiveFam:
    print(family, PositiveFam[family])


# compute alpha according to the Smith-EyreWalker 2002 method for each family

# generate dict of dicts with polymorphism amd divergence counts for each gene
PolymDivCounts = {}
for family in Families:
    for gene in Families[family]:
        if gene in MK:
            if family in PolymDivCounts:
                PolymDivCounts[family][gene] = list(MK[gene])
            else:
                PolymDivCounts[family] = {}
                PolymDivCounts[family][gene] = list(MK[gene])
print('recorded polym and diverg counts')

FamAlpha = {}
for family in PolymDivCounts:
    # compute alpha for each family
    alpha = ComputeAlphaSEW2002(PolymDivCounts[family], 1)
    FamAlpha[family] = alpha
    
# create a list of [alpha, family] for each family
Names = [[val, key] for key, val in FamAlpha.items()]
# sort according to mean from highest to lowest mean
Names.sort()
Names.reverse()
    
# create parallel lists of means, SEN amd family names, sorted accoring to mean values
Alpha, FamNames = [], []
for i in range(len(Names)):
    Alpa.append(Names[i][0])
    FamNames.append(Names[i][1])
print('created alpha and family names sorted according to alpha')


# create figure
fig = plt.figure(1, figsize = (5, 2))
# add a plot to figure (1 row, 1 column, 1 plot)
ax = fig.add_subplot(1, 1, 1)  

# set width of bar
width = 0.2

# plot nucleotide divergence
ax.bar([i / 10 for i in range(0, 38, 2)], Alpha, width, color = '#3182bd', edgecolor = 'black', linewidth = 1)

# write Y axis label
ax.set_ylabel('Alpha', size = 10, ha = 'center', fontname = 'Arial')

# set y limits
#plt.ylim([0, 0.20])

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

fig.savefig('testfile.pdf', bbox_inches = 'tight')
