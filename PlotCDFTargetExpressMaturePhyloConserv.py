# -*- coding: utf-8 -*-
"""
Created on Thu May 26 22:04:23 2016

@author: Richard
"""

# use this script to plot the cumulative density function of the expression of target genes
# for different level of mature miRNA phylogenetic conservation



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
from crm_expression import *

# get the set of valid transcripts
valid_transcripts = get_valid_transcripts('../Genome_Files/unique_transcripts.txt')
print('made list of valid transcripts')

# get the target genes for each miRNA family
# make a dictionary of {seed : set(target genes)}
TargetGenes = {}
infile = open('Cremanei_miRNA_sites.txt')
infile.readline()
for line in infile:
    line = line.rstrip()
    if line != '':
        line = line.split('\t')
        # get gene and seed, convert seed to DNA
        gene, seed = line[0], line[1].replace('U', 'T')
        # add only valid transcripts
        if gene in valid_transcripts:
            if seed in TargetGenes:
                TargetGenes[seed].add(gene)
            else:
                TargetGenes[seed] = set()
                TargetGenes[seed].add(gene)
infile.close()
print('got target genes')


# determine the conservation level of each mature miR
# create a dict {seed: conservation level}
PhyloCons = {}
infile = open('CRM_miRNAsCoordinatesFinal.txt')
infile.readline()
for line in infile:
    line = line.rstrip()
    if line != '':
        line = line.split('\t')
        seed = line[7]
        conservation = line[-1]
        PhyloCons[seed] = conservation
infile.close()
print('got conservation')

# create a dict of {gene: expression leve;}
GeneExpression = expression_developmental_stages('../Genome_Files/WBPaper00041190.cre.mr.csv', '../Genome_Files/c_remanei.PRJNA53967.WS248.geneIDs.txt')
print('computed expression level of target genes')


# create lists of expression of mature miRNAs with different level of conservation
conserved, restricted, specific = [], [], []
# loop over seed in target genes dict
for seed in TargetGenes:
    # check conservation
    if PhyloCons[seed] == 'Caeno':
        # add expression of target genes to list conserved
        for gene in TargetGenes[seed]:
            if gene in GeneExpression:
                conserved.append(GeneExpression[gene])
    elif PhyloCons[seed] == 'CrmCla':
        # add expression of target genes to list restricted
        for gene in TargetGenes[seed]:
            if gene in GeneExpression:
                restricted.append(GeneExpression[gene])
    elif PhyloCons[seed] == 'Crm':
        # add expression of target genes to list specific
        for gene in TargetGenes[seed]:
            if gene in GeneExpression:
                specific.append(GeneExpression[gene])
print('sorted expression according to marure miR phylogenetic conservation')
print('conserved', len(conserved))
print('restricted', len(restricted))
print('specific', len(specific))


# sort expression levels
conserved.sort()
restricted.sort()
specific.sort()
print('values are sorted')


# create figure
fig = plt.figure(1, figsize = (3.5, 2.5))
# add a plot to figure (1 row, 1 column, 1 plot)
ax = fig.add_subplot(1, 1, 1)  


# plot expression of conserved miR target genes
graph1 = ax.step(conserved, np.linspace(0, 1, len(conserved), endpoint=False), linewidth = 1.2, color = '#7fc97f', alpha = 0.7)
# plot expression of restricted miR target genes
graph2 = ax.step(restricted, np.linspace(0, 1, len(restricted), endpoint=False), linewidth = 1.2, color = '#beaed4', alpha = 0.7)
# plot expression of specific miR target genes
graph3 = ax.step(specific, np.linspace(0, 1, len(specific), endpoint=False), linewidth = 1.2, color = '#fdc086', alpha = 0.7)
print('plotted CDF')

# add label for the Y axis
ax.set_ylabel('Proportion of genes', size = 10, ha = 'center', fontname = 'Arial')
# set x axis label
ax.set_xlabel('Target gene expression level', size = 10, ha = 'center', fontname = 'Arial')

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
lns = graph1+graph2+graph3
# get labels
labs = ['Conserved', 'Restricted', 'Specific']
# plot legend
ax.legend(lns, labs, loc=2, fontsize = 8, frameon = False)

fig.savefig('testfile.pdf', bbox_inches = 'tight')



