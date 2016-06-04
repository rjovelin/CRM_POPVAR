# -*- coding: utf-8 -*-
"""
Created on Sat Jun  4 10:42:53 2016

@author: Richard
"""

# use this script to compare gene expression between chemoreceptor genes and non-chemoreceptor genes

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
from diversity_chemo_partitions import *
from genomic_coordinates import *
from protein_divergence import *
from crm_expression import *


# get gene expression (average across developmental time points) {CRE_gene_ID: mean expression}
Expression = expression_developmental_stages('../Genome_Files/WBPaper00041190.cre.mr.csv', '../Genome_Files/c_remanei.PRJNA53967.WS248.geneIDs.txt')
print('parsed expression level')
 
# get the set of chemoreceptors from the iprscan outputfile
chemo = get_chemoreceptors('../Genome_Files//PX356_protein_seq.tsv') 
print('parsed chemo genes')

# create a set of valid transcripts
transcripts = get_valid_transcripts('../Genome_Files/unique_transcripts.txt')
print('got list of valid transcripts')

# create a set of valid chemoreceptors
GPCRs = set(gene for gene in chemo if gene in transcripts)
print('GPCRS', len(GPCRs))
print('made set of valid chemo genes')

NonGPCRs = set(gene for gene in transcripts if gene not in chemo)
print('NonGPCRs', len(NonGPCRs))
print('made set of valid non-chemo genes')

# create lists of expression for chemo and nonchemo genes
# make lists of dN values for chemo and non-chemo genes
chemoExp, NCExp = [], []

# populate lists
for gene in Expression:
    if gene in GPCRs:
        chemoExp.append(Expression[gene])
    elif gene in NonGPCRs:
        NCExp.append(Expression[gene])
print('populated lists with expression values')
    
a = [chemoExp, NCExp]
b = ['Chemo', 'NC']
for i in range(len(a)):
    print('N', b[i], len(a[i]))
for i in range(len(a)):
    print('max', b[i], max(a[i]))

# compare expression for chemo and non-chemo using Wilcoxon sum rank tests
P = stats.ranksums(chemoExp, NCExp)[1]
print('P', P)

# create a list of means
Means = [np.mean(chemoExp), np.mean(NCExp)]

# create a lit of SEM
SEM = [np.std(chemoExp) / math.sqrt(len(chemoExp)),
       np.std(NCExp) / math.sqrt(len(NCExp))]

# create figure
fig = plt.figure(1, figsize = (1, 2))
# add a plot to figure (1 row, 1 column, 1 plot)
ax = fig.add_subplot(1, 1, 1)  

# set width of bar
width = 0.2
# set colors
colorscheme = ['#8856a7', '#9ebcda']

# plot nucleotide divergence
ax.bar([0, 0.2], Means, width, yerr = SEM, color = colorscheme, 
                edgecolor = 'black', linewidth = 1,
                error_kw=dict(elinewidth=1, ecolor='black', markeredgewidth = 1))

ax.set_ylabel('Expression level', size = 10, ha = 'center', fontname = 'Arial')

# set y limits
plt.ylim([0, 10])

# determine tick position on x axis
xpos =  [0.1, 0.3]
xtext = ['GPCR', 'NC']
# set up tick positions and labels
plt.xticks(xpos, xtext, fontsize = 10, fontname = 'Arial')

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

# add margin on the x-axis
plt.margins(0.05)

# I already determined that all site categories are significantly different
# using Wilcoxon rank sum tests, so we need now to add letters to show significance
# P_rep, P_syn and P_omega < 0.001 ---> P = ***
P = '***'

# annotate figure to add significance
# add bracket
ax.annotate("", xy=(0.1, 8.2), xycoords='data',
            xytext=(0.3, 8.2), textcoords='data',
            arrowprops=dict(arrowstyle="-", ec='#aaaaaa', connectionstyle="bar,fraction=0.2", linewidth = 1))
# add stars for significance
ax.text(0.2, 8.8, P, horizontalalignment='center',
        verticalalignment='center', color = 'grey', fontname = 'Arial', size = 6)

fig.savefig('ExpressionChemoNonChemo.pdf', bbox_inches = 'tight')

