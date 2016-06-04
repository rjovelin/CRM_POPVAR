# -*- coding: utf-8 -*-
"""
Created on Fri Jun  3 16:39:34 2016

@author: RJovelin
"""

# use this script to compare divergence between chemoreceptor genes and non-chemoreceptor genes

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

# create lists of divergence for chemo and nonchemo genes
# make lists of dN values for chemo and non-chemo genes
chemodN, NCdN = [], []
#make lists of dS values for chemo and non-chemo genes
chemodS, NCdS = [], []
# make lists of omega values for chemo and non-chemo genes
chemoomega, NComega = [], [] 

# populate lists
for gene in ProtDiverg:
    if gene in GPCRs:
        chemodN.append(ProtDiverg[gene][0])
        chemodS.append(ProtDiverg[gene][1])
        if ProtDiverg[gene][-1] != 'NA':
            chemoomega.append(ProtDiverg[gene][-1])
    elif gene in NonGPCRs:
        NCdN.append(ProtDiverg[gene][0])
        NCdS.append(ProtDiverg[gene][1])
        if ProtDiverg[gene][-1] != 'NA':
            NComega.append(ProtDiverg[gene][-1])
print('populated lists with divergence values')
    

a = [chemodN, NCdN, chemodS, NCdS, chemoomega, NComega]
b = ['dN-Chemo', 'dN-NC', 'dS-Chemo', 'dS-NC', 'omega-Chemo', 'omega-NC']
for i in range(len(a)):
    print('N', b[i], len(a[i]))
print('\n')
for i in range(len(a)):
    print('max', b[i], max(a[i]))

# compare divergence for chemo and non-chemo using wilcoxon sum rank tests
P_rep = stats.ranksums(chemodN, NCdN)[1]
P_syn = stats.ranksums(chemodS, NCdS)[1]
P_omega = stats.ranksums(chemoomega, NComega)[1]

print('P_rep', P_rep)
print('P_syn', P_syn)
print('P_omega', P_omega)

# create a list of means
Means = [np.mean(chemodN), np.mean(NCdN),
         np.mean(chemodS), np.mean(NCdS),
         np.mean(chemoomega), np.mean(NComega)]

# create a lit of SEM
SEM = [np.std(chemodN) / math.sqrt(len(chemodN)),
       np.std(NCdN) / math.sqrt(len(NCdN)),
       np.std(chemodS) / math.sqrt(len(chemodS)),
       np.std(NCdS) / math.sqrt(len(NCdS)),
       np.std(chemoomega) / math.sqrt(len(chemoomega)),
       np.std(NComega) / math.sqrt(len(NComega))]

# create figure
fig = plt.figure(1, figsize = (3, 2))
# add a plot to figure (1 row, 1 column, 1 plot)
ax = fig.add_subplot(1, 1, 1)  

# set width of bar
width = 0.2
# set colors
colorscheme = ['#2ca25f', '#99d8c9','#2ca25f', '#99d8c9', '#2ca25f', '#99d8c9']

# plot nucleotide divergence
ax.bar([0, 0.2, 0.5, 0.7, 1, 1.2], Means, width, yerr = SEM, color = colorscheme, 
                edgecolor = 'black', linewidth = 1,
                error_kw=dict(elinewidth=1, ecolor='black', markeredgewidth = 1))

ax.set_ylabel('Nucleotide divergence', size = 10, ha = 'center', fontname = 'Arial')

# set y limits
plt.ylim([0, 0.31])

# determine tick position on x axis
xpos =  [0.2, 0.7, 1.2]
xtext = ['dN', 'dS', 'omega']
# set up tick positions and labels
plt.xticks(xpos, xtext, fontsize = 10, fontname = 'Arial')

# set x axis label
ax.set_xlabel('Sites in coding sequences', size = 10, ha = 'center', fontname = 'Arial')

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

# create legend
ChemoGene = mpatches.Patch(facecolor = '#2ca25f' , edgecolor = 'black', linewidth = 1, label= 'GPCR')
NCGene = mpatches.Patch(facecolor = '#99d8c9', edgecolor = 'black', linewidth = 1, label = 'NC')
plt.legend(handles=[ChemoGene, NCGene], loc = 2, fontsize = 8, frameon = False)

# I already determined that all site categories are significantly different
# using Wilcoxon rank sum tests, so we need now to add letters to show significance
# P_rep, P_syn and P_omega < 0.001 ---> P = ***
P = '***'

# annotate figure to add significance
# add bracket
ax.annotate("", xy=(0.1, 0.08), xycoords='data',
            xytext=(0.3, 0.08), textcoords='data',
            arrowprops=dict(arrowstyle="-", ec='#aaaaaa', connectionstyle="bar,fraction=0.2", linewidth = 1))
# add stars for significance
ax.text(0.2, 0.10, P, horizontalalignment='center',
        verticalalignment='center', color = 'grey', fontname = 'Arial', size = 6)

ax.annotate("", xy=(0.6, 0.27), xycoords='data',
            xytext=(0.8, 0.27), textcoords='data',
            arrowprops=dict(arrowstyle="-", ec='#aaaaaa', connectionstyle="bar,fraction=0.2", linewidth = 1))
# add stars for significance
ax.text(0.7, 0.29, P, horizontalalignment='center',
        verticalalignment='center', color = 'grey', fontname = 'Arial', size = 6)

ax.annotate("", xy=(1.1, 0.23), xycoords='data',
            xytext=(1.3, 0.23), textcoords='data',
            arrowprops=dict(arrowstyle="-", ec='#aaaaaa', connectionstyle="bar,fraction=0.2", linewidth = 1))
# add stars for significance
ax.text(1.2, 0.25, P, horizontalalignment='center',
        verticalalignment='center', color = 'grey', fontname = 'Arial', size = 6)

fig.savefig('ChemoNonChemoDivergence.pdf', bbox_inches = 'tight')

