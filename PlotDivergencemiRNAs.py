# -*- coding: utf-8 -*-
"""
Created on Tue Jun 28 21:28:04 2016

@author: Richard
"""


# use this script to plot divergence at synonymous, replacement sites and mirna hairpin and mature seq

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
from manipulate_sequences import *
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
        ProtDiverg[gene] = [dN, dS]
print('parse divergence table')

# create a set of valid transcripts
transcripts = get_valid_transcripts('../Genome_Files/unique_transcripts.txt')
print('got list of valid transcripts')

# remove genes that are not valid
to_remove = [i for i in ProtDiverg if i not in transcripts]
if len(to_remove) != 0:
    for gene in to_remove:
        del ProtDiverg[gene]
print('removed {0} genes from divergence dict'.format(len(to_remove)))

# create lists with dN and dS
dN = [ProtDiverg[gene][0] for gene in ProtDiverg]
dS = [ProtDiverg[gene][1] for gene in ProtDiverg]

# create lists of divergence for mature and hairpin
hairpin, mature = [], []

# extract divergence values from files
infile = open('CrmClamiRNAHairpinDivergence.txt')
infile.readline()
for line in infile:
    if line.rstrip() != '':
        line = line.rstrip().split('\t')
        hairpin.append(float(line[2]))
infile.close()

infile = open('CrmClamiRNAMatureDivergence.txt')
infile.readline()
for line in infile:
    if line.rstrip() != '':
        line = line.rstrip().split('\t')
        mature.append(float(line[2]))
infile.close()


# create lists of data and data type
data = [dS, dN, hairpin, mature]
names = ['Syn.', 'Rep.', 'miRNA', 'miR']
for i in range(len(data)):
    print('N', names[i], len(data[i]))
print('\n')
for i in range(len(data)):
    print('max', names[i], max(data[i]))

# compare data divergence pairwise
for i in range(0, len(data) -1):
    for j in range(i+1, len(data)):
        P = stats.ranksums(data[i], data[j])
        print(names[i] + '-' + names[j] + '\t' + str(P))

# create a list of means
Means = [np.mean(dS), np.mean(dN), np.mean(hairpin), np.mean(mature)]

# create a lit of SEM
SEM = [np.std(dS) / math.sqrt(len(dS)),
       np.std(dN) / math.sqrt(len(dN)),
       np.std(hairpin) / math.sqrt(len(hairpin)),
       np.std(mature) / math.sqrt(len(mature))]
       

# create figure
fig = plt.figure(1, figsize = (3, 2))
# add a plot to figure (1 row, 1 column, 1 plot)
ax = fig.add_subplot(1, 1, 1)  

# set width of bar
width = 0.2
# set colors
colorscheme = ['#33a02c', '#b2df8a', '#1f78b4', '#a6cee3']

# plot nucleotide divergence
ax.bar([0, 0.2, 0.4, 0.6], Means, width, yerr = SEM, color = colorscheme, 
                edgecolor = 'black', linewidth = 1,
                error_kw=dict(elinewidth=1, ecolor='black', markeredgewidth = 1))

ax.set_ylabel('Nucleotide divergence', size = 10, ha = 'center', fontname = 'Arial')

# set y limits
plt.ylim([0, 0.25])

# determine tick position on x axis
xpos =  [0.1, 0.3, 0.5, 0.7]
# set up tick positions and labels
plt.xticks(xpos, names, fontsize = 10, fontname = 'Arial')

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

# I already determined that only dN and miRNA are not significantly different from each other
# using Wilcoxon rank sum tests, so we need now to add letters to show significance

# annotate figure to add significance
# get the x and y coordinates
ypos = [0.247, 0.075, 0.06, 0.04]
xpos =  [0.1, 0.3, 0.5, 0.7]
diff = ['A', 'B', 'B', 'C']
for i in range(len(diff)):
    ax.text(xpos[i], ypos[i], diff[i], horizontalalignment='center',
            verticalalignment='center', color = 'black', fontname = 'Arial', size = 10)

fig.savefig('DivergencemiRNAs.pdf', bbox_inches = 'tight')
