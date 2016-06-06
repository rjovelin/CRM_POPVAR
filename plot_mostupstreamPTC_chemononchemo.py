# -*- coding: utf-8 -*-
"""
Created on Sun Jun  5 23:29:57 2016

@author: Richard
"""

# use this script to plot the distribution of the 5' most upstream PTC allele
# in chemo and non chemo genes 


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
import random
# import custom modules
from chemoreceptors import *
from manipulate_sequences import *
from mk_test import *
from premature_stops import *


# create a set of valid transcripts
transcripts = get_valid_transcripts('../Genome_Files/unique_transcripts.txt')
print('got list of valid transcripts')

# get the set of chemoreceptors from the iprscan outputfile
chemo = get_chemoreceptors('../Genome_Files//PX356_protein_seq.tsv') 
print('parsed chemo genes')

# create a set of valid chemoreceptors
GPCRs = set(gene for gene in chemo if gene in transcripts)
print('GPCRS', len(GPCRs))
print('made set of valid chemo genes')

# create a set of non-chemo gene
NonGPCRs = set(gene for gene in transcripts if gene not in chemo)
print('NonGPCRs', len(NonGPCRs))
print('made set of valid non-chemo genes')

# get the relative position of the 5' most PTC allele in chemo and non-chemo genes
ChemoUpstream = upstream_PTC = position_first_PTC('../Genome_Files/CDS_SNP_DIVERG.txt',
                                                  GPCRs, '../Genome_Files/transcripts_indels_CDS.txt',
                                                  '../Genome_Files/noamb_PX356_all_CDS.fasta')
print('grabbed positions of PTC alleles for chemo genes')

NCUpstream = upstream_PTC = position_first_PTC('../Genome_Files/CDS_SNP_DIVERG.txt',
                                                  NonGPCRs, '../Genome_Files/transcripts_indels_CDS.txt',
                                                  '../Genome_Files/noamb_PX356_all_CDS.fasta')
print('grabbed positions of PTC alleles for non-chemo genes')

# express relative position in % of the CDS length
for i in range(len(ChemoUpstream)):
    ChemoUpstream[i] = ChemoUpstream[i] * 100
for i in range(len(NCUpstream)):
    NCUpstream[i] = NCUpstream[i] * 100

# make a histogram
ChemoHist = np.histogram(ChemoUpstream, range(0, 110, 10))
NCHist = np.histogram(NCUpstream, range(0, 110, 10))

# transform the gene counts to proportions
ChemoFreq = [i / sum(ChemoHist[0]) for i in ChemoHist[0]]
NCFreq = [i / sum(NCHist[0]) for i in NCHist[0]]

print(len(NCFreq))
print(NCFreq)

# test uniformity of the distribution of the 5' most upstream PTC
chemotest = stats.chisquare(ChemoHist[0])
NCtest = stats.chisquare(NCHist[0])
print('chemo Pval uniformity', chemotest[1])
print('non-chemo Pval uniformity', NCtest[1])

# create figure
fig = plt.figure(1, figsize = (4, 2))
# add a plot to figure (1 row, 1 column, 1 plot)
ax = fig.add_subplot(1, 1, 1)  

# set width of bar
width = 0.1

# plot positions chemo
graph1 = ax.bar([i / 10 for i in range(10)], ChemoFreq, width, color = '#de2d26', edgecolor = 'black', linewidth = 1, alpha = 0.5)
# plot positions non-chemo
graph2 = ax.bar([i / 10 for i in range(10)], NCFreq, width, color = '#3182bd', edgecolor = 'grey', linewidth = 1, alpha = 0.5)

ax.set_ylabel('Proportion of genes with a PTC', size = 10, ha = 'center', fontname = 'Arial')

# determine tick position on x axis
xpos =  [0, 0.2, 0.4, 0.6, 0.8, 1, 1.2, 1.4, 1,6, 1.8, 2]
xtext = [0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
xtext = list(map(lambda x : str(x), xtext))
# set up tick positions and labels
plt.xticks(xpos, xtext, fontsize = 10, fontname = 'Arial')

# set x axis label
ax.set_xlabel('Decile of CDS length', size = 10, ha = 'center', fontname = 'Arial')

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
chemoptc = mpatches.Patch(facecolor = '#de2d26' , edgecolor = 'grey', linewidth = 1, label= 'GPCR')
NCptc = mpatches.Patch(facecolor = '#3182bd', edgecolor = 'grey', linewidth = 1, label = 'NC')
plt.legend(handles=[chemoptc, NCptc], loc = 1, fontsize = 8, frameon = False)

# add margin on the x-axis
plt.margins(0.05)

fig.savefig('testfile.pdf', bbox_inches = 'tight')