# -*- coding: utf-8 -*-
"""
Created on Wed Jun  8 12:27:25 2016

@author: Richard
"""


# use this script to compare diversity between chemoreceptor genes and non-chemoreceptor genes

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
from divergence import *


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


# compute theta at synonymous sites
chemo_syn = compute_theta_diversity('../Genome_Files/CDS_SNP_DIVERG.txt', GPCRs, 'SYN', 10)
# make a list of theta values
chemoSyn = []
for gene in chemo_syn:
    chemoSyn.append(chemo_syn[gene])

nonchemo_syn = compute_theta_diversity('../Genome_Files/CDS_SNP_DIVERG.txt', NonGPCRs, 'SYN', 10)
# make a list of theta values
NCSyn = []
for gene in nonchemo_syn:
    NCSyn.append(nonchemo_syn[gene])

print('done computing theta at synonymous sites')
print('Sites\tmean\tmin\tmax')
print(' chemo SYN', np.mean(chemoSyn), min(chemoSyn), max(chemoSyn))
print('non-chemo SYN', np.mean(NCSyn), min(NCSyn), max(NCSyn))


# compute theta at nonsynonymous sites
chemo_rep = compute_theta_diversity('../Genome_Files/CDS_SNP_DIVERG.txt', GPCRs, 'REP', 10)
# make a list of theta values
chemoRep = []
for gene in chemo_rep:
    chemoRep.append(chemo_rep[gene])

# compute theta at nonsynonymous sites
nonchemo_rep = compute_theta_diversity('../Genome_Files/CDS_SNP_DIVERG.txt', NonGPCRs, 'REP', 10)
# make a list of theta values
NCRep = []
for gene in nonchemo_rep:
    NCRep.append(nonchemo_rep[gene])

print('done computing theta at replacement sites')
print('Sites\tmean\tmin\tmax')
print('chemo REP', np.mean(chemoRep), min(chemoRep), max(chemoRep))
print('non-chemo REP', np.mean(NCRep), min(NCRep), max(NCRep))

# compute the ratio of diversity
chemoRatio = []
for gene in chemo_rep:
    if gene in chemo_syn and chemo_syn[gene] != 0:
        chemoRatio.append(chemo_rep[gene] / chemo_syn[gene])
NCRatio = []
for gene in nonchemo_rep:
    if gene in nonchemo_syn and nonchemo_syn[gene] != 0:
        NCRatio.append(nonchemo_rep[gene] / nonchemo_syn[gene])
        
# compare divergence for chemo and non-chemo using wilcoxon sum rank tests
P_rep = stats.ranksums(chemoRep, NCRep)[1]
P_syn = stats.ranksums(chemoSyn, NCSyn)[1]
P_omega = stats.ranksums(chemoRatio, NCRatio)[1]

print('P_rep', P_rep)
print('P_syn', P_syn)
print('P_omega', P_omega)

# create a list of means
Means = [np.mean(chemoRep), np.mean(NCRep),
         np.mean(chemoSyn), np.mean(NCSyn),
         np.mean(chemoRatio), np.mean(NCRatio)]

# create a lit of SEM
SEM = [np.std(chemoRep) / math.sqrt(len(chemoRep)),
       np.std(NCRep) / math.sqrt(len(NCRep)),
       np.std(chemoSyn) / math.sqrt(len(chemoSyn)),
       np.std(NCSyn) / math.sqrt(len(NCSyn)),
       np.std(chemoRatio) / math.sqrt(len(chemoRatio)),
       np.std(NCRatio) / math.sqrt(len(NCRatio))]

# create figure
fig = plt.figure(1, figsize = (3, 2))
# add a plot to figure (1 row, 1 column, 1 plot)
ax = fig.add_subplot(1, 1, 1)  


# set width of bar
width = 0.2
# set colors
colorscheme = ['#8856a7', '#9ebcda','#8856a7', '#9ebcda', '#8856a7', '#9ebcda']

# plot nucleotide divergence
ax.bar([0, 0.2, 0.5, 0.7, 1, 1.2], Means, width, yerr = SEM, color = colorscheme, 
                edgecolor = 'black', linewidth = 1,
                error_kw=dict(elinewidth=1, ecolor='black', markeredgewidth = 1))

ax.set_ylabel('Nucleotide divergence', size = 10, ha = 'center', fontname = 'Arial')

# set y limits
plt.ylim([0, 0.36])

# determine tick position on x axis
xpos =  [0.2, 0.7, 1.2]
xtext = ['$' + chr(952) +'_{rep}$', '$' + chr(952) + '_{syn}$', '$' + chr(952) + '_{syn}$' + '/' + '$' + chr(952) + '_{syn}$']
# set up tick positions and labels
#plt.xticks(xpos, xtext, fontsize = 10, fontname = 'Arial')
plt.xticks(xpos, xtext, fontsize = 10)
# set x axis label
#ax.set_xlabel('Sites in coding sequences', size = 10, ha = 'center', fontname = 'Arial')

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
ChemoGene = mpatches.Patch(facecolor = '#8856a7' , edgecolor = 'black', linewidth = 1, label= 'GPCR')
NCGene = mpatches.Patch(facecolor = '#9ebcda', edgecolor = 'black', linewidth = 1, label = 'NC')
plt.legend(handles=[ChemoGene, NCGene], loc = 2, fontsize = 8, frameon = False)

# I already determined that all site categories are significantly different
# using Wilcoxon rank sum tests, so we need now to add letters to show significance
# P_rep, P_syn and P_omega < 0.001 ---> P = ***
P = '***'

# annotate figure to add significance
# add bracket
ax.annotate("", xy=(0.6, 0.08), xycoords='data',
            xytext=(0.8, 0.08), textcoords='data',
            arrowprops=dict(arrowstyle="-", ec='#aaaaaa', connectionstyle="bar,fraction=0.2", linewidth = 1))
# add stars for significance
ax.text(0.7, 0.10, P, horizontalalignment='center',
        verticalalignment='center', color = 'grey', fontname = 'Arial', size = 6)

ax.annotate("", xy=(1.1, 0.32), xycoords='data',
            xytext=(1.3, 0.32), textcoords='data',
            arrowprops=dict(arrowstyle="-", ec='#aaaaaa', connectionstyle="bar,fraction=0.2", linewidth = 1))
# add stars for significance
ax.text(1.2, 0.35, P, horizontalalignment='center',
        verticalalignment='center', color = 'grey', fontname = 'Arial', size = 6)

fig.savefig('testfile.pdf', bbox_inches = 'tight')

