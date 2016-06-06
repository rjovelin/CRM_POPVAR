# -*- coding: utf-8 -*-
"""
Created on Sun Jun  5 15:24:32 2016

@author: Richard
"""


# use this script to plot the Pn/Ps ratio and Dn/Ds ratio for chemo and non chemo genes
# and the distribution of these 2 ratios per gene


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

# count divergence and polymorphisms in Ontario strains
# eliminate singleton polymorphisms
# dict in form {gene : [PN, PS, DN, DS]}        
MK = count_polym_diverg('../Genome_Files/CDS_SNP_DIVERG.txt', 'KSR_PX', 'count', 0, 1)
print('counted number of syn and rep polymorphisms and divergences')

# remove genes that are not valid transcripts
to_remove = [gene for gene in MK if gene not in transcripts]
for gene in to_remove:
    del MK[gene]
print('removed {0} genes from MK analysis'.format(len(to_remove)))

# generate dicts with polymorphism amd divergence counts for chemo and nonchemo genes
ChemoPolymDivCounts, NCPolymDivCounts = {}, {}
for gene in MK:
    if gene in GPCRs:
        ChemoPolymDivCounts[gene] = list(MK[gene])
    elif gene in NonGPCRs:
        NCPolymDivCounts[gene] = list(MK[gene])
print('sorted snps and div counts for chemo and non-chemo')

# create lists of Pn/Ps for chemo and non chemo genes
chemoPolymRatio = sum([ChemoPolymDivCounts[gene][0] for gene in ChemoPolymDivCounts]) / sum([ChemoPolymDivCounts[gene][1] for gene in ChemoPolymDivCounts])
NCPolymRatio = sum([NCPolymDivCounts[gene][0] for gene in NCPolymDivCounts]) / sum([NCPolymDivCounts[gene][1] for gene in NCPolymDivCounts])
# create lists of Dn/Ds for chemo and non chemo genes
chemoDivRatio = sum([ChemoPolymDivCounts[gene][2] for gene in ChemoPolymDivCounts]) / sum([ChemoPolymDivCounts[gene][3] for gene in ChemoPolymDivCounts])
NCDivRatio = sum([NCPolymDivCounts[gene][2] for gene in NCPolymDivCounts]) / sum([NCPolymDivCounts[gene][3] for gene in NCPolymDivCounts])
print('computed polym and divergence ratios')
print('chemo polym ratio', chemoPolymRatio)
print('nonchemo polym ratio', NCPolymRatio)
print('chemo div ratio', chemoDivRatio)
print('nonchemo div ratio', NCDivRatio)

# create lists of ratio for individial genes if Ps and Ds are > 0
chemoGenePolym = [ChemoPolymDivCounts[gene][0] / ChemoPolymDivCounts[gene][1] for gene in ChemoPolymDivCounts if ChemoPolymDivCounts[gene][1] > 0]
NCGenePolym = [NCPolymDivCounts[gene][0] / NCPolymDivCounts[gene][1] for gene in NCPolymDivCounts if NCPolymDivCounts[gene][1] > 0]
chemoGeneDiv = [ChemoPolymDivCounts[gene][2] / ChemoPolymDivCounts[gene][3] for gene in ChemoPolymDivCounts if ChemoPolymDivCounts[gene][3] > 0]
NCGeneDiv = [NCPolymDivCounts[gene][2] / NCPolymDivCounts[gene][3] for gene in NCPolymDivCounts if NCPolymDivCounts[gene][3] > 0]
print('computed polyms and divergence ratios for each gene')


# create figure
fig = plt.figure(1, figsize = (3, 2))

# create subplot in figure
# add a plot to figure (N row, N column, plot N)
ax1 = fig.add_subplot(1, 2, 1)
# set width of bar
width = 0.2
# set colors
colorpolym = ['#fcae91', '#bdd7e7']
# plot ratio polym
graph1 = ax1.bar([0, 0.2], [chemoPolymRatio, NCPolymRatio], width, color = colorpolym,
                  edgecolor = 'black', linewidth = 1)
ax1.set_ylabel('Pn/Ps ratio', size = 10, ha = 'center', fontname = 'Arial')
# set y limits
#plt.ylim([0, 0.31])

# add a light grey horizontal grid to the plot, semi-transparent, 
ax1.yaxis.grid(True, linestyle='--', which='major', color='lightgrey', alpha=0.5, linewidth = 0.5)  
# hide these grids behind plot objects
ax1.set_axisbelow(True)

# do not show lines around figure, keep bottow line  
ax1.spines["top"].set_visible(False)    
ax1.spines["bottom"].set_visible(True)    
ax1.spines["right"].set_visible(False)    
ax1.spines["left"].set_visible(False)      
# offset the spines
for spine in ax1.spines.values():
  spine.set_position(('outward', 5))
  
# add a light grey horizontal grid to the plot, semi-transparent, 
ax1.yaxis.grid(True, linestyle='--', which='major', color='lightgrey', alpha=0.5, linewidth = 0.5)  
# hide these grids behind plot objects
ax1.set_axisbelow(True)

# do not show ticks on 1st graph
ax1.tick_params(
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
ax1.tick_params(
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

for label in ax1.get_yticklabels():
    label.set_fontname('Arial')

print('plotted 1st subplot')

# set up tick positions and labels
xpos =  [0.1, 0.3]
xtext = ['GPCRs', 'NC']
# set up tick positions and labels
plt.xticks(xpos, xtext, rotation = 30, ha = 'right', fontsize = 8, fontname = 'Arial')


# create subplot in figure
# add a plot to figure (N row, N column, plot N)
ax2 = fig.add_subplot(1, 2, 2)

colordiv = ['#cb181d', '#2171b5']     
# plot divergence ratio
graph2 = ax2.bar([0.1, 0.3], [chemoDivRatio, NCDivRatio], width, color = colordiv,
                 edgecolor = 'black', linewidth = 1)

# set up tick positions and labels
xpos =  [0.1, 0.3]
xtext = ['GPCRs', 'NC']
# set up tick positions and labels
plt.xticks(xpos, xtext, rotation = 30, ha = 'right', fontsize = 8, fontname = 'Arial')

# set y axis label
ax2.set_ylabel('Dn/Ds ratio', size = 10, ha = 'center', fontname = 'Arial')

# do not show lines around figure, keep bottow line  
ax2.spines["top"].set_visible(False)    
ax2.spines["bottom"].set_visible(False)    
ax2.spines["right"].set_visible(False)    
ax2.spines["left"].set_visible(False)      

# do not show ticks on 1st graph
ax2.tick_params(
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
ax2.tick_params(
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

for label in ax2.get_yticklabels():
    label.set_fontname('Arial')


## add the second subplot
#ax3 = fig.add_subplot(1, 2, 2)
#
## create a histogram for the distribution of ratios
#Bins = 20
#
## make histograms to get the counts in each MAF bins
#chemoPolymHist, a = np.histogram(chemoGenePolym, bins = Bins)
#NCPolymHist, b = np.histogram(NCGenePolym, bins = Bins)
#chemoDivHist, c = np.histogram(chemoGeneDiv, bins = Bins)
#NCDivHist, d = np.histogram(NCGeneDiv, bins = Bins)
#print('generated histograms')
#
### get the proportions of genes with given ratio
##chemoPolymFreq = [i / sum(chemoPolymHist[0]) for i in chemoPolymHist[0]]
##NCPolymFreq = [i / sum(NCPolymHist[0]) for i in NCPolymHist[0]]
##chemoDivFreq = [i / sum(chemoDivHist[0]) for i in chemoDivHist[0]]
##NCDivFreq = [i / sum(NCDivHist[0]) for i in NCDivHist[0]]
##print('computed proportions of genes in each bin')
#
## get the proportions of genes with given ratio
#chemoPolymFreq = [i / sum(chemoPolymHist) for i in chemoPolymHist]
#NCPolymFreq = [i / sum(NCPolymHist) for i in NCPolymHist]
#chemoDivFreq = [i / sum(chemoDivHist) for i in chemoDivHist]
#NCDivFreq = [i / sum(NCDivHist) for i in NCDivHist]
#print('computed proportions of genes in each bin')
#
## plot each histogram as a line
#
#
##data=np.array(np.random.rand(1000))
##y,binEdges=np.histogram(data,bins=100)
##bincenters = 0.5*(binEdges[1:]+binEdges[:-1])
##p.plot(bincenters,y,'-')
##p.show()
#
## plot polym ratio for chemo
#bincenters = 0.5*(a[1:] + a[:-1])
#graph3 = ax3.plot(bincenters, chemoPolymFreq, linewidth = 1.2, color = '#fcae91', alpha = 7)
#
## plot polym ratio for non chemo
#bincenters = 0.5*(b[1:] + b[:-1])
#graph4 = ax3.plot(bincenters, NCPolymFreq, linewidth = 1.2, color = '#bdd7e7', alpha = 7)
#
## plot div ratio for chemo
#bincenters = 0.5*(c[1:] + c[:-1])
#graph5 = ax3.plot(bincenters, chemoDivFreq, linewidth = 1.2, color = '#cb181d', alpha = 7)
#
## plot div ratio for non chemo
#bincenters = 0.5*(d[1:] + d[:-1])
#graph6 = ax3.plot(bincenters, NCDivFreq, linewidth = 1.2, color = '#2171b5', alpha = 7)
#
#print('added lines to 2nd subplot')
#
#
## do not show lines around figure, keep bottow line  
#ax3.spines["top"].set_visible(False)    
#ax3.spines["bottom"].set_visible(True)    
#ax3.spines["right"].set_visible(False)    
#ax3.spines["left"].set_visible(False)      
## offset the spines
#for spine in ax3.spines.values():
#  spine.set_position(('outward', 5))
#  
## add a light grey horizontal grid to the plot, semi-transparent, 
#ax3.yaxis.grid(True, linestyle='--', which='major', color='lightgrey', alpha=0.5, linewidth = 0.5)  
## hide these grids behind plot objects
#ax3.set_axisbelow(True)
#
## do not show ticks on 1st graph
#ax3.tick_params(
#    axis='x',       # changes apply to the x-axis and y-axis (other option : x, y)
#    which='both',      # both major and minor ticks are affected
#    bottom='on',      # ticks along the bottom edge are off
#    top='off',         # ticks along the top edge are off
#    right = 'off',
#    left = 'off',          
#    labelbottom='on', # labels along the bottom edge are off 
#    colors = 'black',
#    labelsize = 10,
#    direction = 'out') # ticks are outside the frame when bottom = 'on
#
## do not show ticks
#ax3.tick_params(
#    axis='y',       # changes apply to the x-axis and y-axis (other option : x, y)
#    which='both',      # both major and minor ticks are affected
#    bottom='off',      # ticks along the bottom edge are off
#    top='off',         # ticks along the top edge are off
#    right = 'off',
#    left = 'off',          
#    labelbottom='off', # labels along the bottom edge are off 
#    colors = 'black',
#    labelsize = 10,
#    direction = 'out') # ticks are outside the frame when bottom = 'on
#
#for label in ax3.get_yticklabels():
#    label.set_fontname('Arial')
#
#
## set y axis label
#ax3.set_ylabel('Proportion of genes', size = 10, ha = 'center', fontname = 'Arial')
#
#print('plotted 2nd subplot')

# add margin on the x-axis
plt.margins(0.05)

## add lines
#lns = graph3+graph4+graph5+graph6
## get labels
#labs = ['GPCR Pn/Ps', 'NC Pn/Ps', 'GPCR Dn/Ds', 'NC Dn/Ds']
## plot legend
#ax3.legend(lns, labs, loc=2, fontsize = 8, frameon = False)





# make sure subplots do not overlap
plt.tight_layout()

# save figure
fig.savefig('testfile.pdf', bbox_inches = 'tight')



