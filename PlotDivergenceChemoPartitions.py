# -*- coding: utf-8 -*-
"""
Created on Mon May 30 20:23:16 2016

@author: Richard
"""

# use this script to plot a bar graph comparing divergence for TransMembrane and Extramemebrane domains


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



# set minimum number of sites (number of cosons = N sites / 3)
MinimumSites = 30

# set threshold for dN, dS and omega, remove genes instead of removing values
# genes are removed if any divergence estimate is greater than following thresholds
MaxdN = 2
MaxdS = 1.5
MaxOmega = 5

# use this function to get the number of sites without gaps in partitions used to compute divergence
def GetNumbSitesNoGaps(folder):
    '''
    (str) -> dict
    Return a dict with remanei gene name as key and the number of sites
    without gaps in commun between remanei and latens used to cpmute divergence
    in a chemo partition
    ''' 
    # create a dict of gene: N sites pairs
    L = {}
    # determine suffix of gene name
    if folder == './Partitions/Membrane/':
        suffix = '_TM.txt'
    elif folder == './Partitions/Extra_membrane/':
        suffix = '_ExtraTM.txt'
    elif folder == './Partitions/Inside/':
        suffix = '_inside.txt'
    elif folder == './Partitions/Outside/':
        suffix = '_outside.txt'
    
    # loop over file in folder
    for filename in os.listdir(folder):
        # grab the alignment file
        if suffix in filename:
            gene = filename[:filename.index(suffix)]
            infile = open(folder + filename)
            # skip seq count
            infile.readline()
            # skip crm gene name
            infile.readline()
            # grab crem seq without gaps
            crmseq = infile.readline().rstrip().replace('-', '')
            # skip cla gene name
            infile.readline()
            # grab cla seq without gaps
            claseq = infile.readline().rstrip().replace('-', '')
            infile.close()
            L[gene] = min([len(crmseq), len(claseq)])
    return L


# makde dictionaries with {gene: length partition}
LenTM = GetNumbSitesNoGaps('./Partitions/Membrane/')
LenExtra = GetNumbSitesNoGaps('./Partitions/Extra_membrane/')
print('generated dicts with Num sites')

# get the list of codeml output file in each parition
TM_out = [filename for filename in os.listdir('./Partitions/Membrane/') if '.out' in filename]
Extra_out = [filename for filename in os.listdir('./Partitions/Extra_membrane/') if '.out' in filename]
print('generated file lists')


# create dictionnaries to store the divergence of each partition {gene1 :[dN, dS, omega], gene2 : [dN, dS, omega]}
TM_diverg, Extra_diverg = {}, {}

# loop over codeml output files
for filename in TM_out:
    # parse codeml outputfile, return a dict {gene :[ dN, dS. omega]}
    gene_diverg = parse_codeml_output('./Partitions/Membrane/' + filename)
    # populate divergence dict
    for gene in gene_diverg:
        TM_diverg[gene] = list(gene_diverg[gene])
print('got divergence in Transmembrane')

for filename in Extra_out:
    # parse codeml outputfile, return a dict {gene :[ dN, dS. omega]}
    gene_diverg = parse_codeml_output('./Partitions/Extra_membrane/' + filename)
    # populate divergence dict
    for gene in gene_diverg:
        Extra_diverg[gene] = list(gene_diverg[gene])
print('got divergence in Extramembrane')


# create a dict to store divergence for Transmebrane and Extramembrane domains
# {gene: [TM_dN', 'TM_dS', 'TM_dN/dS', 'ExtraTM_dN', 'ExtraTM_dS', 'ExtraTM_dN/dS']}
ChemoDiv = {}

# do not consider partitions of genes that do not have TM and extra-TM domains
# loop over genes in TM_diverg
for gene in TM_diverg:
    if gene in Extra_diverg:
        # check if partitions have minimum number of codons
        if LenTM[gene] >= MinimumSites and LenExtra[gene] >= MinimumSites:
            ChemoDiv[gene] = [TM_diverg[gene][0], TM_diverg[gene][1], TM_diverg[gene][2]]
            # gene has extra-transmembrane domain, add divergence to list
            for i in range(len(Extra_diverg[gene])):
                ChemoDiv[gene].append(Extra_diverg[gene][i])

# make lists of dN values for TM and ExtraMb domains
dNTM, dNEX = [], []
#make lists of dS values for TM and ExtraMb domains
dSTM, dSEX = [], []
# make lists of omega values for TM and ExtraMb domains
omegaTM, omegaEX = [], [] 

# populate lists
for gene in ChemoDiv:
    # check if divergence is defined
    if ChemoDiv[gene][0] != 'NA' and ChemoDiv[gene][3] != 'NA':
        # check if dN is greater than max val
        if ChemoDiv[gene][0] <= MaxdN and ChemoDiv[gene][3] <= MaxdN:
            dNTM.append(ChemoDiv[gene][0])
            dNEX.append(ChemoDiv[gene][3])
    if ChemoDiv[gene][1] != 'NA' and ChemoDiv[gene][4] != 'NA':
        if ChemoDiv[gene][1] <= MaxdS and ChemoDiv[gene][4] <= MaxdS:
            dSTM.append(ChemoDiv[gene][1])
            dSEX.append(ChemoDiv[gene][4])
    if ChemoDiv[gene][2] != 'NA' and ChemoDiv[gene][5] != 'NA':
        if ChemoDiv[gene][2] <= MaxOmega and ChemoDiv[gene][5] <= MaxOmega:
            omegaTM.append(ChemoDiv[gene][2])
            omegaEX.append(ChemoDiv[gene][5])

a = [dNTM, dNEX, dSTM, dSEX, omegaTM, omegaEX]
b = ['dNTM', 'dNEX', 'dSTM', 'dSEX', 'omegaTM', 'omegaEX']

for i in range(len(a)):
    print('N', b[i], len(a[i]))
print('\n')
for i in range(len(a)):
    print('max', b[i], max(a[i]))

# compare dN for TM and EX using paired tests
P_rep = stats.wilcoxon(dNTM, dNEX)[1]
# compare dS for TM and EX using paired tests
P_syn = stats.wilcoxon(dSTM, dSEX)[1]
# compare omega for TM and EX using paired test
P_omega = stats.wilcoxon(omegaTM, omegaEX)[1]

print('P_rep', P_rep)
print('P_syn', P_syn)
print('P_omega', P_omega)

# create a list of means
Means = [np.mean(dNTM), np.mean(dNEX),
         np.mean(dSTM), np.mean(dSEX),
         np.mean(omegaTM), np.mean(omegaEX)]

# create a lit of SEM
SEM = [np.std(dNTM) / math.sqrt(len(dNTM)),
       np.std(dNEX) / math.sqrt(len(dNEX)),
       np.std(dSTM) / math.sqrt(len(dSTM)),
       np.std(dSEX) / math.sqrt(len(dSEX)),
       np.std(omegaTM) / math.sqrt(len(omegaTM)),
       np.std(omegaEX) / math.sqrt(len(omegaEX))]

# create figure
fig = plt.figure(1, figsize = (3, 2))
# add a plot to figure (1 row, 1 column, 1 plot)
ax = fig.add_subplot(1, 1, 1)  

# set width of bar
width = 0.2
# set colors
colorscheme = ['#31a354', '#e5f5e0','#31a354', '#e5f5e0', '#31a354', '#e5f5e0']


# plot nucleotide divergence
ax.bar([0, 0.2, 0.5, 0.7, 1, 1.2], Means, width, yerr = SEM, color = colorscheme, 
                edgecolor = 'black', linewidth = 1,
                error_kw=dict(elinewidth=1, ecolor='black', markeredgewidth = 1))

ax.set_ylabel('Nucleotide divergence', size = 10, ha = 'center', fontname = 'Arial')

# set y limits
plt.ylim([0, 0.45])

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
TransMb = mpatches.Patch(facecolor = '#31a354' , edgecolor = 'black', linewidth = 1, label= 'Transmembrane')
ExtraMb = mpatches.Patch(facecolor = '#e5f5e0', edgecolor = 'black', linewidth = 1, label = 'Extra-membrane')
plt.legend(handles=[TransMb, ExtraMb], loc = 2, fontsize = 8, frameon = False)

# I already determined that all site categories are significantly different
# using Wilcoxon rank sum tests, so we need now to add letters to show significance
# P_rep, P_syn and P_omega < 0.001 ---> P = ***
P = '***'

# annotate figure to add significance
# add bracket
ax.annotate("", xy=(0.1, 0.10), xycoords='data',
            xytext=(0.3, 0.10), textcoords='data',
            arrowprops=dict(arrowstyle="-", ec='#aaaaaa', connectionstyle="bar,fraction=0.2", linewidth = 1))
# add stars for significance
ax.text(0.2, 0.12, P, horizontalalignment='center',
        verticalalignment='center', color = 'grey', fontname = 'Arial', size = 6)

ax.annotate("", xy=(0.6, 0.30), xycoords='data',
            xytext=(0.8, 0.30), textcoords='data',
            arrowprops=dict(arrowstyle="-", ec='#aaaaaa', connectionstyle="bar,fraction=0.2", linewidth = 1))
# add stars for significance
ax.text(0.7, 0.32, P, horizontalalignment='center',
        verticalalignment='center', color = 'grey', fontname = 'Arial', size = 6)

ax.annotate("", xy=(1.1, 0.23), xycoords='data',
            xytext=(1.3, 0.23), textcoords='data',
            arrowprops=dict(arrowstyle="-", ec='#aaaaaa', connectionstyle="bar,fraction=0.2", linewidth = 1))
# add stars for significance
ax.text(1.2, 0.25, P, horizontalalignment='center',
        verticalalignment='center', color = 'grey', fontname = 'Arial', size = 6)

fig.savefig('DivergenceChemoPartitions.pdf', bbox_inches = 'tight')

