# -*- coding: utf-8 -*-
"""
Created on Mon Jun 29 19:27:18 2015

@author: Richard
"""


# use this script to plot a bar graph comparing theta for TransMembrane and Extramemebrane domains
# for different site categories


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
from divergence import *

# set up minimum number of sites in partition
MinimumSites = 30

# get the set of chemoreceptors from the iprscan outputfile
chemo = get_chemoreceptors('../Genome_Files//PX356_protein_seq.tsv') 
print('parsed chemo genes')

# create a set of valid transcripts
transcripts = get_valid_transcripts('../Genome_Files/unique_transcripts.txt')
print('got list of valid transcripts')

# create a set of valid chemoreceptors
GPCRs = set(gene for gene in chemo if gene in transcripts)
print('made list of valid chemo genes')

# create a dict proba_domain {gene: {codon_index: [list of probabilities]}}
proba_domain = {}

# create list of phobius output files
files = [filename for filename in os.listdir('./Chemo_genes/') if '_proba.txt' in filename]

# loop over genes in GPCRs
for gene in GPCRs:
    # check has predictions
    for filename in files:
        if gene == filename[:filename.index('_proba')]:
            # initialize dict
            proba_domain[gene] = {}
            # get the dict of probabilities
            probabilities = parse_phobius_output(gene + '_proba.txt', './Chemo_genes/')
            # populate dict
            for key, val in probabilities.items():
                proba_domain[gene][key] = val
print('extracted proba for chemo domains')



# count SNPs, degenerate sites for partitions for nonsynonymous sites
TM_rep, TMrepsample, EX_rep, EXrepsample = count_SNPs_degenerate_sites_chemo_paritions('../Genome_Files/CDS_SNP_DIVERG.txt', proba_domain, 'REP',10, 0.95)
# compute theta at replacement sites for partitions
TM_theta_rep, EX_theta_rep = compute_theta_chemo_partitions(TM_rep, TMrepsample, EX_rep, EXrepsample, 'REP', transcripts, MinimumSites)
print('computed theta for replacement sites')

# count SNPs, degenerate sites for partitions for synonymous sites
TM_syn, TMsynsample, EX_syn, EXsynsample = count_SNPs_degenerate_sites_chemo_paritions('../Genome_Files/CDS_SNP_DIVERG.txt', proba_domain, 'SYN',10, 0.95)
# compute theta at synonymous sites for partitions
TM_theta_syn, EX_theta_syn = compute_theta_chemo_partitions(TM_syn, TMsynsample, EX_syn, EXsynsample, 'SYN', transcripts, MinimumSites)
print('computed theta for synonymous sites')

# create lists for theta in different partitions keeping the same gene between lists
theta_rep_TM = []
theta_rep_EX = []
for gene in TM_theta_rep:
    theta_rep_TM.append(TM_theta_rep[gene])
    theta_rep_EX.append(EX_theta_rep[gene])
    
theta_syn_TM = []
theta_syn_EX = []
for gene in TM_theta_syn:
    theta_syn_TM.append(TM_theta_syn[gene])
    theta_syn_EX.append(EX_theta_syn[gene])
    
# create dicts {gene: theta_rep / theta_syn} for each partition
TM_omega , EX_omega = {}, {}
# loop over genes in TM_theta_rep
for gene in TM_theta_rep:
    # check if gene TM_theta_syn
    if gene in TM_theta_syn:
        # check that theta syn different than 0
        if TM_theta_syn[gene] != 0:
            TM_omega[gene] = TM_theta_rep[gene] / TM_theta_syn[gene]
# loop over genes in EX_theta_rep
for gene in EX_theta_rep:
    # check if gene in EX_theta_syn
    if gene in EX_theta_syn:
        # check if theta syn is defined
        if EX_theta_syn[gene] != 0:
            EX_omega[gene] = EX_theta_rep[gene] / EX_theta_syn[gene]

# create lists of theta ratio for each partition, keeping the same gene order
theta_omega_TM , theta_omega_EX = [], []
for gene in TM_omega:
    # check that gene in EX_omega
    if gene in EX_omega:
        theta_omega_TM.append(TM_omega[gene])
        theta_omega_EX.append(EX_omega[gene])
print('computed theta for ratio')

a = [theta_rep_TM, theta_rep_EX, theta_syn_TM, theta_syn_EX, theta_omega_TM, theta_omega_EX]
b = ['theta_rep_tm', 'theta_rp_ex', 'theta_syn_tm', 'theta_syn_ex', 'theta_omega_tm', 'theta_omega_ex']
for i in range(len(a)):
    print('N', b[i], len(a[i]))
print('\n')
for i in range(len(a)):
    print('max', b[i], max(a[i]))

# compare theta at replacement sites for TM and EX using paired tests
P_rep = stats.wilcoxon(theta_rep_TM, theta_rep_EX)[1]
# compare theta at synonymous sites using paired tests
P_syn = stats.wilcoxon(theta_syn_TM, theta_syn_EX)[1]
# compare ratio theta rep / theta syn using paired test
P_omega = stats.wilcoxon(theta_omega_TM, theta_omega_EX)[1]

print('P_rep', P_rep)
print('P_syn', P_syn)
print('P_omega', P_omega)

# create a list of means
Means = [np.mean(theta_rep_TM), np.mean(theta_rep_EX),
         np.mean(theta_syn_TM), np.mean(theta_syn_EX),
         np.mean(theta_omega_TM), np.mean(theta_omega_EX)]

# create a lit of SEM
SEM = [np.std(theta_rep_TM) / math.sqrt(len(theta_rep_TM)),
       np.std(theta_rep_EX) / math.sqrt(len(theta_rep_EX)),
       np.std(theta_syn_TM) / math.sqrt(len(theta_syn_TM)),
       np.std(theta_syn_EX) / math.sqrt(len(theta_syn_EX)),
       np.std(theta_omega_TM) / math.sqrt(len(theta_omega_TM)),
       np.std(theta_omega_EX) / math.sqrt(len(theta_omega_EX))]

# create figure
fig = plt.figure(1, figsize = (3, 2))
# add a plot to figure (1 row, 1 column, 1 plot)
ax = fig.add_subplot(1, 1, 1)  

# set width of bar
width = 0.2
# set colors
colorscheme = ['#f03b20', '#ffeda0','#f03b20', '#ffeda0', '#f03b20', '#ffeda0']

# plot nucleotide divergence
ax.bar([0, 0.2, 0.5, 0.7, 1, 1.2], Means, width, yerr = SEM, color = colorscheme, 
                edgecolor = 'black', linewidth = 1,
                error_kw=dict(elinewidth=1, ecolor='black', markeredgewidth = 1))

ax.set_ylabel('Nucleotide diversity', size = 10, ha = 'center', fontname = 'Arial')

# set y limits
plt.ylim([0, 0.20])

# determine tick position on x axis
xpos =  [0.2, 0.7, 1.2]
# use greek letters and subscripts
xtext = ['$' + chr(952) +'_{rep}$', '$' + chr(952) + '_{syn}$', '$' + chr(952) + '_{syn}$' + '/' + '$' + chr(952) + '_{syn}$']

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

# create legend
TransMb = mpatches.Patch(facecolor = '#f03b20' , edgecolor = 'black', linewidth = 1, label= 'Transmembrane')
ExtraMb = mpatches.Patch(facecolor = '#ffeda0', edgecolor = 'black', linewidth = 1, label = 'Extra-membrane')
plt.legend(handles=[TransMb, ExtraMb], loc = 2, fontsize = 8, frameon = False)

# I already determined that all site categories are significantly different
# using Wilcoxon rank sum tests, so we need now to add letters to show significance
# P_rep, P_syn and P_omega < 0.001 ---> P = ***
P = '***'

# annotate figure to add significance
# add bracket
ax.annotate("", xy=(0.1, 0.02), xycoords='data',
            xytext=(0.3, 0.02), textcoords='data',
            arrowprops=dict(arrowstyle="-", ec='#aaaaaa', connectionstyle="bar,fraction=0.2", linewidth = 1))
# add stars for significance
ax.text(0.2, 0.03, P, horizontalalignment='center',
        verticalalignment='center', color = 'grey', fontname = 'Arial', size = 6)

ax.annotate("", xy=(1.1, 0.19), xycoords='data',
            xytext=(1.3, 0.19), textcoords='data',
            arrowprops=dict(arrowstyle="-", ec='#aaaaaa', connectionstyle="bar,fraction=0.2", linewidth = 1))
# add stars for significance
ax.text(1.2, 0.20, P, horizontalalignment='center',
        verticalalignment='center', color = 'grey', fontname = 'Arial', size = 6)

fig.savefig('DiversityChemoPartitions.pdf', bbox_inches = 'tight')

