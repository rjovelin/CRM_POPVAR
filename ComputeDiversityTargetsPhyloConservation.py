# -*- coding: utf-8 -*-
"""
Created on Tue May 17 13:43:10 2016

@author: RJovelin
"""


# use this script to generate a bar graph comparing diversity at mirna target sites
# with different levels of phylogenetic conservation

# use Agg backend on server without X server
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from manipulate_sequences import *
from sliding_windows import *
from sites_with_coverage import *
from miRNA_target import *
from divergence import *
from genomic_coordinates import *
import numpy as np
from scipy import stats
import math
import sys


# convert genome fasta to dict
genome = convert_fasta('../Genome_Files/noamb_356_v1_4.txt')
print('converted genome fasta to dict')
# get the allele counts for all sites with coverage, exclude sites with sample size < 10
chromo_sites = get_non_coding_snps('../SNP_files/', 10) 
print('got allele counts at all sites')
# get mirna target coordinates 
# {gene: [[chromo, start, end, orientation, seed, N_mirnas, site_type, conservation, utr]]}
target_coord = get_miRNA_target_coord('Cremanei_miRNA_sites.txt')
print('got miRNA target coordinates')
# get the set of valid transcripts
valid_transcripts = get_valid_transcripts('../Genome_Files/unique_transcripts.txt')
print('made list of valid transcripts')
# remove genes that are not not valid (ie to keep a single transcript per gene)
to_remove = [gene for gene in target_coord if gene not in valid_transcripts]
for gene in to_remove:
    del target_coord[gene]
print('deleted {0} non-valid genes'.format(len(to_remove)))
print('filtered non-valid transcripts') 


# compute diversity for all mirna targets, sort according to phylogenetic conservation
crmclacel_targets, crmcla_targets, crm_targets = [], [], []

# loop over gene in target cood
for gene in target_coord:
    # loop over targets in gene
    for i in range(len(target_coord[gene])):
        # get chromo
        chromo = target_coord[gene][i][0]
        # get start, end position
        start = target_coord[gene][i][1]
        end = target_coord[gene][i][2]
        # get conservation
        conservation = target_coord[gene][i][7]
        # check that chromo is in chromo_sites
        if chromo in chromo_sites:
            # compute theta, accepting 0 missing sites
            theta = compute_theta_non_coding(chromo_sites, chromo, start, end, 0)
            # check that theta is defined
            if theta != 'NA':
                # check conservation
                if conservation == 'crm':
                    # add theta to specific targets
                    crm_targets.append(theta)
                elif conservation == 'crm:cla':
                    # add theta to the crm-cla conservaed targets
                    crmcla_targets.append(theta)
                elif conservation == 'crm:cla:cel':
                    # add theta to the crm-la-cel conserved targets
                    crmclacel_targets.append(theta)
   
# compare diversity between target sites, make a list of theta lists
alldata = [crmclacel_targets, crmcla_targets, crm_targets]
# make a list of conservation level
site_types = ['Conserved', 'Restricted', 'Specific']


# create figure
fig = plt.figure(1, figsize = (2,2))
# add a plot to figure (1 row, 1 column, 1st plot)
ax = fig.add_subplot(1, 1, 1)

width = 0.8
ind = np.arange(len(alldata))
Means = [np.mean(i) for i in alldata]
SEM = []
for i in alldata:
    SEM.append(np.std(i) / math.sqrt(len(i)))
# use a bar plot
ax.bar(ind, Means, width, yerr = SEM,
       color = ['#2ca25f', '#99d8c9', '#e5f5f9'], linewidth = 1.5,
       error_kw=dict(elinewidth=1.5, ecolor='black', markeredgewidth = 1.5))               
    
# restrict the x and y axis to the range of data
ax.set_ylim([0, 0.030])

# remove lines around the frame
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['left'].set_visible(False)
ax.spines['bottom'].set_visible(True)
# offset the spines
for spine in ax.spines.values():
    spine.set_position(('outward', 5))

# set up tick positions and labels
# set x tick positions
xtickpos = [i + width/2 for i in range(len(alldata))]
ax.set_xticks(xtickpos)

# do not show ticks
plt.tick_params(
    axis='both',       # changes apply to the x-axis and y-axis (other option : x, y)
    which='both',      # both major and minor ticks are affected
    bottom='off',      # ticks along the bottom edge are off
    top='off',         # ticks along the top edge are off
    right = 'off',
    left = 'off',          
    labelbottom= 'on', # labels along the bottom edge are off 
    colors = 'black',
    labelsize = 10,
    direction = 'out')

plt.xticks(xtickpos, site_types, rotation = 20, ha = 'center', size = 10, fontname = 'Arial')

# Set the tick labels font name
for label in ax.get_yticklabels():
    label.set_fontname('Arial')

# write label for x and y axis
ax.set_ylabel('Diversity at target sites', color = 'black',  size = 10, ha = 'center', fontname = 'Arial')
ax.set_xlabel('Phylogenetic conservation', color = 'black', size = 10, ha = 'center', fontname = 'Arial')

# add margins
ax.margins(0.05)

# add a light grey horizontal grid to the plot, semi-transparent, 
ax.yaxis.grid(True, linestyle='--', which='major', color='lightgrey', alpha=0.5)  
# hide these grids behind plot objects
ax.set_axisbelow(True)

# compare mean differences among conservation levels
for i in range(len(alldata) -1):
    for j in range(i+1, len(alldata)):
        # get the P value of Wilcoxon rank sum test
        Pval = stats.ranksums(alldata[i], alldata[j])[1]
        # get stars for significance
        if Pval > 0.05:
            P = 'N.S.'
        elif Pval < 0.05 and Pval > 0.01:
            P = '*'
        elif Pval < 0.01 and Pval > 0.001:
            P = '**'
        elif Pval < 0.001:
            P = '***'
        print(site_types[i], site_types[j], Pval)    
        
# I already determined that all site categories are significantly different
# using Wilcoxon rank sum tests, so we need now to add letters to show significance

# annotate figure to add significance
# get the x and y coordinates
y_pos = [0.0075, 0.0125, 0.028]
x_pos = [i + width/2 for i in range(len(site_types))]
diff = ['A', 'B', 'C']

for i in range(len(diff)):
    ax.text(x_pos[i], y_pos[i], diff[i], horizontalalignment='center',
            verticalalignment='center', color = 'black', fontname = 'Arial', size = 10)

# save figure
fig.savefig('DiversityTargetsConservation.pdf', bbox_inches = 'tight')

