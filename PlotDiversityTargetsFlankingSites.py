# -*- coding: utf-8 -*-
"""
Created on Fri May 13 13:55:12 2016

@author: RJovelin
"""


# use this script to generate a graph comparing diversity at mirna target sites with adjacent windows

# load matplotlib, change backend when X server is not available
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
# load custome modules
from manipulate_sequences import *
from sliding_windows import *
from sites_with_coverage import *
from miRNA_target import *
from divergence import *
from genomic_coordinates import *
# load scipy and numpy
import numpy as np
from scipy import stats
import math

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


# compute diversity for all mirna targets
targets_theta = [] 
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
        # get utr type
        utr = target_coord[gene][i][8]
        # check that chromo is in chromo_sites
        if chromo in chromo_sites:
            # compute theta, accepting 0 missing sites
            theta = compute_theta_non_coding(chromo_sites, chromo, start, end, 0)
            # check that theta is defined
            if theta != 'NA':
                # add to list 
                targets_theta.append(theta)
print('computed theta for miRNA targets')

# create a dict for upstream and downstream flanking sites
# {window_position : [list of thetas across all regions at that position]}
upstream_pos, downstream_pos = {}, {}

# compute theta for each non-overlapping windows upstream of target sites
# loop over targets
# compute theta for each upstream region
for gene in target_coord:
    # loop over target on gene
    for i in range(len(target_coord[gene])):
        # get target coordinates
        chromo = target_coord[gene][i][0]
        start = target_coord[gene][i][1]
        end = target_coord[gene][i][2]
        orientation = target_coord[gene][i][3]
        # check that chromo in chromo-sites
        if chromo in chromo_sites:
            # determine the length of the target sites
            target_size = end - start
            upstream_size = target_size * 3
            # compute upstream coordinates
            # take 3 windows on each side of the target
            if orientation == '+':
                up_end = start
                if start - upstream_size < 0:
                    up_start = 0
                else:
                    up_start = start - upstream_size
            elif orientation == '-':
                up_start = end
                if up_start + upstream_size < len(genome[chromo]):
                    up_end = up_start + upstream_size
                else:
                    up_end = up_start + upstream_size
            # compute theta for each non-verlapping window
            theta_windows = sequence_sliding_window(chromo, up_start, up_end, orientation, True, target_size, target_size, chromo_sites, 0)
            # add theta at each position in the upstream dict
            for j in theta_windows:
                # check if j in dict
                if j in upstream_pos:
                    # add theta if theta is defined
                    if len(theta_windows[j]) != 0:
                        upstream_pos[j].append(theta_windows[j][0])
                else:
                    # initiate list value and add theta if theta is defined
                    upstream_pos[j] = []
                    if len(theta_windows[j]) != 0:
                        upstream_pos[j].append(theta_windows[j][0])
print('computed sliding windows for upstream')                           
            
# loop over targets
for gene in target_coord:
    # loop over target on gene
    for i in range(len(target_coord[gene])):
        # get target coordinates
        chromo = target_coord[gene][i][0]
        start = target_coord[gene][i][1]
        end = target_coord[gene][i][2]
        orientation = target_coord[gene][i][3]
        # check that chromo in chromo-sites
        if chromo in chromo_sites:
            # determine the length of the target sites
            target_size = end - start
            downstream_size = target_size * 3
            # compute downstream coordinates
            # take 3 windows on each side of the target
            if orientation == '+':
                down_start = end
                if down_start + downstream_size < len(genome[chromo]):
                    down_end = down_start + downstream_size
                else:
                    down_end = len(genome[chromo])
            elif orientation == '-':
                down_end = start
                if start - downstream_size < 0:
                    down_start = 0
                else:
                    down_start = start - downstream_size
            theta_windows = sequence_sliding_window(chromo, down_start, down_end, orientation, False, target_size, target_size, chromo_sites, 0)
            # add theta at each position in the downstream dict
            for j in theta_windows:
                # check if j is key
                if j in downstream_pos:
                    # add theta if theta is defined
                    if len(theta_windows[j]) != 0:
                        downstream_pos[j].append(theta_windows[j][0])
                else:
                    # initiate list value and add theta if theta is defined
                    downstream_pos[j] = []
                    if len(theta_windows[j]) != 0:
                        downstream_pos[j].append(theta_windows[j][0])
print('computed sliding windows for downstream')


# create a list of lists of diversity values ordered according to position around target sites
alldata = []

# lower index in upstream_pos corresponds to end the upstream region
# need to invert positions
# create a list of keys in upstream, reverse sorted
up_pos = [i for i in upstream_pos]
up_pos.sort()
up_pos.reverse()
down_pos = [i for i in downstream_pos]
down_pos.sort()
for i in up_pos:
    alldata.append(upstream_pos[i])
alldata.append(targets_theta)
for i in down_pos:
    alldata.append(downstream_pos[i])
  
   
# create figure
fig = plt.figure(1, figsize = (4, 3))
# add a plot to figure (1 row, 1 column, 1st plot)
ax = fig.add_subplot(1, 1, 1)

width = 0.8
ind = np.arange(len(alldata))
Means = [np.mean(i) for i in alldata]
SEM = []
for i in alldata:
    SEM.append(np.std(i) / math.sqrt(len(i)))

# use a bar plot
graph = ax.bar(ind, Means, width, yerr = SEM,
               color = ['#99d8c9', '#99d8c9', '#99d8c9', '#2ca25f','#99d8c9', '#99d8c9', '#99d8c9'],
               linewidth = 1.5, error_kw=dict(elinewidth=1.5, ecolor='black', markeredgewidth = 1.5))               

ax.margins(0.05)

xvals = [i + 0.5 for i in range(len(alldata))]

# Set a buffer around the edge of the x-axis
#plt.xlim([-0.5, len(alldata) + 0.5])
# zom in by setting up y limits
plt.ylim([0.020, 0.025])


# do not show ticks
plt.tick_params(
    axis='both',       # changes apply to the x-axis and y-axis (other option : x, y)
    which='both',      # both major and minor ticks are affected
    bottom='off',      # ticks along the bottom edge are off
    top='off',         # ticks along the top edge are off
    right = 'off',
    left = 'off',          
    labelbottom= 'off', # labels along the bottom edge are off 
    colors = 'black',
    labelsize = 10,
    direction = 'out')

   
# Set the tick labels font name
for label in ax.get_yticklabels():
    label.set_fontname('Arial')

# write label for x and y axis
ax.set_ylabel('Nucleotide polymorphism', color = 'black',  size = 10, ha = 'center', fontname = 'Arial')
ax.set_xlabel('\nUpstream\t\tTargets\t\tDownstream', color = 'black', size = 10, ha = 'center', fontname = 'Arial')

# add labels to x-ticks, rotate and align right
#ax.set_xticklabels(site_types, ha = 'center', size = 10, fontname = 'Arial')

# remove lines around the frame
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['left'].set_visible(False)
ax.spines['bottom'].set_visible(True)

# offset the spines
for spine in ax.spines.values():
  spine.set_position(('outward', 5))

# add a light grey horizontal grid to the plot, semi-transparent, 
ax.yaxis.grid(True, linestyle='--', which='major', color='lightgrey', alpha=0.5)  
# hide these grids behind plot objects
ax.set_axisbelow(True)

# annotate figure to add significance
# get the x and y coordinates
y_min, y_max = 0.020, 0.025

# compare upstream and downstream bins to targets
for i in range(len(list(ind))):
    # do not compare targets with targets
    if i != 3:
        # get the P value of Wilcoxon rank sum test
        Pval = stats.ranksums(alldata[i], alldata[3])[1]
        # get stars for significance
        if Pval > 0.05:
            P = 'N.S.'
        elif Pval < 0.05 and Pval > 0.01:
            P = '*'
        elif Pval < 0.01 and Pval > 0.001:
            P = '**'
        elif Pval < 0.001:
            P = '***'
    
        # add stars for significance
        if P == 'N.S.':
            ax.text(i + width/2, y_max + abs(y_max - y_min)*0.03, P, horizontalalignment='center',
                    verticalalignment='center', color = 'grey', fontname = 'Arial', size = 6)
        else:
            ax.text(i + width/2, y_max + abs(y_max - y_min)*0.01, P, horizontalalignment='center',
                    verticalalignment='center', color = 'grey', fontname = 'Arial')

# save figure
fig.savefig('DiversityTargetsFlankingSites.pdf', bbox_inches = 'tight')
    
    
