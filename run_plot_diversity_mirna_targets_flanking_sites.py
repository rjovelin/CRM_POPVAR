# -*- coding: utf-8 -*-
"""
Created on Fri May 13 13:55:12 2016

@author: RJovelin
"""


# use this script to generate a graph comparing diversity at mirna target sites with adjacent windows

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


for i in alldata:
    print(np.mean(i))
    
    
##################################
    
# create figure
fig = plt.figure(1, figsize = (6,4))
# add a plot to figure (1 row, 1 column, 1st plot)
ax = fig.add_subplot(1, 1, 1)

width = 0.3
ind = np.arange(len(alldata))
Means = [np.mean(i) for i in alldata]
SEM = []
for i in alldata:
    SEM.append(np.std(i) / math.sqrt(len(i)))

# use a boxplot
graph = ax.bar(ind, Means, yerr = SEM,
               color = ['#99d8c9', '#99d8c9', '#99d8c9', '#2ca25f','#99d8c9', '#99d8c9', '#99d8c9'],
               linewidth = 2)

ax.margins(0.05)

xvals = [i + 0.5 for i in range(len(alldata) + 1)]
# Set a buffer around the edge of the x-axis
plt.xlim([min(xvals)- 0.5, max(xvals)+ 0.5])


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

   
# Set the tick labels font name
for label in ax.get_yticklabels():
    label.set_fontname('Arial')

# write label for x and y axis
ax.set_ylabel('Nucleotide polymorphism\n', color = 'black',  size = 10, ha = 'center', fontname = 'Arial')
ax.set_xlabel('Upstream\tTargets\tDownstream', color = 'black', size = 10, ha = 'center', fontname = 'Arial')

# add labels to x-ticks, rotate and align right
#ax.set_xticklabels(site_types, ha = 'center', size = 10, fontname = 'Arial')

# remove lines around the frame
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['left'].set_visible(False)
ax.spines['bottom'].set_visible(False)

# add a light grey horizontal grid to the plot, semi-transparent, 
ax.yaxis.grid(True, linestyle='--', which='major', color='lightgrey', alpha=0.5)  
# hide these grids behind plot objects
ax.set_axisbelow(True)

# give more space to the y labels
fig.subplots_adjust(left=0.2)

        
# save figure
fig.savefig('testfile.pdf', bbox_inches = 'tight')
    
    
    
    
    
    
    
    
    
    
    
    
####################################



## create lists of lists to store the mean theta at each position 
## [sample size, mean, stderror]
#    
## lower index in upstream_pos corresponds to end the upstream region
## need to invert positions
## create a list of keys in upstream
#up_pos = [i for i in upstream_pos]
## sort keys
#up_pos.sort()
## reverse sort keys
#up_pos.reverse()
#    
## loop over reverse sorted keys in up_pos
#for i in up_pos:
#    # compute the mean theta at that position
#    mean_theta = np.mean(upstream_pos[i])
#    # compute standard error
#    stderror = np.std(upstream_pos[i]) / math.sqrt(len(upstream_pos[i]))
#    summary_file.write('\t'.join(['upstream_' + str(i), str(len(upstream_pos[i])), str(mean_theta), str(stderror)]) + '\n')
#    
#   
#    
## create a list of keys in downstream
#down_pos = [i for i in downstream_pos]
## sort list
#down_pos.sort()
#
## loop over sorted keys in down_pos
#for i in down_pos:
#    # compute the mean theta at that position
#    mean_theta = np.mean(downstream_pos[i])
#    # compyte the standard error
#    stderror = np.std(downstream_pos[i]) / math.sqrt(len(downstream_pos[i]))
#    summary_file.write('\t'.join(['downstream_' + str(i), str(len(downstream_pos[i])), str(mean_theta), str(stderror)]) + '\n')
#
## save theta value for targets and flaning sites to a single file
#newfile = open('theta_target_flanking_sites.txt', 'w')
#for i in up_pos:
#    for theta in upstream_pos[i]:
#        newfile.write('upstream_' + str(i) + '\t' + str(theta) + '\n')
#for theta in targets_theta:
#    newfile.write('targets' + '\t' + str(theta) + '\n')
#for i in down_pos:
#    for theta in downstream_pos[i]:
#        newfile.write('downstream_' + str(i) + '\t' + str(theta) + '\n')
#newfile.close()
#
#print('saved theta values to file')
#
## test mean differences between theta at target sites and flanking sites
#
#summary_file.write('\n')
#summary_file.write('Mean differences between all targets and flanking sites\n')
#summary_file.write('-' * 56 + '\n')
#summary_file.write('sites' + '\t' + 'wilcoxon' + 'P-val' + '\n')
#for i in up_pos:
#    wilcoxon, p = stats.ranksums(upstream_pos[i], targets_theta)
#    summary_file.write('\t'.join(['upstream_' + str(i) + '_vs_targets', str(wilcoxon), str(p)]) + '\n')
#for i in down_pos:
#    wilcoxon, p = stats.ranksums(downstream_pos[i], targets_theta)
#    summary_file.write('\t'.join(['downstream_' + str(i) + '_vs_targets', str(wilcoxon), str(p)]) + '\n')
#
#print('tested mean differences between targets and flanking sites')
#
#summary_file.close()
#
#
#
## create list of data
#alldata = [SYN_theta, REP_theta, mirna_theta, mature_theta, caeno_mir_theta, crmcla_mir_theta, crm_mir_theta]
## make a list of datatype
## create a list of corrsponding site type
#site_types = ['Syn', 'Rep', 'miRNA', 'miR', 'Caeno', 'Crm,Cla', 'Crm']
#
## create figure
#fig = plt.figure(1, figsize = (6,4))
## add a plot to figure (1 row, 1 column, 1st plot)
#ax = fig.add_subplot(1, 1, 1)
#
## use a boxplot
#bp = ax.boxplot(alldata, showmeans = False, showfliers = False, widths = 0.7, labels = site_types, patch_artist = True) 
# 
## create a list of colors (seee http://colorbrewer2.org/)
#color_scheme = ['#6e016b', '#88419d', '#8c6bb1', '#8c96c6','#9ebcda', '#bfd3e6', '#edf8fb']
## color boxes for the different sites
#i = 0    
## change box, whisker color to black
#for box in bp['boxes']:
#    # change line color
#    box.set(color = 'black', linewidth = 1.5)
#    box.set(facecolor = color_scheme[i])
#    # upate color
#    i += 1
## change whisker color ro black
#for wk in bp['whiskers']:
#    wk.set(color = 'black', linestyle = '-', linewidth = 1.5)
## change color of the caps
#for cap in bp['caps']:
#    cap.set(color = 'black', linewidth = 1.5)
## change the color and line width of the medians
#for median in bp['medians']:
#    median.set(color = 'black', linewidth = 1.5)
## change the mean marker and marker if showmean is True
##for mean in bp['means']:
##    mean.set(marker = 'o', markeredgecolor = 'black', markerfacecolor = 'black', markersize = 4)
#    
## restrict the x and y axis to the range of data
#ax.set_ylim([0, 0.15])
## create a list with range of x-axis values
#xvals = [i + 0.5 for i in range(len(site_types) + 1)]
## Set a buffer around the edge of the x-axis
#plt.xlim([min(xvals)- 0.5, max(xvals)+ 0.5])
#
#
## do not show ticks
#plt.tick_params(
#    axis='both',       # changes apply to the x-axis and y-axis (other option : x, y)
#    which='both',      # both major and minor ticks are affected
#    bottom='off',      # ticks along the bottom edge are off
#    top='off',         # ticks along the top edge are off
#    right = 'off',
#    left = 'off',          
#    labelbottom= 'on', # labels along the bottom edge are off 
#    colors = 'black',
#    labelsize = 10,
#    direction = 'out')
#
#   
## Set the tick labels font name
#for label in ax.get_yticklabels():
#    label.set_fontname('Arial')
#
## write label for x and y axis
#ax.set_ylabel('Nucleotide polymorphism\n', color = 'black',  size = 10, ha = 'center', fontname = 'Arial')
#ax.set_xlabel('Genomic features', color = 'black', size = 10, ha = 'center', fontname = 'Arial')
#
## add labels to x-ticks, rotate and align right
#ax.set_xticklabels(site_types, ha = 'center', size = 10, fontname = 'Arial')
#
## remove lines around the frame
#ax.spines['top'].set_visible(False)
#ax.spines['right'].set_visible(False)
#ax.spines['left'].set_visible(False)
#ax.spines['bottom'].set_visible(False)
#
## add a light grey horizontal grid to the plot, semi-transparent, 
#ax.yaxis.grid(True, linestyle='--', which='major', color='lightgrey', alpha=0.5)  
## hide these grids behind plot objects
#ax.set_axisbelow(True)
#
## give more space to the y labels
#fig.subplots_adjust(left=0.2)
#
    
