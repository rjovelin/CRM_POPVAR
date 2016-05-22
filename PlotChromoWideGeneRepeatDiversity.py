# -*- coding: utf-8 -*-
"""
Created on Thu Aug 20 13:17:35 2015

@author: RJovelin
"""

# use this script to plot nucleotide diversity with repeat or gene density for 
# linkage groups LG1, LG2 or LG4

# usage python3 PlotChromoWideGeneRepeatDiversity.py [options]
#- [repeats/genes] : plot gene or repeat density along with nucleotide diversity
#- [LG1/LG2/LG4] : chromo to plot
#- [PB, noPB] : if diversity is computed using Ontario strains only (noPB) ot with all strains (PB)

# use Agg backend on server without X server
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
# import modules
import numpy as np
from scipy import stats
import math
import sys
# import custom modules
from manipulate_sequences import *
#from piRNAs import *
from miRNA_target import *
from repeats_TEs import *
from sliding_windows import *
from sites_with_coverage import *
from divergence import *


# get density to plot with nucleotide diversity
density = sys.argv[1]
# get chromo 
LG = sys.argv[2]
# get corresponding chromos
if LG == 'LG1':
    chromosome = 'linkage_group_1'
elif LG == 'LG2':
    chromosome = 'linkage_group_2'
elif LG == 'LG4':
    chromosome = 'linkage_group_4'
# choose strains to compute diversity
strains = sys.argv[3]

# convert genome fasta to dict
genome = convert_fasta('../Genome_Files/noamb_356_v1_4.txt')
print('genome converted to fasta dict')

# get the gene coordinates {gene1 : [chromo, start, end, sense]}
genes_coord = get_genes_coordinates('../Genome_Files/356_10172014.gff3')
print('got gene coordinates')

# get the set of valid transcripts
valid_transcripts = get_valid_transcripts('../Genome_Files/unique_transcripts.txt')
print('got the set of valid transcripts')

# get a dict with {chromo: [start position (0-based) of the valid transcripts]}
genes_start = {}
for gene in genes_coord:
    if gene in valid_transcripts:
        chromo = genes_coord[gene][0]
        start = genes_coord[gene][1] -1
        # check if chromo in genes_start
        if chromo in genes_start:
            genes_start[chromo].append(start)
        else:
            genes_start[chromo] = [start]
# sort positions
for chromo in genes_start:
    genes_start[chromo].sort()
print('got the genes\'s first positions')

# get the coordinates of the repeats
# {repeat : [[chromo, start, end, orientation], [chromo, start, end, orientation]]}
repeats_coord = get_repeats_coord('../Genome_Files/356_v1_4.fasta.out', False)
print('got the repeat coordinates')

# create a dict with {chromo: [list of repeat start positions 0-based]}
repeats_start = {}
# loop over repeat family
for fam in repeats_coord:
    # loop over instances of repeat
    for i in range(len(repeats_coord[fam])):
        # get chromo
        chromo = repeats_coord[fam][i][0]
        # get start position
        start = repeats_coord[fam][i][1]
        # check if chromo in dict
        if chromo in repeats_start:
            # add start position to list
            repeats_start[chromo].append(start)
        else:
            # initiate list
            repeats_start[chromo] = [start]
# sort positions
for chromo in repeats_start:
    repeats_start[chromo].sort()
print('got the repeats\'s first positions')


# create list with count of gene o repeat in 50000 bp windows
range_counts = [0] * (len(genome[chromosome]) // 50000)
# check if gene or repeat density is recorded
if density == 'repeats':
    for start in repeats_start[chromosome]:
        which_range = start // 50000
        if which_range == len(range_counts):
            which_range -= 1
        # count repeats
        range_counts[which_range] += 1
elif density == 'genes':
    for start in genes_start[chromosome]:
        which_range = start // 50000
        if which_range == len(range_counts):
            which_range -= 1
        # counts genes
        range_counts[which_range] += 1

# create a list with the position of each window interval
positions = [i for i in range(len(range_counts))]      

print('position', len(positions))
print(LG, len(genome[chromosome]))
print('interval length', len(genome[chromosome]) // 50000)

# create figure
fig = plt.figure(1, figsize = (4, 1))
# add a plot to figure (1 row, 1 column, 1 plot)
ax1 = fig.add_subplot(1, 1, 1)  

# plot the repeat of gene density per window
ax1.plot(positions, range_counts, linewidth = 1, color = '#b2df8a')

# restrict the x and y axis to the range of data
#ax.set_xlim([0, len(Pos)])
if density == 'genes':
    ax1.set_ylim([0,25])
elif density == 'repeats':
    ax1.set_ylim([0, 200])
            

# set y axis label
if density == 'genes':
    YaxisText = 'Gene density'    
elif density == 'repeats':
    YaxisText = 'Repeat density'
ax1.set_ylabel(YaxisText, size = 10, ha = 'center', fontname = 'Arial')
 
# add labels to x-ticks
#ax.set_xticklabels([list of values], rotation = 30, ha = 'right', size = 10, fontname = 'Arial', family = 'sans-serif')

#plt.yticks(fontsize = 10, fontname = 'Arial')

# determine tick position on x axis
xpos =  [j for j in range(0, len(positions) + 50, 50)]
# convert interval windows numbers to genomic positions
xtext = list(map(lambda x : (x * 50000) / 1000000, xpos))
xtext = list(map(lambda x : str(x), xtext))
# set up tick positions and labels
plt.xticks(xpos, xtext, fontsize = 10, fontname = 'Arial')

# set x axis label
ax1.set_xlabel('Position along linkage group (Mb)', size = 10, ha = 'center', fontname = 'Arial')

# remove top axes and right axes ticks
ax1.get_xaxis().tick_bottom()
ax1.get_yaxis().tick_left()

 
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

# do not show ticks
plt.tick_params(
    axis='both',       # changes apply to the x-axis and y-axis (other option : x, y)
    which='both',      # both major and minor ticks are affected
    bottom='on',      # ticks along the bottom edge are off
    top='off',         # ticks along the top edge are off
    right = 'off',
    left = 'off',          
    labelbottom='on', # labels along the bottom edge are off 
    colors = 'black',
    labelsize = 10,
    direction = 'out') # ticks are outside the frame when bottom = 'on
      

# compute theta per 50 Kb window
# get the allele counts for all sites with coverage, keep all sites 
if strains == 'noPB':
    chromo_sites = get_non_coding_snps('../SNP_files/', 0) 
elif strains == 'PB':
    # compute theta per 50 Kb window in all strains
    chromo_sites = get_all_strains_snps('../PB_ON_SNP_files/', 0)
print('got allele counts')


# create a list of starting positions of each window on chromo
LG_positions = [i for i in range(0, len(genome[chromosome]), 50000)]
print('chromo positions', len(LG_positions))

# create a dict to store theta for each widow interval
diversity = {}
# loop over start position in LG:
for i in LG_positions:
    # compute theta per window
    # use a large threshold > window length to accept any number of missing site  
    theta = compute_theta_non_coding(chromo_sites, chromosome, i, i + 50000, 100000)
    diversity[i] = theta
print('computed diversity')

# create a list of positions for which diversity is calculated, sort positions
divpos = [i for i in diversity]
divpos.sort()
# create parallel list with theta values
polymorphism = [diversity[i] for i in divpos]
print('positions diversity windows', len(divpos))

# add another graph on top of previous one
ax2 = ax1.twinx()

# plot theta per window
ax1.plot(divpos, polymorphism, linewidth = 1, color = '#1f78b4')

## restrict the x and y axis to the range of data
##ax.set_xlim([0, len(Pos)])
#if density == 'genes':
#    ax1.set_ylim([0,25])
#elif density == 'repeats':
#    ax1.set_ylim([0, 200])
            

# set y axis label
ax2.set_ylabel('Nucleotide polymorphism', size = 10, ha = 'center', fontname = 'Arial')
 
# add labels to x-ticks
#ax.set_xticklabels([list of values], rotation = 30, ha = 'right', size = 10, fontname = 'Arial', family = 'sans-serif')

#plt.yticks(fontsize = 10, fontname = 'Arial')

## determine tick position on x axis
#xpos =  [j for j in range(0, len(positions) + 50, 50)]
## convert interval windows numbers to genomic positions
#xtext = list(map(lambda x : (x * 50000) / 1000000, xpos))
#xtext = list(map(lambda x : str(x), xtext))
## set up tick positions and labels
#plt.xticks(xpos, xtext, fontsize = 10, fontname = 'Arial')
#
## set x axis label
#ax1.set_xlabel('Position along linkage group (Mb)', size = 10, ha = 'center', fontname = 'Arial')
#
## remove top axes and right axes ticks
#ax1.get_xaxis().tick_bottom()
#ax1.get_yaxis().tick_left()
#
# 
## do not show lines around figure, keep bottow line  
#ax1.spines["top"].set_visible(False)    
#ax1.spines["bottom"].set_visible(True)    
#ax1.spines["right"].set_visible(False)    
#ax1.spines["left"].set_visible(False)      
## offset the spines
#for spine in ax1.spines.values():
#  spine.set_position(('outward', 5))
#  
## add a light grey horizontal grid to the plot, semi-transparent, 
#ax1.yaxis.grid(True, linestyle='--', which='major', color='lightgrey', alpha=0.5, linewidth = 0.5)  
## hide these grids behind plot objects
#ax1.set_axisbelow(True)
#
## do not show ticks
#plt.tick_params(
#    axis='both',       # changes apply to the x-axis and y-axis (other option : x, y)
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
#


# save figure
fig.savefig('testfile.pdf', bbox_inches = 'tight')
