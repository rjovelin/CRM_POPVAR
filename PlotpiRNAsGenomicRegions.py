# -*- coding: utf-8 -*-
"""
Created on Mon Aug 10 09:59:13 2015

@author: RJovelin
"""


# use this script to plot the number of piRNAs located in different regions


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
import sys
import json
# import custom modules
from manipulate_sequences import *
from miRNA_target import *
from genomic_coordinates import *
from piRNAs import *
from get_coding_sequences import *


# convert genome fasta to dict
genome = convert_fasta('../Genome_Files/noamb_356_v1_4.txt')
# get piRNA coordinates
pirna_coord = get_pirna_loci('PX356_piRNA_coord.txt') 
print('got piRNA coord')

# get CDS_coord {TS1: [chromo, sense, [(s1, end1), (s2, end2)]]}
CDS_coord = get_CDS_positions('../Genome_Files/356_10172014.gff3')
# create new dict in the form {chromo: {set of positions}}
CDS_pos = {}
for gene in CDS_coord:
    # get chromo
    chromo = CDS_coord[gene][0]
    # loop over cds coordinates 
    for i in range(len(CDS_coord[gene][2])):
        # get start and end positions in a list
        # convert to 0-based index
        start  = CDS_coord[gene][2][i][0] - 1
        end = CDS_coord[gene][2][i][1]
        # check if chromo in CDS_pos
        if chromo in CDS_pos:
            for j in range(start, end):
                CDS_pos[chromo].add(j)
        else:
            CDS_pos[chromo] = set()
            for j in range(start, end):
                CDS_pos[chromo].add(j)
print('got CDS coord')

# get intron coordinates
# create dict
intron_pos = {}
# loop over genes in CDS_pos
for gene in CDS_coord:
    # get chromo
    chromo = CDS_coord[gene][0]
    # loop over cds coordinates for that gene
    for i in range(0, len(CDS_coord[gene][2]) -1 ):
        # convert intron coordinates to 0-based index
        start = CDS_coord[gene][2][i][1]
        end = CDS_coord[gene][2][i+1][0]
        # check if chromo in intron pos
        if chromo in intron_pos:
            for j in range(start, end):
                intron_pos[chromo].add(j)
        else:
            intron_pos[chromo] = set()
            for j in range(start, end):
                intron_pos[chromo].add(j)
print('got intron coord')


# get UTR coord {TS1 : [chromo, start, end, orientation]}
# load UTR coordinates from json file
infile = open('../miRNA_Target_sites/CremUTRCoordsNo.json')
UTR_coord = json.load(infile)
infile.close()
print('UTR coords', len(UTR_coord))
print('got UTR coordinates from json file')

# create a dict UTR_pos
UTR_pos = {}
# loop over genes in three_prime
for gene in UTR_coord:
    # get chromo
    chromo = UTR_coord[gene][0]
    # convert to 0-based
    start = UTR_coord[gene][1] -1
    end = UTR_coord[gene][2]
    # check if chromo in UTR_pos
    if chromo in UTR_pos:
        for j in range(start, end):
            UTR_pos[chromo].add(j)
    else:
        UTR_pos[chromo] = set()
        for j in range(start, end):
            UTR_pos[chromo].add(j)
print('got UTR positons')

# get the intergenic positons

# make a dict with positions in intergenic regions
# {chromo: {set of positions}}
intergenic_pos = keep_intergenic_positions(genome, CDS_pos, UTR_pos, intron_pos)
print('got intergenic positions')    

# find pirna locations
locations = find_pirna_locations(genome, pirna_coord, CDS_pos, UTR_pos, intron_pos, intergenic_pos)
print('got the piRNA locations')

# get the counts for each region
counts = [[val, key] for key, val in locations.items()]
# sort list according to count
counts.sort()
# get the pirnas numbers
nums, regions = [], []
for i in counts:
    nums.append(i[0])
    regions.append(i[1])
print(regions)
print(nums)
print(['intergenic', 'intron', 'UTR', 'overlapping', 'CDS'])



# create figure
fig = plt.figure(1, figsize = (2.5, 1.5))
# add a plot to figure (1 row, 1 column, 1 plot)
ax = fig.add_subplot(1, 1, 1)  

# set width of bar
width = 0.4
# set colors
colorscheme = ['#b30000', '#e34a33', '#fc8d59', '#fdcc8a', '#fef0d9']

# plot number of piRNAs in each region
ax.bar([0, 0.8, 1.2, 1.6, 2], nums, width, color = colorscheme, edgecolor = 'black', linewidth = 1)

# write label
ax.set_ylabel('Number of piRNAs', size = 10, ha = 'center', fontname = 'Arial')
# set x axis label
ax.set_xlabel('Genomic regions', size = 10, ha = 'center', fontname = 'Arial')

# determine tick position on x axis
xpos =  [0, 0.8, 1.2, 1.6, 2]
# set up tick positions and labels
plt.xticks(xpos, regions, rotation = 20, fontsize = 10, fontname = 'Arial')

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
intergenic = mpatches.Patch(facecolor = '#b30000', edgecolor = 'black', linewidth = 1, label= 'intergenic')
intron = mpatches.Patch(facecolor = '#e34a33', edgecolor = 'black', linewidth = 1, label = 'intron')
UTR = mpatches.Patch(facecolor = '#fc8d59', edgecolor = 'black', linewidth = 1, label = 'UTR')
overlapping = mpatches.Patch(facecolor = '#fdcc8a', edgecolor = 'black', linewidth = 1, label = 'overlapping')
CDS = mpatches.Patch(facecolor = '#fef0d9', edgecolor = 'black', linewidth = 1, label = 'CDS')

plt.legend(handles=[intergenic, intron, UTR, overlapping, CDS], loc = 1, fontsize = 8, frameon = False)

# add margin on the x-axis
plt.margins(0.05)

fig.savefig('testfile.pdf', bbox_inches = 'tight')


