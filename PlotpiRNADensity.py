# -*- coding: utf-8 -*-
"""
Created on Thu May 26 15:24:34 2016

@author: RJovelin
"""

# use this script to plot piRNA density for linkage groups LG4

# use Agg backend on server without X server
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import rc
rc('mathtext', default='regular')
# import modules
import numpy as np
from scipy import stats
import math
import sys
# import custom modules
from manipulate_sequences import *
from piRNAs import *


# look for clusters on linkage group 4
chromo = 'linkage_group_4'

# get the start positions of piRNAs
pirna_start = chromo_pirna_start('PX356_piRNA_coord.txt')
# sort positions
for i in pirna_start:
    pirna_start[i].sort()

print('sorted the pirna start positions')
    
# convert genome fasta to dict
genome = convert_fasta('../Genome_Files/noamb_356_v1_4.txt')
print('genome converted to fasta dict')
print(chromo, len(genome[chromo]))
print((len(genome[chromo]) // 100000))

# create list with count of gene o repeat in 100000 bp windows
range_counts = [0] * (len(genome[chromo]) // 100000)
for start in pirna_start[chromo]:
    which_range = start // 100000
    if which_range == len(range_counts):
        which_range -= 1
    # count repeats
    range_counts[which_range] += 1
    

# create a list with the position of each window interval
positions = [i for i in range(len(range_counts))]      
print('position', len(positions))
print('interval length', len(genome[chromo]) // 100000)

# create figure
fig = plt.figure(1, figsize = (4, 1))
# add a plot to figure (1 row, 1 column, 1 plot)
ax = fig.add_subplot(1, 1, 1)  

# plot the repeat of gene density per window
ax.plot(positions, range_counts, linewidth = 1.2, color = 'black')

ax.set_ylabel('# piRNAs /100 Kb', size = 10, ha = 'center', fontname = 'Arial')
 
# set x axis label
ax.set_xlabel('Position along linkage group (Mb)', size = 10, ha = 'center', fontname = 'Arial')

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

# determine tick position on x axis
xpos =  [j for j in range(0, len(positions), 10)]
# convert interval windows numbers to genomic positions
xtext = list(map(lambda x : (x * 100000) / 1000000, xpos))
xtext = list(map(lambda x : str(int(x)), xtext))
# set up tick positions and labels
plt.xticks(xpos, xtext, fontsize = 10, fontname = 'Arial')

fig.savefig('piRNADensityLG4.pdf', bbox_inches = 'tight')