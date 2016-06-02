# -*- coding: utf-8 -*-
"""
Created on Thu Jun  2 11:01:47 2016

@author: RJovelin
"""

# use this script to generate a figure with histograms of protein distance
# between pairs of chemoreceptor within each family separately


# use Agg backend on server without X server
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib import rc
rc('mathtext', default='regular')
# import modules
import os
import numpy as np
# import custom modules
from chemoreceptors import *
from manipulate_sequences import *
from tcoffee_alignment import *
from divergence import *


# make a list of folders
folders = [folder for folder in os.listdir('./Pairwise_Chemos/') if '_family' in folder]
# sort folders
folders.sort()

# create a dictionnary 
Distances = {}

# loop over directories
for folder in folders:
    # create a list of filenames of aligned sequences
    files = [filename for filename in os.listdir('./Pairwise_Chemos/' + folder) if '.txt' in filename]
    # create a list to store the distances
    distances = []
    # loop over files, convert to fasta
    for filename in files:
        proteins = convert_fasta('./Pairwise_Chemos/' + folder + '/' + filename)
        # get the sequences
        genes = [i for i in proteins]
        seq1 = proteins[genes[0]]
        seq2 = proteins[genes[1]]
        # get p-distances
        p_distance = pairwise_distance(seq1, seq2, 'protein') 
        if p_distance != 'NA':
            # multiply by 100 to get %
            p_distance *= 100
            distances.append(p_distance)
    # create a histogram
    hist = np.histogram(distances, range(0, 101, 10))
    # grab family name
    family = folder[:folder.index('_family')]
    # populate dict with pairs counts
    Distances[family] = list(hist[0])

# create figure
fig = plt.figure(1, figsize = (6, 6))
## adjust white space between subplots
#fig.subplots_adjust(left=0.2, wspace=0.6)

# make a list of families
families = [i for i in Distances]
families.sort()


# create a function to format the subplots
def CreateAx(Columns, Rows, Position, Xscale, Data, figure, Title, width, XLabel = False, YLabel = False, BottomLine = False):
    # create subplot in figure
    # add a plot to figure (N row, N column, plot N)
    ax = figure.add_subplot(Rows, Columns, Position)
    ax.bar(Xscale, Data, width, color = 'black', edgecolor = 'black', linewidth = 1)
    ax.set_title(Title, size = 10)
    
    if YLabel == True:
        # set y axis label
        ax.set_ylabel('Protein pairs', size = 7, ha = 'center', fontname = 'Helvetica', family = 'sans-serif')
    if XLabel == True:
        # set x axis label
        ax.set_xlabel('Protein distance', size = 7, ha = 'center', fontname = 'Helvetica', family = 'sans-serif')

    # do not show lines around figure, keep bottow line  
    ax.spines["top"].set_visible(False)    
    ax.spines["right"].set_visible(False)    
    ax.spines["left"].set_visible(False)      
    
    if BottomLine == True:
        ax.spines["bottom"].set_visible(True)
    else:
        ax.spines["bottom"].set_visible(False)
    
        # add a light grey horizontal grid to the plot, semi-transparent, 
    ax.yaxis.grid(True, linestyle=':', which='major', color='lightgrey', alpha=0.5, linewidth = 0.5)  
    # hide these grids behind plot objects
    ax.set_axisbelow(True)    
        
    if BottomLine == True:
        # do not show ticks
        plt.tick_params(
            axis='both',       # changes apply to the x-axis and y-axis (other option : x, y)
            which='both',      # both major and minor ticks are affected
            bottom='on',      # ticks along the bottom edge are off
            top='off',         # ticks along the top edge are off
            right = 'off',
            left = 'off',          
            labelbottom='on', # labels along the bottom edge are on
            colors = 'black',
            labelsize = 7,
            direction = 'out') # ticks are outside the frame when bottom = 'on'  
    else:
        # do not show ticks
        plt.tick_params(
            axis='both',       # changes apply to the x-axis and y-axis (other option : x, y)
            which='both',      # both major and minor ticks are affected
            bottom='off',      # ticks along the bottom edge are off
            top='off',         # ticks along the top edge are off
            right = 'off',
            left = 'off',          
            labelbottom='off', # labels along the bottom edge are on
            colors = 'black',
            labelsize = 7,
            direction = 'out') # ticks are outside the frame when bottom = 'on'  

    if XLabel == True:
        # determine tick position on x axis
        xpos =  [i /10  for i in range(11)]
        Dist = ['0', '', '', '', '', '0.5', '', '', '', '', '1']
        # set up tick positions and labels
        plt.xticks(xpos, Dist, rotation = 0, fontsize = 7, ha = 'center', fontname = 'Helvetica')

    if YLabel == True:
        plt.yticks(fontsize = 8)
  
    return ax

# set width for all subplots    
width = 0.1

ax1 = CreateAx(4, 5, 1, [i / 10 for i in range(10)], Distances[families[0]], fig, families[0], width, YLabel = True)
ax2 = CreateAx(4, 5, 2, [i / 10 for i in range(10)], Distances[families[1]], fig, families[1], width)
ax3 = CreateAx(4, 5, 3, [i / 10 for i in range(10)], Distances[families[2]], fig, families[2], width)
ax4 = CreateAx(4, 5, 4, [i / 10 for i in range(10)], Distances[families[3]], fig, families[3], width)
ax5 = CreateAx(4, 5, 5, [i / 10 for i in range(10)], Distances[families[4]], fig, families[4], width, YLabel = True)
ax6 = CreateAx(4, 5, 6, [i / 10 for i in range(10)], Distances[families[5]], fig, families[5], width)
ax7 = CreateAx(4, 5, 7, [i / 10 for i in range(10)], Distances[families[6]], fig, families[6], width)
ax8 = CreateAx(4, 5, 8, [i / 10 for i in range(10)], Distances[families[7]], fig, families[7], width)
ax9 = CreateAx(4, 5, 9, [i / 10 for i in range(10)], Distances[families[8]], fig, families[8], width, YLabel = True)
ax10 = CreateAx(4, 5, 10, [i / 10 for i in range(10)], Distances[families[9]], fig, families[9], width)
ax11 = CreateAx(4, 5, 11, [i / 10 for i in range(10)], Distances[families[10]], fig, families[10], width)
ax12 = CreateAx(4, 5, 12, [i / 10 for i in range(10)], Distances[families[11]], fig, families[11], width)
ax13 = CreateAx(4, 5, 13, [i / 10 for i in range(10)], Distances[families[12]], fig, families[12], width, YLabel = True)
ax14 = CreateAx(4, 5, 14, [i / 10 for i in range(10)], Distances[families[13]], fig, families[13], width)
ax15 = CreateAx(4, 5, 15, [i / 10 for i in range(10)], Distances[families[14]], fig, families[14], width)
ax16 = CreateAx(4, 5, 16, [i / 10 for i in range(10)], Distances[families[15]], fig, families[15], width)
ax17 = CreateAx(4, 5, 17, [i / 10 for i in range(10)], Distances[families[16]], fig, families[16], width, XLabel = True, YLabel = True, BottomLine = True)
ax18 = CreateAx(4, 5, 18, [i / 10 for i in range(10)], Distances[families[17]], fig, families[17], width, XLabel = True, BottomLine = True)
ax19 = CreateAx(4, 5, 19, [i / 10 for i in range(10)], Distances[families[18]], fig, families[18], width, XLabel = True, BottomLine = True)

# make sure subplots do not overlap
plt.tight_layout()

# save figure
fig.savefig('ChemoProtDivergenceFamilies.pdf', bbox_inches = 'tight')

