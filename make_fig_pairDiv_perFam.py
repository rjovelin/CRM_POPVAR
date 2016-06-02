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
print('generated a list of directories')

# multiply by 100 to get %
for i in range(len(FamDist)):
    FamDist[i] = FamDist[i] * 100

# create a dictionnary 
Distances = {}

# loop over directories
for folder in folders:
    # create a list of filenames of aligned sequences
    files = [filename for filename in os.listdir('./Pairwise_Chemos/' + folder) if '.txt' in filename]
    print(folder, len(files))    
    # create a list to store the distances
    distances = []
    # loop over files, convert to fasta
    for filename in files:
        proteins = convert_fasta('./Pairwise_Chemos/' + folder + './' + filename)
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
    Distance[family] = hist[0]


# print pairs counts    
print('Distances', len(Distances))
for family in Distances:
    print(family, sum(Distances[family]))
    

# create figure
fig = plt.figure(1, figsize = (6,2.5))
# add a plot to figure (1 row, 1 column, 1 plot)
ax = fig.add_subplot(1, 1, 1)  

# make a list of families
families = [i for i in Distances]


#############################


# set width of bar
width = 0.1


# plot SNP proportions SYN
graph1 = ax.bar([0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9], Distances[families[0]], width, color = 'black', edgecolor = 'black', linewidth = 1)
# plot SNP proportions REP
graph2 = ax.bar([1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2], Distances[families[1]], width, color = 'black', edgecolor = 'black', linewidth = 1)
# plot SNP proportions miRNAs
graph3 = ax.bar([2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3, 3.1], Distances[families[2]], width, color = 'black', edgecolor = 'black', linewidth = 1)









################################

# set y axis label
ax.set_ylabel('Number of protein pairs', size = 10, ha = 'center', fontname = 'Helvetica', family = 'sans-serif')

# set x axis label
ax.set_xlabel('Protein distance', size = 10, ha = 'center', fontname = 'Helvetica', family = 'sans-serif')

## remove top axes and right axes ticks
#ax.get_xaxis().tick_bottom()
#ax.get_yaxis().tick_left()

# do not show lines around figure, keep bottow line  
ax.spines["top"].set_visible(False)    
ax.spines["bottom"].set_visible(True)    
ax.spines["right"].set_visible(False)    
ax.spines["left"].set_visible(False)      
# offset the spines
for spine in ax.spines.values():
  spine.set_position(('outward', 5))  
  
# add a light grey horizontal grid to the plot, semi-transparent, 
ax.yaxis.grid(True, linestyle='--', which='major', color='lightgrey', alpha=0.5)  
# hide these grids behind plot objects
ax.set_axisbelow(True)
  
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
    labelsize = 10,
    direction = 'out') # ticks are outside the frame when bottom = 'on'  



# determine tick position on x axis
xpos =  [i for i in range(0, 110, 10)]
Dist = ['0', '0.1', '0.2', '0.3', '0.4', '0.5', '0.6', '0.7', '0.8', '0.9', '1']
# set up tick positions and labels
plt.xticks(xpos, Dist, rotation = 0, fontsize = 10, ha = 'center', fontname = 'Helvetica')

## add labels to x-ticks, rotate and align right
#ax.set_xticklabels(range(0, 110, 10), rotation = 0, ha = 'center', size = 10, fontname = 'Helvetica', family = 'sans-serif')
#
#plt.yticks(fontsize = 10)
#plt.xticks(range(0, 110, 10))

# add margins on right and left of gragh
plt.margins(0.05)


# save figure
fig.savefig('testfile.pdf', bbox_inches = 'tight')



    
    
    
    
