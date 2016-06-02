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

# loop over directories
for folder in folders:
    print(folder)
    # create a list of filenames of aligned sequences
    files = [filename for filename in os.listdir('./Pairwise_Chemos/' + folder) if '.txt' in filename]
    print(len(files))    
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
            # add distance to list
            FamDist.append(p_distance)
            
# multiply by 100 to get %
for i in range(len(FamDist)):
    FamDist[i] = FamDist[i] * 100




######################### EDIT THIS ASCRIP TO make hist for each famile
    
# loop over directories
for folder in folders:
    print(folder)
    # create a list of filenames of aligned sequences
    files = [filename for filename in os.listdir('./Pairwise_Chemos/' + folder) if '.txt' in filename]
    print(len(files))    
    
    # create a list to store the distances
    distances = []
    # loop over files, convert to fasta
    for filename in files:
        proteins = convert_fasta(filename)
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
    # get the folder name
    folder_name = folder[:folder.index('_')] + '_hist'
    # create new file in parent directory
    newfile = open('../' + folder_name + '.txt', 'w')
    for i in range(len(hist[0])):
        newfile.write(str(hist[1][i]) + ':' + str(hist[1][i] + 10) + '\t' + str(hist[0][i]) + '\n')
    # close file
    newfile.close()
    # go back to parent directory
    os.chdir('../')







    
# create figure
fig = plt.figure(1, figsize = (4.3,2.56))
# add a plot to figure (1 row, 1 column, 1 plot)
ax = fig.add_subplot(1, 1, 1)  

# create histogram
ax.hist(FamDist, range(0, 110, 10), color = '#e34a33', edgecolor = '#e34a33')

## add title
#ax.set_title('Stop codon mutations along coding sequences\n', size = 10, ha = 'center', fontname = 'Helvetica', family = 'sans-serif')

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



    
    
    
    
