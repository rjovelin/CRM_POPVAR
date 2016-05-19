# -*- coding: utf-8 -*-
"""
Created on Thu May 19 16:28:06 2016

@author: RJovelin
"""

# use this script to plot the correlation between nucleotide diversity in mature sequences
# and expression level of mature sequences









Name	PX356_scaffold	start	end	orientation	Hairpin	Mature	Seed	hairpin_arm	Emb _average	A _average	B _average	C _average	D _average	E _average	F _average	G _average	remanei_male _average	Emb _stdev	A _stdev	B _stdev	C _stdev	D _stdev	E _stdev	F _stdev	G _stdev	remanei_male _stdev	family_conservation











Skip to content
This repository
Search
Pull requests
Issues
Gist
 @rjovelin
 Unwatch 1
  Star 0
  Fork 0 rjovelin/MitoVariants
 Code  Issues 0  Pull requests 0  Wiki  Pulse  Graphs  Settings
Branch: master Find file Copy path
MitoVariants/PlotCorrelationReadDepthHeteroplasmies.py
70061e6  on Mar 24
@rjovelin rjovelin Initial commit
1 contributor
RawBlameHistory     215 lines (162 sloc)  6 KB
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 14 11:47:55 2016
@author: RJovelin
"""

# use this script to generate a scatter plot between the number of heteroplasmies
# and the average or median read depth per individual

# usage PlotCorrelationReadDepthHeteroplasmies.py [options]
# - ReadDepthFile: file with average and median read depth per cancer and per dataset
# - HeteroSummaryFile: file with mitoseek output
# - tumor
# - [WGS/RNAseq]: use WGS or RNAseq data set
# - [median/mean]: use median or mean read depth for each individual
# - output file: save figure file

import os
import sys
from scipy import stats
import numpy as np
import matplotlib.pyplot as plt

# get file with median and read depth per individual and cancer and data type
ReadDepthFile = sys.argv[1]
# get the summary file with mitoseek output
HeteroSummaryFile = sys.argv[2]
# get tumor
tumor = sys.argv[3]
# get datatype
datatype = sys.argv[4]
# use median or mean read depth
read_depth = sys.argv[5]
# save figure to outputfile
outputfile = sys.argv[6]

# parse file to extract read depth values
# open file for reading
infile = open(ReadDepthFile, 'r')
# skip header
infile.readline()
# create a dict {participant: read_depth}
ReadDepth = {}
# loop over file
for line in infile:
    line = line.rstrip()
    if line != '':
        line = line.split('\t')
        # get participant
        participant = line[0]
        # check if correct tumor and data set 
        if tumor == line[3] and datatype == line[-1]:
            # check data type
            if read_depth == 'mean':
                data = float(line[1])
            elif read_depth == 'median':
                data = float(line[2])
            # populate dict
            ReadDepth[participant] = data
infile.close()

# create a dict {participant: set(variable positions)}
Mutations = {}
# open heteroSummary file for reading
infile = open(HeteroSummaryFile, 'r')
infile.readline()
for line in infile:
    line = line.rstrip()
    if line != '':
        line = line.split('\t')
        participant = line[1]
        position = int(line[0]) - 1
        if participant in Mutations:
            Mutations[participant].add(position)
        else:
            Mutations[participant] = set()
            Mutations[participant].add(position)
infile.close()

# create a dict {participant: [# mutations, read_depth]}
data = {}
for participant in Mutations:
    if participant in ReadDepth:
        data[participant] = [len(Mutations[participant]), ReadDepth[participant]]

# create 2 parallel lists
mutations, reads = [], []
for participant in data:
    if data[participant][0] < 1000:
        mutations.append(data[participant][0])
        reads.append(data[participant][1])

# make a scatter plot mutations X reads

# create figure
fig = plt.figure(1, figsize = (4.3,2.56))

# add a plot to figure (1 row, 1 column, 1 plot)
ax = fig.add_subplot(1, 1, 1)    

# write label for y axis
ytext = ax.set_ylabel('Number of variable positions', color = 'black', size = 10, ha = 'center', fontname = 'Arial', family = 'sans-serif')
xtext = ax.set_xlabel('{0} read depth for {1}'.format(read_depth.capitalize(), tumor), color = 'black', size = 10, ha = 'center', fontname = 'Arial', family = 'sans-serif')

# add a light grey horizontal grid to the plot, semi-transparent, 
ax.yaxis.grid(True, linestyle='-', which='major', color='lightgrey', alpha=0.5)  
# hide these grids behind plot objects
ax.set_axisbelow(True)

# remove top axes and right axes ticks
ax.get_xaxis().tick_bottom()
ax.get_yaxis().tick_left()

ax.scatter(reads, mutations, edgecolor = 'red', facecolor = 'r', alpha = 0.8)

# calc the trendline
z = np.polyfit(reads, mutations, 1)
p = np.poly1d(z)
ax.plot(reads,p(reads), linewidth = 1, color = 'black', linestyle = '-', alpha = 0.5)


# compute Spearman's rank correlation
rho, pval = stats.spearmanr(mutations, reads)
if pval > 0.05:
    P = 'NS'
elif pval < 0.05 and pval > 0.01:
    P = '< 0.05'
elif pval < 0.01 and pval > 0.001:
    P = '< 0.01'
elif pval < 0.001:
    P = '< 0.001'
if P == 'NS':
    correlation = 'Spearman\'s r = {0}, P > 0.05'.format(round(rho, 3))
else:
    correlation = 'Spearman\'s r = {0}, P '.format(round(rho, 3)) + P


# annotate text
y_min, y_max = min(mutations), max(mutations) 
x_min, x_max = min(reads), max(reads)

ax.text(x_max, y_max, correlation, horizontalalignment='right',
        verticalalignment='center', fontsize = 8, fontname = 'Arial')
        
# do not show lines around figure  
ax.spines["top"].set_visible(False)    
ax.spines["bottom"].set_visible(False)    
ax.spines["right"].set_visible(False)    
ax.spines["left"].set_visible(False)  


## create a list with range of x-axis values
#xvals = [i + 0.5 for i in range(len(names) + 1)]
## Set a buffer around the edge of the x-axis
#plt.xlim([min(xvals)- 0.5, max(xvals)+ 0.5])

# do not show ticks
plt.tick_params(
    axis='y',       # changes apply to the x-axis and y-axis (other option : x, y)
    which='both',      # both major and minor ticks are affected
    bottom='off',      # ticks along the bottom edge are off
    top='off',         # ticks along the top edge are off
    right = 'off',
    left = 'off',          
    labelbottom='off', # labels along the bottom edge are off 
    colors = 'black',
    labelsize = 10)
      

# do not show ticks
plt.tick_params(
    axis='x',       # changes apply to the x-axis and y-axis (other option : x, y)
    which='both',      # both major and minor ticks are affected
    bottom='on',      # ticks along the bottom edge are off
    top='off',         # ticks along the top edge are off
    right = 'off',
    left = 'off',          
    labelbottom='on', # labels along the bottom edge are on 
    colors = 'lightgrey',
    labelcolor = 'black',
    labelsize = 10,
    direction = 'out')


# Set the tick labels font name
for label in ax.get_yticklabels():
    label.set_fontname('Arial')
    
# add title
plt.title('Heteroplasmies as a function of the {0} read depth\n'.format(read_depth), size = 12, fontname = 'Arial')  


















# save figure
fig.savefig(outputfile, bbox_inches = 'tight')
    
    
 