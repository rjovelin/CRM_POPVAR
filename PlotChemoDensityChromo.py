# -*- coding: utf-8 -*-
"""
Created on Fri Jun  3 10:41:19 2016

@author: RJovelin
"""


# use this script to plot chemoreceptor gene density for linkage groups 

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
from chemoreceptors import *
from genomic_coordinates import *


# get the set of chemoreceptors from the iprscan outputfile
chemo = get_chemoreceptors('../Genome_Files//PX356_protein_seq.tsv') 
print('parsed chemo genes')

# create a set of valid transcripts
transcripts = get_valid_transcripts('../Genome_Files/unique_transcripts.txt')
print('got list of valid transcripts')

# create a set of valid chemoreceptors
GPCRs = set(gene for gene in chemo if gene in transcripts)
print('GPCRS', len(GPCRs))
print('made list of valid chemo genes')

# get transcript coordinates {transcript : [chromo, start, end, sense]}
Coords = get_genes_coordinates('../Genome_Files/356_10172014.gff3')
print('got gene coordinates')

# count the number of chemo genes per chromo
ChemoCounts = {}
for gene in Coords:
    # check if gene is GPCR:
    if gene in GPCRs:
        # grab chromo
        chromo = Coords[gene][0]
        if chromo in ChemoCounts:
            ChemoCounts[chromo] += 1
        else:
            ChemoCounts[chromo] = 1
print('counted the number of chemo genes per chromo')

# get frequencies of chemo genes on each LG
Freq = {}
for chromo in ChemoCounts:
    Freq[chromo] = round((ChemoCounts[chromo] / len(GPCRs)) * 100, 4)

# make a list of [count, chromo]
Counts = [[val, key] for key, val in ChemoCounts.items()]
# sort according to count
Counts.sort()

# Only consider LG with > 5% of all chemo genes
# count LG with a low propotion of chemo genes
lowfreq = 0
# create a set to store the LG with high chemo proportions
HighFreq = set()
for chromo in Freq:
    if Freq[chromo] < 8:
        lowfreq += 1
    else:
        HighFreq.add(chromo)
print('low chemo proportion', lowfreq)
for chromo in HighFreq:
    print('high chemo proportion', chromo, Freq[chromo])

# get the length of each chromo with high chem proportion 
# convert genome fasta to dict
genome = convert_fasta('../Genome_Files/noamb_356_v1_4.txt')
print('genome converted to fasta dict')

for chromo in genome:
    if chromo in HighFreq:
        print(chromo, str(len(genome[chromo])) + ' bp', str(len(genome[chromo]) // 1000000) + ' Mb', str(len(genome[chromo]) // 100000) + ' per 100Kb window')

# look for clusters on LG with high frequency chemo

# get the start positions of chemo genes on high proportion LG
chemo_start = {}
for gene in Coords:
    if gene in GPCRs:
        chromo = Coords[gene][0]
        if chromo in HighFreq:
            # check orientation
            if Coords[gene][-1] == '+':
                start = Coords[gene][1]
            elif Coords[gene][-1] == '-':
                start = Coords[gene][2]
            if chromo in HighFreq:
                if chromo in chemo_start:
                    chemo_start[chromo].append(start)
                else:
                    chemo_start[chromo] = [start]
print('got chemo genes start positions')

# sort positions
for chromo in chemo_start:
    chemo_start[chromo].sort()
print('sorted the chemo start positions')
    
# create a function to count the number of chemo per window interval
def CountChemoWindow(chemo_start, genome, chromo, window):
    '''
    (dict, str, int) -> list
    Take the dictionary with chemo gene start positions per chromo, 
    the dict with genome sequence, a chromosome and a window interval in bp
    and return a list with the number of chemo genes on chromo per window interval'
    '''
    # create list with count of gene o repeat in 100000 bp windows
    range_counts = [0] * (len(genome[chromo]) // window)
    for start in chemo_start[chromo]:
        which_range = start // window
        if which_range == len(range_counts):
            which_range -= 1
        # count repeats
        range_counts[which_range] += 1
    return range_counts
    
    
# set up interval length in bp
Interval = 100000
print('Interval window: {0} bp'.format(Interval))    

# get the count of chemo gene per window
WindowCount = {}
for chromo in chemo_start:
    range_counts = CountChemoWindow(chemo_start, genome, chromo, Interval)     
    WindowCount[chromo] = range_counts
print('got chemo count per window')    

# create a list with the position of each window interval
Positions = {}
for chromo in WindowCount:
    pos = [i for i in range(len(WindowCount[chromo]))] 
    Positions[chromo] = pos     
    print('position', chromo, len(pos))
    print('interval length', len(genome[chromo]) // Interval)
print('got positions of window intervals')


# create figure
fig = plt.figure(1, figsize = (4, 1))
# add a plot to figure (1 row, 1 column, 1 plot)
ax = fig.add_subplot(1, 1, 1)  

# find the longtest chromo
chromoLength = [[len(genome[chromo]), chromo] for chromo in HighFreq]
chromoLength.sort()
size = [chromoLength[i][0] for i in range(len(chromoLength))]
LG = [chromoLength[i][1] for i in range(len(chromoLength))]
longest, maxlength = chromoLength[-1][-1], chromoLength[-1][0] 

# create a list of colors
#colorscheme = ['#a6cee3','#1f78b4','#b2df8a','#33a02c']
colorscheme = ['#1b9e77','#d95f02','#7570b3']

Graph = {}
# loop over chromo, from chromo with lowest to highest count
for i in range(len(LG)):
    print(LG[i])
    # plot the repeat of gene density per window
    graph = ax.plot(Positions[LG[i]], WindowCount[LG[i]], linewidth = 1.2, color = colorscheme[i], alpha = 0.7)
    Graph[LG[i]] = graph
    
ax.set_ylabel('GPCRs /100 Kb', size = 10, ha = 'center', fontname = 'Arial')
 
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

print('longest chromo', longest, maxlength)

# determine tick position on x axis
xpos =  [j for j in range(0, len(Positions[longest]), 10)]
# convert interval windows numbers to genomic positions
xtext = list(map(lambda x : (x * Interval) / 1000000, xpos))
Xtext = []
for i in xtext:
    if i % 2 == 0:
        Xtext.append(str(int(i)))
    else:
        Xtext.append('')

# set up tick positions and labels
plt.xticks(xpos, Xtext, rotation = 0, fontsize = 10, fontname = 'Arial')

# add lines
lns = Graph[LG[0]]
for chromo in LG[1:]:
    lns += Graph[chromo]
# get labels
labs = []
for chromo in LG:
    assert chromo.count('_') == 2
    lg = chromo[chromo.index('_', chromo.index('_')+1, -1)+1:]
    labs.append('LG' + lg)
# plot legend
ax.legend(lns, labs, loc=1, fontsize = 8, frameon = False)

fig.savefig('ChemoDensityLG.pdf', bbox_inches = 'tight')