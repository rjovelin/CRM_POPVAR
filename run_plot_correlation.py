# -*- coding: utf-8 -*-
"""
Created on Thu May 19 16:28:06 2016

@author: RJovelin
"""

# use this script to plot the correlation between nucleotide diversity in mature sequences
# and expression level of mature sequences





# use Agg backend on server without X server
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
# load modules
import numpy as np
import math
from scipy import stats
import sys
# load custom modules
from manipulate_sequences import *
from divergence import *
from miRNA_target import *
from sites_with_coverage import *



# get option to use expression in males or across development
expression = sys.argv[1]


# create a dict to store the expression level of mature miRNAs
mature_expression = {}
# get mature expression from file
infile = open('CRM_miRNAsCoordinatesFinal.txt')
infile.readline()
for line in infile:
    line = line.rstrip()
    if line != '':
        line = line.split('\t')
        mirna = line[0]
        # check that mirna has expression
        if len(line) >= 17:
            # get expression
            if expression == 'male':
                # take expression in remanei males
                level = float(line[17])
            elif expression == 'development':
                # take mean of expression at different time points
                level = np.mean(list(map(lambda x: float(x), line[9:17])))
            # populate dict
            mature_expression[mirna] = level
        else:
            print('no expression for {0}'.format(mirna))
print('recorded miR expression')


############################


# create a set of valid transcripts (1 transcript mapped to 1 gene)
valid_transcripts = get_valid_transcripts('../Genome_Files/unique_transcripts.txt')

# get the allele counts for all sites with coverage
chromo_sites = get_non_coding_snps('../SNP_files/', 10)
print('got allele counts at all sites')

# create a dict with coordinates of mature sequences
miR_coord = {}
infile = open('CRM_MatureCoordinatesFinal.txt')
infile.readline()
for line in infile:
    line = line.rstrip()
    if line != '':
        line = line.split('\t')
        name, chromo, start, end, orientation = line[0], line[1], int(line[2]) -1, int(line[3]), line[4]
        miR_coord[name] = [chromo, start, end, orientation]
infile.close()

# compute theta for each mirna
# create dict to store the mirna: diversity pair
diversity = {}
# loop over mirna name in {name : [chromo, start, end, orientation, conservation]}
for name in miR_coord:
    chromo, start, end = miR_coord[name][0], miR_coord[name][1], miR_coord[name][2]
    conservation = famCons[name]
    theta = compute_theta_non_coding(chromo_sites, chromo, start, end, 2)
    # check if theta is defined
    if theta != 'NA':
        diversity[name] = theta        

# remove mirnas from expression for which couldn't be computed
to_remove = []
for mirna in mature_expression:
    if mirna not in diversity:
        to_remove.append(mirna)
print('remove {0} mirnas with expression but no diversity'.format(len(to_remove)))
if len(to_remove) != 0:
    for mirna in to_remove:
        del mature_expression[mirna]
# check that dicts have same keys
assert mature_expression.keys() == diversity.keys(), 'expression and diversity should have same mirna keys'

# create parallel lists of expression and diversity
mature_theta = []
expression_level = []
for mirna in mature_expression:
    expression_level.append(mature_expression[mirna])
    mature_theta.append(diversity[mirna])
print('made lists of expression and diversity values')


# make a scatter plot diversity X expression

# create figure
fig = plt.figure(1, figsize = (4.3,2.56))
# add a plot to figure (1 row, 1 column, 1 plot)
ax = fig.add_subplot(1, 1, 1)    

# write label for y axis
ytext = ax.set_ylabel('Mature miR polymorphism', color = 'black', size = 10, ha = 'center', fontname = 'Arial', family = 'sans-serif')
if expression == 'male':
    xtext = ax.set_xlabel('Mature miR expression level in males', color = 'black', size = 10, ha = 'center', fontname = 'Arial', family = 'sans-serif')
elif expression == 'development':
    xtext = ax.set_xlabel('Mature miR expression level across development', color = 'black', size = 10, ha = 'center', fontname = 'Arial', family = 'sans-serif')

# add a light grey horizontal grid to the plot, semi-transparent, 
ax.yaxis.grid(True, linestyle='-', which='major', color='lightgrey', alpha=0.5)  
# hide these grids behind plot objects
ax.set_axisbelow(True)

# remove top axes and right axes ticks
ax.get_xaxis().tick_bottom()
ax.get_yaxis().tick_left()

# add data to plot
ax.scatter(expression_level, mature_theta, edgecolor = 'red', facecolor = 'r', alpha = 0.8)

# calc the trendline
z = np.polyfit(expression_level, mature_theta, 1)
p = np.poly1d(z)
ax.plot(expression_level,p(expression_level), linewidth = 1, color = 'black', linestyle = '-', alpha = 0.5)


# compute Spearman's rank correlation
rho, pval = stats.spearmanr(expression_level, mature_theta)
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
y_min, y_max = min(mature_theta), max(mature_theta) 
x_min, x_max = min(expression_level), max(expression_level)

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
#plt.title('Heteroplasmies as a function of the {0} read depth\n'.format(read_depth), size = 12, fontname = 'Arial')  

# save figure
fig.savefig('testfile.pdf', bbox_inches = 'tight')
    
    
 