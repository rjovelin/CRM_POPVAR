# -*- coding: utf-8 -*-
"""
Created on Thu May 19 16:28:06 2016

@author: RJovelin
"""

# use this script to plot mature miR polymorphism as a function of mature expression level


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

# expression is the average expression level across different time points
# including larval and embryonic development

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
            # take mean of expression at different time points
            level = np.mean(list(map(lambda x: float(x), line[9:17])))
            # populate dict
            mature_expression[mirna] = level
        else:
            print('no expression for {0}'.format(mirna))
print('recorded miR expression')

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
# removes mirnas with no expression
to_remove = []
for mirna in diversity:
    if mirna not in mature_expression:
        to_remove.append(mirna)
print('remove {0} mirnas with diversity but no expression'.format(len(to_remove)))        
if len(to_remove) != 0:
    for mirna in to_remove:
        del diversity[mirna]
# check that dicts have same keys
assert mature_expression.keys() == diversity.keys(), 'expression and diversity should have same mirna keys'

# create list of expression values
expression_level = [mature_expression[mirna] for mirna in mature_expression]
print('made list of expression values')

# compute quartiles of expression values
Q1 = np.percentile(expression_level, 25)
Q2 = np.percentile(expression_level, 50)
Q3 = np.percentile(expression_level, 75)
print('computed quartiles of expression distribution')
print(Q1, Q2, Q3)

# partition thetas according to expression quartiles
# create lists to store theta for different levels of expression
lowExp, moderateExp, mediumExp, highExp = [], [], [], []
# loop over mirna in theta dict
for mirna in diversity:
    # check expression value
    if mature_expression[mirna] < Q1:
        lowExp.append(diversity[mirna])
    elif mature_expression[mirna] >= Q1 and mature_expression[mirna] < Q2:
        moderateExp.append(diversity[mirna])
    elif mature_expression[mirna] >= Q2 and mature_expression[mirna] < Q3:
        mediumExp.append(diversity[mirna])
    elif mature_expression[mirna] >= Q3:
        highExp.append(diversity[mirna])

# verify that mirnas belong to all levels of expression
assert len(lowExp) != 0
assert len(moderateExp) != 0
assert len(mediumExp) != 0
assert len(highExp) != 0
print('low', len(lowExp))
print('moderate', len(moderateExp))
print('medium', len(mediumExp))
print('high', len(highExp))

# create list of data
alldata = [lowExp, moderateExp, mediumExp, highExp]
# create a list of expression levels
ExpLevels = ['Low', 'Moderate', 'Medium', 'High']

# create figure
fig = plt.figure(1, figsize = (2,2))
# add a plot to figure (1 row, 1 column, 1st plot)
ax = fig.add_subplot(1, 1, 1)

width = 0.8
ind = np.arange(len(alldata))
Means = [np.mean(i) for i in alldata]
SEM = []
for i in alldata:
    SEM.append(np.std(i) / math.sqrt(len(i)))
# use a bar plot
ax.bar(ind, Means, width, yerr = SEM,
       color = ['#fef0d9','#fdcc8a','#fc8d59','#d7301f'],
       linewidth = 1.5, error_kw=dict(elinewidth=1.5, ecolor='black', markeredgewidth = 1.5))               
    
# restrict the x and y axis to the range of data
ax.set_ylim([0, 0.014])

# remove lines around the frame
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['left'].set_visible(False)
ax.spines['bottom'].set_visible(True)
# offset the spines
for spine in ax.spines.values():
    spine.set_position(('outward', 5))

# set x tick positions
xtickpos = [i + width/2 for i in range(len(alldata))]
ax.set_xticks(xtickpos)

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

# set up xtick labels
plt.xticks(xtickpos, ExpLevels, rotation = 20, ha = 'right', size = 10, fontname = 'Arial')

# Set the tick labels font name
for label in ax.get_yticklabels():
    label.set_fontname('Arial')

# write label for x and y axis
ax.set_ylabel('Nucleotide polymorphism', color = 'black',  size = 10, ha = 'center', fontname = 'Arial')
ax.set_xlabel('miR expression level', color = 'black', size = 10, ha = 'center', fontname = 'Arial')

# add margins
ax.margins(0.05)

# add a light grey horizontal grid to the plot, semi-transparent, 
ax.yaxis.grid(True, linestyle='--', which='major', color='lightgrey', alpha=0.5)  
# hide these grids behind plot objects
ax.set_axisbelow(True)

# compare mean differences among conservation levels
for i in range(len(alldata) -1):
    for j in range(i+1, len(alldata)):
        # get the P value of Wilcoxon rank sum test
        Pval = stats.ranksums(alldata[i], alldata[j])[1]
        # get stars for significance
        if Pval > 0.05:
            P = 'N.S.'
        elif Pval < 0.05 and Pval > 0.01:
            P = '*'
        elif Pval < 0.01 and Pval > 0.001:
            P = '**'
        elif Pval < 0.001:
            P = '***'
        print(ExpLevels[i], ExpLevels[j], Pval)    
        
# I already determined that all site categories are significantly different
# using Wilcoxon rank sum tests, so we need now to add letters to show significance

# annotate figure to add significance
# get the x and y coordinates
y_pos = [0.0142, 0.0123, 0.0045, 0.0025]
x_pos = [i + width/2 for i in range(len(site_types))]
diff = ['A', 'A,B', 'B,C', 'C']

for i in range(len(diff)):
    ax.text(x_pos[i], y_pos[i], diff[i], horizontalalignment='center',
            verticalalignment='center', color = 'black', fontname = 'Arial', size = 10)

# save figure
fig.savefig('DiversityMatureExpression.pdf', bbox_inches = 'tight')
    
    
 