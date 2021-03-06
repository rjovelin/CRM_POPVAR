# -*- coding: utf-8 -*-
"""
Created on Wed May 18 16:20:05 2016

@author: RJovelin
"""

# use this script to make a box plot comparing theta for mature miRNAs with different levels of conservation
# usage PlotDiversityConservationMature.py [options]
# [box/bar] : type of graphich to plot


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

# get graph type from command
graphtype = sys.argv[1]

# get the allele counts for all sites with coverage
chromo_sites = get_non_coding_snps('../SNP_files/', 10)
print('got allele counts at all sites')

# make a dict with family level conservation for all miRNAs {name : conservarion}
famCons = {}
infile = open('CRM_miRNAsCoordinatesFinal.txt')
infile.readline()
for line in infile:
    line = line.rstrip()
    if line != '':
        line = line.split('\t')
        name, conservation = line[0], line[-1]
        famCons[name] = conservation
infile.close()
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

# partition thetas for miRs according to level of conservation
# create lists to store theta at miRNA for different level of family conservation
caeno_mir_theta, crmcla_mir_theta, crm_mir_theta = [], [], []
# loop over mirna name in {name : [chromo, start, end, orientation, conservation]}
for name in miR_coord:
    chromo, start, end = miR_coord[name][0], miR_coord[name][1], miR_coord[name][2]
    conservation = famCons[name]
    theta = compute_theta_non_coding(chromo_sites, chromo, start, end, 2)
    # check if theta is defined
    if theta != 'NA':
        # check conservation level
        if conservation == 'Crm':
            crm_mir_theta.append(theta)
        elif conservation == 'CrmCla':
            crmcla_mir_theta.append(theta)
        elif conservation == 'Caeno':
            caeno_mir_theta.append(theta)
# verify that mirnas belong to all 3 levels of conservation
assert len(caeno_mir_theta) != 0
assert len(crmcla_mir_theta) != 0
assert len(crm_mir_theta) != 0
print(len(caeno_mir_theta))
print(len(crmcla_mir_theta))
print(len(crm_mir_theta))


# create list of data
alldata = [caeno_mir_theta, crmcla_mir_theta, crm_mir_theta]
# make a list of datatype
# create a list of corrsponding site type
site_types = ['Conserved', 'Restricted', 'Specific']

# create figure
fig = plt.figure(1, figsize = (2,2))
# add a plot to figure (1 row, 1 column, 1st plot)
ax = fig.add_subplot(1, 1, 1)

# check graphich type in option
if graphtype == 'box':
    width = 0.8
    # use a boxplot
    bp = ax.boxplot(alldata, showmeans = False, showfliers = False, widths = width, patch_artist = True) 
    
    # create a list of colors (seee http://colorbrewer2.org/)
    color_scheme = ['#8856a7', '#9ebcda', '#e0ecf4']

    # color boxes for the different sites
    i = 0    
    # change box, whisker color to black
    for box in bp['boxes']:
        # change line color
        box.set(color = 'black', linewidth = 1.5)
        box.set(facecolor = color_scheme[i])
        # upate color
        i += 1
    # change whisker color ro black
    for wk in bp['whiskers']:
        wk.set(color = 'black', linestyle = '-', linewidth = 1.5)
    # change color of the caps
    for cap in bp['caps']:
        cap.set(color = 'black', linewidth = 1.5)
    # change the color and line width of the medians
    for median in bp['medians']:
        median.set(color = 'black', linewidth = 1.5)
    
    # restrict the x and y axis to the range of data
    ax.set_ylim([0, 0.07])


elif graphtype == 'bar':
    width = 0.8
    ind = np.arange(len(alldata))
    Means = [np.mean(i) for i in alldata]
    SEM = []
    for i in alldata:
        SEM.append(np.std(i) / math.sqrt(len(i)))
    # use a bar plot
    ax.bar(ind, Means, width, yerr = SEM,
           color = ['#8856a7', '#9ebcda', '#e0ecf4'],
           
           linewidth = 1.5, error_kw=dict(elinewidth=1.5, ecolor='black', markeredgewidth = 1.5))               
    
    # restrict the x and y axis to the range of data
    ax.set_ylim([0, 0.025])

# remove lines around the frame
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['left'].set_visible(False)
if graphtype == 'box':
    ax.spines['bottom'].set_visible(False)
elif graphtype == 'bar':
    ax.spines['bottom'].set_visible(True)
    # offset the spines
    for spine in ax.spines.values():
        spine.set_position(('outward', 5))

# set up tick positions and labels
if graphtype == 'box':
    # set x tick positions
    xtickpos = [i + 0.5 for i in range(len(alldata))]
elif graphtype == 'bar':
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

if graphtype == 'box':
    plt.xticks(xtickpos, site_types, rotation = 20, ha = 'center', size = 10, fontname = 'Arial')
elif graphtype == 'bar':
    plt.xticks(xtickpos, site_types, rotation = 20, ha = 'right', size = 10, fontname = 'Arial')

# Set the tick labels font name
for label in ax.get_yticklabels():
    label.set_fontname('Arial')

# write label for x and y axis
ax.set_ylabel('Nucleotide polymorphism', color = 'black',  size = 10, ha = 'center', fontname = 'Arial')
ax.set_xlabel('Phylogenetic conservation', color = 'black', size = 10, ha = 'center', fontname = 'Arial')

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
        print(site_types[i], site_types[j], Pval)    
        
# I already determined that all site categories are significantly different
# using Wilcoxon rank sum tests, so we need now to add letters to show significance

# annotate figure to add significance
# get the x and y coordinates
if graphtype == 'bar':
    y_pos = [0.004, 0.0125, 0.026]
    x_pos = [i + width/2 for i in range(len(site_types))]
elif graphtype == 'box':
    y_pos = [0.035, 0.035, 0.065]
    x_pos = [i + 1 for i in range(len(site_types))]
diff = ['A', 'B', 'C']

for i in range(len(diff)):
    ax.text(x_pos[i], y_pos[i], diff[i], horizontalalignment='center',
            verticalalignment='center', color = 'black', fontname = 'Arial', size = 10)

if graphtype == 'box':
    # save figure
    fig.savefig('DiversityMatureConservation_box.pdf', bbox_inches = 'tight')
elif graphtype == 'bar':
    # save figure
    fig.savefig('DiversityMatureConservation_bar.pdf', bbox_inches = 'tight')

