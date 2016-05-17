# -*- coding: utf-8 -*-
"""
Created on Thu May 12 10:52:40 2016

@author: RJovelin
"""


# use this script to make a box plot comparing theta for miRNA loci and different site categories
# use Agg backend on server without X server
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
# load modules
import numpy as np
import math
from scipy import stats
# load custom modules
from manipulate_sequences import *
from divergence import *
from miRNA_target import *
from sites_with_coverage import *

# create a set of valid transcripts (1 transcript mapped to 1 gene)
valid_transcripts = get_valid_transcripts('../Genome_Files/unique_transcripts.txt')

# compute theta at synonymous sites
theta_syn = compute_theta_diversity('../Genome_Files/CDS_SNP_DIVERG.txt', valid_transcripts, 'SYN', 10)
# make a list of theta values
SYN_theta = []
for gene in theta_syn:
    SYN_theta.append(theta_syn[gene])
print('done computing theta at synonymous sites')
print('Sites\tmean\tmin\tmax')
print('SYN', np.mean(SYN_theta), min(SYN_theta), max(SYN_theta))

# compute theta at nonsynonymous sites
theta_rep = compute_theta_diversity('../Genome_Files/CDS_SNP_DIVERG.txt', valid_transcripts, 'REP', 10)
# make a list of theta values
REP_theta = []
for gene in theta_rep:
    REP_theta.append(theta_rep[gene])
print('done computing theta at replacement sites')
print('Sites\tmean\tmin\tmax')
print('REP', np.mean(REP_theta), min(REP_theta), max(REP_theta))

# get the allele counts for all sites with coverage
chromo_sites = get_non_coding_snps('../SNP_files/', 10)
print('got allele counts at all sites')

# get miRNA coordinates {chromo: [[start, end, orientation]]}
mirna_coord = get_mirna_loci('CRM_miRNAsCoordinatesFinal.txt')
print('got miRNA coordinates')

# get mature coordinates {chromo: [[start, end orientation]]}
mature_coord = get_mirna_loci('CRM_MatureCoordinatesFinal.txt')
print('got mature miR coordinates')

# make a dict with family level conservation for all miRNAs {name : [chromo, start, end, orientation, conservation]}
famCons = {}
infile = open('CRM_miRNAsCoordinatesFinal.txt')
infile.readline()
for line in infile:
    line = line.rstrip()
    if line != '':
        line = line.split('\t')
        name, chromo, start, end, orientation, conservation = line[0], line[1], int(line[2]) -1, int(line[3]), line[4], line[-1]
        famCons[name] = [chromo, start, end, orientation, conservation]
infile.close()

# compute theta for mirnas
# create a list to store theta at miRNA loci
mirna_theta = []
# loop over chromo in mirna_coord
for chromo in mirna_coord:
    # loop over all mirnas on that chromo
    for i in range(len(mirna_coord[chromo])):
        # get start
        start = mirna_coord[chromo][i][0]
        end = mirna_coord[chromo][i][1]
        # compute theta, accepting maximum of 2 missing sites
        theta = compute_theta_non_coding(chromo_sites, chromo, start, end, 2)
        # check that theta is defined
        if theta != 'NA':
            mirna_theta.append(theta)
print('done computing theta at miRNAs')
print('Sites\tmean\tmin\tmax')
print('miRNAs', np.mean(mirna_theta), min(mirna_theta), max(mirna_theta))

# compute theta at mature miRs
# create a list to store theta at miR
mature_theta = []
# loop over chromo in mature coord
for chromo in mature_coord:
    # loop over all mature on that chromo
    for i in range(len(mature_coord[chromo])):
        # get start
        start = mature_coord[chromo][i][0]
        end = mature_coord[chromo][i][1]
        # compute theta, accpeting only 2 missing sites
        theta = compute_theta_non_coding(chromo_sites, chromo, start, end, 2)
        # check if theta is defined
        if theta != 'NA':
            mature_theta.append(theta)
print('done computing theta at mature miRs')
print('Sites\tmean\tmin\tmax')
print('mature', np.mean(mature_theta), min(mature_theta), max(mature_theta))


# compute thetas for miRNAs with different level of conservation
# create lists to store theta at miRNA for different level of family conservation
caeno_mirna_theta, crmcla_mirna_theta, crm_mirna_theta = [], [], []
# loop over mirna name in {name : [chromo, start, end, orientation, conservation]}
for name in famCons:
    chromo, start, end, conservation = famCons[name][0], famCons[name][1], famCons[name][2], famCons[name][-1]
    theta = compute_theta_non_coding(chromo_sites, chromo, start, end, 2)
    # check if theta is defined
    if theta != 'NA':
        # check conservation level
        if famCons[name][-1] == 'Crm':
            crm_mirna_theta.append(theta)
        elif famCons[name][-1] == 'CrmCla':
            crmcla_mirna_theta.append(theta)
        elif famCons[name][-1] == 'Caeno':
            caeno_mirna_theta.append(theta)



# create list of data
alldata = [SYN_theta, REP_theta, mirna_theta, mature_theta, caeno_mirna_theta, crmcla_mirna_theta, crm_mirna_theta]
# make a list of datatype
# create a list of corrsponding site type
site_types = ['Syn', 'Rep', 'miRNA', 'miR', 'Caeno', 'Crm,Cla', 'Crm']

# create figure
fig = plt.figure(1, figsize = (6,4))
# add a plot to figure (1 row, 1 column, 1st plot)
ax = fig.add_subplot(1, 1, 1)

# use a boxplot
bp = ax.boxplot(alldata, showmeans = False, showfliers = False, widths = 0.7, labels = site_types, patch_artist = True) 
 
# create a list of colors (seee http://colorbrewer2.org/)
color_scheme = ['#6e016b', '#88419d', '#8c6bb1', '#8c96c6','#9ebcda', '#bfd3e6', '#edf8fb']
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
# change the mean marker and marker if showmean is True
#for mean in bp['means']:
#    mean.set(marker = 'o', markeredgecolor = 'black', markerfacecolor = 'black', markersize = 4)
    
# restrict the x and y axis to the range of data
ax.set_ylim([0, 0.15])
# create a list with range of x-axis values
xvals = [i + 0.5 for i in range(len(site_types) + 1)]
# Set a buffer around the edge of the x-axis
plt.xlim([min(xvals)- 0.5, max(xvals)+ 0.5])


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

   
# Set the tick labels font name
for label in ax.get_yticklabels():
    label.set_fontname('Arial')

# write label for x and y axis
ax.set_ylabel('Nucleotide polymorphism\n', color = 'black',  size = 10, ha = 'center', fontname = 'Arial')
ax.set_xlabel('Site categories\n', color = 'black', size = 10, ha = 'center', fontname = 'Arial')

# add labels to x-ticks, rotate and align right
ax.set_xticklabels(site_types, ha = 'center', size = 10, fontname = 'Arial')

# remove lines around the frame
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['left'].set_visible(False)
ax.spines['bottom'].set_visible(False)

# add a light grey horizontal grid to the plot, semi-transparent, 
ax.yaxis.grid(True, linestyle='--', which='major', color='lightgrey', alpha=0.5)  
# hide these grids behind plot objects
ax.set_axisbelow(True)

# give more space to the y labels
fig.subplots_adjust(left=0.2)



#################

 
## annotate figure to add significance
## get the x and y coordinates
#yvalues = []
#for i in range(len(all_data)):
#    yvalues.extend(all_data[i])
## get max and min y values
#y_min, y_max = min(yvalues), max(yvalues)
#
 
 
 
 
 
 
#for i in range(len(theta_sites)):
#    newfile.write('\t'.join([site_types[i], str(len(theta_sites[i])), str(np.mean(theta_sites[i])), str(np.std(theta_sites[i]) / math.sqrt(len(theta_sites[i])))])+'\n')
#
#newfile.write('\n')
#
#newfile.write('Wilcoxon rank sum test of mean differences\n')
#newfile.write('-' * 43 + '\n')
#newfile.write('sites' + '\t' + 'wilcoxon' + '\t' + 'P' + '\n')
#
## loop over list
#for i in range(0, len(theta_sites) - 1):
#    for j in range(i+1, len(theta_sites)):
#        wilcoxon, p = stats.ranksums(theta_sites[i], theta_sites[j])
#        newfile.write('\t'.join([site_types[i] + '_vs_' + site_types[j], str(wilcoxon), str(p)]) + '\n')
## close file after writing
#newfile.close()
#
## open file to dump all theta values
#newfile = open('theta_values_miRNAs.txt', 'w')
#for i in REP_theta:
#    newfile.write('REP' + '\t' + str(i) + '\n')
#for i in SYN_theta:
#    newfile.write('SYN' + '\t' + str(i) + '\n')
#for i in mirna_theta:
#    newfile.write('miRNA' + '\t' + str(i) + '\n')
#for i in mature_theta:
#    newfile.write('mature' + '\t' + str(i) + '\n')
## close file after writing
#newfile.close()
#
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
## compare read depth between WGS and RNAseq for each cancer
#for i in range(0, len(all_data), 2):
#    # get the P value of Wilcoxon rank sum test
#    Pval = stats.ranksums(all_data[i], all_data[i+1])[1]
#    # get stars for significance
#    if Pval > 0.05:
#        P = 'N.S.'
#    elif Pval < 0.05 and Pval > 0.01:
#        P = '*'
#    elif Pval < 0.01 and Pval > 0.001:
#        P = '**'
#    elif Pval < 0.001:
#        P = '***'
#    
#    # add bracket
#    ax.annotate("", xy=(i+1, y_max + 200), xycoords='data', xytext=(i+2, y_max + 200),
#                textcoords='data', arrowprops=dict(arrowstyle="-", ec='#aaaaaa',
#                connectionstyle="bar,fraction=0.2", linewidth = 0.75))
#    # add stars for significance
#    if P == 'N.S.':
#        ax.text(i + 1.5, y_max + abs(y_max - y_min)*0.08, P, horizontalalignment='center',
#                verticalalignment='center', color = 'grey', fontname = 'Arial', size = 6)
#    else:
#        ax.text(i + 1.5, y_max + abs(y_max - y_min)*0.05, P, horizontalalignment='center',
#                verticalalignment='center', color = 'grey', fontname = 'Arial')
#

##################




# save figure
fig.savefig('testfile.pdf', bbox_inches = 'tight')
    

