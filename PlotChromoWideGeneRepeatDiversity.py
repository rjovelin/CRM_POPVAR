# -*- coding: utf-8 -*-
"""
Created on Thu Aug 20 13:17:35 2015

@author: RJovelin
"""

# use this script to plot nucleotide diversity with repeat or gene density for 
# linkage groups LG1, LG2 and LG4

# usage python3 PlotChromoWideGeneRepeatDiversity.py [options]
#- [repeats/genes] : plot gene or repeat density along with nucleotide diversity


# use Agg backend on server without X server
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
# import modules
import numpy as np
from scipy import stats
import math
import sys
# import custom modules
from manipulate_sequences import *
#from piRNAs import *
from miRNA_target import *
from repeats_TEs import *
from sliding_windows import *
from sites_with_coverage import *
from divergence import *


# convert genome fasta to dict
genome = convert_fasta('../Genome_Files/noamb_356_v1_4.txt')
print('genome converted to fasta dict')

# get the gene coordinates {gene1 : [chromo, start, end, sense]}
genes_coord = get_genes_coordinates('../Genome_Files/356_10172014.gff3')
print('got gene coordinates')

# get the set of valid transcripts
valid_transcripts = get_valid_transcripts('../Genome_Files/unique_transcripts.txt')
print('got the set of valid transcripts')

# get a dict with {chromo: [start position (0-based) of the valid transcripts]}
genes_start = {}
for gene in genes_coord:
    if gene in valid_transcripts:
        chromo = genes_coord[gene][0]
        start = genes_coord[gene][1] -1
        # check if chromo in genes_start
        if chromo in genes_start:
            genes_start[chromo].append(start)
        else:
            genes_start[chromo] = [start]
# sort positions
for chromo in genes_start:
    genes_start[chromo].sort()
print('got the genes\'s first positions')

# get the coordinates of the repeats
# {repeat : [[chromo, start, end, orientation], [chromo, start, end, orientation]]}
repeats_coord = get_repeats_coord('../Genome_Files/356_v1_4.fasta.out', False)
print('got the repeat coordinates')

# create a dict with {chromo: [list of repeat start positions 0-based]}
repeats_start = {}
# loop over repeat family
for fam in repeats_coord:
    # loop over instances of repeat
    for i in range(len(repeats_coord[fam])):
        # get chromo
        chromo = repeats_coord[fam][i][0]
        # get start position
        start = repeats_coord[fam][i][1]
        # check if chromo in dict
        if chromo in repeats_start:
            # add start position to list
            repeats_start[chromo].append(start)
        else:
            # initiate list
            repeats_start[chromo] = [start]
# sort positions
for chromo in repeats_start:
    repeats_start[chromo].sort()
print('got the repeats\'s first positions')




########################## REMOVE function

# use this function to generate save the number of pirnas per window
def cluster_pirnas(pirna_start, genome, chromo, window_size, outputfile):
    '''
    (dict, dict, str, int, file) -> file
    Take a dictionnary with chromo: list of start position pairs, a dict
    with chromo : sequence pairs, the focal chromosome, the size of the window,
    and save the number of pirnas per window in outputfile
    '''
    
    # make list of size window that contains only 0s:
    # each value in the list is the count of position for the range [0 - window[ etc
    range_counts = [0] * (len(genome[chromo]) // window_size)
    
    # loop over starting positions
    for start in pirna_start[chromo]:
        # determine the index in the list range_count where the position should be added
        which_range = start // window_size
        if which_range == len(range_counts):
            which_range -= 1
        # count pirnas
        range_counts[which_range] += 1
        
    # open file for writing
    newfile = open(outputfile, 'w')
    # write header
    newfile.write('Range' + '\t' + 'Lower_point' + '\t' + 'Midpoint' + '\t' + 'Higher_point' + '\t' + 'Count' + '\t' + '\n')
        
    # loop over indices of list
    for i in range(len(range_counts)):
        newfile.write('[' + str(i * window_size) + '-' + str((i * window_size) + window_size -1) + ']' + '\t')
        newfile.write(str(i * window_size) + '\t')
        newfile.write(str(int(((i * window_size) + window_size) / 2)) + '\t')
        newfile.write(str((i * window_size) + window_size) + '\t')
        newfile.write(str(range_counts[i]) + '\n')
        
    newfile.close()












################################
















# create files with gene densities per 50 Kb window 
cluster_pirnas(genes_start, genome, 'linkage_group_1', 50000, 'LG1_gene_density_50Kb.txt')
cluster_pirnas(genes_start, genome, 'linkage_group_2', 50000, 'LG2_gene_density_50Kb.txt')
cluster_pirnas(genes_start, genome, 'linkage_group_4', 50000, 'LG4_gene_density_50Kb.txt')
print('saved gene densities to files')

# create files with repeat densities per 50 Kb window
cluster_pirnas(repeats_start, genome, 'linkage_group_1', 50000, 'LG1_repeat_density_50Kb.txt')
cluster_pirnas(repeats_start, genome, 'linkage_group_2', 50000, 'LG2_repeat_density_50Kb.txt')
cluster_pirnas(repeats_start, genome, 'linkage_group_4', 50000, 'LG4_repeat_density_50Kb.txt')
print('saved repeat densities to files')

# compute theta per 50 Kb window
# get the allele counts for all sites with coverage, keep all sites 
chromo_sites = get_non_coding_snps('../SNP_files/', 0) 
print('got allele counts at all sites')

# create a list of starting positions of each window on chromo
LG1_positions = [i for i in range(0, len(genome['linkage_group_1']), 50000)]
LG2_positions = [i for i in range(0, len(genome['linkage_group_2']), 50000)]
LG4_positions = [i for i in range(0, len(genome['linkage_group_4']), 50000)]

# open file for writing
newfile = open('LG1_theta_50Kb.txt', 'w')
# loop over start position in LG1:
for i in LG1_positions:
    # compute theta per window
    # use a large threshold > window length to accept any number of missing site  
    theta = compute_theta_non_coding(chromo_sites, 'linkage_group_1', i, i + 50000, 100000)
    newfile.write(str(i) + '\t' + str(theta) + '\n')
# close file
newfile.close()
print('theta per windows saved to file for LG1')

# open file for writing
newfile = open('LG2_theta_50Kb.txt', 'w')
# loop over start positions in LG2
for i in LG2_positions:
    # compute theta per window
    # use a large threshold > window length to accept any number of missing site 
    theta = compute_theta_non_coding(chromo_sites, 'linkage_group_2', i, i + 50000, 100000)
    newfile.write(str(i) + '\t' + str(theta) + '\n')
# close file
newfile.close()
print('theta per windows saved to file for LG2')

# open file for writing
newfile = open('LG4_theta_50Kb.txt', 'w')
# loop over start positions in LG4:
for i in LG4_positions:
    # compute theta per window
    # use a large threshold > window length to accept any number of missing site  
    theta = compute_theta_non_coding(chromo_sites, 'linkage_group_4', i, i + 50000, 100000)
    newfile.write(str(i) + '\t' + str(theta) + '\n')
# close file
newfile.close()
print('theta per windows saved to file for LG4')


###################################################################

# -*- coding: utf-8 -*-
"""
Created on Fri Sep 11 21:10:55 2015

@author: Richard
"""



from accessories import *
from piRNAs import *
from miRNA_target import *
from repeats_TEs import *
from sliding_windows import *
from sites_with_coverage import *
from divergence import *
import numpy as np
from scipy import stats
import math


# convert genome fasta to dict
genome = convert_fasta('noamb_356_v1_4.txt')
print('genome converted to fasta dict')

# compute theta per 50 Kb window in all strains
# get the allele counts for all sites with coverage, keep all sites 
chromo_sites = get_all_strains_snps('../PB_ON_SNP_files/', 0)
print('got allele counts at all sites')

# create a list of starting positions of each window on chromo
LG1_positions = [i for i in range(0, len(genome['linkage_group_1']), 50000)]
LG2_positions = [i for i in range(0, len(genome['linkage_group_2']), 50000)]
LG4_positions = [i for i in range(0, len(genome['linkage_group_4']), 50000)]

# open file for writing
newfile = open('LG1_theta_PB+ON_50Kb.txt', 'w')
# loop over start position in LG1:
for i in LG1_positions:
    # compute theta per window
    # use a large threshold > window length to accept any number of missing site  
    theta = compute_theta_non_coding(chromo_sites, 'linkage_group_1', i, i + 50000, 100000)
    newfile.write(str(i) + '\t' + str(theta) + '\n')
# close file
newfile.close()
print('theta per windows saved to file for LG1')

# open file for writing
newfile = open('LG2_theta_PB+ON_50Kb.txt', 'w')
# loop over start positions in LG2
for i in LG2_positions:
    # compute theta per window
    # use a large threshold > window length to accept any number of missing site 
    theta = compute_theta_non_coding(chromo_sites, 'linkage_group_2', i, i + 50000, 100000)
    newfile.write(str(i) + '\t' + str(theta) + '\n')
# close file
newfile.close()
print('theta per windows saved to file for LG2')

# open file for writing
newfile = open('LG4_theta_PB+ON_50Kb.txt', 'w')
# loop over start positions in LG4:
for i in LG4_positions:
    # compute theta per window
    # use a large threshold > window length to accept any number of missing site  
    theta = compute_theta_non_coding(chromo_sites, 'linkage_group_4', i, i + 50000, 100000)
    newfile.write(str(i) + '\t' + str(theta) + '\n')
# close file
newfile.close()
print('theta per windows saved to file for LG4')




