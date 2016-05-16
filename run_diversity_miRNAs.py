# -*- coding: utf-8 -*-
"""
Created on Thu May 12 10:52:40 2016

@author: RJovelin
"""


# use this script to compute theta at miRNA loci
# compare theta with synonymous, replacement sites

from accessories import *
from divergence import *
from miRNA_target import *
import numpy as np
from scipy import stats
import math


# compute theta at synonymous sites
theta_syn = compute_theta_diversity('../Genome_Files/CDS_SNP_DIVERG.txt', '../Genome_Files/unique_transcripts.txt', 'SYN', 10)
                             
# make a list of theta values
SYN_theta = []
for gene in theta_syn:
    SYN_theta.append(theta_syn[gene])
print('done computing theta at synonymous sites')
print('Sites\tmean\tmin\tmax')
print('SYN', np.mean(SYN_theta), min(SYN_theta), max(SYN_theta))

# compute theta at nonsynonymous sites
theta_rep = compute_theta_diversity('../Genome_Files/CDS_SNP_DIVERG.txt', '../Genome_Files/unique_transcripts.txt', 'REP', 10)

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
mature_coord = get_mirna_loci('../miRNA_Target_sites/CRM_MatureCoordinatesFinal.txt')
print('got mature miR coordinates')

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


# create a list of theta
theta_sites = [SYN_theta, REP_theta, mirna_theta, mature_theta]
# create a list of corrsponding site type
site_types = ['SYN', 'REP', 'miRNA', 'mature_miR']

# open file to store the results of the analusis
newfile = open('diversity_miRNAs.txt', 'w')
# write header
newfile.write('sites' + '\t' + 'N' + '\t' + 'mean_theta' + '\t' + 'SEM' + '\n')
for i in range(len(theta_sites)):
    newfile.write('\t'.join([site_types[i], str(len(theta_sites[i])), str(np.mean(theta_sites[i])), str(np.std(theta_sites[i]) / math.sqrt(len(theta_sites[i])))])+'\n')

newfile.write('\n')

newfile.write('Wilcoxon rank sum test of mean differences\n')
newfile.write('-' * 43 + '\n')
newfile.write('sites' + '\t' + 'wilcoxon' + '\t' + 'P' + '\n')

# loop over list
for i in range(0, len(theta_sites) - 1):
    for j in range(i+1, len(theta_sites)):
        wilcoxon, p = stats.ranksums(theta_sites[i], theta_sites[j])
        newfile.write('\t'.join([site_types[i] + '_vs_' + site_types[j], str(wilcoxon), str(p)]) + '\n')
# close file after writing
newfile.close()

# open file to dump all theta values
newfile = open('theta_values_miRNAs.txt', 'w')
for i in REP_theta:
    newfile.write('REP' + '\t' + str(i) + '\n')
for i in SYN_theta:
    newfile.write('SYN' + '\t' + str(i) + '\n')
for i in mirna_theta:
    newfile.write('miRNA' + '\t' + str(i) + '\n')
for i in mature_theta:
    newfile.write('mature' + '\t' + str(i) + '\n')
# close file after writing
newfile.close()

