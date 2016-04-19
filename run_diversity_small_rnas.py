# -*- coding: utf-8 -*-
"""
Created on Wed Aug 12 22:05:38 2015

@author: Richard
"""
from accessories import *
from divergence import *
from miRNA_target import *
import numpy as np
from scipy import stats
import math


# compute theta at synonymous sites
theta_syn = compute_theta_diversity('../CREM_CLA_protein_divergence/CDS_SNP_DIVERG.txt',
                                    '../CREM_CLA_protein_divergence/unique_transcripts.txt',
                                    'SYN', 10)
                                   
# make a list of theta values
SYN_theta = []
for gene in theta_syn:
    SYN_theta.append(theta_syn[gene])

print('done computing theta at synonymous sites')

# compute theta at nonsynonymous sites
theta_rep = compute_theta_diversity('../CREM_CLA_protein_divergence/CDS_SNP_DIVERG.txt',
                                    '../CREM_CLA_protein_divergence/unique_transcripts.txt',
                                    'REP', 10)


# make a list of theta values
REP_theta = []
for gene in theta_rep:
    REP_theta.append(theta_rep[gene])
    
print('done computing theta at replacement sites')

# get the allele counts for all sites with coverage
chromo_sites = get_non_coding_snps('../SNP_files/', 10)

print('got allele counts at all sites')

# get miRNA coordinates {chromo: [[start, end, orientation]]}
mirna_coord = get_mirna_loci('../miRNA_Target_sites/crm_miRBase21_premina_coordinates.txt')

print('got miRNA coordinates')

# get mature coordinates {chromo: [[start, end orientation]]}
mature_coord = get_mature_loci('../miRNA_Target_sites/crm_mature_miRBase_genomic_coordinates.txt')
 
print('got mature miR coordinates')

# get piRNA coordinates {chromo: [[start, end, orienation]]}
pirna_coord = get_pirna_loci('PX356_piRNA_coord.txt')

print('got piRNA coordinates')

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

# compute theta at pirnas
# create a list to store theta at pirnasd
pirna_theta = []
# loop over chromo
for chromo in pirna_coord:
    # loop over all prina on that chromo:
    for i in range(len(pirna_coord[chromo])):
        # get start, end
        start = pirna_coord[chromo][i][0]
        end = pirna_coord[chromo][i][1]
        # compute theta, accepting only 2 missing sites
        theta = compute_theta_non_coding(chromo_sites, chromo, start, end, 2)
        # check if theta is defined
        if theta != 'NA':
            pirna_theta.append(theta)
            
print('done computing theta at piRNAs')

# create a list of theta
theta_sites = [SYN_theta, REP_theta, mirna_theta, mature_theta, pirna_theta]
# create a list of corrsponding site type
site_types = ['SYN', 'REP', 'miRNA', 'mature_miR', 'piRNA']

# open file to store the results of the analusis
newfile = open('diversity_small_RNAs.txt', 'w')
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
newfile = open('theta_values_small_RNAs.txt', 'w')
for i in REP_theta:
    newfile.write('REP' + '\t' + str(i) + '\n')
for i in SYN_theta:
    newfile.write('SYN' + '\t' + str(i) + '\n')
for i in mirna_theta:
    newfile.write('miRNA' + '\t' + str(i) + '\n')
for i in mature_theta:
    newfile.write('mature' + '\t' + str(i) + '\n')
for i in pirna_theta:
    newfile.write('piRNA' + '\t' + str(i) + '\n')
    
# close file after writing
newfile.close()


    



