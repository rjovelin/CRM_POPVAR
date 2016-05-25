# -*- coding: utf-8 -*-
"""
Created on Wed Aug 19 09:21:34 2015

@author: RJovelin
"""

import math
import numpy as np
from premature_stops import *
from scipy import stats
from accessories import *
from miRNA_target import *
from piRNAs import *
from sites_with_coverage import *
from Cel_UTR import *
from cel_UTR_length import *
from randomize_SNPs import *


Gstr = lambda x : str(x)


# get the allele counts for all sites with coverage, exclude sites with sample size < 10
chromo_sites = get_non_coding_snps('../SNP_files/', 10)

print('got allele counts in genome')

# get CDS_coord {TS1: [chromo, sense, [(s1, end1), (s2, end2)]]}
CDS_coord = get_CDS_positions('../CREM_CLA_protein_divergence/356_10172014.gff3')
# create new dict in the form {chromo: {set of positions}}
CDS_pos = {}
for gene in CDS_coord:
    # get chromo
    chromo = CDS_coord[gene][0]
    # loop over cds coordinates 
    for i in range(len(CDS_coord[gene][2])):
        # get start and end positions in a list
        # convert to 0-based index
        start  = CDS_coord[gene][2][i][0] - 1
        end = CDS_coord[gene][2][i][1]
        # check if chromo in CDS_pos
        if chromo not in CDS_pos:
            CDS_pos[chromo] = set()
        for j in range(start, end):
            CDS_pos[chromo].add(j)
        
print('got CDS coord')

# get the coordinates of the piRNA loci
pirnas_coord = get_pirna_loci('PX356_piRNA_coord.txt')

# get the allele counts for piRNA sites
pirna_sites = get_feature_sites(chromo_sites, pirnas_coord)

# compute MAF for piRNA sites, (sites with sample size < 10 are already excluded)
MAF_pirna = MAF_non_coding(pirna_sites)

print('MAF for piRNA sites done')

# get all the piRNAs positions in the genome
pirna_pos = get_small_rna_sites(pirnas_coord)

print('got piRNA positions')

# get SNPs flanking piRNAs within 500 bp of the piRNAs
pirna_flanking_snps = get_small_rna_flanking_SNPs(chromo_sites, 'PX356_piRNA_coord.txt', 'piRNA', 500)

# remove positions falling in coding sequences
for chromo in CDS_pos:
    # check if chromo in flanking sites
    if chromo in pirna_flanking_snps:
        # loop over CDS positions
        for i in CDS_pos[chromo]:
            # check if site in flanking
            if i in pirna_flanking_snps[chromo]:
                # remove site
                del pirna_flanking_snps[chromo][i]

print('got SNPs flanking piRNAs')

pirna_snps = 0
for chromo in pirna_flanking_snps:
    pirna_snps += len(pirna_flanking_snps[chromo])
print('SNPs within 500 bp of piRNAs: ', pirna_snps)

# remove chromos if no SNPs are on chromo
to_remove = []
for chromo in pirna_flanking_snps:
    if len(pirna_flanking_snps[chromo]) == 0:
        to_remove.append(chromo)
for chromo in to_remove:
    del pirna_flanking_snps[chromo]
    
print('removed {0} chromo with no SNPs after filtering'.format(len(to_remove)))


# resample SNPs and compute MAF
pirna_resampled_MAF = SNP_MAF_randomization(pirna_flanking_snps, 5000, 1000)

# get the proportions of SNPs in each MAF bin
pirna_MAF_proportions = get_MAF_distribution_from_replicates(pirna_resampled_MAF)

# get the SNP proportions for the observed SNPs
empirical_pirna_MAF = SNP_proportions_MAF_bin(MAF_pirna)

# creat a list of keys, being the MAF lower bound in the dicts with MAF proportions
maf_limit = [i for i in empirical_pirna_MAF]
# sort list
maf_limit.sort()

# open file to save results of Z-test
newfile = open('summary_MAF_resampled_SNPs_piwiRNAs.txt', 'w')

newfile.write('Z-test of SNP proportions in MAF bins\n')
newfile.write('random resampling of 5000 SNPs within 500 bp of piRNAs, exluding piRNA and CDS sites with 1000 replicates\n')

newfile.write('\t'.join(['MAF', 'mean_sample', 'margin', 'stdev_sample', 'l95', 'h95', 'observed', 'z-score', 'alpha_0.05', 'alpha_0.01', 'alpha_0.001']) + '\n')

# loop over sorted keys
for i in maf_limit:
    # compute mean
    average = np.mean(pirna_MAF_proportions[i])
    # get standard deviation
    stdev = np.std(pirna_MAF_proportions[i])
    # compute 95% CI
    stderror = stdev / math.sqrt(len(pirna_MAF_proportions[i]))
    # compute the margin error (critical value = 1.96 for 95% CI)
    margin = 1.96 * stderror
    lCI = average - margin
    hCI = average + margin
    observed = empirical_pirna_MAF[i]
    z_score = (observed - average) / stdev
    # critical values for 1-sample 2-tailed z-test: 0.05: +- 1.96, 0.01: +- 2.58, 0.001: +-3.27
    # Ho : obsvered == mean, H1: observed != mean
    if z_score < -1.96 or z_score > 1.96:
        P5 = '*'
    elif -1.96 <= z_score <= 1.96:
        P5 = 'NS'
    if z_score < -2.58 or z_score > 2.58:
        P1 = '*'
    elif -2.58 <= z_score <= 2.58:
        P1 = 'NS'
    if z_score < -3.27 or z_score > 3.27:
        P01 = '*'
    elif -3.27 <= z_score <= 3.27:
        P01 = 'NS'
    newfile.write('\t'.join([str(i), str(average), str(margin), str(stdev), str(lCI), str(hCI), str(observed), str(z_score), str(P5), str(P1), str(P01)]) + '\n')
    
# close file after writing
newfile.close()
    



