# -*- coding: utf-8 -*-
"""
Created on Tue Aug 18 22:50:48 2015

@author: Richard
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
from get_coding_sequences import *


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

# get the coordinates of the miRNA loci
mirnas_coord = get_mirna_loci('crm_miRBase21_premina_coordinates.txt')

# get the allele counts for miRNA sites
mirna_sites = get_feature_sites(chromo_sites, mirnas_coord)

# compute MAF for mirna sites (sites with sample size < 10 are already excluded)
MAF_mirna = MAF_non_coding(mirna_sites)

print('MAF for miRNA sites done')

# get SNPs flanking miRNAs within 500 bp of the miRNAs
mirna_flanking_snps = get_small_rna_flanking_SNPs(chromo_sites, 'crm_miRBase21_premina_coordinates.txt',
                                                  'miRNA', 500)

# remove positions falling in coding sequences
for chromo in CDS_pos:
    # check if chromo in flanking sites
    if chromo in mirna_flanking_snps:
        # loop over CDS positions
        for i in CDS_pos[chromo]:
            # check if site in flanking
            if i in mirna_flanking_snps[chromo]:
                # remove site
                del mirna_flanking_snps[chromo][i]
                
print('got SNPs flanking miRNAs')

mirna_snps = 0
for chromo in mirna_flanking_snps:
    mirna_snps += len(mirna_flanking_snps[chromo])
print('SNPs within 500 bp of miRNAs: ', mirna_snps)

# resample SNPs and compute MAF
mirna_resampled_MAF = SNP_MAF_randomization(mirna_flanking_snps, 5000, 1000)

# get the proportions of SNPs in each MAF bin
mirna_MAF_proportions = get_MAF_distribution_from_replicates(mirna_resampled_MAF)

# get the SNP proportions for the observed SNPs
empirical_mirna_MAF = SNP_proportions_MAF_bin(MAF_mirna)

# creat a list of keys, being the MAF lower bound in the dicts with MAF proportions
maf_limit = [i for i in empirical_mirna_MAF]
# sort list
maf_limit.sort()

# open file to save results of Z-test
newfile = open('summary_MAF_resampled_SNPs_miRNAs.txt', 'w')

newfile.write('Z-test of SNP proportions in MAF bins\n')
newfile.write('random resampling of 5000 SNPs within 500 bp of miRNAs, exluding miRNA and CDS sites with 1000 replicates\n')

newfile.write('\t'.join(['MAF', 'mean_sample', 'margin', 'stdev_sample', 'l95', 'h95', 'observed', 'z-score', 'alpha_0.05', 'alpha_0.01', 'alpha_0.001']) + '\n')

# loop over sorted keys
for i in maf_limit:
    # compute mean
    average = np.mean(mirna_MAF_proportions[i])
    # get standard deviation
    stdev = np.std(mirna_MAF_proportions[i])
    # compute 95% CI
    stderror = stdev / math.sqrt(len(mirna_MAF_proportions[i]))
    # compute the margin error (critical value = 1.96 for 95% CI)
    margin = 1.96 * stderror
    lCI = average - margin
    hCI = average + margin
    observed = empirical_mirna_MAF[i]
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
    



