# -*- coding: utf-8 -*-
"""
Created on Fri Sep  4 16:33:29 2015

@author: RJovelin
"""


from Cel_UTR import *
from accessories import *
from piRNAs import *
from miRfam_targetscan_input import *
from parse_targetscan_output import *
import os



# get the allele counts for all sites with coverage, exclude sites with sample size < 10
chromo_sites = get_non_coding_snps('../SNP_files/', 10)
print('got allele counts in genome')

# get the remanei targets with high DAF
crm_targets = find_miRNA_targets_with_high_DAF(chromo_sites, '../CREM_CLA_protein_divergence/noamb_356_v1_4.txt',
                                               'Cremanei_Clatens_miRNA_sites.txt', 'Crm_Cla_UTR_sequences/',
                                               '../CREM_CLA_protein_divergence/unique_transcripts.txt', 0.8, 1)

# get the 99th percentile of UTR length in elegans
UTR_length = celegans_three_prime_UTR_length('c_elegans.PRJNA13758.WS248.annotations.gff3')
threshold = get_percentile(UTR_length, 99)

# get the latens miRNA targets with high DAF
cla_targets = find_non_remanei_targets_in_latens_with_high_DAF(chromo_sites, 'Cremanei_miRBase21_mature.txt',
                                                               '../CREM_CLA_protein_divergence/noamb_356_v1_4.txt', '../CREM_CLA_protein_divergence/356_10172014.gff3',
                                                               threshold, 'Clatens_specific_miRNA_targets_coordinates.txt',
                                                               'Crm_Cla_UTR_sequences/', '../CREM_CLA_protein_divergence/unique_transcripts.txt', 0.8, 1)
                                                               
# consider only heptamer seed matches: remove 7mer-1a and convert 8mer to 7mer-m8
crm_targets = merge_8mer_with_7merm8(crm_targets)
cla_targets = miRNA_target.merge_8mer_with_7merm8(cla_targets)

                                                               
# make a list of seeds
# create a dict of seed sequence : list of mirnas sharing the seed
seed_mature_pairs = seed_mirnas('Cremanei_miRBase21_mature.txt') 
# make a list of remanei seeds
seeds = [i for i in seed_mature_pairs]

# sort sites into site loss, gain and site to site conversion
loss, gain, site_to_site =  infer_target_gain_loss(crm_targets, cla_targets, seeds)
                                                

# continue here                
