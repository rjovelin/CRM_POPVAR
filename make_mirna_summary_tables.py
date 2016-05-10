# -*- coding: utf-8 -*-
"""
Created on Thu Jun  4 16:55:43 2015

@author: Richard
"""

# import script
import os
from miRNA_target import *
from genomic_oordinates import *
from manipulate_sequences import *
from parse_targetscan_output import *


# use this script to generate the summary tables of the remanei miRNA target
# sites for analyzing divergence and diversity

# compute threshold based on the distribution of elegans UTR length
UTR_length = celegans_three_prime_UTR_length('c_elegans.PRJNA13758.WS248.annotations.gff3')
threshold = get_percentile(UTR_length, 99)

# make summary table of the target sites from the analysis of the remanei - latens orthologs
make_table_cremanei_clatens_sites('../CREM_CLA_protein_divergence/356_10172014.gff3', 'Crm_Cla_UTR_seq_targetscan.txt',
                                  'Crm_Cla_predicted_sites_targetscan.txt', 'Cremanei_miRBase21_mature.txt', '../CREM_CLA_protein_divergence/noamb_356_v1_4.txt',
                                  threshold, 'Cremanei_Clatens_miRNA_sites.txt')

# parse targetscan outputs
crm_sites = parse_targetscan_output('Crm_UTR_seq_targetscan.txt', 'Cremanei_predicted_sites_targetscan.txt', 'Cremanei_miRBase21_mature.txt', 'remanei')
crm_cla_conserved_sites = parse_targetscan_output('Crm_Cla_UTR_seq_targetscan.txt', 'Crm_Cla_predicted_sites_targetscan.txt', 'Cremanei_miRBase21_mature.txt', 'latens')
crm_cel_conserved_sites = parse_targetscan_output('Crm_Cel_UTR_seq_targetscan.txt', 'Crm_Cel_predicted_sites_targetscan.txt', 'Cremanei_miRBase21_mature.txt', 'elegans')

# make summary table for all the remanei target sites
summary_table_cremanei_target_sites(crm_sites, crm_cla_conserved_sites, crm_cel_conserved_sites, '../CREM_CLA_protein_divergence/noamb_356_v1_4.txt',
                                    '../CREM_CLA_protein_divergence/356_10172014.gff3', 'Cremanei_miRBase21_mature.txt', threshold, 'Cremanei_miRNA_sites.txt')


