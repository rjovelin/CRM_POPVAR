# -*- coding: utf-8 -*-
"""
Created on Thu Jun  4 16:55:43 2015

@author: Richard
"""


import os
from miRNA_target import *
from genomic_coordinates import *
from manipulate_sequences import *
from parse_targetscan_output import *

# use this script to generate the summary tables of the remanei miRNA target
# sites for analyzing divergence and diversity

# compute threshold based on the distribution of elegans UTR length
UTR_length = celegans_three_prime_UTR_length('../Genome_Files/c_elegans.PRJNA13758.WS248.annotations.gff3')
threshold = get_percentile(UTR_length, 99)
print('threshold', threshold)

# generate a dict of seeds: mirna pairs
seeds = seed_mirnas('Cremanei_mature.fasta')
print('seeds', len(seeds))

# generate a dict of site_info : [seeds targeting that site]
targets = parse_crm_cla_sites('Crm_Cla_UTR_seq_targetscan.txt', 'Crm_Cla_predicted_sites_targetscan.txt', seeds)
print('targets', len(targets))

# get the coordinates of the UTR/downstream seq
# {TS1 : [chromo, start, end, orientation]}
UTR_coord = get_three_prime_UTR_positions('../Genome_Files/356_10172014.gff3', '../Genome_Files/noamb_356_v1_4.txt', threshold)
print('UTR', len(UTR_coord))

# generate the summary file used to polarise changes in remanei target sites
# and to compute divergence at target sites between remanei and latens
make_table_cremanei_clatens_sites('../Genome_Files/356_10172014.gff3', '../Genome_Files/noamb_356_v1_4.txt', threshold, seeds, targets, 'Cremanei_Clatens_miRNA_sites.txt')
print('done with crm cla sites table')

# parse targetscan outputs into dicts of site info: [seeds]
crm_sites = parse_targetscan_output('Crm_UTR_seq_targetscan.txt', 'Cremanei_predicted_sites_targetscan.txt', seeds, 'remanei')
crm_cla_conserved_sites = parse_targetscan_output('Crm_Cla_UTR_seq_targetscan.txt', 'Crm_Cla_predicted_sites_targetscan.txt', seeds, 'latens')
crm_cel_conserved_sites = parse_targetscan_output('Crm_Cel_UTR_seq_targetscan.txt', 'Crm_Cel_predicted_sites_targetscan.txt', seeds, 'elegans')
print('done parsing file')

# make summary table for all the remanei target sites used to compute diversity at remanei target sites
summary_table_cremanei_target_sites(crm_sites, crm_cla_conserved_sites, crm_cel_conserved_sites, '../Genome_Files/noamb_356_v1_4.txt',
                                    '../Genome_Files/356_10172014.gff3', seeds, threshold, 'Cremanei_miRNA_sites.txt')
print('done writing summary file')
