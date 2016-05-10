# -*- coding: utf-8 -*-
"""
Created on Mon Aug  3 02:03:56 2015

@author: Richard
"""

from piRNAs import *
from accessories import *


# define a small function to convert a num to string
Gstr = lambda x : str(x)

# tolerate 0 mismatches between pirna and targets

# find targets in coding sequences
CDS_targets = find_pirna_CDS_targets('PX356_piRNA_coord.txt', '../CREM_CLA_protein_divergence/noamb_PX356_all_CDS.fasta',
                                     '../CREM_CLA_protein_divergence/noamb_356_v1_4.txt', 0, '../CREM_CLA_protein_divergence/unique_transcripts.txt')
print('found ', len(CDS_targets), ' targets in CDS with ', 0, ' mismatches')
# check if targets are present in CDS
if len(CDS_targets) != 0:
    # open file for writing
    newfile = open('piRNA_CDS_targets_' + str(0) + '_mismatch.txt', 'w')
    # write header
    newfile.write('# start and end position are 1-based indices relative to the CDS sequence and not to the genome sequence\n')
    newfile.write('\t'.join(['gene', 'start', 'end', 'pirna_name', 'target_site_sequence']) + '\n')
    # loop over gene in dict
    for gene in CDS_targets:
        # loop over target in CDS
        for target_site in CDS_targets[gene]:
            # convert start and position to 1-based index
            target_site[0] += 1
            # convert all items in target site to strings
            target_site = list(map(Gstr, target_site))
            # write to file
            newfile.write(gene + '\t')
            newfile.write('\t'.join(target_site) + '\n')
    # close after writing
    newfile.close()
    
    
    
# find targets in repeat sequences

TE_targets = find_pirna_TE_targets('PX356_piRNA_coord.txt', '356_v1_4.fasta.out',
                                   '../CREM_CLA_protein_divergence/noamb_356_v1_4.txt', 0)
print(len(TE_targets), ' repeat type are targeted by piRNAs')
total = 0
if len(TE_targets) != 0:
    for repname in TE_targets:
        for target in TE_targets[repname]:
            total += 1
print('found ', total, ' targets in TEs with ', 0, ' mismatches')
# check if targets are present in TE
if len(TE_targets) != 0:
    # open file for writing
    newfile = open('piRNA_TE_targets_' + str(0) + '_mismatch.txt', 'w')
    # write header
    newfile.write('# start and end position are 1-based indices relative to the genome sequence\n')
    newfile.write('\t'.join(['repeat', 'family', 'chromo', 'start', 'end', 'sense', 'piRNA_name', 'target_site_sequence']) + '\n')    
    # loop over repeat in dict
    for repname in TE_targets:
        # loop over target in repat
        for target_site in TE_targets[repname]:
            # convert start to 1-based index
            target_site[2] += 1
            # convert all items in target site to string
            target_site = list(map(Gstr, target_site))
            # write to file
            newfile.write(repname + '\t')
            newfile.write('\t'.join(target_site) + '\n')
    # close file after writing
    newfile.close()
 

    
