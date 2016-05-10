# -*- coding: utf-8 -*-
"""
Created on Wed Sep 16 10:01:56 2015

@author: RJovelin
"""

from accessories import *
from sites_with_coverage import *
import os

# get the allele counts for all sites with coverage, exclude sites with sample size < 10
# [chromo: {site : [ref_allele, alt_allele, count_ref, count alt]}]
chromo_sites = get_non_coding_snps('../SNP_files/', 10)
print('got allele counts in genome')


# make a list of positions for scaffold_144
positions = [i for i in chromo_sites['linkage_group_1']]
# sort list
positions.sort()


# make input file
newfile = open('testfileSweeD.SF', 'w')
# write header
newfile.write('position\tx\tn\tfolded' + '\n')
for i in positions:
    ref_count  = chromo_sites['linkage_group_1'][i][2]
    alt_count = chromo_sites['linkage_group_1'][i][3]
    if ref_count != 0 and alt_count != 0:
        # site is polymorphic
        # find minor allele
        ref = chromo_sites['linkage_group_1'][i][0]
        alt = chromo_sites['linkage_group_1'][i][1]
        if alt_count + ref_count >= 10:
            if alt_count <= ref_count:
                # alt is minor
                x = alt_count
            elif ref_count < alt_count:
                # ref is minor
                x = ref_count
            n = ref_count + alt_count
            newfile.write('\t'.join([str(i+1), str(x), str(n), '0']) + '\n')

newfile.close()

# get the gridsize paramater
genome = convert_fasta('../CREM_CLA_protein_divergence/noamb_356_v1_4.txt')
GRIDSIZE = int(len(genome['linkage_group_1']) / 50)
print(GRIDSIZE)

# run SweepFinder
os.system('SweeD -name testfileSweeDoutscaffold144 -input testfileSweeD.SF -grid ' + str(GRIDSIZE))

