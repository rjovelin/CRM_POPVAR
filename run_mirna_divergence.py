# -*- coding: utf-8 -*-
"""
Created on Tue Jun 28 13:50:20 2016

@author: RJovelin
"""

# use this script to compute divergence between remanei and latens orthologs


# parse ortholog file to get the remanei and latens orthologs


# create a dict of latens and remanei mirnas orthologs
orthos = {}
infile = open('CremClamiRNAOrthologs.txt')
for line in infile:
    if line.rstrip() != '':
        line = line.rstrip().split('\t')
        orthos[line[0]] = line[1]
infile.close()








# create a dict with latens hairpin

# create a dict with latens mature sequences

# create a dict with remanei hairpin

# create a dict with remanei mature

# align premirnas

# align matures mirnas

# create a file with alignments

# manually inspect alignments

# create a table with divergence value between hairpin and mirnas






# plot divergence for premirnas and for mature mirnas



# use paml to compute divergence for comparison with dN and dS?



# or can I compare Jukes-Cantor with ML estimates?