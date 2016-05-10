# -*- coding: utf-8 -*-
"""
Created on Thu May 28 21:23:49 2015

@author: Richard
"""

from manipulate_sequences import *
from miRNA_target import *



# use this script to generate a fasta file with cremanei mature sequences
# and the targetscan miRNA family input file


# generate a fasta file with mature mirna sequences
infile = open('CRM_miRNAsCoordinatesFinal.txt')
infile.readline()
mature = {}
for line in infile:
    line = line.rstrip()
    if line != '':
        line = line.split('\t')
        # replace T by U
        mature[line[0]] = line[6].upper().replace('T', 'U')
infile.close()

newfile = open('Cremanei_mature.fasta', 'w')
for mir in mature:
    newfile.write('>' + mir + '\n')
    newfile.write(mature[mir] + '\n')
newfile.close()

# generate the targetscan mirna family imput file
# create a dict of seed sequence and list of mirna pairs
seeds = seed_mirnas('Cremanei_mature.fasta')

# open file for writing
newfile = open('Cremanei_miRfam_info.txt', 'w')
# loop over the seed sequences
for motif in seeds:
    # record only 1 mirna per family, grab first mirna
    # write mirna, seed, species ID
    newfile.write(seeds[motif][0] + '\t' + motif + '\t' + '31234' + '\n')
# close file after writing
newfile.close()

