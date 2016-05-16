# -*- coding: utf-8 -*-
"""
Created on Mon May 16 13:19:28 2016

@author: RJovelin
"""

# use this script to generate files with 7bp seeds for different species

from manipulate_sequences import *


# convert fasta files to dict
CBNmirna = convert_fasta('CBN_mirbase_miRs.txt')
CBRmirna = convert_fasta('CBR_mirbase_miRs.txt')
CELmirna = convert_fasta('CEL_mirbase_miRs.txt')

# open file with list of mature miRs in latens
CLA = set()
infile = open('CLA_new_miRs.txt')
for line in infile:
    line = line.rstrip()
    if line != '':
        seed = line.upper()[1:8]
        CLA.add(seed)
infile.close()

# create sets to store seeds for the other species
CBN, CBR, CEL = set(), set(), set()
# loop over mature sequences, convert U to T, and upper case and extract seed to sets
for mirna in CBNmirna:
    seed = CBNmirna[mirna].upper().replace('U', 'T')[1:8]
    CBN.add(seed)
for mirna in CBRmirna:
    seed = CBRmirna[mirna].upper().replace('U', 'T')[1:8]
    CBR.add(seed)
for mirna in CELmirna:
    seed = CELmirna[mirna].upper().replace('U', 'T')[1:8]
    CEL.add(seed)
    
# save seeds to files
newfile = open('CLA_seeds.txt', 'w')
for seed in CLA:
    newfile.write(seed + '\n')
newfile.close()

newfile = open('CBN_seeds.txt', 'w')
for seed in CBN:
    newfile.write(seed + '\n')
newfile.close()

newfile = open('CBR_seeds.txt', 'w')
for seed in CBR:
    newfile.write(seed + '\n')
newfile.close()

newfile = open('CEL_seeds.txt', 'w')
for seed in CEL:
    newfile.write(seed + '\n')
newfile.close()
 
    