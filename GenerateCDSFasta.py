# -*- coding: utf-8 -*-
"""
Created on Wed Apr 29 13:21:35 2015

@author: Richard
"""

# use this script to generate FASTA files with a single string sequence for all scaffolds/LGs

# place this script in the same directory with input files

from manipulate_sequences import *
from get_coding_sequences import *
from divergence import *
import os

# convert genomes to fasta
remanei = convert_fasta('noamb_356_v1_4.fasta')
latens = convert_fasta('noamb_534_v1.fasta')
print('converted fasta files to dicts')

# save genomes to files
remfile = open('noamb_356_v1_4.txt', 'w')
for chromo in remanei:
    remfile.write('>' + chromo + '\n')
    remfile.write(remanei[chromo] + '\n')
remfile.close()

latfile = open('noamb_534_v1.txt', 'w')
for chromo in latens:
    latfile.write('>' + chromo + '\n')
    latfile.write(latens[chromo] + '\n')
latfile.close()
print('wrote fasta files')

# save the CDS sequences of each transcript to file
grab_CDS_sequences('356_10172014.gff3', 'noamb_356_v1_4.txt','noamb_PX356_all_CDS.fasta')
grab_CDS_sequences('534_10172014.gff3', 'noamb_534_v1.txt', 'noamb_PX534_all_CDS.fasta')
print('generated CDS fasta files')



