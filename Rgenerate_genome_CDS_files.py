# -*- coding: utf-8 -*-
"""
Created on Wed Apr 29 13:21:35 2015

@author: Richard
"""

from accessories import *
from get_coding_sequences import *
from crem_cla_divergence import *
import os


# convert genomes to fasta
remanei = convert_fasta('../../noamb_356_v1_4.fasta')
latens = convert_fasta('../../noamb_534_v1.fasta')
print('fasta_conversion done')

# save genomes to files
remfile = open('../Utilitary_files/noamb_356_v1_4.txt', 'w')
for chromo in remanei:
    remfile.write('>' + chromo + '\n')
    remfile.write(remanei[chromo] + '\n')
remfile.close()

latfile = open('../Utilitary_files/noamb_534_v1.txt', 'w')
for chromo in latens:
    latfile.write('>' + chromo + '\n')
    latfile.write(latens[chromo] + '\n')
latfile.close()

print('fasta to file done')


# save the CDS sequences of each transcript to file
grab_CDS_sequences('../Utilitary_files/356_10172014.gff3', '../Utilitary_files/noamb_356_v1_4.txt','../Utilitary_files/noamb_PX356_all_CDS.fasta')
grab_CDS_sequences('../Utilitary_files/534_10172014.gff3', '../Utilitary_files/noamb_534_v1.txt', '../Utilitary_files/noamb_PX534_all_CDS.fasta')

print('CDS fasta files done')



