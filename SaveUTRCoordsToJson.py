# -*- coding: utf-8 -*-
"""
Created on Wed May 25 00:40:50 2016

@author: Richard
"""

# use this script to save the coordinates of remanei UTR sequences as a json file


import json
from genomic_coordinates import *
from manipulate_sequences import *
from miRNA_target import *

# compute threshold based on the distribution of elegans UTR length
UTR_length = celegans_three_prime_UTR_length('../Genome_Files/c_elegans.PRJNA13758.WS248.annotations.gff3')
threshold = get_percentile(UTR_length, 99)
# get UTR coord {TS1 : [chromo, start, end, orientation]}
three_prime = get_three_prime_UTR_positions('../Genome_Files/356_10172014.gff3', '../Genome_Files/noamb_356_v1_4.txt', threshold)
print('got UTR coordinates')

# save data as json file
datafile = open('CremUTRCoordsNo.json', 'w')
json.dump(three_prime, datafile, sort_keys=True, indent=4)
datafile.close()
