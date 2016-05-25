# -*- coding: utf-8 -*-
"""
Created on Wed May 25 11:05:02 2016

@author: RJovelin
"""

# use this script to save the miRNA targets coordinates to json file

# import modules
import json
from miRNA_target import *

# get the coordinates of all target sites
target_coord = get_miRNA_target_loci('Cremanei_miRNA_sites.txt', '../Genome_Files/unique_transcripts.txt', 'all')
print('got miRNA target site coordinates')

newfile = open('AllmiRNATargetsCoords.json', 'w')
json.dump(target_coord, newfile, sort_keys=True, indent=4)
newfile.close()
