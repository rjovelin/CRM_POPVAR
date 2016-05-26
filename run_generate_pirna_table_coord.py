# -*- coding: utf-8 -*-
"""
Created on Sat Jul 11 09:10:53 2015

@author: Richard
"""

# use this script to generate a file with the 


from piRNAs import *
from manipulate_sequences import *

# map the piRNAs and their upstream regions to the remanei genome 
pirna_coord = map_remanei_pirnas('../Genome_Files/noamb_356_v1_4.txt', 'remanei_piRNA_upstream70.fa')

# write piRNAs coordinates to file
write_pirna_coordinates_to_file(pirna_coord, 'PX356_piRNA_coord.txt')