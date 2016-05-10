# -*- coding: utf-8 -*-
"""
Created on Sat Jul 11 09:10:53 2015

@author: Richard
"""

from piRNAs import *
from accessories import *

pirna_coord = map_remanei_pirnas('../CREM_CLA_protein_divergence/noamb_356_v1_4.txt', 'remanei_piRNA_upstream70.fa')

write_pirna_coordinates_to_file(pirna_coord, 'PX356_piRNA_coord.txt')