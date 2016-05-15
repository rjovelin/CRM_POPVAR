# -*- coding: utf-8 -*-
"""
Created on Sun May 15 12:17:34 2016

@author: Richard
"""

# use this script to get the SNPs at noncoding sites in ONtario strains(KSR + PX),
# removing positions with sample size < 10 and save data to json format for fast retrieval



import json
from sites_with_coverage import *


# create a dictionary with all sites with coverage in the genome
# consider only site with a minimum sample size of 10
chromo_sites = get_non_coding_snps('../SNP_files/', 10)
print('got sites with coverage')
print(len(chromo_sites))

# save data as json file
datafile = open('NonCodingSNPsOntario.json', 'w')
json.dump(chromo_sites, datafile, sort_keys=True, indent=4)
datafile.close()

