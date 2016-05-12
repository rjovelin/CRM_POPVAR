# -*- coding: utf-8 -*-
"""
Created on Thu Sep 10 14:04:22 2015

@author: RJovelin
"""

import os
from parse_vcf_comp_files import *
from manipulate_sequences import *
from repeats_TEs import *


# use this script to generate SNP files that includes Ontario strains (PX+KSR) and PB strains
# use these files to compute diversity at non-coding sites

# path to the vcf comp files
vcf_dir = '/disk3/wei.wang/c_rem_analysis/split_strain_vcf_ran/vcf/'

# create a dictionnary to store the SNP files
os.mkdir('../PB_ON_SNP_files/')

# path to the snp_files
snp_dir = '../PB_ON_SNP_files/'

# create a dict of repeat name coordinate
repeat_coord = get_repeats_coord('356_v1_4.fasta.out', False)
# create a dict of chromo : set of repeat indices
repeat_pos = get_repeat_positions(repeat_coord)

# make a list of filename
files = [filename for filename in os.listdir(vcf_dir) if 'all' in filename and 'filtered.com' in filename]
    
# loop over filename
for i in files:
    print(i)
    # get the full path of the file
    filename = vcf_dir + i
    # get chromo and dict with SNP counts
    # remove SNPs in regions with repeats
    chromo, snps = parse_vcf_comp_PB_ON(filename, repeat_pos, True)
    # check that chromo and snps are defined (ie. that comp file is not empty)
    if chromo != '' and len(snps) != 0:
        # get name of outputfile
        outputfile = snp_dir + 'PBON_' + chromo + '_SNPs.txt'        
        # save SNP info to outputfile
        generate_PB_ON_snp_file(chromo, snps, outputfile)   
    else:
        print(i, 'empty')





