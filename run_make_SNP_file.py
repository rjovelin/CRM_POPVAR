# -*- coding: utf-8 -*-
"""
Created on Thu Jul 16 00:16:27 2015

@author: Richard
"""

import os
from parse_vcf_comp_files import *
from manipulate_sequences import *
from repeats_TEs import *

# path to the vcf comp files
vcf_dir = '/disk3/wei.wang/c_rem_analysis/split_strain_vcf_combKSRPX/vcf/'

# create a dictionnary to store the SNP files
os.mkdir('../SNP_files/')

# path to the snp_files
snp_dir = '../SNP_files/'

# create a dict of repeat name coordinate
repeat_coord = get_repeats_coord('../piRNAs/356_v1_4.fasta.out', False)
# create a dict of chromo : set of repeat indices
repeat_pos = get_repeat_positions(repeat_coord)


# check that all chromo in the GFF file are in the vcf comp file
chromo_gff, chromo_dir = check_chromo(vcf_dir, '356_10172014.gff3')
print(len(chromo_gff), len(chromo_dir))

if chromo_gff == chromo_dir:
    print('same number of chromos')
    # proceed to generate the summary SNP file
    
    # make a list of filename
    files = [filename for filename in os.listdir(vcf_dir) if 'ksrpx' in filename and 'filtered.com' in filename]
    
    # loop over filename
    for i in files:
        print(i)
        # get the full path of the file
        filename = vcf_dir + i
        # get chromo and dict with SNP counts
        chromo, snps = from_vcf_comp_to_dict(filename, repeat_pos, True)
        # check that chromo and snps are defined (ie. that comp file is not empty)
        if chromo != '' and snps != 0:
            # get name of outputfile
            outputfile = snp_dir + 'PXKSR_' + chromo + 'SNPs.txt'        
            # save SNP info to outputfile
            generate_snp_file(chromo, snps, outputfile)
    
elif chromo_gff != chromo_dir:
    print('not the same number of chromos, maybe something is wrong')





