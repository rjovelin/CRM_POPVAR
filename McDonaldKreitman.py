# -*- coding: utf-8 -*-
"""
Created on Tue Jun  9 14:25:11 2015

@author: Richard
"""


#!/usr/bin/env python3
from manipulate_sequences import *
from scipy import stats


# define global variable with genetic code
genetic_code = {'TTT': 'F', 'CTT': 'L', 'ATT': 'I', 'GTT': 'V',
                   'TTC': 'F', 'CTC': 'L', 'ATC': 'I', 'GTC': 'V',
                   'TTA': 'L', 'CTA': 'L', 'ATA': 'I', 'GTA': 'V',
                   'TTG': 'L', 'CTG': 'L', 'ATG': 'M', 'GTG': 'V',
                   'TCT': 'S', 'CCT': 'P', 'ACT': 'T', 'GCT': 'A',
                   'TCC': 'S', 'CCC': 'P', 'ACC': 'T', 'GCC': 'A',
                   'TCA': 'S', 'CCA': 'P', 'ACA': 'T', 'GCA': 'A',
                   'TCG': 'S', 'CCG': 'P', 'ACG': 'T', 'GCG': 'A',
                   'TAT': 'Y', 'CAT': 'H', 'AAT': 'N', 'GAT': 'D',
                   'TAC': 'Y', 'CAC': 'H', 'AAC': 'N', 'GAC': 'D',
                   'TAA': '*', 'CAA': 'Q', 'AAA': 'K', 'GAA': 'E',
                   'TAG': '*', 'CAG': 'Q', 'AAG': 'K', 'GAG': 'E',
                   'TGT': 'C', 'CGT': 'R', 'AGT': 'S', 'GGT': 'G',
                   'TGC': 'C', 'CGC': 'R', 'AGC': 'S', 'GGC': 'G',
                   'TGA': '*', 'CGA': 'R', 'AGA': 'R', 'GGA': 'G',
                   'TGG': 'W', 'CGG': 'R', 'AGG': 'R', 'GGG': 'G'}


# use this function to count the number of differences between 2 codons
def diff_codon(codon1, codon2):
    '''
    (str, str) -> int
    Take 2 codons names codon1 and codon2 and return the number of nucleotide
    differences 0, 1 , 2 or 3
    '''
    
    # set up counter
    diff = 0
    
    # loop over codon1, compare each position to codon2
    for i in range(len(codon1)):
        if codon1[i] != codon2[i]:
            diff += 1
            
    return diff
        




# use this function to count the number of replacement and synonymous changes
def count_polym_diverg(snp_file, strains, rare_sites, cutoff, raw_count):
    '''
    (file, str, str, float, int) -> dict 
    Take the file with the SNPs in coding sequences, a given focal group of strains
    and return a dictionnary with gene as key and a list of containing PN, PS,
    DN, DS as value. If rare_sites is freq then polymorphic sites with a
    frequency < cutoff are ignored, and if rare_sites is count then sites that
    are < than raw_count in the sample are ignored
    '''
    
    # PN: number of replacement polymorphisms
    # PS: number of synonymous polymorphisms
    # DN: number of fixed replacements
    # DS: number of fixed synonymous changes    
    
    # create a set of stop codons
    stop_codons = {'TAG', 'TGA', 'TAA'} 
    
    
    # initialize dict
    SNPs = {}
    # initialize gene : empty list pairs
    # open file for reading
    infile = open(snp_file, 'r')
    # skip header
    infile.readline()
    # loop over file, get the gene name as key and initialize list
    # list contains PN, PS, DN, DS
    for line in infile:
        if line.rstrip() != '':
            line = line.rstrip().split()
            gene = line[2]
            if gene not in SNPs:
                # {gene : [PN, PS, DN, DS]}
                SNPs[gene] = [0, 0, 0, 0]
    # close file after reading
    infile.close()
    
    # open file for reading
    infile = open(snp_file, 'r')
    # skip header
    infile.readline()
    # loop over file
    for line in infile:
        line = line.rstrip()
        if line != '':
            line = line.split()
            # gett dict key
            gene = line[2]
            # record only sites for which ancestral state is defined
            if line[19] in {'A', 'C', 'T', 'G'}:
                # record only sites labeled snp or no_snp and valid snp type
                if line[6] in {'snp', 'no_snp'} and line[9] in {'NA', 'COD', 'SYN', 'REP'}:
                    # get reference codon
                    ref_codon = line[3]
                    # get alternative codon
                    alt_codon = line[8]
                    # get latens codon
                    cla_codon = line[18]
                    # get reference allele
                    ref = line[5]
                    # get alternative allele
                    alt = line[7]
                    # get latens base
                    cla_base = line[19]
                    # get SNP type
                    snp = line[9]
                    # do not consider stop codons
                    if ref_codon not in stop_codons and alt_codon not in stop_codons and cla_codon not in stop_codons:
                        # do not consider codons with 'Ns'
                        if 'N' not in ref_codon and 'N' not in alt_codon and 'N' not in cla_codon:
                            # check focal group, consider only KSR strains or KRS and PX combined
                            if strains == 'KSR':
                                # get reference allele count
                                ref_count = int(line[10])
                                # get alternative allele count
                                alt_count = int(line[14])
                            elif strains == 'KSR_PX':
                                # get reference allele count
                                ref_count = int(line[13])
                                # get alternative allele count
                                alt_count = int(line[17])
                            # Consider sites that have a sample size >= 10
                            if ref_count + alt_count >= 10:
                                # determine if site is polymorphic or fixed
                                # fixed diff if alternative allele fixed and different from latens 
                                if ref_count == 0 and alt_count != 0 and cla_base != alt:
                                    # fixed difference between latens and remanei
                                    # check the number of substitutions between codons
                                    if diff_codon(cla_codon, alt_codon) <= 1:
                                        # do not consider codons with more than 1 substitutions
                                        # check if change is synonymous or nonsynonymous
                                        if genetic_code[cla_codon] != genetic_code[alt_codon]:
                                            # fixed nonsynonymous change
                                            SNPs[gene][2] += 1
                                        elif genetic_code[cla_codon] == genetic_code[alt_codon]:
                                            # fixed synonymous change
                                            SNPs[gene][3] += 1
                                elif ref_count != 0 and alt_count == 0 and cla_base != ref:
                                    # fixed difference between latens and remanei
                                    # check the number of difference between latens and remanei
                                    if diff_codon(cla_codon, ref_codon) <= 1:
                                        # do not consider codons with more than 1 substitutions
                                        # check if change is synonymous or nonsynonymous
                                        if genetic_code[cla_codon] != genetic_code[ref_codon]:
                                            # fixed nonsynonymous change
                                            SNPs[gene][2] += 1
                                        elif genetic_code[cla_codon] == genetic_code[ref_codon]:
                                            # fixed synonymous change
                                            SNPs[gene][3] += 1
                                elif ref_count != 0 and alt_count != 0:
                                    # site is polymorphic
                                    # check the number of differences between ref_codon and alt_codon within remanei
                                    if diff_codon(ref_codon, alt_codon) <= 1:
                                        # do not consider codons with more than 1 substitutions
                                        # check if cutoff or raw_count applies
                                        if rare_sites == 'freq':
                                            # use frequency cutoff
                                            # compare the snp frequency to cutoff
                                            if ref_count >= alt_count:
                                                freq = alt_count / ref_count
                                            elif ref_count < alt_count:
                                                freq = ref_count / alt_count
                                            if freq >= cutoff:
                                                # determine the type of snp
                                                if snp == 'SYN' and genetic_code[ref_codon] == genetic_code[alt_codon]:
                                                    # mutation is synonymous
                                                    SNPs[gene][1] += 1
                                                elif snp == 'REP' and genetic_code[ref_codon] != genetic_code[alt_codon]:
                                                    # mutation is nonsynonymous
                                                    SNPs[gene][0] += 1
                                        elif rare_sites == 'count':
                                            # use allele count to filter sites
                                            # determine the allele with lower count
                                            if ref_count >= alt_count:
                                                # compre the count of alt_codon to raw_count
                                                if alt_count >= raw_count:
                                                    # alt_count is greater than minimum required threshold
                                                    # dtermine the type of mutation
                                                    if snp == 'SYN' and genetic_code[ref_codon] == genetic_code[alt_codon]:
                                                        # mutation is synonymous
                                                        SNPs[gene][1] += 1
                                                    elif snp == 'REP' and genetic_code[ref_codon] != genetic_code[alt_codon]:
                                                        # mutation is nonsynonymous
                                                        SNPs[gene][0] += 1
                                            elif ref_count < alt_count:
                                                # compare the count of ref-codon to raw_count
                                                if ref_count >= raw_count:
                                                    # ref_count is greater than minimum required thereshold
                                                    # determine the type of mutation
                                                    if snp == 'SYN' and genetic_code[ref_codon] == genetic_code[alt_codon]:
                                                        # mutation is synonymous
                                                        SNPs[gene][1] += 1
                                                    elif snp == 'REP' and genetic_code[ref_codon] != genetic_code[alt_codon]:
                                                        # mutation is nonsynonymous
                                                        SNPs[gene][0] += 1
                                                    
    # close file after readling
    infile.close()
    
    # remove genes with no polymorhisms or no divergence
    # because MK test cannot be computed for these genes
    
    to_remove = []
    for gene in SNPs:
        if SNPs[gene][0] == 0 and SNPs[gene][1] == 0:
            to_remove.append(gene)
        elif SNPs[gene][2] == 0 and SNPs[gene][3] == 0:
            to_remove.append(gene)
        
    for gene in to_remove:
      del SNPs[gene]  
    
    # return dict
    return SNPs
    
    

def MK_test(SNPs, test_mode):
    '''
    (dict, str) -> dict
    Take a dict of gene : [PN, PS, DN, DS] pairs and a string fisher or G_test
    and a return  a new dict with gene : [PN, PS, DN, DS, p-val] pairs 
    with PN and DN being respectively replacement polymorphisms and divergence
    and PS and DS being respectively synonymous polymorphisms and divergence 
    and p-val being the p-value of the contingency test using either Fisher's
    two-sided exact test or the G-test with Yate's correction
    '''
    
    # create new dict
    MK = {}    
    
    # loop over genes in dict
    for gene in SNPs:
        # initialize list with PN, PS
        polym = [SNPs[gene][0], SNPs[gene][1]]
        # initialize list with DN, DS
        diverg = [SNPs[gene][2], SNPs[gene][3]]
        # perform the MK test according to fisher 2-tailed or G-test
        if test_mode == 'fisher':
            # get the p-value
            P = stats.fisher_exact([polym, diverg])[1]
        elif test_mode == 'G_test':
            P = stats.chi2_contingency([polym, diverg], lambda_ = 'log-likelihood')[1]
        # add p-val to list
        MK[gene] = list(SNPs[gene])
        MK[gene].append(P)
        
    # return modified dict
    return MK
            
    
    