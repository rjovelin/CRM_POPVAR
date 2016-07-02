# -*- coding: utf-8 -*-
"""
Created on Tue Jun  9 14:25:11 2015

@author: Richard
"""


from manipulate_sequences import *
from scipy import stats
import numpy as np
import random


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
            # get dict key
            gene = line[2]
            # record only sites for which ancestral state is defined
            if line[19] in {'A', 'C', 'T', 'G'}:
                # record only sites labeled snp or no_snp and valid snp type
                if line[6] in {'snp', 'no_snp'} and line[9] in {'NA', 'COD', 'SYN', 'REP'}:
                    # get reference and alternative codons and latens codon
                    ref_codon, alt_codon, cla_codon = line[3], line[8], line[18]
                    # get reference, alternative and latens alleles
                    ref, alt, cla_base = line[5], line[7], line[19]
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
                                # do not consider codons with more than 1 substitutions
                                if diff_codon(cla_codon, alt_codon) <= 1:
                                    # determine if site is polymorphic or fixed
                                    # fixed diff if alternative allele fixed and different from latens 
                                    if ref_count == 0 and alt_count != 0 and cla_base != alt:
                                        # fixed difference between latens and remanei
                                        # check if change is synonymous or nonsynonymous
                                        if genetic_code[cla_codon] != genetic_code[alt_codon]:
                                            # fixed nonsynonymous change
                                            SNPs[gene][2] += 1
                                        elif genetic_code[cla_codon] == genetic_code[alt_codon]:
                                            # fixed synonymous change
                                            SNPs[gene][3] += 1
                                    elif ref_count != 0 and alt_count == 0 and cla_base != ref:
                                        # fixed difference between latens and remanei
                                        # do not consider codons with more than 1 substitutions
                                        # check if change is synonymous or nonsynonymous
                                        if genetic_code[cla_codon] != genetic_code[ref_codon]:
                                            # fixed nonsynonymous change
                                            SNPs[gene][2] += 1
                                        elif genetic_code[cla_codon] == genetic_code[ref_codon]:
                                            # fixed synonymous change
                                            SNPs[gene][3] += 1
                                    elif ref_count != 0 and alt_count != 0 and (ref == cla_base or alt == cla_base):
                                        # site is polymorphic
                                        # set up boolean to identify polymorphic site after filtering based on MAF or raw count
                                        PolymorphicSite = False
                                        # check if cutoff or raw_count applies
                                        if rare_sites == 'freq':
                                            # use frequency cutoff
                                            # compare the snp MAF frequency to cutoff
                                            if ref_count >= alt_count:
                                                freq = alt_count / (ref_count + alt_count)
                                            elif ref_count < alt_count:
                                                freq = ref_count / (ref_count + alt_count)
                                            if freq >= cutoff:
                                                # record polymorphic site
                                                PolymorphicSite = True
                                        elif rare_sites == 'count':
                                            # use allele count to filter sites
                                            # check that allele with lowest count > raw_count threshold
                                            if ref_count >= alt_count and alt_count > raw_count:
                                                # alt_count is greater than minimum required threshold
                                                # record polymorphic site
                                                PolymorphicSite = True
                                            elif ref_count < alt_count and ref_count > raw_count:
                                                # ref_count is greater than minimum required thereshold
                                                # record polymorphic site
                                                PolymorphicSite = True
                                        # determine the type of mutation if polymorphic site is to be recorded
                                        if PolymorphicSite == True:
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
    

# use this function to compute the MK test     
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
        
    return MK
            
    
# use this function determine genes significantly under positive or negative selection
def NaturalSelection(PolymDivCounts, Significant):
    '''
    (dict, list) -> (list, list)
    Take the dictionary with polymorphism and divergence counts and a list of
    genes with significant MK test and return a tuple with lists of genes
    that are respectively under negative and positive selection
    '''
    negative, positive = [], []
    for gene in PolymDivCounts:
        # check if gene has significant departure from neutrality
        if gene in Significant:
            # check if PS is defined
            if PolymDivCounts[gene][1] != 0:
                # compute ratio PN / PS
                polymorphism = PolymDivCounts[gene][0] / PolymDivCounts[gene][1]
            elif PolymDivCounts[gene][1] == 0:
                # add small number to PS
                polymorphism = PolymDivCounts[gene][0] / (PolymDivCounts[gene][1] + 0.1)
            # check if DS is defined
            if PolymDivCounts[gene][3] != 0:
                # compute DN / DS
                divergence = PolymDivCounts[gene][2] / PolymDivCounts[gene][3]
            elif PolymDivCounts[gene][3] == 0:
                # add small number to DS
                divergence = PolymDivCounts[gene][2] / (PolymDivCounts[gene][3] + 0.1)
            # check if divergence is lower r higher than polymorphism
            if divergence > polymorphism:
                # positive selection
                positive.append(gene)
            elif divergence < polymorphism:
                negative.append(gene)
    return negative, positive    
    
    
    
# use this function to compute alpha according to Smith-EyreWalker 2002    
def ComputeAlphaSEW2002(PolymDivCounts, MinimumPS):
     '''
     (dict, int) -> float
     Take a dictionary with polymorphim and divergence counts at synonymous 
     and replacement sites and the minimum number of synonymous polymorphisms
     and compute alpha, the average proportion of amino-acid substitutions
     driven by positive selection according to the method of Smith-Eyre-Walker 2002
     '''
     
     # alpha = 1 - ((MeanDS / MeanDN) * (Mean(PN / (PS + 1))))
     # MeanDS: average number of fixed differences at synonymous sites
     # MeanDN: average number of of fixed differences at nonsynonymous sites
     # PN: number of nonymous replacement polymorphisms for a given gene
     # PS: number of synonymous polymorphisms or a given gene
     # PolymDivCounts  if in the form {gene : [PN, PS, DN, DS]}
     
     # create lists with divergence counts
     DN = [PolymDivCounts[gene][2] for gene in PolymDivCounts]
     DS = [PolymDivCounts[gene][3] for gene in PolymDivCounts]
     
     # create list with polymorphism ratio
     Polym = [(PolymDivCounts[gene][0] / (PolymDivCounts[gene][1] + 1)) for gene in PolymDivCounts if PolymDivCounts[gene][1] >= MinimumPS]
     
     alpha = 1 - ((np.mean(DS) / np.mean(DN)) * np.mean(Polym))
     return alpha     
     

# use this function to bootstrap genes to compoute a distribution of alpha values
def BootstrapAlphaSEW2002(PolymDivCounts, MinimumPS, replicates, Ngenes):
    '''
    (dict, int, int, int) -> list
    Take a dictionary with polymorphim and divergence counts at synonymous 
    and replacement sites, the minimum number of synonymous polymorphisms,
    the number of bootstrap replicates, the number of genes to draw with replacement,
    and return a list with distribution of alpha (Smith-Eyre-Walker 2002)
    for each replicate
    '''
    
    # create a list to store alpha values
    AlphaDistribution = []
    
    # create a list of genes to draw genes at random    
    GeneNames = [gene for gene in PolymDivCounts]
        
    # loop over number of replicates
    while replicates != 0:
        # create a dictionary with {gene : [PN, PS, DN, DS]} from genes draw at random
        RandomDraw = {}
        i = Ngenes
        # draw Ngenes from PolymDivCounts
        while i != 0:
            # draw the index of gene at random (randint include last number)
            # draw genes with replacement
            position = random.randint(0, len(GeneNames)-1)
            gene = GeneNames[position]
            # populate dict
            RandomDraw[gene] = list(PolymDivCounts[gene])
            # update counter
            i -= 1
        # compute alpha
        alpha = ComputeAlphaSEW2002(RandomDraw, MinimumPS)
        AlphaDistribution.append(alpha)
        # update counter
        replicates -= 1
        
    return AlphaDistribution


# use this function to count the number of fixed diffs and polymorphisms at 4-fold degenerate sites
def CountPolymDivergFourFold(snp_file, rare_alleles, threshold):
    '''
    (file, bool, int) -> dict 
    Take the file with the SNPs in coding sequences, whether rare alleles 
    should be ignored or not, and the number of mutations to ignore per site      
    (ie, singletons) and return a dict with D and P the number of fixed 
    differences and polymorphisms at 4-fold degenerate sites for each gene
    '''
    
    # create a dict with 4 fold degenerate codons
    fourfold = {'S': ['TCA', 'TCT', 'TCG', 'TCC'] , 'L': ['CTA', 'CTT', 'CTG', 'CTC'],
                'P': ['CCA', 'CCT', 'CCG', 'CCC'], 'R': ['CGA', 'CGT', 'CGG', 'CGC'],
                'T': ['ACA', 'ACT', 'ACG', 'ACC'], 'A': ['GCA', 'GCT', 'GCG', 'GCC'],
                'G': ['GGA', 'GGT', 'GGG', 'GGC'], 'V': ['GTA', 'GTT', 'GTG', 'GTC']} 

    # create a dict to record D and P {gene: [D, P]}
    # D: number of fixed replacements at old sites
    # P: number of polymorphism at 4-fold sites    
    
    # initialize dict {gene: [0, 0]}
    SNPs = {}
    infile = open(snp_file, 'r')
    infile.readline()
    # loop over file, get the gene name as key and initialize list
    for line in infile:
        if line.rstrip() != '':
            line = line.rstrip().split()
            gene = line[2]
            if gene not in SNPs:
                SNPs[gene] = [0, 0]
    # close file after reading
    infile.close()
    
    # count D and P
    infile = open(snp_file, 'r')
    infile.readline()
    for line in infile:
        line = line.rstrip()
        if line != '':
            line = line.split()
            # get dict key
            gene = line[2]
            # check that ancestral state is defined, and that sites are labeled snp or no snp and valid snp
            if line[19] in {'A', 'C', 'T', 'G'} and line[6] in {'snp', 'no_snp'} and line[9] in {'NA', 'COD', 'SYN', 'REP'}:
                # get reference and alternative codons and latens codon
                ref_codon, alt_codon, cla_codon = line[3], line[8], line[18]
                # get reference, alternative and latens alleles
                ref, alt, cla_base = line[5], line[7], line[19]
                # get ref and alt counts for KSR + PX strains
                ref_count, alt_count = int(line[13]), int(line[17])
                # check that ref codon, alt codon and cla-codon translate to AAs with 4-fold degenerate sites
                if cds_translate(ref_codon) in fourfold and cds_translate(alt_codon) in fourfold and cds_translate(cla_codon) in fourfold:
                    # check that mutation is synonymous
                    if cds_translate(ref_codon) == cds_translate(alt_codon) and cds_translate(ref_codon) == cds_translate(cla_codon) and cds_translate(alt_codon) == cds_translate(cla_codon):
                        # check that codons correspond to the 4-fold degenerate codons
                        if ref_codon in fourfold[cds_translate(ref_codon)] and alt_codon in fourfold[cds_translate(alt_codon)] and cla_codon in fourfold[cds_translate(cla_codon)]:
                            # consider only sites with sample size >= 10
                            if ref_count + alt_count >= 10:
                                # determine if site is fixed or polymorphic
                                if ref_count == 0 and alt_count != 0 and cla_base != alt:
                                    # fixed difference between latens and remanei
                                    SNPs[gene][0] += 1
                                elif ref_count != 0 and alt_count == 0 and cla_base != ref:
                                    # fixed difference between latens and remanei
                                    SNPs[gene][0] += 1
                                elif ref_count != 0 and alt_count != 0 and (ref == cla_base or alt == cla_base):
                                    # site is polymorphic
                                    # check if some polymorphic sites need to be ignored
                                    if rare_alleles == True:
                                        # use allele count to filter sites
                                        # check that allele with lowest count > raw_count threshold
                                        if min(ref_count, alt_count) > threshold:
                                            # record polymorphic site
                                            SNPs[gene][1] += 1
                                    elif rare_alleles == False:
                                        # record polymorphic sites
                                        SNPs[gene][1] += 1
    # close file after readling
    infile.close()
    
    # remove genes with no polymorhisms and no divergence
    to_remove = [i for i in SNPs if SNPs[i] == [0, 0]]
    for i in to_remove:
      del SNPs[i]  
    # return dict
    return SNPs
    
    
    
 # use this function to count the number of fixed differences and polymorphisms in miRNAs   
def CountPolymDivergmiRNAs(hairpins, hairpin_coord, CrmGenome, chromo_sites, rare_alleles, threshold):
    '''
    (dict, dict, dict, dict, bool, int)
    Take a dict with aligned sequences between remanei and latens mirna orthologs,
    a dict with mirna coordinates, a dict with genome sequence, a dict with SNP data    
    a boolean indicating whether rare alleles should be ignored or not, and the
    number of mutations to ignore per site (eg, singletons) and return a dict
    with D and P the number of fixed differences and polymorphisms in each mirna
    '''

    # hairpins is a dict {crm_mirna_name: {crm_mirna_name: aligned_seq, cla_mirna_name: aligned_seq}}
    # hairpin_coord is a dict with mirna coordinates {name: [chromo, start, end, orientation]}
    # CrmGenome is a dict with genome {chromo: sequence}
    # chromo-sites is a dict with SNP data {chromo: {site : [ref_allele, alt_allele, count_ref, count alt]}]}


    # create a dict to store the fixed diffs and polyms {mirna: [D, P]}
    hairpin_diffs = {}

    # loop over aligned hairpins
    for mirna in hairpins:
        # get chromo, start, end and orientation
        chromo, start, end, orientation = hairpin_coord[mirna][0], hairpin_coord[mirna][1], hairpin_coord[mirna][2], hairpin_coord[mirna][3] 
        # extract sequence from genome
        sequence = CrmGenome[chromo][start: end]
        if orientation == '-':
            sequence = reverse_complement(sequence)
        # get mirna positions on chromo
        positions = [i for i in range(start, end)]    
        # get the position in decreasing order if orientation is -    
        if orientation == '-':
            positions.reverse()
        # get the remanei and latens mirna name, and correspsonding sequences
        for name in hairpins[mirna]:
            if name.startswith('crm'):
                crmmirna = name
                assert crmirna == mirna, 'mirna names do not match'
                crmseq = hairpins[mirna][crmmirna]
            elif name.startswith('cla'):
                clamirna = name
                claseq = hairpins[mirna][clamirna]
        # set up a gap counter
        gaps = 0
        # loop over crm sequence
        # keep track of index to look for SNP data when gaps are present
        for i in range(len(crmseq)):
            # get ancestral allele and remanei allele
            ancestral, crmallele = claseq[i], crmseq[i]
            if crmallele != '-':
                # get the index to look in positions
                j = i - gaps
                # check that nucleotide in crmseq correspond to sequence from genome
                assert sequence[j] == crmallele, 'nucleotides do not match between extracted sequence and aligned sequence'
                # check that index in positions is correct
                # check that ref allele in SNP dictionary is correct
                if orientation == '+':
                    assert CrmGenome[chromo][positions[j]] == crmallele, 'no match with nucleotide extracted with list index on +'
                    if positions[j] in chromo_sites[chromo]:
                        assert chromo_sites[chromo][positions[j]][0] == crmallele, 'no match with ref allele in SNP dict in +'
                    else:
                        print(j, positions[j], orientation)
                elif orientation == '-':
                    assert seq_complement(CrmGenome[chromo][positions[j]]) == crmallele, 'no match with nucleotide extracted with list index on -'
                    if positions[j] in chromo_sites[chromo]:
                        assert seq_complement(chromo_sites[chromo][positions[j]][0]) == crmallele, 'no match with ref allele in SNP dict in -'
                    else:
                        print(j, positions[j], orientation)
            else:
                gaps += 1
            # determine if site a fized difference or a polymorphism
            # check that ancestral allele is valid base
            if ancestral in 'ATCG' and crmallele in 'ATCG':
                # check that site has coverage            
                if positions[j] in chromo_sites[chromo]:
                    # fixed diff if ref_count != 0 and alt_count = 0 and ref != ancestral 
                    # fixed diff is ref_count = 0 and alt_count != 0 and alt != ancestral
                    # polymorphism if (ref_count != 0 and alt_count != 0) and (ref = amcestral or alt = ancestral)
                    # get ref and alt counts
                    ref_count, alt_count = chromo_sites[chromo][positions[j]][2], chromo_sites[chromo][positions[j]][3]
                    # get reference and alternative alleles
                    ref, alt = chromo_sites[chromo][positions[j]][0], chromo_sites[chromo][positions[j]][1]
                    # consider positions with sample size > 10 (positions are already filtered in chromo_sites)
                    if ref_count + alt_count >= 10:
                        if ref_count != 0 and alt_count == 0 and ref != ancestral:
                            # fixed difference, populate dict
                            if mirna in hairpin_diffs:
                                hairpin_diffs[mirna][0] += 1
                            else:
                                hairpin_diffs[mirna] = [1, 0]
                        elif ref_count == 0 and alt_count != 0 and alt != ancestral:
                            # fixed difference, populate dict
                            if mirna in hairpin_diffs:
                                hairpin_diffs[mirna][0] += 1
                            else:
                                hairpin_diffs[mirna] = [1, 0]
                        elif (ref_count != 0 and alt_count != 0) and (ref == ancestral or alt == ancestral):
                            # check if some polymorphic sites need to be ignored
                            if rare_alleles == True:
                                # use allele count to filter sites
                                # check that allele with lowest count > raw_count threshold
                                if min(ref_count, alt_count) > threshold:
                                    # polymorphism, populate dict
                                    if mirna in hairpin_diffs:
                                        hairpin_diffs[mirna][1] += 1
                                    else:
                                        hairpin_diffs[mirna] = [0, 1]
   
    return hairpin_diffs
    
