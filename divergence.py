# -*- coding: utf-8 -*-
"""
Created on Fri Jun 26 17:16:38 2015

@author: Richard
"""


from mk_test import *
from genomic_coordinates import *
from manipulate_sequences import *
from sites_with_coverage import *
import numpy as np
import math




# use this function to compute theta when all sites have the same sample size
def watterson_theta(k, N):
    '''
    (int, int) -> float
    Take the number of polymorphic sites and the sample size and return
    Watterson's theta estimator of polymorphism (per sequence)
    Precondition: all sites have the same sample size
    '''
    
    # compute a (1/i from i = 1 to i = N-1)
    a = 0
    for i in range(1, N):
        a += 1/i
    # compute theta
    theta = k/a
    
    return theta
    

# use this function to count the number of SNPs and sites for 
# SNPs in coding sequence, irrespective of changes in protein sequence
def count_coding_sites_for_diversity(snp_file, threshold):
    '''
    (dict, int) -> (dict, dict)
    Take the file with SNPs in the CDS, a threshold indicating the minimum
    acceptable sample size and return a tuple with one dict with gene as key
    and and list with the number of of SNPs and the number of sites, 
    and one dict with gene as key and a list of sample size for polymorphic sites
    '''
           
    # create a set of stop codons
    stop_codons = {'TAG', 'TGA', 'TAA'} 
    
    # create dict {gene: [snps, sites]}
    SNPs = {}
 
    # create a dict {gene : [sample size for polymorphic sites]}
    # record the sample size of the polymorphic sites
    sample_size = {}

    # set up gene variable to check which gene is being recorded
    gene = ''    
    
    # open file for reading
    infile = open(snp_file, 'r')
    # skip header
    infile.readline()
    # loop over file
    for line in infile:
        line = line.rstrip()
        if line != '':
            line = line.split()
            # check which gene is being read
            if line[2] != gene:
                # new gene, assign newgene to variable, and initialize dicts
                gene = line[2]
                SNPs[gene] = []
                sample_size[gene] = []
            # record only sites labeled snp or no_snp and valid snp type
            if line[6] in {'snp', 'no_snp'} and line[9] in {'NA', 'COD', 'SYN', 'REP'}:
                # get reference codon
                ref_codon = line[3]
                # get alternative codon
                alt_codon = line[8]
                # get reference allele
                ref = line[5]
                # get alternative allele
                alt = line[7]
                # get allele counts
                ref_count = int(line[13])
                alt_count = int(line[17])
                # check if sample size is greater than threshold
                if ref_count + alt_count >= threshold:
                    # sample size is large enough, count site and check if polymorphic
                    # do not consider stop codons
                    if ref_codon not in stop_codons and alt_codon not in stop_codons:
                        # do not consider codons with 'Ns'
                        if 'N' not in ref_codon and 'N' not in alt_codon:
                            # check the number of differences between ref_codon and alt_codon within remanei
                            if diff_codon(ref_codon, alt_codon) <= 1:
                                # do not consider codons with more than 1 substitutions
                                # update site counter
                                SNPs[gene][1] += 1
                                # check if site is polymorphic
                                if ref_count != 0 and alt_count != 0:
                                    # site is polymorphic
                                    # check that ref alleles and alt alleles are different
                                    assert ref != alt, 'ref and alt counts > 0 but alleles are the same'
                                    # update SNP counter
                                    SNPs[gene][0] += 1
                                    # record the sample size of the polymorphic site
                                    sample_size[gene].append(ref_count + alt_count)
                                    
                        
    # close file after readling
    infile.close()
    
    # return dict
    return SNPs, sample_size


# use this function to count the number of SNPs and 0-fold, 2-fold and 4-fold degenerate sites    
def count_rep_synonymous_sites_for_diversity(snp_file, sites, threshold):
    '''
    (file, str, int) -> (dict, dict)
    Take the file with SNPs in the CDS, the type of changes to record 
    (REP or SYN), a threshold indicating the mininum acceptable sample size
    and return a tuple with one dict with gene as key a list with the number of SNPs, 
    0-fold, 2-fold and 4-fold degenerate sites in that gene, and one dict with
    gene as key and a list of sample size for polymorphic sites
    '''
    
    # create dict {gene: [snp, 0-fold, 2-fold, 4-fold]}
    SNPs = {}
    
    # create a dict to store the sample size of polymorphic sites
    sample_size = {}
    
    # create a set of stop codons
    stop_codons = {'TAG', 'TGA', 'TAA'} 
       
    # set up gene variable to check which gene is being recorded
    gene = ''    
    
    # open file for reading
    infile = open(snp_file, 'r')
    # skip header
    infile.readline()
    # loop over file
    for line in infile:
        line = line.rstrip()
        if line != '':
            line = line.split()
            # check which gene is being read
            if line[2] != gene:
                # new gene, assign newgene to variable, and initialize dict
                gene = line[2]
                SNPs[gene] = [0, 0, 0, 0]
                sample_size[gene] = []
                
            # record only sites labeled snp or no_snp and valid snp type
            if line[6] in {'snp', 'no_snp'} and line[9] in {'NA', 'COD', 'SYN', 'REP'}:
                # get reference codon
                ref_codon = line[3]
                # get alternative codon
                alt_codon = line[8]
                # get reference allele
                ref = line[5]
                # get alternative allele
                alt = line[7]
                # get allele counts
                ref_count = int(line[13])
                alt_count = int(line[17])
                # check is sample size is greater than threshold
                if ref_count + alt_count >= threshold:
                    # sample size is large enough, record sites and check polymorphism
                    # do not consider stop codons
                    if ref_codon not in stop_codons and alt_codon not in stop_codons:
                        # do not consider codons with 'Ns'
                        if 'N' not in ref_codon and 'N' not in alt_codon:
                            # check the number of differences between ref_codon and alt_codon within remanei
                            if diff_codon(ref_codon, alt_codon) <= 1:
                                # do not consider codons with more than 1 substitutions
                                # get site degeneracy
                                degeneracy = count_degenerate_sites_in_codon(ref_codon, line[4])
                                # update the degenerate site counter
                                for i in range(len(degeneracy)):
                                    # add the value to indices 1 to 3 in list
                                    SNPs[gene][i+1] += degeneracy[i]
                                # check if site is polymorphic
                                if ref_count != 0 and alt_count != 0:
                                    # site is polymorphic
                                    # check that ref alleles and alt alleles are different
                                    assert ref != alt, 'ref and alt counts > 0 but alleles are the same'
                                    # determine the type of snp
                                    if sites == 'SYN' and genetic_code[ref_codon] == genetic_code[alt_codon]:
                                        # mutation is synonymous
                                        # update snp counter
                                        SNPs[gene][0] += 1
                                        # record sample size of polymorphic site
                                        sample_size[gene].append(ref_count + alt_count)
                                    elif sites == 'REP' and genetic_code[ref_codon] != genetic_code[alt_codon]:
                                        # mutation is nonsynonymous
                                        # update snp counter
                                        SNPs[gene][0] += 1
                                        # record sample size of polymorphic site
                                        sample_size[gene].append(ref_count + alt_count)
                                        
                                               
    # close file
    infile.close()
    
    # return dict
    return SNPs, sample_size


# use this function to count the number of synonymous and nonsynonymous sites in the reference
# based on the number of 0-fold, 2-fold and 4-fold degenerate sites
def count_nonsynonymous_synonymous_sites(SNPs):
    '''
    (dict) -> dict
    Take a dict with genes as key and a list with the number of SNPs, 0-fold,
    2-fold and 4-fold degenerate sites as value, and return a new dictionnary
    with gene as key and a list containing the number of
    nonsynonymous and synonymous sites in the gene's sequence
    '''
    
    # SNPs dict is {gene: [snps, nondegenerate, 2-fold, 4-fold]}
      
    # create a new dict {gene: [nonsyn, syn]}
    nonsyn_syn_sites = {}
    
    # rules to compute synonymous and nonsynonymous sites
    # number of synonymous sites is number of 4-fold degenerate sites + one third of
    # the number of 2-fold degenerate sites
    # N_syn = N_4-fold + 1/3 * N_2-fold
    
    # number of nonsynonymous sites is number of nondegenerate sites + two third
    # of the number of 2-fold degenerate sites
    # N_nonsyn = N_nondegenerate + 2/3 * N_2-fold
    
    # loop over genes
    for gene in SNPs:
        # compute the number of nonsynonymous sites
        N_nonsyn = SNPs[gene][1] + 2/3 * SNPs[gene][2]
        # compute the number of synonymous sites
        N_syn = SNPs[gene][3] + 1/3 * SNPs[gene][2]
        # populate dict
        nonsyn_syn_sites[gene] = [N_nonsyn, N_syn]
        
    return nonsyn_syn_sites



# use this function to count the degenerate sites in a codon
def count_degenerate_sites_in_codon(codon, position):
    '''
    (str, str) -> list
    Take a codon, and a string specifying the position of a given site in codon
    and return a 3-item list with counts of the number of nondegenerate, 2-fold
    and 4-fold degenerate sites in the codon
    '''
    
    # verify that codon has correct case
    codon = codon.upper()
    
    # initialise list [0-fold, 2-fold, 4-fold]
    sites = [0, 0, 0]
    
    # check all codons and count the number of degenerate, 2-fold and 4-fold sites
    # look up 2-fold degenerate codons
    if codon in {'AAA', 'AAC', 'AAT', 'AGC', 'AGT', 'ATA', 'ATC', 'ATG', 'ATT',
                     'CAA', 'CAC', 'CAG', 'CAT', 'GAA', 'GAC', 'GAG', 'GAT', 'TTC',
                     'TTT', 'TGT', 'TGC', 'TAA', 'TAC', 'TAG', 'TAT'}:
        # check position:
        if position == '3':
            # add 2-fold
            sites[1] += 1
        else:
            # non-degenerate
            sites[0] += 1
    # look up nondegenerate codons
    elif codon in {'TGG', 'AAG'}:
        # nondegenerate
        sites[0] += 1
    # look up 4-fold degenerate codons
    elif codon in {'ACA', 'ACC', 'ACG', 'ACT', 'CCA', 'CCC', 'CCG',
                   'CCT', 'CGC', 'CGT', 'CTC', 'CTT', 'GCA', 'GCC', 'GCG', 'GCT',
                   'GGA', 'GGC', 'GGG', 'GGT', 'GTA', 'GTC', 'GTG', 'GTT', 
                   'TCA', 'TCC', 'TCG', 'TCT'}:
        # check position
        if position == '3':
            # 4-fold
            sites[2] += 1
        else:
            # nondegenerate
            sites[0] += 1
    # look up the special codons
    elif codon in {'TTA', 'TTG', 'AGA', 'AGG'}:
        # check position
        if position == '3' or position == '1':
            # 2-fold
            sites[1] += 1
        elif position == '2':
            # nondegenerate
            sites[0] += 1
    elif codon in {'CTA', 'CTG', 'CGA', 'CGG'}:
        # check position
        if position == '3':
            # 4-fold
            sites[2] += 1
        elif position == '1':
            # 2-fold
            sites[1] += 1
        elif position == '2':
            # nondegenerate
            sites[0] += 1
           
    return sites


# use this function to compute theta for each gene   
def compute_theta_diversity(snp_file, valid_transcripts, sites, threshold):
    '''
    (file, set, str, int) -> dict
    Take the file with SNPs in CDS, the set of valid transcripts (mapping a 
    single transcript to a parent gene) and the type of sites to consider
    (REP, SYN or coding), a threshold indicating the minimum sample size acceptable
    and return a dict with theta per site for replacement, synonymous or all
    sites in CDS depending on the sites parameter
    '''
    
    # create a dict {gene: theta} to store theta per site for each gene
    theta_sites = {}
        
    # check site type
    if sites == 'coding':
        # count SNPs and number of sites for each gene
        # get the sample size of polymorphich sites
        SNPs, sample_size = count_coding_sites_for_diversity(snp_file, threshold)
        
        # loop over genes:
        for gene in SNPs:
            # check that gene has sites
            if SNPs[gene][1] != 0:
                # compute theta per site
                theta = SNPs[gene][0] / SNPs[gene][1]
                # check if theta > 0
                if theta != 0:
                    # need to correct for varying sample size
                    theta = theta / math.log(np.mean(sample_size[gene]) -1)
                theta_sites[gene] = theta

       
    elif sites == 'REP' or sites == 'SYN':
        # count SNPs and 0-fold, 2-fold and 4-fold genes for each gene
        # get the sample size of polymorphic sites
        SNPs, sample_size = count_rep_synonymous_sites_for_diversity(snp_file, sites, threshold)
        # compute the number of replacement and synonymous sites per gene
        N_nonsyn_syn = count_nonsynonymous_synonymous_sites(SNPs)
        # loop over genes in SNPs:
        for gene in SNPs:
            # check site type
            if sites == 'REP':
                # check that gene has nonsynonymous sites
                if N_nonsyn_syn[gene][0] != 0:
                    # compute theta per site
                    theta = SNPs[gene][0] / N_nonsyn_syn[gene][0]
                    # check if theta > 0
                    if theta != 0:
                        # need to correct for varying sample size
                        theta = theta / math.log(np.mean(sample_size[gene]) -1)
                    theta_sites[gene] = theta
            elif sites == 'SYN':
                # check that gene has synonymopus sites
                if N_nonsyn_syn[gene][1] != 0:
                    # compute theta per site
                    theta = SNPs[gene][0] / N_nonsyn_syn[gene][1]
                    # check if theta > 0
                    if theta != 0:
                        # need to correct for varying sample size
                        theta = theta / math.log(np.mean(sample_size[gene]) - 1)
                    theta_sites[gene] = theta
                
    # remove genes from theta_sites that are not in transcripts
    # to avoid counting SNPs lutiple times for transcripts of the same parent gene
    # create a lit of gene to remove
    to_remove = [gene for gene in theta_sites if gene not in valid_transcripts]
    # loop over genes to remove, delete from theta_sites
    for gene in to_remove:
        del theta_sites[gene]
        
    return theta_sites
    

# use this function to compute the p-distance between 2 aligned sequences
def pairwise_distance(seq1, seq2, bioseq):
    '''
    (str, str, str) -> float
    Return the p-distance between seq1 and seq2: the fraction of mismatched nuceotides
    if bioseq is 'DNA' or the or mismatched amino acids if bioseq is 'protein', 
    between seq1 and seq2. Precondition: seq1 and se2 must have the same length

    >>> pairwise_distance('CAGAATGCCT', 'CAGAAGGCCT')
    0.1
    >>> pairwise_distance('CAGAA-GTAA', 'CAGATTGTTA')
    0.222222
    >>> pairwise_distance('CAG-AAGTAA', 'CAGATNGTTA')
    0.25
    '''

    # test prerequisite of equal length
    assert len(seq1) == len(seq2), 'The sequences have different lengths, check alignment and/or missing end gaps'

    # accept - or ~ as gaps
    SEQ1 = seq1.replace('~', '-')
    SEQ2 = seq2.replace('~', '-')
    # make sure the sequences have the same case 
    SEQ1 = seq1.upper()
    SEQ2 = seq2.upper()
    
    # set up counters
    D = 0
    undefined = 0
    
    valid_nucleotides = {'A', 'T', 'C', 'G'}    
    valid_AA = {'A', 'R', 'N', 'D', 'C', 'E', 'Q', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V'}    
    
    # check bioseq
    if bioseq == 'DNA':
        # count differences between seq1 and seq2,
        # gapped positions are excluded, undefined positions are excluded
        # loop over sequence 1    
        for i in range(len(SEQ1)):
            # check if site is valid
            if SEQ1[i] not in valid_nucleotides or SEQ2[i] not in valid_nucleotides:
                undefined += 1
            elif SEQ1[i] in valid_nucleotides and SEQ2[i] in valid_nucleotides:
                # check if sites are different
                if SEQ1[i] != SEQ2[i]:
                    # if different, update distance counter
                    D += 1
    elif bioseq == 'protein':
        # count differences between seq1 and seq2
        # gapped positions are excluded, undefined positions are excluded
        # loop over sequence 1
        for i in range(len(SEQ1)):
            # check if site is valid
            if SEQ1[i] not in valid_AA or SEQ2[i] not in valid_AA:
                undefined += 1
            elif SEQ1[i] in valid_AA and SEQ2[i] in valid_AA:
                # check if sites are different
                if SEQ1[i] != SEQ2[i]:
                    # if different, update distance counter
                    D += 1
                    
    # distance = N differences / length of sequence (without undefined sites)    
    L = len(SEQ1) - undefined
    if L != 0:
        p = D/ L
        return round(p, 6)
    else:
        return 'NA'


# use this function to compute the Jukes-Cantor distance between 2 sequences
def compute_K_JC(seq1, seq2):
    '''
    (str, str) -> float
    Return the Jukes-Cantor distance between sequences seq1 and seq
    Prerequisite: seq1 and seq2 have same length
    >>> compute_K_JC('CAGAAGGCCTCGCCGAATATGACACTGTTGCGTAGCATGACCAATGAGACTCTACAAA-AGAGTGAATGGATCTGACTACGCTCAGTGGAATACCCGGCGAGTGCCTTCT',
    'CAGAAGGCCTCGCCGAATATGACACTGTTGCGTAGCATGACCAATGAGACTCTCCAAAGAGAGTGAATGGATCTGACTACGCTCAGTGGAATACCCGGCGAGTGCCTTCT')
    0.009231
    >>> compute_K_JC('AGCCGACGGAACGAGTAAATCTCATCCTAATCTGGTTG-ACACAACA-CAAGGAAGTGCTACCGATTTGGCTTGGGATTGACTTGTGGAAATGGCT',
    'AGCCGACGGAACGAGTAAATCTCATCCTAATCTGGTAGCACACAACAACAAGGAAGTGCTACCGATTTGGCTTGGGATTGACTTGTGGAATTGGCG')
    0.032614
    >>> compute_K_JC('TCAACGATGTGCTGCCTAACAGCGGATTCCAGTAACTCTGAGCATATACATGAAGTTTCAAAACAATGTAAATGCCAGTCGTTGCTGGAATTCGATAAGCACTCGTGAC',
    'TCAACGCTGTGCTGCCTAACAGCGGATTCCAGTAACTCTGAGCATATACATGAAGTTTCAAAACAATGTAAATGCCAGTCGTTGCTGGAATTCGATAAGCACTCGTGAC')
    0.009231
    '''

    # test prerequisite of equal length
    assert len(seq1) == len(seq2), 'The sequences have different lengths, check alignment and/or missing end gaps'

    # accept - or ~ as gaps
    SEQ1 = seq1.replace('~', '-').upper()
    SEQ2 = seq2.replace('~', '-').upper()
    
    D = 0
    gap = 0
    N = 0

    # count differences between seq1 and seq2
    # gapped positions and Ns are excluded
    for i in range(len(SEQ1)):
        if SEQ1[i] == '-' or SEQ2[i] == '-':
            gap += 1
        elif SEQ1[i] == 'N' or SEQ2[i] == 'N':
            N += 1
        elif SEQ1[i] != SEQ2[i]:
            D += 1
    L = len(SEQ1) - gap - N
    
    if L != 0:
        p = D/ L
        K = (-3/4) * math.log(1-(4/3 * p))
        return round(abs(K), 6)
    else:
        return 'NA'
            
# use this function to count the number of differences between 2 sequences
def match_diff(seq1, seq2):
    '''
    (str, str) -> int
    Take 2 sequences and return the number of differences between them
    Precondition: seq1 and seq2 have same length    
    '''
    # set up counter
    diff = 0
    # loop over seq1, compare each position to seq2
    for i in range(len(seq1)):
        if seq1[i].upper() != seq2[i].upper():
            diff += 1
            
    return diff


# use this function to find the indices where 2 sequences have a difference
def find_mismatch_positions(seq1, seq2):
    '''
    (str, str) -> list
    Take 2 sequences and return a list of indices where the 2 sequences differ
    Preconditions: the 2 sequences are aligned and/or have the same length
    '''
    
    # create a list to store the indices
    pos = []
    # loop over seq1, compare each position to seq2
    for i in range(len(seq1)):
        # compare seq1 and se2
        if seq1.upper()[i] != seq2.upper()[i]:
            pos.append(i)
            
    return pos

# use this function to find the indices of a SNP
def find_SNP_positions(chromo_sites, chromo, positions):
    '''
    (dict, str, list) -> list
    Take the dictionary with site positions : allele counts, the chromo where
    the sequence of interest is located and a sorted list with positions of the
    sequence of interest on chromo and return a list with the indices of the
    SNPs in the sequence
    Precondition: positions are 0-based, the sequence has coverage and
    minumum sample size for all its sites and the list of its positions is
    already ordered such that the sequence is orientated 5'-3'
    (positions are in decreasing order is the sequence orientation is '-')    
    '''
    
    # create a list to store the indisces of the snps
    SNPs = []    
    # loop over positions
    for i in range(len(positions)):
        # get the ref and alt allele counts
        ref_count = chromo_sites[chromo][positions[i]][2]
        alt_count = chromo_sites[chromo][positions[i]][3]
        # check if site is polymorphic
        if ref_count != 0 and alt_count != 0:
            # site is polymorphic, add i to list
            SNPs.append(i)
    return SNPs



# use this function to compute theta per site for a sequence of interest
# does not consider synonymous or replacement sites
def compute_theta_non_coding(chromo_sites, chromo, start, end, threshold):
    '''
    (dict, str, int, int, int) -> float
    Take the dictionary with allele counts in the genome, a chromo where the 
    sequence of interest is located, the coordinates of the sequence, and the
    maximum number of accptable missing sites and return theta per site
    adjusted for varying sample size among sites
    Precondition: all positions in the dict and start and end are 0-based    
    '''
    
    # get the positions of the sequence of interest
    positions = [i for i in range(start, end)]
    
    # count the number of polymorphic sites
    diff = 0
    
    # create a list to store the sample size of the polymorphic sites    
    sample_size = []
    
    # count the number of sites
    N = 0
    
    # loop over indices
    for i in positions:
        # check if site in chromo_sites
        if i in chromo_sites[chromo]:
            # adjust site counter
            N += 1
            # get reference, alternative allele and their counts
            ref, alt = chromo_sites[chromo][i][0], chromo_sites[chromo][i][1]
            ref_count, alt_count = chromo_sites[chromo][i][2], chromo_sites[chromo][i][3]
            # check if site is polymorphic
            if ref_count != 0 and alt_count != 0:
                # site is polymorphich
                # double check
                assert ref != alt, 'counts of both allele different than 0 but ref and alt are the same'
                # adjust polym site counter
                diff += 1
                # add sample size of polymorphich site to list
                sample_size.append(ref_count + alt_count)
    # check that a maximum of threshold site are missing
    if N >= len(positions) - threshold:
        # compute theta per site
        theta = diff / N
        if theta != 0:
            # apply a correction for varying sample size
            theta = theta / math.log(np.mean(sample_size) -1)
                
        return theta
    
    else:
        return 'NA'


    