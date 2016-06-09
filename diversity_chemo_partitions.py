# -*- coding: utf-8 -*-
"""
Created on Mon Jun 29 12:37:23 2015

@author: Richard
"""

from chemoreceptors import *
from divergence import *

    

# use this function to count SNPs, 0-fold, 2-fold and 4-fold sites
# in chemoreceptor partitions
def count_SNPs_degenerate_sites_chemo_paritions(snp_file, proba_domain, sites, threshold, cutoff):
    '''
    (file, dict, str, int, float) -> (dict, dict)
    Take the file with SNPs in CDS, a dictionnary with chemoreceptor gene as key
    and an inner dictionnary of codon index: list of probabilities pairs, 
    the type of sites (REP, SYN or coding) to consider, the minimum number of
    individuals with coverage at that site, the probability of identifying a domain, 
    and return a tuple of dictionniaries: one with each chemoreceptor gene as
    key and a list with the number of SNPs and the number of degenerate sites
    for partitions corresponding to transmembrane sites, one with the same count for
    extra-transmembrane sites and one dictionary with gene as key and the sample
    size of each polymorphic transmembrane site, and one dictionary with sample
    size of each extra-membrane polymorphic sites
    '''
    
    # open file for reading    
    infile = open(snp_file, 'r')
    # skip header
    infile.readline()
    
    # initiate gene variable
    gene = ''
    
    # initiate dictionnary for transmembrane sites
    # {gene : [snp, 0-fold, 2-fold, 4-fold]}
    TM = {}
    for gene in proba_domain:
        TM[gene] = [0, 0, 0, 0]
    # initiate dictionnary for extra-transmembrane sites
    # {gene : [snp, 0-fold, 2-fold, 4-fold]}    
    EX = {}
    for gene in proba_domain:
        EX[gene] = [0, 0, 0, 0]

    # create a dict to store the sample size of polymorphic sites
    TM_sample, EX_sample = {}, {}
    
    # create a set of stop codons
    stop_codons = {'TAG', 'TGA', 'TAA'} 
        

    # loop over file
    for line in infile:
        if line.rstrip() != '':
            line = line.rstrip().split()
            # check which gene is being read
            if line[2] != gene:
                # new gene, update gene variable
                gene = line[2]
                # initialize sample size dict
                TM_sample[gene], EX_sample[gene] = [], []
                # initiate codon index
                if line[4] == '1':
                    j = 0
            elif line[2] == gene:
                # still same gene, update codon index
                if line[4] == '1':
                    j += 3
            # check that gene is chemoreceptor
            if gene in proba_domain:
                # record only sites snpn, no_snp and with valid snp type
                if line[6] in {'snp', 'no_snp'} and line[9] in {'NA', 'COD', 'SYN', 'REP'}:
                    # get reference and alternative codons
                    ref_codon, alt_codon = line[3], line[8]
                    # get reference and alternative alleles
                    ref, alt = line[5], line[7]
                    # get reference anc alternative allele counts
                    ref_count, alt_count = int(line[13]), int(line[17]) 
                    # check that sample size is greater than threshold
                    if ref_count + alt_count >= threshold:
                        # do not consider stop codons
                        if ref_codon not in stop_codons and alt_codon not in stop_codons:
                            # do not consider codons with N
                            if 'N' not in ref_codon and 'N' not in alt_codon:
                                # do not consider codons with more than 1 differences
                                if diff_codon(ref_codon, alt_codon) <= 1:
                                    # get site degeneracy
                                    degeneracy = count_degenerate_sites_in_codon(ref_codon, line[4])
                                    # verify that ref codon correspond to the correct AA in domain dict
                                    assert cds_translate(ref_codon) == proba_domain[gene][j][0], 'the amino acid don\'t match'
                                    # check the domain probability of site 
                                    if proba_domain[gene][j][3] >= cutoff:
                                        # codon in TM, update the degenerate site counter
                                        for i in range(len(degeneracy)):
                                            # add the value to indices 1 to 3 in list
                                            TM[gene][i+1] += degeneracy[i]
                                    elif proba_domain[gene][j][1] >= cutoff or proba_domain[gene][j][2] >= cutoff or proba_domain[gene][j][4] >= cutoff:
                                        # codon in EX, update the degenerate site counter
                                        for i in range(len(degeneracy)):
                                            # add the value to indices 1 to 3 in list
                                            EX[gene][i+1] += degeneracy[i]
                                    # check if site is polymorphic
                                    if ref_count != 0 and alt_count != 0:
                                        # site is polymorphic, # verify SNP by comparing ref and alt
                                        assert ref != alt,'ref and alt counts > 0 but alleles are the same'
                                        # site is polymorphic
                                        # check the domain probability of site 
                                        if proba_domain[gene][j][3] >= cutoff:
                                            # site in TM
                                            # check what type of SNP
                                            if sites == 'SYN' and cds_translate(ref_codon) == cds_translate(alt_codon):
                                                # snp is synonymous, add snp to TM parition
                                                TM[gene][0] += 1
                                                # record sample size of polymorphic site
                                                TM_sample[gene].append(ref_count + alt_count)
                                            elif sites == 'REP' and cds_translate(ref_codon) != cds_translate(alt_codon):
                                                # snp is nonsynonymous, add snp to TM parition
                                                TM[gene][0] += 1
                                                # record sample size of polymorphic site
                                                TM_sample[gene].append(ref_count + alt_count)
                                        elif proba_domain[gene][j][1] >= cutoff or proba_domain[gene][j][2] >= cutoff or proba_domain[gene][j][4] >= cutoff:
                                            # codon in EX
                                            # check what type of SNP
                                            if sites == 'SYN' and cds_translate(ref_codon) == cds_translate(alt_codon):
                                                # snp is synonymous, add snp to EX partition
                                                EX[gene][0] += 1
                                                # record sample size of polymorphic site
                                                EX_sample[gene].append(ref_count + alt_count)
                                            elif sites == 'REP' and cds_translate(ref_codon) != cds_translate(alt_codon):
                                                # snp is nonsynonymous, add snp to EX partition
                                                EX[gene][0] += 1
                                                # record sample size of polymorphic site
                                                EX_sample[gene].append(ref_count + alt_count)
                                            
    # close file
    infile.close()
    
    # remove genes that do not have both partitions
    to_remove = set()
    for gene in TM:
        if sum(TM[gene]) == 0 or sum(EX[gene]) == 0:
            to_remove.add(gene)
    # loop over genes and delete from dicts
    for gene in to_remove:
        del TM[gene]
        del EX[gene]
        
    return TM, TM_sample, EX, EX_sample
    

# use this function to compute theta for TM and EX paritions in chemoreceptors
def compute_theta_chemo_partitions(TM, TM_sample, EX, EX_sample, sites, valid_transcripts, site_threshold):
    '''
    (dict, dict, dict, dict, int, set, int) -> dict, dict
    Take the dictionaries with number of SNPs and number of degenerate sites for 
    TM and EX partitions, the dictionaries with sample size for for each gene,
    the type of sites to consider (SYN and REP) and the set of valid transcripts,
    and the minimum number of sites in a given's partition and compute theta
    corrected for differences in sample size for each partition of each gene
    Precondition: compute theta for Ontario sample (KSR + PX)
    '''
    
    # create a dict {gene : theta pairs} for TN and EX partitions
    TM_theta, EX_theta = {}, {}
   
    # compute the number of replacement and synonymous sites per gene
    TMSynNonSyn = count_nonsynonymous_synonymous_sites(TM)
    EXSynNonSyn = count_nonsynonymous_synonymous_sites(EX)
    
    # loop over genes in TM:
    for gene in TM:
        # only consider genes found in EX
        if gene in EX:
            # check site type, and that genes has site type
            if sites == 'REP':
                TMsites = TMSynNonSyn[gene][0]
                EXsites = EXSynNonSyn[gene][0]            
            elif sites == 'SYN':
                TMsites = TMSynNonSyn[gene][1]
                EXsites = EXSynNonSyn[gene][1]
            # check that gene has nonsynonymous sites    
            if TMsites != 0 and TMsites > site_threshold:
                # compute theta per site
                theta = TM[gene][0] / TMsites
                # check if theta > 0
                if theta != 0:
                    # need to correct for varying sample size
                    theta = theta / math.log(np.mean(TM_sample[gene]) -1)
                TM_theta[gene] = theta
            if EXsites != 0 and EXsites > site_threshold:
                # compute theta per site
                theta = EX[gene][0] / TMsites
                # check if theta > 0
                if theta != 0:
                    # need to correct for varying sample size
                    theta = theta / math.log(np.mean(EX_sample[gene]) -1)
                EX_theta[gene] = theta
      
    # remove genes not in valid_transcripts to avoid counting SNPs multiple times for transcripts of the same parent gene
    to_remove = [gene for gene in TM_theta if gene not in valid_transcripts]
    # loop over genes to remove, delete from theta_sites
    for gene in to_remove:
        del TM_theta[gene]
    to_remove = [gene for gene in EX_theta if gene not in valid_transcripts]        
    for gene in to_remove:
        del EX_theta[gene]
        
    # remove genes that do not have both partitions
    to_remove = set()
    for gene in TM_theta:
        if gene not in EX_theta:
            to_remove.add(gene)
    for gene in EX_theta:
        if gene not in TM_theta:
            to_remove.add(gene)
    for gene in to_remove:
        if gene in TM_theta:
            del TM_theta[gene]
        if gene in EX_theta:
            del EX_theta[gene]
        
    return TM_theta, EX_theta