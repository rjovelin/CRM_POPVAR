# -*- coding: utf-8 -*-
"""
Created on Mon Jun 29 12:37:23 2015

@author: Richard
"""

from chemoreceptors import *
from divergence import *

    






# use this function to count SNPs, 0-fold, 2-fold and 4-fold sites
# in chemoreceptor partitions
def count_SNPs_degenerate_sites_chemo_paritions(snp_file, proba_domain, sites, cutoff):
    '''
    (file, dict, str) -> (dict, dict)
    Take the file with SNPs in CDS, a dictionnary with chemoreceptor gene as key
    and an inner dictionnary of codon index: list of probabilities pairs, 
    the type of sites (REP, SYN or coding) to consider and return a tuple
    of dictionniaries with each chemoreceptor gene as key and a list with
    the number of SNPs and the number of degenerate sites for partitions
    corresponding to transmembrane sites and extra-transmembrane sites
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
                    # get refeence codon
                    ref_codon = line[3]
                    # get alternatiove codon
                    alt_codon = line[8]
                    # get reference allele
                    ref = line[5]
                    # get alternative allele
                    alt = line[7]
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
                                # get reference allele count
                                ref_count = int(line[13])
                                # get alternative allele count
                                alt_count = int(line[17])
                                # check the domain probability of site 
                                if proba_domain[gene][j][3] >= cutoff:
                                    # codon in TM
                                    # update the degenerate site counter
                                    for i in range(len(degeneracy)):
                                        # add the value to indices 1 to 3 in list
                                        TM[gene][i+1] += degeneracy[i]
                                elif proba_domain[gene][j][1] >= cutoff or proba_domain[gene][j][2] >= cutoff or proba_domain[gene][j][4] >= cutoff:
                                    # codon in EX
                                    # update the degenerate site counter
                                    for i in range(len(degeneracy)):
                                        # add the value to indices 1 to 3 in list
                                        EX[gene][i+1] += degeneracy[i]
                                # check if site is polymorphic
                                if ref_count != 0 and alt_count != 0:
                                    # site is polymorphic
                                    # verify SNP by comparing ref and alt
                                    if ref != alt:
                                        # site is polymorphic
                                        # check the domain probability of site 
                                        if proba_domain[gene][j][3] >= cutoff:
                                            # site in TM
                                            # check what type of SNP
                                            if sites == 'SYN' and cds_translate(ref_codon) == cds_translate(alt_codon):
                                                # snp is synonymous, add snp to TM parition
                                                TM[gene][0] += 1
                                            elif sites == 'REP' and cds_translate(ref_codon) != cds_translate(alt_codon):
                                                # snp is nonsynonymous, add snp to TM parition
                                                TM[gene][0] += 1
                                            elif sites == 'coding':
                                                # record all SNPs in CDS
                                                TM[gene][0] += 1
                                        elif proba_domain[gene][j][1] >= cutoff or proba_domain[gene][j][2] >= cutoff or proba_domain[gene][j][4] >= cutoff:
                                            # codon in EX
                                            # check what type of SNP
                                            if sites == 'SYN' and cds_translate(ref_codon) == cds_translate(alt_codon):
                                                # snp is synonymous, add snp to EX partition
                                                EX[gene][0] += 1
                                            elif sites == 'REP' and cds_translate(ref_codon) != cds_translate(alt_codon):
                                                # snp is nonsynonymous, add snp to EX partition
                                                EX[gene][0] += 1
                                            elif sites == 'coding':
                                                # record all SNPs in CDS
                                                EX[gene][0] += 1
                                            
                                            
    # close file
    infile.close()
    
    # remove genes that do not have both partitions
    to_remove = []
    for gene in TM:
        if sum(TM[gene]) == 0 or sum(EX[gene]) == 0:
            to_remove.append(gene)
    # loop over genes and delete from dicts
    for gene in to_remove:
        del TM[gene]
        del EX[gene]
        
    return TM, EX
    





# use this function to compute theta for partition of chemoreceptor genes
def from_counts_to_diversity_partition(partition, sites, site_threshold):
    '''
    (dict, str, int) -> dict
    Take a dictionnary with counts of SNPs and degenerate sites for a given
    chemoreceptor gene partition, a string specifying the type
    of site to consider (REP, SYN, coding), the mininum number of sites to
    accept to compute theta, and return and a dicttionnary with gene : theta
    pairs for that partition
    Precondition: compute theta for Ontario sample
    '''
    
    # N_syn = N_4-fold + 1/3 * N_2-fold
    # N_nonsyn = N_nondegenerate + 2/3 * N_2-fold
    
    # initiate dicts
    partition_theta ={}

    # loop over genes in partition
    for gene in partition:
        # check the site type
        if sites == 'REP':
            # count the number of nonsynonymous sites
            N_sites = partition[gene][1] + 2/3 * partition[gene][2]
        elif sites == 'SYN':
            # count the number of synonymous sites
            N_sites = partition[gene][3] + 1/3 * partition[gene][2]
        elif sites == 'coding':
            # count the number of sites
            N_sites = sum(partition[gene][1:])
        # check sites greater than threshold
        if N_sites >= site_threshold:
            # compute theta per site
            theta = watterson_theta(partition[gene][0], 49) / N_sites
            # populate dict
            partition_theta[gene] = theta
    
    return partition_theta


# use this function to compute theta for TM and EX paritions in chemoreceptors
def compute_theta_chemo_partitions(TM, EX, sites, site_threshold):
    '''
    (dict, dict, str) -> dict, dict
    Take dictionnaries with counts of SNPs and degenerate sites for transmbrane
    (TM) and extra-transmembrane (X) partitions, a string specifying the type
    of site to consider (REP, SYN, coding), the mininum number of sites to
    accept to compute theta and return two dictionnaries with gene : theta pairs
    for TM and EX partitions
    Precondition: compute theta for Ontario sample
    '''
    
    # compute theta for transmebrane partition
    TM_theta = from_counts_to_diversity_partition(TM, sites, site_threshold)
    # compute theta for extra-transmebrane partition
    EX_theta = from_counts_to_diversity_partition(EX, sites, site_threshold)
    
    # remove genes that are not in both dictionnaries
    # because paired tests are used to compare TM and EX partitions
    to_remove = []
    for gene in TM_theta:
        if gene not in EX_theta:
            to_remove.append(gene)
    for gene in to_remove:
        del TM_theta[gene]
    to_remove = []
    for gene in EX_theta:
        if gene not in TM_theta:
            to_remove.append(gene)
    for gene in to_remove:
        del EX_theta[gene]
    
    
    return TM_theta, EX_theta
            
    