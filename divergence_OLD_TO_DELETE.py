# -*- coding: utf-8 -*-
"""
Created on Fri Jun 26 17:16:38 2015

@author: Richard
"""


from MK_test import *
from accessories import *



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
        
        
# use this function to compute theta        
def watterson_theta(k, N):
    '''
    (int, int) -> float
    Take the number of polymorphic sites and the sample size and return
    Watterson's theta estimator of polymorphism (per sequence)
    '''
    
    # compute a (1/i from i = 1 to i = N-1)
    a = 0
    for i in range(1, N):
        a += 1/i
    # compute theta
    theta = k/a
    
    return theta
    
    
    
# use this function to count the number of sites to analyze per gene
def count_sites(snp_file):
    '''
    (file) -> dict
    Take the file with SNPs in the CDS and return a dict with gene as key and
    the number of sites invariant and polymorphic as value couting only valid sites
    '''
    
    # set gene variable
    gene = ''
    
    # create dict: {gene : number of sites}
    sites = {}
    
    # create a set of stop codons
    stop_codons = {'TAG', 'TGA', 'TAA'} 
    
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
                # new gene, assign new gene to gene and initialize dict
                gene = line[2]
                sites[gene] = 0
            # record only sites labeled snp or no_snp and valid snp type
            if line[6] in {'snp', 'no_snp'} and line[9] in {'NA', 'COD', 'SYN', 'REP'}:
                # get reference codon
                ref_codon = line[3]
                # get alternative codon
                alt_codon = line[8]
                # do not consider stop codons
                if ref_codon not in stop_codons and alt_codon not in stop_codons:
                    # do not consider codons with 'Ns'
                    if 'N' not in ref_codon and 'N' not in alt_codon:
                        # do not consider codons with more than 1 differeneces
                        if diff_codon(ref_codon, alt_codon) <= 1:
                            # count the number of sites
                            sites[gene] += 1
                        
    # close file after reading
    infile.close()
    
    return coding
    

# use this function to determine the number of nondegenerate, 2-fold and 4-fold degenerate sites
def count_degenerate_sites(snp_file):
    '''
    (file) -> dict
    Take the file with SNPs in CDS and return a dictionnary with gene as key
    and a list with counts of non-degenerate, 2-fold and 4-fold degenerate sites
    as value
    '''

    # create dict {gene: [nondegenerate, 2-fold, 4-fold]}
    coding = {}    
    
    # initialize gene : empty list pairs
    # open file for reading
    infile = open(snp_file, 'r')
    # skip header
    infile.readline()
    # loop over file, get the gene name as key and initialize list
    for line in infile:
        if line.rstrip() != '':
            line = line.rstrip().split()
            gene = line[2]
            if gene not in coding:
                coding[gene] = [0, 0, 0]
    # close file after reading
    infile.close()
    
    # create a set of stop codons
    stop_codons = {'TAG', 'TGA', 'TAA'} 
    
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
            # record only sites labeled snp or no_snp and valid snp type
            if line[6] in {'snp', 'no_snp'} and line[9] in {'NA', 'COD', 'SYN', 'REP'}:
                # get reference codon
                ref_codon = line[3]
                # get alternative codon
                alt_codon = line[8]
                # do not consider stop codons
                if ref_codon not in stop_codons and alt_codon not in stop_codons:
                    # do not consider codons with 'Ns'
                    if 'N' not in ref_codon and 'N' not in alt_codon:
                        # check that codons differ by at most 1 change
                        if diff_codon(ref_codon, alt_codon) <= 1:
                            # check all codons and count the number of degenerate, 2-fold and 4-fold sites
                    
                            # look up 2-fold degenerate codons
                            if ref_codon in {'AAA', 'AAC', 'AAT', 'AGC', 'AGT', 'ATA', 'ATC', 'ATG', 'ATT',
                                             'CAA', 'CAC', 'CAG', 'CAT', 'GAA', 'GAC', 'GAG', 'GAT', 'TTC',
                                             'TTT', 'TGT', 'TGC', 'TAA', 'TAC', 'TAG', 'TAT'}:
                                # check position:
                                if line[4] == '3':
                                    # add 2-fold 
                                    coding[gene][1] += 1
                                else:
                                    # non-degenerate
                                    coding[gene][0] += 1
                       
                            # look up nondegenerate codons
                            elif ref_codon in {'TGG', 'AAG'}:
                                # nondegenerate
                                coding[gene][0] += 1
                            
                            # look up 4-fold degenerate codons
                            elif ref_codon in {'ACA', 'ACC', 'ACG', 'ACT', 'CCA', 'CCC', 'CCG',
                                               'CCT', 'CGC', 'CGT', 'CTC', 'CTT', 'GCA', 'GCC', 'GCG', 'GCT',
                                               'GGA', 'GGC', 'GGG', 'GGT', 'GTA', 'GTC', 'GTG', 'GTT',
                                               'TCA', 'TCC', 'TCG', 'TCT'}:
                                # check position
                                if line[4] == '3':
                                    # 4-fold
                                    coding[gene][2] += 1
                                else:
                                    # nondegenerate
                                    coding[gene][0] += 1
                        
                            # look up the special codons                        
                            elif ref_codon in {'TTA', 'TTG', 'AGA', 'AGG'}:
                                # check position
                                if line[4] == '3' or line[4] == '1':
                                    # 2-fold
                                    coding[gene][1] += 1
                                elif line[4] == '2':
                                    # nondegenerate
                                    coding[gene][0] += 1
                            elif ref_codon in {'CTA', 'CTG', 'CGA', 'CGG'}:
                                # check position
                                if line[4] == '3':
                                    # 4-fold
                                    coding[gene][2] += 1
                                elif line[4] == '1':
                                    # 2-fold
                                    coding[gene][1] += 1
                                elif line[4] == '2':
                                    # nondegenerate
                                    coding[gene][0] += 1
                   
    infile.close()
    
    
    # remove genes with empty lists (ie. genes not detected)
    to_remove = []
    for gene in coding:
        if coding[gene] == [0,0,0]:
            to_remove.append(gene)
    for gene in to_remove:
        del coding[gene]
    
    return coding


# use this function to count the number of synonymous and nonsynonymous sites in the reference
def count_nonsynonymous_synonymous_sites(coding):
    '''
    (dict) -> dict
    Take a dict with genes as key and a list with the number of non-degenerate
    sites, 2-fold and 4-fold degenerate sites in gene as value, and return
    a new dictionnary with gene as key and a list containing the number of
    nonsynonymous and synonymous sites in the gene's sequence
    '''
    
    # coding dict is {gene: [nondegenerate, 2-fold, 4-fold]}
      
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
    for gene in coding:
        # compute the number of nonsynonymous sites
        N_nonsyn = coding[gene][0] + 2/3 * coding[gene][1]
        # compute the number of synonymous sites
        N_syn = coding[gene][2] + 1/3 * coding[gene][1]
        # populate dict
        nonsyn_syn_sites[gene] = [N_nonsyn, N_syn]
        
    return nonsyn_syn_sites
        

# use this function to count the number of nonsynonymous and synonymous changes    
def count_polymorphisms_coding(snp_file, sites):
    '''
    (file, str) -> dict
    Take the file with SNPs in the CDS, the type of changes to record 
    (REP, SYN or coding) and return a dict with gene as key and the number
    of nonynonymous changes or synonymous changes or total number of changes
    in the CDS
    '''
       
    # create a set of stop codons
    stop_codons = {'TAG', 'TGA', 'TAA'} 
    
    # create dict {gene: snps}
    SNPs = {}
    
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
                SNPs[gene] = 0
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
                # do not consider stop codons
                if ref_codon not in stop_codons and alt_codon not in stop_codons:
                    # do not consider codons with 'Ns'
                    if 'N' not in ref_codon and 'N' not in alt_codon:
                        # get reference allele count
                        ref_count = int(line[13])
                        # get alternative allele count
                        alt_count = int(line[17])
                        # check if site is polymorphic
                        if ref_count != 0 and alt_count != 0:
                            # site is polymorphic
                            # check that ref alleles and alt alleles are different
                            if ref != alt:
                                # check the number of differences between ref_codon and alt_codon within remanei
                                if diff_codon(ref_codon, alt_codon) <= 1:
                                    # do not consider codons with more than 1 substitutions
                                    # determine the type of snp
                                    if sites == 'SYN' and genetic_code[ref_codon] == genetic_code[alt_codon]:
                                        # mutation is synonymous
                                        SNPs[gene] += 1
                                    elif sites == 'REP' and genetic_code[ref_codon] != genetic_code[alt_codon]:
                                        # mutation is nonsynonymous
                                        SNPs[gene] += 1
                                    elif sites == 'coding':
                                        # record all SNPs in the CDS
                                        SNPs[gene] += 1
                        
    # close file after readling
    infile.close()
    
    # return dict
    return SNPs
    
    

   
# use this function to compute theta for each gene   
def compute_theta_diversity(snp_file, unique_transcripts, sites):
    '''
    (file, file, str) -> dict
    Take the file with SNPs in CDS, the file with valid transcripts (mapping a 
    single transcript to a parent gene) and the type of sites to consider
    (REP, SYN or coding) and return a dict with theta per site for replacement,
    synonymous or all sites in CDS depending on the sites parameter
    Precondition: use SNPs in the Ontario strains (N = 49)
    '''
    
    # create a dict {gene: theta} to store theta per site for each gene
    theta_sites = {}
    
    # count the number of SNPs
    N_SNPs = count_polymorphisms_coding(snp_file, sites)
        
    # check the type of sites to consider
    if sites == 'coding':
        # compute theta over the entire CDS
        # count the number of sites
        N_sites = count_sites(snp_file)
        # loop over gene in SNPs
        for gene in N_SNPs:
            # compute theta per sequence
            k = N_SNPs[gene]
            theta = watterson_theta(k, 49)
            # compute theta per site
            theta = theta / N_sites[gene]
            # populate dict
            theta_sites[gene] = theta
        
    else:
        # count the number of degenerate sites
        N_degenerate = count_degenerate_sites(snp_file)
        # count the number of nonsynonymous and synonymous sites
        N_nonsyn_syn = count_nonsynonymous_synonymous_sites(N_degenerate)
        # loop over gene in SNPs:
        for gene in N_SNPs:
            # compute theta per sequence
            k = N_SNPs[gene]
            theta = watterson_theta(k, 49)
            # compute theta per site
            # check which site to consider
            if sites == 'REP':
                # check if gene in N_nonsyn_syn (ie if gene is detected)
                if gene in N_nonsyn_syn:
                    # check if number of nonsynonymous sites is different than 0
                    if N_nonsyn_syn[gene][0] != 0:
                        # compute theta per replacement sites
                        theta = theta / N_nonsyn_syn[gene][0]
                        # populate dict
                        theta_sites[gene] = theta
            elif sites == 'SYN':
                # check if gene in N_nonsyn_syn (ie if gene is detected)
                if gene in N_nonsyn_syn:
                    # check if gene has synonymous sites
                    if N_nonsyn_syn[gene][1] != 0:
                        # compute theta per synonymous site
                        theta = theta / N_nonsyn_syn[gene][1]
                        # populate dict
                        theta_sites[gene] = theta
        
    # get set of valid genes
    transcripts = get_valid_transcripts(unique_transcripts)
    
    # remove genes from theta_sites that are not in transcripts
    # to avoid counting SNPs lutiple times for transcripts of the same parent gene
    # create a lit of gene to remove
    to_remove = [gene for gene in theta_sites if gene not in transcripts]
    # loop over genes to remove, delete from theta_sites
    for gene in to_remove:
        del theta_sites[gene]
        
    return theta_sites
            
    
    
    
    
    
    
    
    
    
    
    