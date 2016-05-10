# -*- coding: utf-8 -*-
"""
Created on Wed Jun 10 11:08:59 2015

@author: Richard
"""

#!/usr/bin/env python3


# import modules
from manipulate_sequences import *
from mk_test import * 
from genomic_coordinates import *

  
# use this function to identify the genes that have premature stop codons caused by SNPs only
def genes_with_premature_stops(snp_file, unique_transcripts, indel_transcripts):
    '''
    (file, file) -> dict
    Take the file of SNPs in CDS, the file with unique transcripts and the file
    with genes with indels and return a dictionnary of transcripts with SNP-caused
    premature stop codons and a list with the number of stop mutations and premature
    stop codons as value
    Note: PTC with > 1 difference from reference allele are considered in gene counts
    '''
    
    # get the set of valid transcripts
    transcripts = get_valid_transcripts(unique_transcripts)
    
    # get the set of transcripts with indels in their CDS
    indels = get_genes_with_indels(indel_transcripts)
    
    # create a set of valid nucleotides
    valid_base = {'A', 'T', 'C', 'G'}    
        
    # create a dict of genes : [number of stop mutations, number of stop codons]
    stops = {}   
    
    # set up boolean to count number of codons with stop mutations
    valid_stop = False    
    
    
    # open file for reading
    infile = open(snp_file, 'r')
    # skip header
    infile.readline()
    # loop over file
    for line in infile:
        line = line.rstrip()
        if line != '':
            line = line.split()
            # check that gene is in valid transcripts
            # and does not have indels
            gene = line[2]
            if gene in transcripts and gene not in indels:
                # check only valid snps
                if line[6] == 'snp' and line[5] != line[7] and line[5] in valid_base and line[7] in valid_base:
                    # get ref codon
                    ref_codon = line[3]
                    # get alternative codon
                    alt_codon = line[8]
                    # verify that reference is not terminal stop codon
                    if genetic_code[ref_codon] != '*':
                        # verify that SNP is in KSR+PX
                        alt_count = int(line[17])
                        if alt_count != 0:
                            # alternative allele present in KSR+PX
                            # check that alternative codon is stop codon
                            if genetic_code[alt_codon] == '*':
                                # mutation results in stop
                                # populate dict
                                if gene in stops:
                                    # update the number of mutation
                                    stops[gene][0] += 1
                                else:
                                    # add 1 to number of mutation, don't count number of codon yet
                                    stops[gene] = [1, 0]
                                # modify boolean value
                                valid_stop = True
                # check if the codon is valid after the codon is evaulated
                if line[4] == '3':
                    # codon has been evaluated, still same gene and same codon
                    if valid_stop == True:
                        # count stop codon, gene is necessarily a key already
                        stops[gene][1] += 1
                        # reset boolean value
                        valid_stop = False
                        
                
    # close file
    infile.close()
    
    return stops


# use this function to count the number of transcripts that have PTC
# regardless of the presence of indels
def count_genes_with_indels_premature_stops(snp_file, unique_transcripts):
    '''
    (file, file) -> set
    Take the file with snp positions in coding sequences and the file with
    valid transcripts and return a set of transcripts that have at least
    1 premature stop codon that could be caused by an indel and or a SNP
    Note: PTC with > 1 difference from reference allele are considered in gene counts
    '''
    
    # get the set of valid transcripts
    transcripts = get_valid_transcripts(unique_transcripts)
    
    # create a set of valid nucleotides
    valid_base = {'A', 'T', 'C', 'G'}    
        
    # create a set of transcripts
    PTC = set()   
    
    # open file for reading
    infile = open(snp_file, 'r')
    # skip header
    infile.readline()
    # loop over file
    for line in infile:
        line = line.rstrip()
        if line != '':
            line = line.split()
            # check that gene is in valid transcripts
            gene = line[2]
            if gene in transcripts:
                # check only valid snps
                if line[6] == 'snp' and line[5] != line[7] and line[5] in valid_base and line[7] in valid_base:
                    # get ref codon
                    ref_codon = line[3]
                    # get alternative codon
                    alt_codon = line[8]
                    # verify that reference is not terminal stop codon
                    if genetic_code[ref_codon] != '*':
                        # verify that allele is in KSR+PX
                        alt_count = int(line[17])
                        if alt_count != 0:
                            # alternative allele present in KSR+PX
                            # check that alternative codon is stop codon
                            if genetic_code[alt_codon] == '*':
                                # mutation results in stop
                                # add gene to set
                                PTC.add(gene)
    # close file
    infile.close()
    
    return PTC


# use this function to count the total number of premature stop mutations
# (relative to the reference)
def count_total_PTC_SNP(snp_file, unique_transcripts):
    '''
    (file, file) -> int
    Take the file of SNPs in the CDS of the remanei genes, a set of transcripts
    from unique parent genes and return the number of SNPs causing premature
    stop codons    
    Note: PTC with > 1 difference from reference allele are considered in gene counts
    '''
    
    # get the set of valid transcripts
    transcripts = get_valid_transcripts(unique_transcripts)    
    
    # set up counter variable
    stops = 0    
    
    # create a set of valid nucleotides
    valid_base = {'A', 'T', 'C', 'G'}    
    
    # open file for reading
    infile = open(snp_file, 'r')
    #skip header
    infile.readline()
    # loop over file
    for line in infile:
        line = line.rstrip()
        if line != '':
            line = line.split()
            # get the transcript name
            gene = line[2]
            # check that transcript is valid
            if gene in transcripts:
                # check only valid snps
                if line[6] == 'snp' and line[5] != line[7] and line[5] in valid_base and line[7] in valid_base:
                    # get the reference codon
                    ref_codon = line[3]
                    # get the alternative codon
                    alt_codon = line[8]
                    # verifiy that reference codon is not the terminal stop codon
                    if genetic_code[ref_codon] != '*':
                        # verify that alternative allele is in the KSR+PX group
                        alt_count = int(line[17])
                        if alt_count != 0:
                            # alternative allele present in the KSR+PX strains
                            # check that alternative codon is stop codon
                            if genetic_code[alt_codon] == '*':
                                # mutation results in stop
                                stops += 1
    # close file
    infile.close()
    
    return stops
                    
                
# use this function to get the MAF in KSR and PX combined    
def MAF_SNP(snp_file, unique_transcripts, indel_transcripts, site_type, threshold):
    '''
    (file, file, file, str, int) -> list
    Take the file of snps, the file with valid transcript names, the file with
    the transcripts having indels, the type of mutation to record (REP, PTC
    or SYN), and a threshold of the minimum sample size to accept and return
    a list with minor allele frequencies in the KSR and PX strains combined
    Precondition: codons with more than 1 SNP are not considered
    '''
    
    # create a list to store the frequencies
    MAF = []    
        
    # get the set of valid transcripts
    transcripts = get_valid_transcripts(unique_transcripts)    
    
    # get the set of transcripts with indels in their CDS
    indels = get_genes_with_indels(indel_transcripts)
    
    # create a set of valid nucleotides
    valid_base = {'A', 'T', 'C', 'G'}    
        
    # open file for reading
    infile = open(snp_file, 'r')
    # skip header
    infile.readline()
    # loop over file
    for line in infile:
        line = line.rstrip()
        if line != '':
            line = line.split()
            # check that gene is in valid transcripts
            gene = line[2]
            if gene in transcripts:
                # check only valid snps
                if line[6] == 'snp' and line[5] != line[7] and line[5] in valid_base and line[7] in valid_base:
                    # get ref codon
                    ref_codon = line[3]
                    # get alternative codon
                    alt_codon = line[8]
                    # do not consider Ns in codons
                    if 'N' not in ref_codon and 'N' not in alt_codon:
                        # get allele counts
                        ref_count = int(line[13])
                        alt_count = int(line[17])
                        # verify that SNP is in KSR+PX
                        if alt_count != 0 and ref_count != 0:
                            # SNP present in KSR+PX
                            # check that sample size is greater than threshold
                            if ref_count + alt_count >= threshold:
                                # check that ref and alt codons have only 1 difference
                                if diff_codon(ref_codon, alt_codon) <= 1:
                                    # don't count reference stop codons
                                    if genetic_code[ref_codon] != '*':
                                        # reference codon is not terminal stop
                                        # check site_type
                                        if site_type == 'PTC':
                                            # check that gene is not in indels
                                            if gene not in indels:
                                                # count premature stops
                                                # check that alternative allele is stop codon
                                                if genetic_code[alt_codon] == '*':
                                                    # mutation results in stop
                                                    # check which allele is minor allele
                                                    if alt_count <= ref_count:
                                                        # alternative alelle is minor
                                                        # add minor frequency to list
                                                        MAF.append(alt_count / (alt_count + ref_count))
                                                    elif alt_count > ref_count:
                                                        # reference allele is minor
                                                        # add minor frequency to list
                                                        MAF.append(ref_count / (alt_count + ref_count))
                                        elif site_type == 'REP':
                                            # count replacement site
                                            # verify that change is non-synonymous but not PTC
                                            if genetic_code[ref_codon] != genetic_code[alt_codon] and genetic_code[alt_codon] != '*':
                                                # site is nonsynonymous
                                                # check which allele is minor 
                                                if alt_count <= ref_count:
                                                    # alternative allele is minor
                                                    # add minor frequency to list
                                                    MAF.append(alt_count / (alt_count + ref_count))
                                                elif alt_count > ref_count:
                                                    # reference codon is minor
                                                    # add minor frequency to list
                                                    MAF.append( ref_count / (alt_count + ref_count))
                                        elif site_type == 'SYN':
                                            # count synonymous sites
                                            # verify that site is synonymous
                                            if genetic_code[ref_codon] == genetic_code[alt_codon]:
                                                # site is synonymous
                                                # check which allele is minor
                                                if alt_count <= ref_count:
                                                    # alternative allele is minor
                                                    # add minor frequency to list
                                                    MAF.append(alt_count / (ref_count + alt_count))
                                                elif alt_count > ref_count:
                                                    # reference allele is minor
                                                    # add minor frequency to list
                                                    MAF.append(ref_count / (ref_count + alt_count))
                                       
    # close file after reading
    infile.close()
    
    return MAF
    
    
# use this function to get the DAF in KSR and PX combined    
def DAF_SNP(snp_file, unique_transcripts, indel_transcripts, site_type, threshold):
    '''
    (file, file, file, str, int) -> list
    Take the file of snps, the file with valid transcript names, the file with
    the transcripts having indels and the type of mutation to record (REP, PTC
    or SYN) a threshold of the minimum sample size to accept and return a list
    with derived allele frequencies in the KSR and PX strains combined
    Precondition: codons with more than 1 SNP are not considered
    '''
    
    # create a list to store the frequencies
    DAF = []    
        
    # get the set of valid transcripts
    transcripts = get_valid_transcripts(unique_transcripts)    
    
    # get the set of transcripts with indels in their CDS
    indels = get_genes_with_indels(indel_transcripts)
    
    # create a set of valid nucleotides
    valid_base = {'A', 'T', 'C', 'G'}    
        
    # open file for reading
    infile = open(snp_file, 'r')
    # skip header
    infile.readline()
    # loop over file
    for line in infile:
        line = line.rstrip()
        if line != '':
            line = line.split()
            # check that gene is in valid transcripts
            gene = line[2]
            if gene in transcripts:
                # check only valid snps,
                if line[6] == 'snp' and line[5] != line[7] and line[5] in valid_base and line[7] in valid_base:
                    # check that ancestral is valid base
                    if line[19] in valid_base:
                        # check that ancestral allele is either reference or alternative allele
                        if line[19] == line[5] or line[19] == line[7]:
                            # direction of mutation can be defined
                            # get ref codon
                            ref_codon = line[3]
                            # get alternative codon
                            alt_codon = line[8]
                            # get allele counts
                            ref_count = int(line[13])
                            alt_count = int(line[17])
                            # get alleles
                            ref_allele = line[5]
                            alt_allele = line[7]
                            cla_allele = line[19]
                            # do not consider codons with N
                            if 'N' not in ref_codon and 'N' not in alt_codon:
                                # check that sample size is greater than threshold
                                if ref_count + alt_count >= threshold:
                                    # verify that SNP is in KSR+PX
                                    if alt_count != 0 and ref_count != 0:
                                        # SNP present in KSR+PX
                                        # check that ref and alt codon have at most 1 difference
                                        if diff_codon(ref_codon, alt_codon) <= 1:
                                            # don't count reference stop codons
                                            if genetic_code[ref_codon] != '*':
                                                # reference codon is not terminal stop
                                                # check site_type
                                                if site_type == 'PTC':
                                                    # check that gene does not have indels
                                                    if gene not in indels:
                                                        # count premature stops
                                                        # check that alternative allele is stop codon
                                                        if genetic_code[alt_codon] == '*':
                                                            # mutation results in stop
                                                            # check which allele is derived
                                                            if cla_allele != alt_allele:
                                                                # alternative allele is derived
                                                                # add DAF to list
                                                                DAF.append(alt_count / (ref_count + alt_count))
                                                            elif cla_allele != ref_allele:
                                                                # reference allele is derived
                                                                # add DAF to list
                                                                DAF.append(ref_count / (ref_count + alt_count))
                                                elif site_type == 'REP':
                                                    # count replacement site
                                                    # verify that change is non-synonymous but not PTC
                                                    if genetic_code[ref_codon] != genetic_code[alt_codon] and genetic_code[alt_codon] != '*':
                                                        # site is nonsynonymous
                                                        # check which allele is derived
                                                        if cla_allele != alt_allele:
                                                            # alternative allele is derived
                                                            # add DAF to list
                                                            DAF.append(alt_count / (alt_count + ref_count))
                                                        elif cla_allele != ref_allele:
                                                            # reference allele is derived
                                                            # add DAF to list
                                                            DAF.append(ref_count / (alt_count + ref_count))
                                                elif site_type == 'SYN':
                                                    # count synonymous sites
                                                    # verify that site is synonymous
                                                    if genetic_code[ref_codon] == genetic_code[alt_codon]:
                                                        # site is synonymous
                                                        # check which allele is derived
                                                        if cla_allele != alt_allele:
                                                            # alternative allele is derived
                                                            # add DAF to list
                                                            DAF.append(alt_count / (ref_count + alt_count))
                                                        elif cla_allele != ref_allele:
                                                            # reference allele is minor
                                                            # add minor frequency to list
                                                            DAF.append(ref_count / (ref_count + alt_count))
                                       
    # close file after reading
    infile.close()
    
    return DAF
    
# use this funtion to find the 5' most upstream PTC SNPs
def position_first_PTC(snp_file, unique_transcripts, indel_transcripts, CDS_fasta):
    '''
    (file, file, file, file) -> list
    Take the file with snps in CDS, the file with unique transcripts and the
    file with transcripts having indels and the remanei fasta CDS sequences and
    return a list with the relative position of the 5' most upstream premature
    in the CDS
    Note: PTC with > 1 difference from reference allele are considered in gene counts    
    '''
    
    # convert fasta file to dict
    CDS = convert_fasta(CDS_fasta)
    
    # get the set of valid transcripts
    transcripts = get_valid_transcripts(unique_transcripts)    
    
    # get the set of transcripts with indels in their CDS
    indels = get_genes_with_indels(indel_transcripts)
    
    # create a set of valid nucleotides
    valid_base = {'A', 'T', 'C', 'G'}    
        
    # create a dict {gene : PTC position}
    PTC = {}
    
    # set up gene variable, gene is updated only when a new gene is reached
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
            # check that a new gene is reached
            if line[2] != gene:
                # a new gene is reached, update gene
                gene = line[2]
                # initialize position variable
                i = 0
                # initialize dict
                PTC[gene] = []
            elif line[2] == gene:
                # still in the same gene
                # update position variable
                i += 1
            # check that gene is in valid transcripts and does not have indels
            if gene in transcripts and gene not in indels:
                # check only valid snps
                if line[6] == 'snp' and line[5] != line[7] and line[5] in valid_base and line[7] in valid_base:
                    # get ref codon
                    ref_codon = line[3]
                    # get alternative codon
                    alt_codon = line[8]
                    # verify that reference is not terminal stop codon
                    if genetic_code[ref_codon] != '*':
                        # verify that SNP is in KSR+PX
                        alt_count = int(line[17])
                        if alt_count != 0:
                            # alternative present in KSR+PX
                            # check that alternative codon is stop codon
                            if genetic_code[alt_codon] == '*':
                                # mutation results in stop
                                # record only the 5' most upstream PTC SNP
                                if PTC[gene] == []:
                                    # PTC SNP not found yet, can record it
                                    # add the position of PTC SNP relative to the CDS length
                                    PTC[gene].append(i / len(CDS[gene]))
                                    
    # close file
    infile.close()
    

    # create list to record the position of the 5' most upstream PTC SNP
    PTC_positions = []
    # loop over genes in PTC
    for gene in PTC:
        # if PTC SNP in gene, add position
        if len(PTC[gene]) != 0:
            PTC_positions.append(PTC[gene][0])
            
    return PTC_positions


# use this function to get the frequency of the 5' most upstream PTC SNP
def freq_first_PTC(snp_file, unique_transcripts, indel_transcripts, CDS_fasta):
    '''
    (file, file, file, file) -> list
    Take the file with snps in CDS, the file with unique transcripts and the
    file with transcripts having indels, the remanei fasta CDS sequences and 
    return a list of lists with the mean frequency and standard error of 
    the PTC allele frequency (the alternative allele) located in each decile
    of the CDS length
    Note: PTC with > 1 difference from reference allele are considered in gene counts
    '''
    
    # create a dict of transcript : PTC freq
    PTC = {}
    
    # convert fasta file to dict
    CDS = convert_fasta(CDS_fasta)
    
    # get the set of valid transcripts
    transcripts = get_valid_transcripts(unique_transcripts)    
    
    # get the set of transcripts with indels in their CDS
    indels = get_genes_with_indels(indel_transcripts)
    
    # create a set of valid nucleotides
    valid_base = {'A', 'T', 'C', 'G'}    
    
    # set up gene variable, gene is updated only when a new gene is reached
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
            # check that a new gene is reached
            if line[2] != gene:
                # a new gene is reached, update gene
                gene = line[2]
                # initialize position variable
                i = 0
                # initialize dict
                PTC[gene] = []
            elif line[2] == gene:
                # still in the same gene
                # update position variable
                i += 1
            # check that gene is in valid transcripts and does not have indels
            if gene in transcripts and gene not in indels:
                # check only valid snps
                if line[6] == 'snp' and line[5] != line[7] and line[5] in valid_base and line[7] in valid_base:
                    # get ref codon
                    ref_codon = line[3]
                    # get alternative codon
                    alt_codon = line[8]
                    # verify that reference is not terminal stop codon
                    if genetic_code[ref_codon] != '*':
                        # verify that alternative is in KSR+PX
                        alt_count = int(line[17])
                        ref_count = int(line[13])
                        if alt_count != 0:
                            # alternative present in KSR+PX
                            # check that alternative codon is stop codon
                            if genetic_code[alt_codon] == '*':
                                # mutation results in stop
                                # record only the 5' most upstream PTC SNP
                                if PTC[gene] == []:
                                    # PTC SNP not found yet, can record it
                                    # add the position of PTC SNP relative to the CDS length
                                    PTC[gene].append(i / len(CDS[gene]))
                                    # record the frequency of the alternative allele
                                    PTC[gene].append(alt_count / (ref_count + alt_count))
                                                                                       
    # close file
    infile.close()
    
    # make a list of empty lists to be updated with the frequency values
    # the index in the outer list corresponds to the position of the PTC SNP
    # in % of the CDS length
    i = 10
    PTC_freq = []
    while i != 0:
        PTC_freq.append([])
        i -= 1
    
    # loop over genes in PTC:
    for gene in PTC:
        # check if gene has PTC
        if len(PTC[gene]) == 2:
            # find the index in the list
            # get the relative position of the PTC SNP in % of the CDS length
            i = PTC[gene][0] * 100        
            which_range = int(i) // 10
            PTC_freq[which_range].append(PTC[gene][1])
            
            
    return PTC_freq
            
    

# use this function to get the MAF of PTC SNPs for X-linked and autosomal genes    
def MAF_PTC_chromo(snp_file, unique_transcripts, indel_transcripts, caeno_gff, chromosome_file, threshold):
    '''
    (file, file, file, file, file, int) -> (list, list)
    Take the file of snps, the file with valid transcript names, the file with
    the transcripts having indels, the remanei GFF file and the file with
    remanei-elegans corresponding chromsomes, a threshold of minimum acceptable 
    sample size and return a tuple with lists of minor allele frequencies in
    the KSR and PX strains combined for X-linked and autosomal genes respectively
    Precondition: codons with more than 1 SNP are not considered
    '''
    
    # create lists to store the frequencies
    MAF_X, MAF_auto = [], []    
        
    # get the set of valid transcripts
    transcripts = get_valid_transcripts(unique_transcripts)    
    
    # get the set of transcripts with indels in their CDS
    indels = get_genes_with_indels(indel_transcripts)
    
    # create a set of valid nucleotides
    valid_base = {'A', 'T', 'C', 'G'}    
    
    # get a set of genes with PTC including genes with indels
    PTC_indel_genes = count_genes_with_indels_premature_stops(snp_file, unique_transcripts)
    
    # get a list of genes with PTC SNPs only, excluding genes with indels
    PTC_genes = [gene for gene in PTC_indel_genes if gene not in indels]    
    
    # get lists of PTC X-linked and autosomal genes
    X_genes, autosomal_genes = X_autosomal_genes(PTC_genes, caeno_gff, chromosome_file)
    
    # open file for reading
    infile = open(snp_file, 'r')
    # skip header
    infile.readline()
    # loop over file
    for line in infile:
        line = line.rstrip()
        if line != '':
            line = line.split()
            # check that gene is in valid transcripts
            gene = line[2]
            if gene in transcripts:
                # check only valid snps
                if line[6] == 'snp' and line[5] != line[7] and line[5] in valid_base and line[7] in valid_base:
                    # get ref codon
                    ref_codon = line[3]
                    # get alternative codon
                    alt_codon = line[8]
                    # get allele counts
                    ref_count = int(line[13])
                    alt_count = int(line[17])
                    # check that sample size is greater than threshold
                    if ref_count + alt_count >= threshold:
                        # do not consider codons with Ns
                        if 'N' not in ref_codon and 'N' not in alt_codon:
                            # check that ref and codons have at most 1 difference
                            if diff_codon(ref_codon, alt_codon) <= 1:
                                # verify that SNP is in KSR+PX
                                if ref_count != 0 and alt_count != 0:
                                    # SNP present in KSR+PX
                                    # don't count reference stop codons
                                    if genetic_code[ref_codon] != '*':
                                        # reference codon is not terminal stop
                                        # check that gene is not in indels
                                        if gene not in indels:
                                            # count premature stops
                                            # check that alternative allele is stop codon
                                            if genetic_code[alt_codon] == '*':
                                                # mutation results in stop
                                                # check which allele is minor allele
                                                if alt_count <= ref_count:
                                                    # alternative alelle is minor
                                                    # add minor frequency to list
                                                    # check linkage
                                                    if gene in X_genes:
                                                        # gene is X-linked
                                                        MAF_X.append(alt_count / (alt_count + ref_count))
                                                    elif gene in autosomal_genes:
                                                        # gene is autosomal
                                                        MAF_auto.append(alt_count / (alt_count + ref_count))
                                                elif alt_count > ref_count:
                                                    # reference allele is minor
                                                    # add minor frequency to list
                                                    # check linkage
                                                    if gene in X_genes:
                                                        # gene is X-linked                                            
                                                        MAF_X.append(ref_count / (alt_count + ref_count))
                                                    elif gene in autosomal_genes:
                                                        # gene is autosomal
                                                        MAF_auto.append(ref_count / (alt_count + ref_count))
                                        
    # close file after reading
    infile.close()
    
    return MAF_X, MAF_auto
    
    
# use this function to get the DAF in KSR and PX combined    
def DAF_PTC_chromo(snp_file, unique_transcripts, indel_transcripts, caeno_gff, chromosome_file, threshold):
    '''
    (file, file, file, file, file, int) -> (list, list)
    Take the file of snps, the file with valid transcript names, the file with
    the transcripts having indels, the remanei GFF file and the file with
    remanei-elegans corresponding chromsomes, a threshold of minimum acceptable
    sample size and return a tuple with lists of derived allele frequencies
    in the KSR and PX strains combined for X-linked and autosomal genes respectively
    '''    
    
    # create lists to store the frequencies
    DAF_X, DAF_auto = [], []    
        
    # get the set of valid transcripts
    transcripts = get_valid_transcripts(unique_transcripts)    
    
    # get the set of transcripts with indels in their CDS
    indels = get_genes_with_indels(indel_transcripts)
    
    # create a set of valid nucleotides
    valid_base = {'A', 'T', 'C', 'G'}    
    
    # get a set of genes with PTC including genes with indels
    PTC_indel_genes = count_genes_with_indels_premature_stops(snp_file, unique_transcripts)
    
    # get a list of genes with PTC SNPs only, excluding genes with indels
    PTC_genes = [gene for gene in PTC_indel_genes if gene not in indels]    
    
    # get lists of PTC X-linked and autosomal genes
    X_genes, autosomal_genes = X_autosomal_genes(PTC_genes, caeno_gff, chromosome_file)    
    
    # open file for reading
    infile = open(snp_file, 'r')
    # skip header
    infile.readline()
    # loop over file
    for line in infile:
        line = line.rstrip()
        if line != '':
            line = line.split()
            # check that gene is in valid transcripts
            gene = line[2]
            if gene in transcripts:
                # check only valid snps
                if line[6] == 'snp' and line[5] != line[7] and line[5] in valid_base and line[7] in valid_base:
                    # check that ancestral is valid base
                    if line[19] in valid_base:
                        # check that ancestral allele is either reference or alternative allele
                        if line[19] == line[5] or line[19] == line[7]:
                            # direction of mutation can be defined
                            # get ref codon
                            ref_codon = line[3]
                            # get alternative codon
                            alt_codon = line[8]
                            # get allele counts
                            ref_count = int(line[13])
                            alt_count = int(line[17])
                            # get alleles
                            ref_allele = line[5]
                            alt_allele = line[7]
                            cla_allele = line[19]
                            # check that sample size is greater than trheshold
                            if ref_count + alt_count >= threshold:
                                # do not consider codons with Ns
                                if 'N' not in ref_codon and 'N' not in alt_codon:
                                    # check that ref and alt codons have at most 1 difference
                                    if diff_codon(ref_codon, alt_codon) <= 1:
                                        # check that SNP is in KSR + PX
                                        # verify that SNP is in KSR+PX
                                        if ref_count != 0 and alt_count != 0:
                                            # SNP present in KSR+PX
                                            # don't count reference stop codons
                                            if genetic_code[ref_codon] != '*':
                                                # reference codon is not terminal stop
                                                # check that gene does not have indels
                                                if gene not in indels:
                                                    # count premature stops
                                                    # check that alternative allele is stop codon
                                                    if genetic_code[alt_codon] == '*':
                                                        # mutation results in stop
                                                        # check which allele is derived
                                                        if cla_allele != alt_allele:
                                                            # alternative allele is derived
                                                            # add DAF to list
                                                            # check linkage
                                                            if gene in X_genes:
                                                                DAF_X.append(alt_count / (ref_count + alt_count))
                                                            elif gene in autosomal_genes:
                                                                DAF_auto.append(alt_count / (ref_count + alt_count))
                                                        elif cla_allele != ref_allele:
                                                            # reference allele is derived
                                                            # add DAF to list
                                                            # check linkage
                                                            if gene in X_genes:
                                                                DAF_X.append(ref_count / (ref_count + alt_count))
                                                            elif gene in autosomal_genes:
                                                                DAF_auto.append(ref_count / (ref_count + alt_count))
                                                
    # close file after reading
    infile.close()
    
    return DAF_X, DAF_auto
