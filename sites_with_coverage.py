# -*- coding: utf-8 -*-
"""
Created on Mon Aug  3 20:30:31 2015

@author: Richard
"""

import os
from miRNA_target import *
from piRNAs import *
from genomic_coordinates import *
from manipulate_sequences import *



# use this function to get allele counts at all sites with coverage in the genome
def get_non_coding_snps(directory, threshold):
    '''
    (str, str) -> dict
    Take the directory where the SNP files, a threshold of minimum sample size to
    accept, and return a dictionary with chromo as key to an inner directory
    of position as key and list with ref allele, alt allele and their counts
    in the PX+KSR sample as value    
    '''
    
    # make a list of files
    files = [i for i in os.listdir(directory) if 'PXKSR' in i]
    
    # create a dict [chromo: {site : [ref_allele, alt_allele, count_ref, count alt]}]
    chromo_sites = {}
    
    # loop of filename in file list
    for filename in files:
        # open file
        infile = open(directory + filename, 'r')
        # skip header
        infile.readline()
        # loop over file
        for line in infile:
            line = line.rstrip()
            if line != '':
                line = line.split()
                # get chromo
                chromo = line[0]
                # convert position to 0-based index
                pos = int(line[1]) -1
                # get reference allele
                ref = line[2]
                alt = line[3]
                # get ref and alt counts in PX+KSR
                ref_count = int(line[7])
                alt_count = int(line[10])
                # populate dict
                # check if sample size is greater than threshold
                if ref_count + alt_count >= threshold:
                    # check if chromo in dict
                    if chromo in chromo_sites:
                        # add the position as key
                        chromo_sites[chromo][pos] = [ref, alt, ref_count, alt_count]
                    else:
                        chromo_sites[chromo] = {}
                        chromo_sites[chromo][pos] = [ref, alt, ref_count, alt_count]
        # close file after reading
        infile.close()
        
    return chromo_sites
  


# use this function to get all the indices in genome corresponding to pirnas or mirnas
def get_small_rna_sites(rna_coord):
    '''
    (dict) -> dict
    Take a dict with coordinates of small RNAs (mirnas or pirnas) and return
    a dictionary with chromo as key and a set of indices corresponding to all
    the positions of piRNAs or miRNAs on each chromo
    Precondition: All position are 0-based indices
    '''
    
    # rna_coord is a dict in the form {chromo: [[start, end, orienation]]}
    # start, end are 0-based
        
    # create a dict to store the pirna indices
    rna_pos = {}

    # loop over chromo in rna_coord
    for chromo in rna_coord:
        # loop over rna on chromo
        for i in range(len(rna_coord[chromo])):
            # get start, end position
            start = rna_coord[chromo][i][0]
            end = rna_coord[chromo][i][1]
            # check if chromo in rna_pos
            if chromo in rna_pos:
                # add all the indices corresponding to each rna
                for j in range(start, end):
                    rna_pos[chromo].add(j)
            else:
                # add chromo as key, and create value set
                rna_pos[chromo] = set()
                # add all the indices corressponding to each pirna
                for j in range(start, end):
                    rna_pos[chromo].add(j)
    return rna_pos
    
    
# use this function to get all the indices in genome corresponding to predicted UTRs
def get_UTR_sites(UTR_coord):
    '''
    (file, file, file, int) -> dict
    Take the dictionary of UTR coordinates and return a dictionary with chromo
    as key and a set of indices correspionding to the positions of remanei UTR sites.
    Precondition: All position are 0-based indices
    '''
        
    # create a dict UTR_pos
    UTR_pos = {}
    # loop over genes in three_prime
    for gene in UTR_coord:
        # get chromo
        chromo = UTR_coord[gene][0]
        # convert to 0-based
        start = UTR_coord[gene][1] -1
        end = UTR_coord[gene][2]
        # check if chromo in UTR_pos
        if chromo in UTR_pos:
            for j in range(start, end):
                UTR_pos[chromo].add(j)
        else:
            UTR_pos[chromo] = set()
            for j in range(start, end):
                UTR_pos[chromo].add(j)
    return UTR_pos


# use this function to get all the indices in genome corresponding to genes
def get_gene_sites(caeno_gff):
    '''
    (file) -> dict
    Take the remanei gff file and return a dictionary with chromo as key and
    a set of indices corresponding to all the positions of remanei genes on each chromo
    Precondition: All position are 0-based indices
    '''
    
    # get coding gene coordinates
    # create ditc {gene1 : [chromo, start, end, sense]}
    genes_coord = {}

    # open annotation file for reading
    infile = open(caeno_gff, 'r')
    for line in infile:
        line = line.rstrip()
        if line != '':
            line = line.split()
            if len(line) > 8:
                if line[2] == 'mRNA':
                    # get coordinates
                    chromo = line[0]
                    # convert to 0-based
                    start = int(line[3]) -1
                    end = int(line[4])
                    sense = line[6]
                    transcript = line[8][line[8].index('ID=')+3 : line[8].index(';')]
                    genes_coord[transcript] = [chromo, start, end, sense]
    infile.close()
    
    # create a dict with genic sites for all chromos
    gene_pos = {}
    # loop over genes:
    for gene in genes_coord:
        # get coordinates
        chromo = genes_coord[gene][0]
        start = genes_coord[gene][1]
        end = genes_coord[gene][2]
        # populate dict
        if chromo in gene_pos:
            for j in range(start, end):
                gene_pos[chromo].add(j)
        else:
            gene_pos[chromo] = set()
            for j in range(start, end):
                gene_pos[chromo].add(j)
                
    return gene_pos
            
             
    
# use this function to get allele counts for intergenic sites only    
def get_intergenic_sites(chromo_sites, gene_sites, pirna_sites, mirna_sites, UTR_sites, remove_pirna_mirna_utr):
    '''
    (dict, dict, dict, dict, dict, bool) -> dict
    Take the dictionnary with allele counts for all sites with coverage in the
    genome, the dictionaries with all sites corresponding to genes, piRNAs, miRNAs
    and predicted UTRs, a boolean to specify whether pirnas, mirna and UTRs
    should be removed (=True) or whether only the sites between annotated genes
    should be removed (=False) and return a modified dictionary with allele
    counts for intergenic sites only
    Precondition: all indices in all dicts are 0-based    
    '''
    
    # remove all sites falling at protein coding gene loci
    # loop over chromo in genic_sites
    for chromo in gene_sites:
        # check that chromo in chromo sites
        if chromo in chromo_sites:
            # loop over gene positions on that chromo
            for i in gene_sites[chromo]:
                # check if site is on chromo in chromo_sites
                if i in chromo_sites[chromo]:
                    del chromo_sites[chromo][i]
                    
    # check if all intergenic sites are kept or of small RNAs and UTR are removed    
    if remove_pirna_mirna_utr == True:
        # do not consider intergenic sites falling at pirna and mirna loci and UTRs
        
        # remove all sites corresponding to piRNAs
        for chromo in pirna_sites:
            # check if chromo in chromo_sites
            if chromo in chromo_sites:
                # loop over pirna positions on that chromo
                for i in pirna_sites[chromo]:
                    # check if site on chromo in chromo-sites
                    if i in chromo_sites[chromo]:
                        # remove site
                        del chromo_sites[chromo][i]
        
        # remove all sites corresponding to mirnas
        for chromo in mirna_sites:
            # check if chromo in chromo sites
            if chromo in chromo_sites:
                # loop over mirna positions on that chromo
                for i in mirna_sites[chromo]:
                    # check if site on chromo in chromo_sites
                    if i in chromo_sites[chromo]:
                        # remove site
                        del chromo_sites[chromo][i]
                        
        # remove sites coreresponding to predicted UTRs
        for chromo in UTR_sites:
            # check if chromo in chromo_sites:
            if chromo in chromo_sites:
                # loop over positions on that chromo
                for i in UTR_sites[chromo]:
                    # check if site on that chromo in chromo_sites
                    if i in chromo_sites[chromo]:
                        # remove site
                        del chromo_sites[chromo][i]
                        
    return chromo_sites


# use this function to get the allele counts for sites in small RNAs (mirnas or pirnas)
def get_feature_sites(chromo_sites, feature_coord):
    '''
    (dict, dict) -> dict
    Take the dictionnary with allele counts for all sites with coverage in the
    genome, and a dictionnary with the coordinates of a given feature
    (ie. miRNA, piRNA, target sites etc) and return a dictionnary
    with allele counts for the feature sites only
    Precondition: feature_coord is in the form {chromo: [[start, end , orientation]]}
    with positions in 0-based indices
    '''    
        
    # create a dict to get the allele counts at feature sites
    feature_sites = {}    
    
    # loop over chromo in feature coord
    for chromo in feature_coord:
        # check that chromo is key in chromo sites
        if chromo in chromo_sites:
            # loop over features on that chromo
            for i in range(len(feature_coord[chromo])):
                # get the positions of the feature
                for j in range(feature_coord[chromo][i][0], feature_coord[chromo][i][1]):
                    # check if position is on chromo
                    if j in chromo_sites[chromo]:
                        # check if chromo in feature sites
                        if chromo in feature_sites:
                            # add pos : list coordinates pair to dict
                            feature_sites[chromo][j] = list(chromo_sites[chromo][j])
                        else:
                            feature_sites[chromo] = {}
                            feature_sites[chromo][j] = list(chromo_sites[chromo][j])
                        
                        
    return feature_sites
    


# use this function to get the MAF in KSR and PX combined     
def MAF_non_coding(feature_sites):
    '''
    (dict) -> list
    Take a dictionnary with allele counts at positions of a given feature
    and return a list with minor allele frequencies of all the SNPs at that feature
    Precondition: do not use this function for synonymous and replacement sites
    '''
  
    # create a dict [chromo: {site : [ref_allele, alt_allele, count_ref, count alt]}]
    
     # create a list to store the frequencies
    MAF = []    
        
    # loop over chromo in dict
    for chromo in feature_sites:
        # loop over the positions on that chromo
        for site in feature_sites[chromo]:
            # check if site is polymorphic
            # get allele counts
            ref_count = feature_sites[chromo][site][2]
            alt_count = feature_sites[chromo][site][3]
            # get alelles
            ref = feature_sites[chromo][site][0]
            alt = feature_sites[chromo][site][1]
            # check if site is polymorphich
            if ref_count != 0 and alt_count != 0:
                # site is polymorphich
                # double check 
                assert ref != alt, 'allele counts are different than 0, but nucleotides are the same'
                # check which allele is minor
                if alt_count <= ref_count:
                    # alternative alelle is minor
                    # add minor frequency to list
                    MAF.append(alt_count / (ref_count + alt_count))
                elif alt_count > ref_count:
                    # reference allele is minor
                    # add minor frequency to list
                    MAF.append(ref_count / (ref_count + alt_count))
                    
    return MAF    
    
    
    
# use this function to get allele counts at all sites with coverage in all strains
def get_all_strains_snps(directory, threshold):
    '''
    (str, str) -> dict
    Take the directory where the SNP files, a threshold of minimum sample size to
    accept, and return a dictionary with chromo as key to an inner directory
    of position as key and list with ref allele, alt allele and their counts
    in the PX+KSR+PB sample as value    
    '''
    
    # make a list of files
    files = [i for i in os.listdir(directory) if 'PBON' in i and 'SNPs' in i]
    
    # create a dict [chromo: {site : [ref_allele, alt_allele, count_ref, count alt]}]
    chromo_sites = {}
    
    # loop of filename in file list
    for filename in files:
        # open file
        infile = open(directory + filename, 'r')
        # skip header
        infile.readline()
        # loop over file
        for line in infile:
            line = line.rstrip()
            if line != '':
                line = line.split()
                # get chromo
                chromo = line[0]
                # convert position to 0-based index
                pos = int(line[1]) -1
                # get reference allele
                ref = line[2]
                alt = line[3]
                # get ref and alt counts in PX+KSR + PB
                ref_count = int(line[7]) + int(line[8])
                alt_count = int(line[11]) + int(line[12]) 
                # populate dict
                # check if sample size is greater than threshold
                if ref_count + alt_count >= threshold:
                    # check if chromo in dict
                    if chromo in chromo_sites:
                        # add the position as key
                        chromo_sites[chromo][pos] = [ref, alt, ref_count, alt_count]
                    else:
                        chromo_sites[chromo] = {}
                        chromo_sites[chromo][pos] = [ref, alt, ref_count, alt_count]
        # close file after reading
        infile.close()
        
    return chromo_sites   
    
# use this function to count the total number os SNPs    
def count_total_snps(chromo_sites):
    '''
    (dict) -> int
    Take the dictionnary with allele counts for all sites with coverage in the
    genome, and return the total number of SNPs
    Precondition: SNPs have been filtered or not according to sample size    
    and sites in repeat regions have been removed
    '''
    
    # create counter variable
    total = 0
    
    # loop ovver chromo
    for chromo in chromo_sites:
        # loop over site on chromo
        for i in chromo_sites[chromo]:
            # get ref and alt counts
            ref_count = chromo_sites[chromo][i][2]
            alt_count = chromo_sites[chromo][i][3]
            ref = chromo_sites[chromo][i][0]
            alt = chromo_sites[chromo][i][1]
            # check if site is polymorphic
            if ref_count != 0 and alt_count != 0:
                # site is polymorphic
                assert ref != alt, 'ref anc alt counts are different but ref and alt alleles are the same'
                total += 1
                
    return total
    
    
# use this function to count the number of shared SNPs between PB abd ON   
def count_shared_snps_PB_ON(directory):
    '''
    (str, int) -> int
    Take the directory where the SNP files, a threshold of minimum sample size to
    accept, and return the number of SNPs shared between PB and KSR+PX strains
    (ie. sites that are SNPs in PB and in ON)    
    '''
    
    # make a list of files
    files = [i for i in os.listdir(directory) if 'PBON' in i and 'SNPs' in i]
    
    # set up counter variable
    shared = 0
    
    # loop of filename in file list
    for filename in files:
        # open file
        infile = open(directory + filename, 'r')
        # skip header
        infile.readline()
        # loop over file
        for line in infile:
            line = line.rstrip()
            if line != '':
                line = line.split()
                # get reference allele
                ref = line[2]
                alt = line[3]
                # get ref and alt counts in PX+KSR
                ref_count_ON = int(line[7])
                alt_count_ON = int(line[11]) 
                # get ref and alt counts in PB
                ref_count_PB = int(line[8])
                alt_count_PB = int(line[12]) 
                # check if site is polymorphic in ON and in PB
                if (ref_count_ON != 0 and alt_count_ON != 0) and (ref_count_PB != 0 and alt_count_PB != 0):
                    # site is polymorphic in PB abd in OB
                    assert ref != alt, 'ref and alt counts are different than 0 but ref is same as alt'
                    # update counter
                    shared += 1
        # close file after reading
        infile.close()
    return shared
                
    