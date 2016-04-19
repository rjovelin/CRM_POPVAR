# -*- coding: utf-8 -*-
"""
Created on Thu May 28 21:23:49 2015

@author: Richard
"""

#!/usr/bin/env python3


from accessories import *


# use this function to slice the 7bp seed from the mature sequence 
def grab_seed(mature_seq):
    '''
    (str) -> str
    Slice the miRNA mature sequence to return the 7bp seed motif
    '''
    seed = mature_seq[1:8]
    return seed


# use this function to generate a dict of seeds and list pf mirnas pairs
def seed_mirnas(mature_fasta):
    '''
    (file) -> dict
    Take a fasta file of miRNA mature sequences and return a dictionnary with
    seed as key and a list of mirnas from the same family (sharing the same seed)
    as value
    '''
    
    # convert fasta file to 
    mirnas = convert_fasta(mature_fasta)
    
    # create a dict of seed : [mir1, mir2]
    seeds = {}
    
    # loop over mature sequences
    for name in mirnas:
        # get seed sequence
        seed_seq = grab_seed(mirnas[name])
        # clean up the mirna name
        mirna_name = name[4:name.index(' ')]
        # populate dict
        if seed_seq in seeds:
            seeds[seed_seq].append(mirna_name)
        else:
            seeds[seed_seq] = [mirna_name]
            
    return seeds
    
    
# use this function to generate the remanei miRNA family input file for TargetScan
def generate_crem_mirfam_targetscan(mature_fasta, outputfile):
    '''
    (file) -> file
    Take a fasta file of miRNA matire sequences and generate the miRNA family
    input file for TargetScan    
    '''
    
    # create a dict of seed sequence and list of mirna pairs
    seeds = seed_mirnas(mature_fasta)
    
    # open file for writing
    newfile = open(outputfile, 'w')
    
    # loop over the seed sequences
    for motif in seeds:
        # record only 1 mirna per family
        # brag first mirna
        # write mirna, seed, species ID
        newfile.write(seeds[motif][0] + '\t' + motif + '\t' + '31234' + '\n')
    
    # close file after writing
    newfile.close()


# use this function to generate a set of seeds for a given species, or group of mature sequences
def get_all_seeds_in_species(mature_fasta):
    '''
    (file) -> set
    Take a fasta file of miRNA mature sequences and return a set of miRNA 
    seed sequences present in a given species
    '''
    
    # create a dict of seed sequence and list of mirna pairs
    seeds = seed_mirnas(mature_fasta)
    
    # create a set of seeds present in species
    seed_species = {i for i in seeds}
    
    return seed_species


# use this function to check if a miRNA family is conserved in a group of miRNAs or species
def is_miRNA_family_conserved(seed_seq, seeds_species):
    '''
    (str, set) -> bool
    Return True if the seed sequence is conserved and present in the set of 
    seeds from species, return False otherwise
    '''
    
    if seed_seq in seeds_species:
        return True
    else:
        return False
    
    
    
    
    
    
