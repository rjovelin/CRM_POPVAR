# -*- coding: utf-8 -*-
"""
Created on Fri Jun 12 16:45:23 2015

@author: Richard
"""

#!/usr/bin/env python3

# import modules
from get_coding_sequences import *
from miRNA_target import *
from accessories import *


# use this function to make a dict with the coordinates of all indels on each chromosome
def get_short_indel_coordinates(short_indel_file):
    '''
    (file) -> dict
    Take the compiled pindel outputfile and return a dictionnary with
    chromosome as key and list of tuples with start and end position of each indel
    '''
    
    # create dict with {chromo: [(start, end), (start, indel)]}
    indel_coord = {}    
    
    # open file for reading
    infile = open(short_indel_file, 'r')
    
    for line in infile:
        if line.rstrip() != '':
            line = line.rstrip().split()
            chromo = line[1]
            # get positions 0-based
            start = int(line[2]) -1
            end = int(line[3])
            if chromo in indel_coord:
                indel_coord[chromo].append((start, end))
            else:
                indel_coord[chromo] = [(start, end)]
    # close file
    infile.close()
    
    return indel_coord 

# use this function to remove indels that overlap with repeats
# before counting indels
def remove_indels_in_repeats(indel_coord, repeat_positions):
    '''
    (dict, dict) -> dict
    Take the dictionary with chromo: list of indel coordinates and a dictionary with chromo:
    set of repeat positions and return a modified dictionary in which indels
    overlapping with repeats have been removed
    '''
    
    # loop over chromo in indel coord
    for chromo in indel_coord:
        # create a list with idnels to remove
        to_remove = []
        # check if chromo in repeats
        if chromo in repeat_positions:
            # loop over indels in chromo
            for i in range(len(indel_coord[chromo])):
                # check if indel overlapps with repeat positions
                # get the indel positions
                indel_pos = set(range(indel_coord[chromo][i][0], indel_coord[chromo][i][1]))
                # check if indel overlapps with repeat
                if indel_pos.intersection(repeat_positions[chromo]):
                    # add indel index to list
                    to_remove.append(i)
            # check if indels need to be removed
            if len(to_remove) != 0:
                # reverse sort list
                to_remove.sort()
                to_remove.reverse()
                # loop over indel indices starting at the end of the list of indel coordinates
                for i in to_remove:
                    indel_coord[chromo].remove(indel_coord[chromo][i])
                    
    return indel_coord
            
                
    
# use this function to count all the indels in a group of strains
# indels may be overlapping
def count_short_indels(indels):
    '''
    (dict) -> int
    Take a dictionnary with chromo as key and a list of indel coordinates as 
    value and return the number of indels    
    '''
    
    # set up counter
    total = 0    
    
    # loop over chromo in dict
    for chromo in indels:
        total += len(indels[chromo])
    
    return total
        
    
# use this function the get the length of indels        
def indel_length(indels):
    '''
    (dict) -> list
    Take a dictionnary of chromo as key and list of indel positions as value
    and return a list of indel size in bp
    '''    
    
    # create a set to store the indel length
    size = set()
    
    # loop over chromo
    for chromo in indels:
        # loop over indels
        for i in range(len(indels[chromo])):
            # compute length, start and end are 0-based
            indel_length = indels[chromo][i][1] - indels[chromo][i][0]
            size.add(indel_length)
            
    return list(size)
    
# use this function to create a dict with unique indels
def combine_indels(indel_coord1, indel_coord2):
    '''
    (dict, dict) -> dict
    Take 2 dicttionnaries with indel coordinates and combine them into a new
    dictionnary with coordinates of unique indels
    '''
    
    # create a new dict with {chromo: {(start, end), (start, end)}}
    new_coord = {}
    
    # loop over chromo of dict1, populate new dict
    for chromo in indel_coord1:
        new_coord[chromo] = set()
        for indel in indel_coord1[chromo]:
            new_coord[chromo].add(indel)
    # loop over chromo in dict2    
    for chromo in indel_coord2:
        # check if chromo is in new dict
        if chromo in new_coord:
            # key in dict, add all indels
            for indel in indel_coord2[chromo]:
                new_coord[chromo].add(indel)
        else:
            # key not in dict, add key and initiate set
            new_coord[chromo] = set()
            for indel in indel_coord2[chromo]:
                new_coord[chromo].add(indel)

    return new_coord


# use this function to find the transcripts with indels in their coding sequences
def identify_genes_with_indels(caeno_gff, short_indel_file):
    '''
    (file, file) -> set
    Take the remanei GFF file and the pindel outputfile with coordinates of
    short_indels and return a set of transcripts with indels in their coding
    sequence
    '''
   
    # get the positions of each transcripts
    # {gene1 : [chromo, start, end, sense]}
    # positions are 1-based
    genes = get_genes_coordinates(caeno_gff)
    
    # get the CDS coordinates
    # {TS1: [chromo, sense, [(s1, end1), (s2, end2)]]}
    # positions are 1-based
    CDS = get_CDS_positions(caeno_gff)

    # get the coordinates of short indels
    # positions are 0-based
    indels = get_short_indel_coordinates(short_indel_file)
    
    # create a set of transcripts with indel in the CDS
    CDS_indels = set()
    
    # loop over transcripts in gene coordinates
    for transcript in genes:
        # set up boolean for each new gene
        indel_found = False        
        # get chromo
        chromo = genes[transcript][0]
        # ask if indels are in chromo
        if chromo in indels:
            # get indices of the indels, make start and end 0-based
            gene_coord = set(range(genes[transcript][1] -1, genes[transcript][2]))
            # check all indels until an indel in CDS is found
            for i in range(len(indels[chromo])):
                # if indel not found, loop over indels
                if indel_found == False:
                    # indel positions are 0-based
                    indel_coord = set(range(indels[chromo][i][0], indels[chromo][i][1]))
                    if len(gene_coord.intersection(indel_coord)) != 0:
                        # indel is within gene
                        for j in range(len(CDS[transcript][2])):
                            # make startand end 0-based
                            CDS_coord = set(range(CDS[transcript][2][j][0] -1, CDS[transcript][2][j][1]))
                            # ask if indel is within CDS
                            if len(CDS_coord.intersection(indel_coord)) != 0:
                                # indel is within CDS, add transcript to set
                                CDS_indels.add(transcript)
                                # update boolean so that next indel is not checked
                                indel_found = True
                                # exit inner loop, no need to check all exons
                                break
                            
    return CDS_indels


# use this function to count the number of genes with indels in CDS and
# the number of indels affecting CDS only for the set of transcripts used in analyses of diversity
def indels_in_CDS(caeno_gff, unique_transcripts, indel_coord):
    '''
    (file, file, dict) -> dict
    Take the remanei GFF file, a dictionnary of chromosome: list of indel positions,
    the file with valid transcript names and a return a dict of transcript: list of
    indel positions for indels affecting the CDS of transcripts    
    '''
    
    # get the set of valid transcripts
    transcripts = get_valid_transcripts(unique_transcripts)

    # get the positions of each transcripts
    # {gene1 : [chromo, start, end, sense]}
    # positions are 1-based
    genes = get_genes_coordinates(caeno_gff)
    
    # get the CDS coordinates
    # {TS1: [chromo, sense, [(s1, end1), (s2, end2)]]}
    # positions are 1-based    
    CDS = get_CDS_positions(caeno_gff)

    # create a dict of transcripts with indel in the CDS
    CDS_indels = {}
    
    # loop over transcripts in gene coordinates
    for transcript in genes:
        # check that transcript is valid
        if transcript in transcripts:
            # get chromo
            chromo = genes[transcript][0]
            # ask if indels are in chromo
            if chromo in indel_coord:
                # make positions 0-based
                gene_pos = set(range(genes[transcript][1] -1, genes[transcript][2]))
                # check all indels to see if they fall in a CDS
                for i in range(len(indel_coord[chromo])):
                    # positions are 0-based
                    indel_pos = set(range(indel_coord[chromo][i][0], indel_coord[chromo][i][1]))
                    if len(gene_pos.intersection(indel_pos)) != 0:
                        # indel is within gene
                        for j in range(len(CDS[transcript][2])):
                            # make positions 0-based
                            CDS_pos = set(range(CDS[transcript][2][j][0] -1, CDS[transcript][2][j][1]))
                            # ask if indel is within CDS
                            if len(CDS_pos.intersection(indel_pos)) != 0:
                                # indel is within CDS, add indel to set
                                if transcript in CDS_indels:
                                    # add 0-based indel positions
                                    CDS_indels[transcript].add((indel_coord[chromo][i][0], indel_coord[chromo][i][1]))
                                else:
                                    CDS_indels[transcript] = set((indel_coord[chromo][i][0], indel_coord[chromo][i][1]))

    return CDS_indels





