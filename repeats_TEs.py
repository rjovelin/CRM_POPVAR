# -*- coding: utf-8 -*-
"""
Created on Fri Jul 31 13:05:22 2015

@author: RJovelin
"""

from manipulate_sequences import *

# use this function to get the coordinates of TEs in remanei
def get_repeats_coord(repeatmasker_output, remove_simple_repeat):
    '''
    (file, bool) -> dict
    Take the RepeatMasker output file and a boolean specifying if simple repeats
    and low complexity regions should be remove (=True) or not (=False) and
    return a dictionary with repeat name as key a list of list coordinates
    for all the given repeats in the genome 
    '''
    
    # get the family: repeat pairs    
    fam = repeat_family(repeatmasker_output)
    
    # create dictionary {repeat : [[chromo, start, end, orientation], [chromo, start, end, orientation]]}
    repeat_coord = {}
    
    # open file for reading
    infile = open(repeatmasker_output, 'r')
    # skip header
    i = 3
    while i != 0:
        infile.readline()
        i -= 1
    # loop over lines in file
    for line in infile:
        # strip line from end of line and from white space at the begining
        line = line.strip()
        if line != '':
            line = line.split()
            # get chromo
            chromo = line[4]
            # get start and end positions, 0-based
            start = int(line[5]) -1
            end = int(line[6])
            # get orientation
            if line[8] == '+':
                orientation = '+'
            elif line[8] == 'C':
                orientation = '-'
            else:
                return -1
            # get repeat, TE name
            repname = line[9]
            # check if repname is key in dict
            if repname in repeat_coord:
                # add coordinates to list
                repeat_coord[repname].append([chromo, start, end, orientation])
            else:
                repeat_coord[repname] = [[chromo, start, end, orientation]]
                
    # close file after reading
    infile.close()
    
    # check if simple repeats need to be removed
    if remove_simple_repeat == True:
        # remove simple repeats and low complexity repeats
        to_remove = []
        for repname in repeat_coord:
            if repname in fam['Simple_repeat']:
                to_remove.append(repname)
            elif repname in fam['Low_complexity']:
                to_remove.append(repname)
        for i in to_remove:
            del repeat_coord[i]
    
    return repeat_coord
            
  
# use this function to get all the indices of repeat positions on each chromo
def get_repeat_positions(repeat_coord):
    '''
    (dict) -> dict
    Take the dictionary of repeat name : list of repeat coordinates
    and return a dictionary with chromo as key and the set of all indices
    corresponding to the repeats on that chromo
    '''
    
    # create a dict {chromo: {set of indices}}
    repeat_pos = {}
    
    # loop over repeat names
    for repname in repeat_coord:
        # loop over the coordinate l;ist for each repeat of that name
        for i in range(len(repeat_coord[repname])):
            # get chromo
            chromo = repeat_coord[repname][i][0]
            # get start and end positions (0-based)
            start = repeat_coord[repname][i][1]
            end = repeat_coord[repname][i][2]
            # check if chromo in repeat_pos
            if chromo not in repeat_pos:
                # initiate set
                repeat_pos[chromo] = set()
            # loop over positions
            for j in range(start, end):
                repeat_pos[chromo].add(j)
                
    return repeat_pos
                
    
  
# use this function to get the family names of each repeat            
def repeat_family(repeatmasker_output):
    '''
    (file) -> dict
    Take the RepeatMasker output file and return a dictionary with the repeat
    class and a set of repeat names from the given class
    '''
    
    # create dictionary {repeat : [[chromo, start, end, orientation], [chromo, start, end, orientation]]}
    repeat_fam = {}
    
    # open file for reading
    infile = open(repeatmasker_output, 'r')
    # skip header
    i = 3
    while i != 0:
        infile.readline()
        i -= 1
    # loop over lines in file
    for line in infile:
        # strip line from end of line and from white space at the begining
        line = line.strip()
        if line != '':
            line = line.split()
            # get repeat, TE name
            repname = line[9]
            # get the repeat class
            repclass = line[10]
            # check if  key in dict
            if repclass in repeat_fam:
                repeat_fam[repclass].append(repname)
            else:
                repeat_fam[repclass] = [repname]
                
    # close file after reading
    infile.close()
    
    return repeat_fam

    
# use this function to get the sequences of all repeats  
def grab_repeat_sequences(repeatmasker_output, genome_fasta):
    '''
    (file, file) -> dict
    Take the RepeatMasker outputfile, the genome fasta file and return a
    dictionary with repeat name as key and its sequence as value
    '''
    
    # convert the genome to dict
    genome = convert_fasta(genome_fasta)
    
    # get the coordinates of the repeats
    repeat_coord = get_repeats_coord(repeatmasker_output)
    
    # create a dict {repeat_name : [seq1, seq2,seq3]}
    repeat_seq = {}
    
    # loop over repeats
    for repeat in repeat_coord:
        # check if the repeat has multiple instances
        for i in range(len(repeat_coord[repeat])):
            # extract the sequence
            chromo = repeat_coord[repeat][i][0]
            # get the start, end positions
            start = repeat_coord[repeat][i][1]
            end = repeat_coord[repeat][i][2]
            # get orientation
            sens = repeat_coord[repeat][i][3]
            seq = genome[chromo][start:end].upper()
            if sens == '-':
                # take reverse complement
                seq = reverse_complement(seq)
            # check if key in dict
            if repeat in repeat_seq:
                repeat_seq[repeat].append(seq)
            else:
                repeat_seq[repeat] = [seq]
                
    return repeat_seq
                
    
    
    