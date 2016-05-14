# -*- coding: utf-8 -*-
"""
Created on Sat May 14 10:48:33 2016

@author: Richard
"""

# use this script to find the genomic coordinates of mature miRNAs and to sabe them to file


from manipulate_sequences import *
from genomic_coordinates import *

# convert genome fasta to dict
genome = convert_fasta('noamb_356_v1_4.txt')

# create a dict with hairpin sequences
hairpin = {}
# create a dict with 0-based hairpin genomic coordinates
mirna_coord = {}
# create a dict with mature sequences
mature = {}

# open hairpin coordinates file for reading, populkate dicts
infile = open('CRM_miRNAsCoordinatesFinal.txt')
infile.readline()
for line in infile:
    line = line.rstrip()
    if line != '':
        line = line.split('\t')
        name = line[0]
        mirna = line[5]
        miR = line[6]
        chromo = line[1]
        start = int(line[2]) - 1
        end = int(line[3])
        orientation = line[4]
        mirna_coord[name] = [chromo, start, end, orientation]
        mature[name] = miR
        hairpin[name] = mirna
infile.close()


# create a dict with mature mirna genomic coords
mature_coord = {}

# loop over mature sequences
for mirna in mature:
    # get the coordinates of the mature seq in hairpin seq
    start = hairpin[mirna].index(mature[mirna])
    end = start + len(mature[mirna])
    # get the coordinates of hairpin on chromo
    seqstart, seqend = mirna_coord[mirna][1], mirna_coord[mirna][2]
    # get orientation of hairpin seq on chromo
    orientation = mirna_coord[mirna][3]
    # get chromo
    chromo = mirna_coord[mirna][0]
    # get the genomic coordinates of the mature seq on chromo
    genomic_start, genomic_end = convert_seq_coord_genome_coord(start, end, seqstart, seqend, orientation, len(genome[chromo]))
    mature_coord[mirna] = [chromo, genomic_start, genomic_end, orientation]

# check that mirna coord is coorect
for mirna in mature_coord:
    seq1 = mature[mirna]
    chromo, start, end = mature_coord[mirna][0], mature_coord[mirna][1], mature_coord[mirna][2]
    orientation = mature_coord[mirna][-1]
    seq2 = genome[chromo][start:end]
    if orientation == '-':
        seq2 = reverse_complement(seq2)
    if seq1 != seq2:
        print(mirna)
    
 
    
 