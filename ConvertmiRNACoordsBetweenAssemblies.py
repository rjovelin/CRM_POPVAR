# -*- coding: utf-8 -*-
"""
Created on Sat May  7 22:16:13 2016

@author: Richard
"""

from Manipulate_Sequences import *
from divergence import *
import os


# add family conservation to mirnas
# get the seeds of distant caeno species
caenoSeeds = set()
infile = open('CBR_seeds.txt')
for line in infile:
    line = line.rstrip()
    if line != '':
        caenoSeeds.add(line)
infile.close()
infile = open('CEL_seeds.txt')
for line in infile:
    line = line.rstrip()
    if line != '':
        caenoSeeds.add(line)
infile.close()
infile = open('CBN_seeds.txt')
for line in infile:
    line = line.rstrip()
    if line != '':
        caenoSeeds.add(line)
infile.close()

# get seeds in latens
latensSeeds = set()
infile = open('CLA_seeds.txt')
for line in infile:
    line = line.rstrip()
    if line != '':
        latensSeeds.add(line)
infile.close()


# create a dict to store information about mirna {name :[list of info]}
mirnas = {}
infile = open('CRM_miRNAs_Coords_rj.txt')
header = infile.readline().rstrip().split('\t')
for line in infile:
    if line.rstrip() != '':
        line = line.rstrip().split('\t')
        mirnas[line[0]] = line[1:]
infile.close()

# make sure strings are upper caps
for i in mirnas:
    mirnas[i][0] = mirnas[i][0].upper()
    mirnas[i][3] = mirnas[i][3].upper()
    mirnas[i][4] = mirnas[i][4].upper()
    
# check family conservation
for mir in mirnas:
    # get the seed sequence
    seed = mirnas[mir][4]
    if seed in caenoSeeds:
        conservation = 'Caeno'
    if seed not in caenoSeeds and seed in latensSeeds:
        conservation = 'CrmCla'
    if seed not in caenoSeeds and seed not in latensSeeds:
        conservation = 'Crm'
    # add conservation as varaible to list
    mirnas[mir].append(conservation)

# add fields to header
header.append('family_conservation')
header.append('strand')

# create a list of mirna names
names = [i for i in mirnas]
names.sort()

# convert genome sequence to dict
genome = convert_fasta('noamb_356_v1_4.txt')

# create a dict with hairpin and sequences (+ and rev complement) extracted from genome
# {name :[hairpin, seq, reverse_complement]}
sequences = {}
for mir in mirnas:
    # get hairpin
    hairpin = mirnas[mir][0]
    # get chromo
    chromo = mirnas[mir][1]
    # get start and end positions 0-based
    start = int(mirnas[mir][2]) - 1
    end = start + len(hairpin)
    # extract forward sequence
    seq = genome[chromo][start: end]
    # extract reverse complement
    revseq = reverse_complement(seq)
    sequences[mir] = [hairpin, seq, revseq]
    
# create a dict to store differences between hairpin and extracted sequences
# {name :[diff_seq, diff_revseq]}
differences = {}
for i in mirnas:
    diff1 = match_diff(sequences[i][0], sequences[i][1])
    diff2 = match_diff(sequences[i][0], sequences[i][2])
    differences[i] = [diff1, diff2]

# create lists of mirna names for divergent and related sequences 
divergent, close, nodiff = [], [], []

try:
    os.listdir('files_to_edit')
except:
    os.mkdir('files_to_edit')

# loop over differences 
for mir in differences:
    # check if differences are large or low
    # get mirna orientation
    if differences[mir][0] <= 3:
        strand = '+'
        assert differences[mir][1] > 3, 'differences are not strand biased'
        if differences[mir][0] > 0:
            close.append(mir)
        elif differences[mir][0] == 0:
            nodiff.append(mir)
    elif differences[mir][1] <= 3:
        strand = '-'
        assert differences[mir][0] > 3, 'differences are not strand biased'
        if differences[mir][1] > 0:
            close.append(mir)
        elif differences[mir][1] == 0:
            nodiff.append(mir)
    elif differences[mir][0] > 3 and differences[mir][1] > 3:
        strand = 'NA'
        divergent.append(mir)
    # add orientation to list
    mirnas[mir].append(strand)
    # get mature miR
    mature = mirnas[mir][3]
    # check if mirna needs manual curation
    curation = False
    if strand in '+-':
        if strand == '+' and differences[mir][0] > 0:
            # check if mature has differences
            if mature not in sequences[mir][1]:
                curation = True
            else:
                # no need for manual curation but replace hairpin with seq
                mirnas[mir][0] = sequences[mir][1]
        elif strand == '-' and differences[mir][1] > 0:
            # check if mature has differences
            if mature not in sequences[mir][2]:
                curation = True
            else:
                # no need for manual curation but replace hairpin with revseq
                mirnas[mir][0] = sequences[mir][2]
    elif strand == 'NA':
        curation = True
    # check if manual curation is needed
    if curation == True:
        # add mature to file and do manual curation
        newfile = open('files_to_edit/' + mir + '.fasta', 'w')
        newfile.write('>crem\n')
        newfile.write(sequences[mir][0] + '\n')
        newfile.write('>seq\n')
        newfile.write(sequences[mir][1] + '\n')
        newfile.write('>revseq\n')
        newfile.write(sequences[mir][2] + '\n')
        newfile.write('>mature\n')
        newfile.write(mature + '\n')
        newfile.close()

# check that all mirnas are dealt with
assert len(mirnas) == len(nodiff) + len(close) + len(divergent), 'mirnas are missing'   
 
# save information to file
newfile = open('CRM_miRNAsCoordinatesEdits.txt', 'w')
newfile.write('\t'.join(header) + '\n')
total = 0
for mir in nodiff:
    newfile.write(mir + '\t' + '\t'.join(mirnas[mir]) + '\n')
    total += 1
for mir in close:
    newfile.write(mir + '\t' + '\t'.join(mirnas[mir]) + '\n')
    total += 1
for mir in divergent:
    newfile.write(mir + '\t' + '\t'.join(mirnas[mir]) + '\n')
    total += 1
newfile.close()

# check that all mirnas are saved to file
assert len(mirnas) == total, 'some mirnas are not saved to file'
