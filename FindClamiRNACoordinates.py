# -*- coding: utf-8 -*-
"""
Created on Thu Jun 23 16:29:59 2016

@author: RJovelin
"""

import os
import sys
from manipulate_sequences import *
from divergence import *

# use this script to extract latens and remanei mirna orthologs
# and to find the latens mirna coordinates

# usage GetClamiRNACoordinates.py [search/found]
# - [search/found]: search to find if mirnas match in genome, and found if mirnas sequences have been edited

# get option from command
step = sys.argv[1]

# check that mutual blast hist retrieved the same best hists
infile = open('../Genome_Files/ClaCrmCla_miRNABestBlast.txt')
total = 0
nonmatching = 0
for line in infile:
    if line.rstrip() != '':
        line = line.rstrip().split('\t')
        # compare set of mirna pairs
        a, b = set(), set()
        a.add(line[0])
        a.add(line[1])
        b.add(line[12])
        b.add(line[13])
        if a != b:
            nonmatching += 1
        total += 1
infile.close()

print('total', total)
print('non-matching', nonmatching)
print('proportion', nonmatching / total * 100)
assert nonmatching == 0, 'some blast hists are not reciprocal'

# remove any miRNAs with multiple hits

# make a dict for cla to crm blast
ClaToCrm = {}
# make a dict for crm to cla blast
CrmToCla = {}
infile = open('../Genome_Files/ClaCrmCla_miRNABestBlast.txt')
for line in infile:
    if line.rstrip() != '':
        line = line.rstrip().split('\t')
        # populate dict
        if line[0] in ClaToCrm:
            ClaToCrm[line[0]].append(line[1])
        else:
            ClaToCrm[line[0]] = [line[1]]
            
        if line[12] in CrmToCla:
            CrmToCla[line[12]].append(line[13])
        else:
            CrmToCla[line[12]] = [line[13]]
infile.close()
# count cla and crm mirnas with more than 1 hist
clamore, crmmore = 0, 0
for mirna in ClaToCrm:
    if len(ClaToCrm[mirna]) > 1:
        clamore += 1
for mirna in CrmToCla:
    if len(CrmToCla[mirna]) > 1:
        crmmore += 1
print('crmmore', crmmore)
print('clamore', clamore)

# make a set of keys to remove
to_delete = set()
removed_cla, removed_crm = 0, 0
for mirna in ClaToCrm:
    if len(ClaToCrm[mirna]) > 1:
        # add cla key to set of keys to remove
        to_delete.add(mirna)
        # also delete the remanei keys from the crm to cla dict
        for i in ClaToCrm[mirna]:
            to_delete.add(i)
for mirna in CrmToCla:
    if len(CrmToCla[mirna]) > 1:
        # add crm key to set of keys to delete
        to_delete.add(mirna)
        # also delete the latens mirnas
        for i in CrmToCla[mirna]:
            to_delete.add(i)
for mirna in to_delete:
    if mirna in ClaToCrm:
        del ClaToCrm[mirna]
        removed_cla += 1
    if mirna in CrmToCla:
        del CrmToCla[mirna]
        removed_crm += 1
print('removed from cla', removed_cla)
print('removed from crm', removed_crm)
print('crm', len(CrmToCla))
print('cla', len(ClaToCrm))

# make a dict of remanei : latens pairs
# make sets of keys and vals for each dicts
crmkeys, crmvals, clakeys, clavals = set(), set(), set(), set()
for i in CrmToCla:
    crmkeys.add(i)
    assert len(CrmToCla[i]) == 1, 'crm should have a single cla match'
    crmvals.add(CrmToCla[i][0])
for i in ClaToCrm:
    clakeys.add(i)
    assert len(ClaToCrm[i]) == 1, 'cla should have a single crm match'
    clavals.add(ClaToCrm[i][0])
# compare sets
assert crmkeys == clavals, 'crm keys different than cla vals'
assert clakeys == crmvals, 'cla keys than crm vals'
# keys of crm to cla are vals of cla to crm, so can use either dict as final dict of remanei latens pairs

# generate a file of remanei and latens mirna orthologs
newfile = open('CremClamiRNAOrthologs.txt', 'w')
for mir in CrmToCla:
    # change names in lower caps
    crmmirna = mir.lower()
    clamirna = CrmToCla[mir][0].lower()
    # change "new" to mir
    if 'new' in crmmirna:
        crmmirna = crmmirna.replace('new', 'mir')
    if 'new' in clamirna:
        clamirna = clamirna.replace('new', 'mir')
    # replace mir -> miR
    if 'mir' in crmmirna:
        assert 'mir' in clamirna, "mir in crm name but not in cla name"
        crmmirna = crmmirna.replace('mir', 'miR')
        clamirna = clamirna.replace('mir', 'miR')
    # replace "-5p" and "-3p" with "_5p" and "_3p" 
    if '-5p' in crmmirna:
        crmmirna = crmmirna.replace('-5p', '_5p')
        if '-5p' in clamirna:
            clamirna = clamirna.replace('-5p', '_5p')
        elif '3p' in clamirna:
            clamirna = clamirna.replace('-3p', '_3p')
    elif '-3p' in crmmirna:
        crmmirna = crmmirna.replace('-3p', '_3p')
        if '-3p' in clamirna:
            clamirna = clamirna.replace('-3p', '_3p')
        elif '-5p' in clamirna:
            clamirna = clamirna.replace('-5p', '_5p')
    # write to file
    newfile.write(crmmirna + '\t' + clamirna + '\n')
newfile.close()


# generate a file with latens coordinates 
mirnas = {}
infile = open('latens_coords.txt')
infile.readline()
for line in infile:
    if line.rstrip() != '':
        line = line.rstrip().split('\t')
        if len(line) == 8 and line[3] != '':
            # get mirna name, hairpin, chromo, hairpin start
            name, hairpin, chromo, start = line[0], line[1].upper(), line[2], int(line[3]) -1
            # get mature, mature coord, and arm
            mature, matcoord, arm = line[4].upper(), line[5], line[7]
            # edit name
            if 'new' in name:
                name = name.replace('new', 'mir')
            if 'mir' in name:
                name = name.replace('mir', 'miR')
            # add arm to name
            name = name + '_' + arm
            # parse mature coords
            matstart = matcoord[matcoord.index(':')+1: matcoord.index('-')]
            matstart = int(matstart) - 1
            end = start + len(hairpin)
            matend = matstart + len(mature)
            # populate dict
            mirnas[name] = [chromo, str(start + 1), str(end), hairpin, str(matstart + 1), str(matend), mature, arm]
infile.close()


newfile = open('Cla_miRNACoordinatesTemporary.txt', 'w')
#write header
header = ['name', 'chromo', 'hairpin_start', 'hairpin_end', 'hairpin', 'mature_start', 'mature_end', 'mature', 'arm']
newfile.write('\t'.join(header) + '\n')

mirNames = [i for i in mirnas]
mirNames.sort()
for mir in mirNames:
    newfile.write(mir + '\t' + '\t'.join(mirnas[mir]) + '\n')
newfile.close()

# convert 
ClaGenome = convert_fasta('../Genome_Files/noamb_534_v1.txt')

# create a dict with hairpin and sequences (+ and rev complement) extracted from genome
# {name :[hairpin, seq, reverse_complement]}
sequences = {}
infile = open('Cla_miRNACoordinatesTemporary.txt')
infile.readline()
for line in infile:
    if line.rstrip() != '':
        line = line.rstrip().split('\t')
        # get mirna name, hairpin, chromo, start
        mirna, hairpin, chromo, start = line[0], line[4],  line[1], int(line[2]) - 1
        # get end position
        end = start + len(hairpin)
        # extract forward sequence
        seq = ClaGenome[chromo][start: end].upper()
        # get reverse sequence
        revseq = reverse_complement(seq)
        sequences[mirna] = [hairpin, seq, revseq]
infile.close()            
    
# create a dict to store differences between hairpin and extracted sequences
# {name :[diff_seq, diff_revseq]}
differences = {}
for i in sequences:
    diff1 = match_diff(sequences[i][0], sequences[i][1])
    diff2 = match_diff(sequences[i][0], sequences[i][2])
    differences[i] = [diff1, diff2]

# create a set of latens mirnas that have remanei orthologs
claorthos = set()
infile = open('CremClamiRNAOrthologs.txt')
for line in infile:
    if line.rstrip() != '':
        line = line.rstrip().split('\t')
        claorthos.add(line[1])
infile.close()


# parse mirna coordinates, no need to convert coordinates
mirnas = {}
infile = open('Cla_miRNACoordinatesTemporary.txt')
infile.readline()
for line in infile:
    if line.rstrip() != '':
        line = line.rstrip().split('\t')
        # get mirna name
        name = line[0]
        mirnas[name] = line[1:]
infile.close()


# check if mirnas sequences need to be edited
if step == 'search':
    # create lists of mirna names for divergent and related sequences 
    divergent, nodiff = [], []
        
    # create a folder with fasta files to edit
    try:
        os.listdir('files_to_edit')
    except:
        os.mkdir('files_to_edit')
    
    missing = set()
    # loop over differences 
    for mir in differences:
        # edit only cla mirna coordinates if they have a remanei ortholog
        if mir in claorthos:
            # check if differences are large or low
            # get mirna orientation
            if differences[mir][0] == 0:
                strand = '+'
                assert differences[mir][1] != 0, 'differences are not strand biased'
                nodiff.append(mir)
            elif differences[mir][1] == 0:
                strand = '-'
                assert differences[mir][0] != 0, 'differences are not strand biased'
                nodiff.append(mir)
            else:
                strand = 'NA'
                divergent.append(mir)
            
            # insert strand after chromo
            mirnas[mir].insert(1, strand)
            # get mature miR
            mature = mirnas[mir][7]
            # edit all mirnas with differences
            if strand == 'NA':
                # add mature to file and do manual curation
                newfile = open('files_to_edit/' + mir + '.fasta', 'w')
                newfile.write('>cla\n')
                newfile.write(sequences[mir][0] + '\n')
                newfile.write('>seq\n')
                newfile.write(sequences[mir][1] + '\n')
                newfile.write('>revseq\n')
                newfile.write(sequences[mir][2] + '\n')
                newfile.write('>mature\n')
                newfile.write(mature + '\n')
                newfile.close()
    
    print('nodiff', len(nodiff))
    print('divergent', len(divergent))
        
    # write mirna information to file for mirnas with 
    mirNames = [i for i in nodiff]
    mirNames.extend(divergent)
        
    newfile = open('Cla_miRNACoordinatesTemporary2.txt', 'w')
    #write header
    header = ['name', 'chromo', 'orientation', 'hairpin_start', 'hairpin_end', 'hairpin', 'mature_start', 'mature_end', 'mature', 'arm']
    newfile.write('\t'.join(header) + '\n')
    for mir in mirNames:
        newfile.write(mir + '\t' + '\t'.join(mirnas[mir]) + '\n')
    newfile.close()


elif step == 'found':
    # create a dict with mirna infos, update coordinates
    mirnas = {}    
    to_remove = []    
    infile = open('Cla_miRNACoordinatesTemporary3.txt')
    infile.readline()
    for line in infile:
        if line.rstrip():
            line = line.rstrip().split('\t')
            # check that all mature are in hairpin
            name, chromo, hairpin, mature = line[0], line[1], line[5], line[8]
            if mature not in hairpin:
                if mature in reverse_complement(hairpin):
                    # update hairpin
                    hairpin = reverse_complement(hairpin)
                else:
                    print('mature does not match', name)
            # get coordinates
            if hairpin in ClaGenome[chromo]:
                genomestart = ClaGenome[chromo].index(hairpin)
                genomeend = genomestart + len(hairpin)
                strand = '+'
            elif reverse_complement(hairpin) in ClaGenome[chromo]:
                genomestart = ClaGenome[chromo].index(reverse_complement(hairpin))
                genomeend = genomestart + len(hairpin)
                strand = '-'
            else:
                print('hairpin does not match', name)
                to_remove.append(name)
            # get 7bp seed
            seed = mature[1:8]
            #populate dict
            if name not in to_remove:
                mirnas[name] = [chromo, strand, str(genomestart + 1), str(genomeend), hairpin, mature, seed, arm] 
                
    infile.close()    
    
    # write coordinates to file
    header = ['name', 'chromo', 'orientation', 'hairpin_start', 'hairpin_end', 'hairpin', 'mature', 'seed', 'arm']    
    newfile = open('Cla_miRNACoordinatesFinal.txt', 'w')
    newfile.write('\t'.join(header) + '\n')
    for name in mirnas:
        newfile.write(name + '\t' + '\t'.join(mirnas[name]) + '\n')
    newfile.close()
    
    # do some QC
    infile = open('Cla_miRNACoordinatesFinal.txt')
    infile.readline()
    for line in infile:
        if line.rstrip() != '':
            line = line.rstrip().split('\t')
            # get mirna name, hairpin seq, mature seq
            name, hairpin, mature = line[0], line[5], line[6]
            # get coordinates
            chromo, strand, start, end = line[1], line[2], int(line[3]) - 1, int(line[4])
            # extract sequence from genome
            seq = ClaGenome[chromo][start: end]
            if strand == '-':
                seq = reverse_complement(seq)
            if seq != hairpin:
                print(name, seq, hairpin, strand)
            if mature not in seq:
                print('mature not in extracted sequence for {0}'.format(name))
            if mature not in hairpin:
                print('mature not in hairpin for {0}'.format(name))
    infile.close()
    

