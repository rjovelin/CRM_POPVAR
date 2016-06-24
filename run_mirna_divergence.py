# -*- coding: utf-8 -*-
"""
Created on Thu Jun 23 16:29:59 2016

@author: RJovelin
"""


from manipulate_sequences import *
from divergence import *

# use this script to compute divergence between remanei and latens miRNAs


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


# remove any miRNAs with multiple hits

# make a dict for cla to crm blast
ClaToCrm = {}
# make a dict for crm to cla blast
CrmToCla = {}

infile = open('../Genome_Files/ClaCrmCla_miRNABestBlast.txt')
total = 0
nonmatching = 0
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


# convert 
ClaGenome = convert_fasta('../Genome_Files/noamb_534_v1.txt')

# get latens coordinates
infile = open('latens_coords.txt')
infile.readline()
# make a dict {name: [hairpin, LG, start]}
clacoords = {}
for line in infile:
    if line.rstrip() != '':
        line = line.rstrip().split('\t')
        # only consider mirnas with hairpins and coordinate/arm information
        if len(line) == 8 and line[3] != '':
            clacoords[line[0]] = [line[1], line[2], int(line[3])-1]
infile.close()        

# check how many mirnas cannot be found in the genome
missing = 0
for mirna in clacoords:
    if clacoords[mirna][0].upper() not in ClaGenome[clacoords[mirna][1]]:
        if reverse_complement(clacoords[mirna][0].upper()) not in ClaGenome[clacoords[mirna][1]]:
            missing += 1

print('missing', missing)
print('total', len(clacoords))


# create a dict with hairpin and sequences (+ and rev complement) extracted from genome
# {name :[hairpin, seq, reverse_complement]}
sequences = {}

infile = open('latens_coords.txt')
infile.readline()
for line in infile:
    if line.rstrip() != '':
        line = line.rstrip().split('\t')
        # only consider mirnas with hairpins and coordinate/arm information
        if len(line) == 8 and line[3] != '':
            # get mirna name
            mirna = line[0]
            # get hairpin
            hairpin = line[1].upper()
            # get chromo
            chromo = line[2]
            # get hairpin start position
            start = int(line[3]) -1
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
Differences = {}
for i in sequences:
    diff1 = match_diff(sequences[i][0], sequences[i][1])
    diff2 = match_diff(sequences[i][0], sequences[i][2])
    Differences[i] = [diff1, diff2]

# create lists of mirna names for divergent and related sequences 
divergent, close, nodiff = [], [], []


a = [i for i in ClaToCrm]
for i in range(len(a)):
    if '5p' in a[i]:
        a[i] = a[i].replace('-5p', '')
    elif '3p' in a[i]:
        a[i] = a[i].replace('-3p', '')
    # replace mir
    a[i] = a[i].lower()
    if 'new' in a[i]:
        a[i] = a[i]. replace('new', 'mir')

print('a', len(a))
a = set(a)
print('a', len(a))        
        
        
differences = {}
for mir in Differences:
    mirna = mir.lower()
    if 'new' in mirna:
        mirna = mirna.replace('new', 'mir')
    differences[mirna] = list(Differences[mir])


# loop over differences 
for mir in differences:
    if mir in a:
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


for mir in divergent:
    print(mir)
    print(sequences[mir][0])
    print(sequences[mir][1])
    print(sequences[mir][2])
    print('*****')
        


print('nodiff', len(nodiff))
print('close', len(close))
print('divergent', len(divergent))
print('total', len(nodiff) + len(divergent) + len(close))























# change names and match names according to the remanei diversity table

# verify that sequences correspond to extracted sequences



# align premirnas

# align matures mirnas


# plot divergence for premirnas and for mature mirnas
# use paml to compute divergence for comparison with dN and dS?
# or can I compare Jukes-Cantor with ML estimates?




#################################

#    # create lists of mirna names for divergent and related sequences 
#    divergent, close, nodiff = [], [], []
#
#    try:
#        os.listdir('files_to_edit')
#    except:
#        os.mkdir('files_to_edit')
#
#    # loop over differences 
#    for mir in differences:
#        # check if differences are large or low
#        # get mirna orientation
#        if differences[mir][0] <= 3:
#            strand = '+'
#            assert differences[mir][1] > 3, 'differences are not strand biased'
#            if differences[mir][0] > 0:
#                close.append(mir)
#            elif differences[mir][0] == 0:
#                nodiff.append(mir)
#        elif differences[mir][1] <= 3:
#            strand = '-'
#            assert differences[mir][0] > 3, 'differences are not strand biased'
#            if differences[mir][1] > 0:
#                close.append(mir)
#            elif differences[mir][1] == 0:
#                nodiff.append(mir)
#        elif differences[mir][0] > 3 and differences[mir][1] > 3:
#            strand = 'NA'
#            divergent.append(mir)
#        # add orientation to list
#        mirnas[mir].append(strand)
#        # get mature miR
#        mature = mirnas[mir][3]
#        # check if mirna needs manual curation
#        curation = False
#        if strand in '+-':
#            if strand == '+' and differences[mir][0] > 0:
#                # check if mature has differences
#                if mature not in sequences[mir][1]:
#                    curation = True
#                else:
#                    # no need for manual curation but replace hairpin with seq
#                    mirnas[mir][0] = sequences[mir][1]
#            elif strand == '-' and differences[mir][1] > 0:
#                # check if mature has differences
#                if mature not in sequences[mir][2]:
#                    curation = True
#                else:
#                    # no need for manual curation but replace hairpin with revseq
#                    mirnas[mir][0] = sequences[mir][2]
#        elif strand == 'NA':
#            curation = True
#        # check if manual curation is needed
#        if curation == True:
#            # add mature to file and do manual curation
#            newfile = open('files_to_edit/' + mir + '.fasta', 'w')
#            newfile.write('>crem\n')
#            newfile.write(sequences[mir][0] + '\n')
#            newfile.write('>seq\n')
#            newfile.write(sequences[mir][1] + '\n')
#            newfile.write('>revseq\n')
#            newfile.write(sequences[mir][2] + '\n')
#            newfile.write('>mature\n')
#            newfile.write(mature + '\n')
#            newfile.close()
#
#    # check that all mirnas are dealt with
#    assert len(mirnas) == len(nodiff) + len(close) + len(divergent), 'mirnas are missing'   
# 
#    # save information to file
#    newfile = open('CRM_miRNAsCoordinatesEdits.txt', 'w')
#    newfile.write('\t'.join(header) + '\n')
#    total = 0
#    for mir in nodiff:
#        newfile.write(mir + '\t' + '\t'.join(mirnas[mir]) + '\n')
#        total += 1
#    for mir in close:
#        newfile.write(mir + '\t' + '\t'.join(mirnas[mir]) + '\n')
#        total += 1
#    for mir in divergent:
#        newfile.write(mir + '\t' + '\t'.join(mirnas[mir]) + '\n')
#        total += 1
#    newfile.close()
#
#    # check that all mirnas are saved to file
#    assert len(mirnas) == total, 'some mirnas are not saved to file'
#
#elif step == 'found':
#    # remove mirnas for which the mature miR doesn't match the hairpin sequence    
#    to_remove = []
#    # convert file to dict, check that all mature seqs are in hairpin
#    infile = open('CRM_miRNAsCoordinatesEdited.txt')
#    header = infile.readline().rstrip().split('\t')
#    mirnas = {}
#    for line in infile:
#        line = line.rstrip()
#        if line != '':
#            line = line.split('\t')
#            if line[4] not in line[1]:
#                # check that mature in reverse complement of hairpin
#                if line[4] in reverse_complement(line[1]):
#                    # change hairpin to its reverse complement
#                    line[1] = reverse_complement(line[1])
#                else:
#                    to_remove.append(line[0])
#            mirnas[line[0]] = line
#    infile.close()
#    
#    # remove mirnas
#    for i in to_remove:
#        del mirnas[i]
#        
#    # convert 6bp seed to 7bp seed
#    for mir in mirnas:
#        mirnas[mir][5] = mirnas[mir][4][1:8]
#        
#    newfile = open('CRM_miRNAsCoordinatesFinal.txt', 'w')   
#    newfile.write('\t'.join(header) + '\n')
#    for i in mirnas:
#        newfile.write('\t'.join(mirnas[i]) + '\n')
#    newfile.close()
#    
#    # do some QC
#    infile = open('CRM_miRNAsCoordinatesFinal.txt')
#    infile.readline()
#    for line in infile:
#        line = line.rstrip()
#        if line != '':
#            line = line.split('\t')
#            # check that mature in hairpin
#            assert line[4] in line[1], 'mature not in hairpin'
#            # check that seed  == 7bp
#            assert len(line[5]) == 7, 'seed is not 6bp'
#            # check that seed is seed
#            assert line[4][1:8] == line[5], 'seed does not match'
#    infile.close()    
#    
#    # add coordinates to each mirnas
#    mirnas = {}
#    infile = open('CRM_miRNAsCoordinatesFinal.txt')
#    header = infile.readline().rstrip().split('\t')
#    for line in infile:
#        line = line.rstrip()
#        if line != '':
#            line = line.split('\t')
#            mir = line[0]
#            mirnas[mir] = line
#    infile.close()
#    
#    # modify header
#    # remove strand
#    header.remove(header[-1])
#    # remove start position
#    header.remove(header[3])
#    # correct typo
#    header[1] = 'Hairpin'    
#    coord = ['start', 'end', 'orientation']    
#    # extract chromo and add at first position in coord list
#    chromo = header.pop(2)
#    # add coordinates after chromo
#    coord.insert(0, chromo)
#    # add coordinates after mirna name
#    j = 1
#    for i in range(len(coord)):
#        header.insert(j+i, coord[i])
#
#    # convert genome sequence to dict
#    genome = convert_fasta('../Genome_Files/noamb_356_v1_4.txt')
#    
#    # loop over mirnas, modifiy conservation and coordinates
#    for mir in mirnas:
#        # remove strand
#        mirnas[mir].remove(mirnas[mir][-1])
#        # remove conservation
#        mirnas[mir].remove(mirnas[mir][-1])
#        # remove start position
#        mirnas[mir].remove(mirnas[mir][3])
#        # get hairpin
#        hairpin = mirnas[mir][1]        
#        # extract chromo
#        chromo = mirnas[mir].pop(2)
#        # find coordinates
#        if hairpin in genome[chromo]:
#            start = genome[chromo].index(hairpin)
#            end = start + len(hairpin)
#            orientation = '+'
#        elif reverse_complement(hairpin) in genome[chromo]:
#            start = genome[chromo].index(reverse_complement(hairpin))
#            end = start + len(hairpin)
#            orientation = '-'
#        else:
#            print(mir)
#        # create a list of coordinates, change coordinates to 1-based
#        coord = [chromo, str(start + 1), str(end), orientation]
#        # add coordinates after mirna name
#        j = 1
#        for i in range(len(coord)):
#            mirnas[mir].insert(j+i, coord[i])
#
#        # check family conservation
#        # get the seed sequence
#        seed = mirnas[mir][7]
#        if seed in caenoSeeds:
#            conservation = 'Caeno'
#        elif seed not in caenoSeeds and seed in latensSeeds:
#            conservation = 'CrmCla'
#        elif seed not in caenoSeeds and seed not in latensSeeds:
#            conservation = 'Crm'
#        # add conservation as varaible to list
#        mirnas[mir].append(conservation)
#    
#    # save data to file
#    # create a list of mirna names
#    names = [i for i in mirnas]
#    names.sort()
#    newfile = open('CRM_miRNAsCoordinatesFinal.txt', 'w')   
#    newfile.write('\t'.join(header) + '\n')
#    for i in names:
#        newfile.write('\t'.join(mirnas[i]) + '\n')
#    newfile.close()