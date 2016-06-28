# -*- coding: utf-8 -*-
"""
Created on Tue Jun 28 13:50:20 2016

@author: RJovelin
"""

# use this script to compute divergence between remanei and latens orthologs

from tcoffee_alignment import *
from manipulate_sequences import *
from divergence import *
import os




# parse ortholog file to get the remanei and latens orthologs
# create a dict of latens and remanei mirnas orthologs
matureorthos = {}
infile = open('CremClamiRNAOrthologs.txt')
for line in infile:
    if line.rstrip() != '':
        line = line.rstrip().split('\t')
        matureorthos[line[0]] = line[1]
infile.close()

# create a dict with orthologs for hairpin sequences
hairpinorthos = {}
for mir in matureorthos:
    # strip name from arm
    crmname = mir[: mir.index('_')]
    claname = matureorthos[mir][:matureorthos[mir].index('_')]
    hairpinorthos[crmname] = claname

# create a dict with chromo : sequence for remanei and for latens genomes
ClaGenome = convert_fasta('../Genome_Files/noamb_534_v1.txt')
CrmGenome = convert_fasta('../Genome_Files/noamb_356_v1_4.txt')

# create a dict with latens hairpin
# create a dict with latens mature
clahairpin, clamature = {}, {}
infile = open('Cla_miRNACoordinatesFinal.txt')
infile.readline()
for line in infile:
    if line.rstrip() != '':
        line = line.rstrip().split('\t')
        # get mirna name, chromo, strand, start and end
        name, chromo, strand, start, end = line[0], line[1], line[2], int(line[3]) - 1, int(line[4])
        # extract hairpin
        hairpin = ClaGenome[chromo][start: end]
        if strand == '-':
            hairpin = reverse_complement(hairpin)
        # get mature sequence
        mature = line[6]
        # check that mature is in hairpin
        assert mature in hairpin, 'mature not in cla hairpin sequence'
        # populate dicts
        clamature[name] = mature
        # remove arm from name
        mirname = name[:name.index('_')]
        clahairpin[mirname] = hairpin
infile.close()

# create a dict with remanei hairpin
# create a dict with remanei mature
crmhairpin, crmmature = {}, {}
infile = open('CRM_miRNAsCoordinatesFinal.txt')
infile.readline()
for line in infile:
    if line.rstrip() != '':
        line = line.rstrip().split('\t')
        # get mirna name, chromo, strand, start, end
        name, chromo, start, end, strand = line[0], line[1], int(line[2]) - 1, int(line[3]), line[4]
        # extract hairpin from genome
        hairpin = CrmGenome[chromo][start: end]
        if strand == '-':
            hairpin = reverse_complement(hairpin)
        # get mature sequence
        mature = line[6]
        # check that mature is in hairpin
        assert mature in hairpin, 'mature not in crm hairpin sequence'
        # populate dicts
        crmmature[name] = mature
        # remove arm from name
        mirname = name[: name.index('_')]
        crmhairpin[mirname] = hairpin
infile.close()

# create directory to save alignment files
try:
    os.listdir('Crm_Cla_miRNA_orthos')
except:
    os.mkdir('Crm_Cla_miRNA_orthos')

# all pairs of hairpin sequences have been manually inspected in a first run    
    
# create files 
for crmmirna in crmhairpin:
    # if mirna has ortholog in latens
    if crmmirna in hairpinorthos:
        # get cla ortholog
        clamirna = hairpinorthos[crmmirna]
        if clamirna in clahairpin:
            # open file for writing
            # get file from mirna name
            newfile = open('Crm_Cla_miRNA_orthos/' + crmmirna + '_hairpin.fas', 'w')
            newfile.write('>' + crmmirna + '\n')
            newfile.write(crmhairpin[crmmirna] + '\n')
            newfile.write('>' + clamirna + '\n')
            if clamirna in ['cla-miR-231', 'cla-miR-50']:
                # take the reverse complement of the cla hairpin
                newfile.write(reverse_complement(clahairpin[clamirna]) + '\n')
            else:
                newfile.write(clahairpin[clamirna] + '\n')
            newfile.close()
            

# all pairs of mature sequences have been manually inspected in a first run

for crmmirna in crmmature:
    # if mirna has ortholog in latens
    if crmmirna in matureorthos:
        # do not align miR-2225_5p, sequences do not align at all
        if crmmirna != 'crm-miR-2225_5p':
            # get cla ortholog
            clamirna = matureorthos[crmmirna]
            if clamirna in clamature:
                # open file for writing
                # get file from mirna name
                newfile = open('Crm_Cla_miRNA_orthos/' + crmmirna + '_mature.fas', 'w')
                newfile.write('>' + crmmirna + '\n')
                newfile.write(crmmature[crmmirna] + '\n')
                # if crm mature is miR-39, align aopposite arms
                if crmmirna == 'crm-miR-39_5p':
                    newfile.write('>cla-miR-39_3p' + '\n')
                    newfile.write(clamature['cla-miR-39_3p'] + '\n')
                elif crmmirna == 'crm-miR-39_3p':
                    newfile.write('>cla-miR-39_5p' + '\n')
                    newfile.write(clamature['cla-miR-39_5p'] + '\n')
                else:
                    # write cla mirna sequence
                    newfile.write('>' + clamirna + '\n')
                    if clamirna in ['cla-miR-231_5p', 'cla-miR-50_3p', 'cla-miR-50_5p']:
                        # take the reverse complement of the cla mature                
                        newfile.write(reverse_complement(clamature[clamirna]) + '\n')
                    else:
                        newfile.write(clamature[clamirna] + '\n')
                newfile.close()
            

# align premirnas and mature sequences
run_tcoffee_noncoding('Crm_Cla_miRNA_orthos/')
# t-coffee outputs are saved in the current directory so move outputfiles
os.system('mv crm-*.* ./Crm_Cla_miRNA_orthos/')
# convert the t-coffee output to fasta format and save in text files
generate_fasta_from_tcoffee('./Crm_Cla_miRNA_orthos/')

# change dir
os.chdir('./Crm_Cla_miRNA_orthos/')

# create a list of hairpin files
hairpinfiles = [i for i in os.listdir() if 'hairpin' in i and '.txt' in i]
# create a list of mature files
maturefiles = [i for i in os.listdir() if 'mature' in i and '.txt' in i]

# create a dict to store the Jukes-Cantor distance for hairpins and for mature
Khairpin, Kmature = {}, {}

# count number of undefined K for hairpin and mature
undefined_hairpin, undefined_mature = 0, 0


# loop over hairpin files
for filename in hairpinfiles:
    # convert file to dict
    fastafile = convert_fasta(filename)
    # get seq names and sequences
    sequences = [[name, seq] for name, seq in fastafile.items()]
    sequences.sort()
    assert sequences[0][0].startswith('cla'), '1st sequence should be cla'
    assert sequences[1][0].startswith('crm'), '2nd sequence should be crm'
    # compute K Jukes-Cantor
    K = compute_K_JC(sequences[0][1], sequences[1][1])
    if K != 'NA':
        # populate dict {crm mirna name: [cla-mirna name, K]}
        Khairpin[sequences[1][0]] = [sequences[0][0], K]
    else:
        print(sequences[1][0])
        undefined_hairpin += 1
        
# loop over mature files
for filename in maturefiles:
    # convert file to dict
    fastafile = convert_fasta(filename)
    # get seq names and sequences
    sequences = [[name, seq] for name, seq in fastafile.items()]
    sequences.sort()
    assert sequences[0][0].startswith('cla'), '1st sequence should be cla'
    assert sequences[1][0].startswith('crm'), '2nd sequence should be crm'
    # compute K Jukes-Cantor
    K = compute_K_JC(sequences[0][1], sequences[1][1])
    if K != 'NA':
        # populate dict {crm mirna name: [cla-mirna name, K]}
        Kmature[sequences[1][0]] = [sequences[0][0], K]
    else:
        print(sequences[1][0])
        undefined_mature += 1
        
print('cannot compute K for {0} hairpins'.format(len(undefined_hairpin)))
print('cannot compute K for {0} matures'.format(len(undefined_mature)))        
        
# move to parent directory
os.chdir('../')

# create a table with divergence value between hairpins
newfile = open('CrmClamiRNAHairpinDivergence.txt', 'w')
newfile.write('Cremanei\tClatens\t\K_JukesCantor\n')
for mir in Khairpin:
    newfile.write(mir + '\t' + Khairpin[mir][0] + '\t' + str(Khairpin[mir][1]) + '\n')
newfile.close()

newfile = open('CrmClamiRNAMatureDivergence.txt', 'w')
newfile.write('Cremanei\tClatens\t\K_JukesCantor\n')
for mir in Kmature:
    newfile.write(mir + '\t' + Kmature[mir][0] + '\t' + str(Kmature[mir][1]) + '\n')
newfile.close()

