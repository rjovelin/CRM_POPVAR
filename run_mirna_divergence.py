# -*- coding: utf-8 -*-
"""
Created on Tue Jun 28 13:50:20 2016

@author: RJovelin
"""

# use this script to compute divergence between remanei and latens orthologs


from manipulate_sequences import *


# parse ortholog file to get the remanei and latens orthologs
# create a dict of latens and remanei mirnas orthologs
orthos = {}
infile = open('CremClamiRNAOrthologs.txt')
for line in infile:
    if line.rstrip() != '':
        line = line.rstrip().split('\t')
        orthos[line[0]] = line[1]
infile.close()

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
        





















# align premirnas

# align matures mirnas

# create a file with alignments

# manually inspect alignments

# create a table with divergence value between hairpin and mirnas






# plot divergence for premirnas and for mature mirnas



# use paml to compute divergence for comparison with dN and dS?



# or can I compare Jukes-Cantor with ML estimates?