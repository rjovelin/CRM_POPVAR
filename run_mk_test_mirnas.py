# -*- coding: utf-8 -*-
"""
Created on Wed Jun 29 16:07:21 2016

@author: RJovelin
"""

# use this script to compute the MK test for miRNAs

import os
from divergence import *
from miRNA_target import *
from sites_with_coverage import *
from mk_test import *
from manipulate_sequences import *
from multiple_testing import *


# convert genome to dict
CrmGenome = convert_fasta('../Genome_Files/noamb_356_v1_4.txt')
print('converted genome to dict')

# create a set of valid transcripts (1 transcript mapped to 1 gene)
valid_transcripts = get_valid_transcripts('../Genome_Files/unique_transcripts.txt')
print('made a set of valid transcripts')

# create a dict with remanei mirna hairpin name and remanei-latens aligned hairpins
hairpins = {}
# get the hairpin alignment files
hairpinfiles = [i for i in os.listdir('Crm_Cla_miRNA_orthos') if 'hairpin' in i and '.txt' in i]
for filename in hairpinfiles:
    # extract crm mir name from file
    crmname = filename[:filename.index('_')]
    # convert file to dict
    fastafile = convert_fasta('Crm_Cla_miRNA_orthos' + '/' + filename)
    # initialize dict
    hairpins[crmname] = {}
    for name in fastafile:
        hairpins[crmname][name] = fastafile[name]
print('got aligned hairpins')
print('aligned hairpins', len(hairpins))


# create a dict with remanei mirna mature name and remanei-latens aligned mature
matures = {}
# get the mature alignment files
maturefiles = [i for i in os.listdir('Crm_Cla_miRNA_orthos') if 'mature' in i and '.txt' in i]
for filename in maturefiles:
    # extract crm mature name from file
    crmname = filename[:filename.index('_mature')]
    # convert file to dict
    fastafile = convert_fasta('Crm_Cla_miRNA_orthos' + '/' + filename)
    # initialize dict
    matures[crmname] = {}
    for name in fastafile:
        matures[crmname][name] = fastafile[name]
print('got aligned matures')
print('aligned matures', len(matures))


# make a dict with family level conservation for all miRNAs {name : conservarion}
famCons = {}
infile = open('CRM_miRNAsCoordinatesFinal.txt')
infile.readline()
for line in infile:
    line = line.rstrip()
    if line != '':
        line = line.split('\t')
        name, conservation = line[0], line[-1]
        famCons[name] = conservation
infile.close()
print('got mirna conservation level')

# create a dict with coordinates of mature sequences
miR_coord = {}
infile = open('CRM_MatureCoordinatesFinal.txt')
infile.readline()
for line in infile:
    line = line.rstrip()
    if line != '':
        line = line.split('\t')
        name, chromo, start, end, orientation = line[0], line[1], int(line[2]) -1, int(line[3]), line[4]
        if orientation == '+':
            assert CrmGenome[chromo][start : end] == line[-1], 'mirna does not match mature on +'
        elif orientation == '-':
            assert reverse_complement(CrmGenome[chromo][start : end]) == line[-1], 'mirna does not match mature on -'
        miR_coord[name] = [chromo, start, end, orientation]
infile.close()
print('got mature coordinates')


# create a dict with coordinates of hairpin sequences
hairpin_coord = {}
infile = open('CRM_miRNAsCoordinatesFinal.txt')
infile.readline()
for line in infile:
    line = line.rstrip()
    if line != '':
        line = line.split('\t')
        name, chromo, start, end, orientation = line[0], line[1], int(line[2]) - 1, int(line[3]), line[4]
        # remove arm from name
        name = name[: name.index('_')]
        if orientation == '+':
            assert CrmGenome[chromo][start : end] == line[5], 'mirna does not match hairpin on +'
        elif orientation == '-':
            assert reverse_complement(CrmGenome[chromo][start : end]) == line[5], 'mirna does not match hairpin on -'
        hairpin_coord[name] = [chromo, start, end, orientation]
infile.close()
print('got hairpin coordinates')


# get the allele counts for all sites with coverage
# {chromo: {site : [ref_allele, alt_allele, count_ref, count alt]}]}
chromo_sites = get_non_coding_snps('../SNP_files/', 10)
print('got allele counts at all sites')

# check that mirna names in aligned sequences dict are in coordinates dict
for mirna in hairpins:
    assert mirna in hairpin_coord, 'mirna with aligned sequences does not have coordinates'
for mirna in matures:
    assert mirna in miR_coord, 'miR with aligned sequences does not have coordinates'        

# create a dict to record the number of fixed diffs and polmorphisms for each mirna hairpin
# {mirna: [D, P]}
hairpin_diffs = CountPolymDivergmiRNAs(hairpins, hairpin_coord, CrmGenome, chromo_sites, True, 1)
print('got the diff counts for hairpins')  
  
# create a dict to record the number of fixed diffs and polmorphisms for each mirna mature
# {mirna: [D, P]}
mature_diffs = CountPolymDivergmiRNAs(matures, miR_coord, CrmGenome, chromo_sites, True, 1)
print('got the diff counts for matures')

# create a dict to record the number of fixed diffs and polymorphisms at 4-fold degenerate sites
fourfold_diffs = CountPolymDivergFourFold('../Genome_Files/CDS_SNP_DIVERG.txt', True, 1)
print('got the diff count at 4-fold degenerate sites')

# remove genes that are not valid
total = 0
to_remove = [i for i in fourfold_diffs if i not in valid_transcripts]
for i in to_remove:
    del fourfold_diffs[i]
    total += 1
print('removed {0} non-valid genes'.format(total))


# compute Mk test according to Lyu et al Plos Genetics 2014:
# New MicroRNAs in Drosophila—Birth, Death and Cycles of Adaptive Evolution
'''
McDonald-Kreitman test to detect positive selection in miRNAs from each age group
based on the polymorphisms within Dmel and divergence between Dmel-Dsim. 
Precursor or mature sequences of each miRNA group were combined and treated as
the functional category, 4-fold degenerate sites in the whole genome were used
as the neutral control.
'''

# compute the MK test for each individual mirnas with 4-fold diffs as neutral control
# make a list with polym and divergence at 4-fold by pooling all the gene-centered counts [D, P]
Neutral = [0, 0]
for gene in fourfold_diffs:
    Neutral[0] += fourfold_diffs[gene][0]
    Neutral[1] += fourfold_diffs[gene][1]
print('pooled diffs for 4-fold sites', Neutral)

# create a dict by combining D and P for mirna sites and for 4-fold
# {mirna: [Pmirna, P4fold, Dmirna, D4fold]}
HairpinCounts, MatureCounts = {}, {}
for mirna in hairpin_diffs:
    HairpinCounts[mirna] = [hairpin_diffs[mirna][1], Neutral[1], hairpin_diffs[mirna][0], Neutral[0]]
for mirna in mature_diffs:
    MatureCounts[mirna] = [mature_diffs[mirna][1], Neutral[1], mature_diffs[mirna][0], Neutral[0]]
print('pooled fixed diffs and polyms for hairpins and matures')

# compute MK for each mirna and generate a dict {mirna: [[Pmirna, P4fold, Dmirna, D4fold, Pval]}
MKhairpin = MK_test(HairpinCounts, 'fisher')
MKmature = MK_test(MatureCounts, 'fisher')
print('computed MK tests for hairpins and matures')

# determine genes that are neutral, under negative and positive selection
# make a list of genes with significant MK test
HairpinSignificant = [mirna for mirna in MKhairpin if MKhairpin[mirna][-1] < 0.05]
HairpinNeutral = [mirna for mirna in MKhairpin if MKhairpin[mirna][-1] >= 0.05]
# determine mirnas with significant MK that are under positive or negative selection
HairpinNegative, HairpinPositive = NaturalSelection(MKhairpin, HairpinSignificant)
MatureSignificant = [mirna for mirna in MKmature if MKmature[mirna][-1] < 0.05]
MatureNeutral = [mirna for mirna in MKmature if MKmature[mirna][-1] >= 0.05]
# determine mirnas with signidicant MK that are under positive or negative selection
MatureNegative, MaturePositive = NaturalSelection(MKmature, MatureSignificant)
print('got lists of mirnas under positive and negative selection')

print('significant', len(HairpinSignificant), len(MatureSignificant))
print('positive', len(HairpinPositive), len(MaturePositive))
print('negative', len(HairpinNegative), len(MatureNegative))
print('neutral', len(HairpinNeutral), len(MatureNeutral))

# apply a Benjamini-Hochberg correction for multiple testing {gene: corrected_Pval}
MKhairpinCorr, MKmatureCorr = {}, {}
for gene in MKhairpin:
    MKhairpinCorr[gene] = MKhairpin[gene][-1]
HairpinPvals = [(p, gene) for gene, p in MKhairpinCorr.items()]
PCorrHairpin = Benjamini_Hochberg_correction(HairpinPvals)
print('corrected P for hairpins')
for gene in MKmature:
    MKmatureCorr[gene] = MKmature[gene][-1]
MaturePvals = [(p, gene) for gene, p in MKmatureCorr.items()]
PCorrMature = Benjamini_Hochberg_correction(MaturePvals)
print('corrected P for matures')


# count the number of genes with significant MK after BJ correction with 10% FDR 
FDR = 0.1
HairpinSignifCorr = [mirna for mirna in PCorrHairpin if PCorrHairpin[mirna] < FDR]
HairpinNeutralCorr = [mirna for mirna in PCorrHairpin if PCorrHairpin[mirna] >= FDR]
# determine mirnas with significant MK that are under positive or negative selection
HairpinNegativeCorr, HairpinPositiveCorr = NaturalSelection(MKhairpin, HairpinSignifCorr)
MatureSignifCorr = [mirna for mirna in PCorrMature if PCorrMature[mirna] < FDR]
MatureNeutralCorr = [mirna for mirna in PCorrMature if PCorrMature[mirna] >= FDR]
MatureNegativeCorr, MaturePositiveCorr = NaturalSelection(MKmature, MatureSignifCorr)

print('significant after correction', len(HairpinSignifCorr), len(MatureSignifCorr))
print('positive after correction', len(HairpinPositiveCorr), len(MaturePositiveCorr))
print('negative after correction', len(HairpinNegativeCorr), len(MatureNegativeCorr))
print('neutral after correction', len(HairpinNeutralCorr), len(MatureNeutralCorr))


for mirna in MKhairpin:
    print(mirna, MKhairpin[mirna])
for mirna in MKmature:
    print(mirna, MKmature[mirna])

# make a summary file with results of the MK test for hairpin
newfile = open('MKtestmiRNAHairpinsNoSingleton.txt', 'w')
newfile.write('\t'.join(['Crm_miRNA', 'Pmirna', 'P4fold', 'Dmirna', 'D4fold', 'MK_P', 'Selection', 'MK_P_Corr', 'Selection_Corr']) + '\n')
for mirna in MKhairpin:
    newfile.write('\t'.join([mirna, str(MKhairpin[mirna][0]), str(MKhairpin[mirna][1]), str(MKhairpin[mirna][2]), str(MKhairpin[mirna][3]), str(MKhairpin[mirna][4])]) + '\t')
    if mirna in HairpinNeutral:
        newfile.write('neutral' + '\t')
    elif mirna in HairpinPositive:
        newfile.write('adaptive' + '\t')
    elif mirna in HairpinNegative:
        newfile.write('negative' + '\t')
    newfile.write(str(PCorrHairpin[mirna]) + '\t')
    if mirna in HairpinNeutralCorr:
        newfile.write('neutral' + '\n')
    elif mirna in HairpinPositiveCorr:
        newfile.write('adaptive' + '\n')
    elif mirna in HairpinNegativeCorr:
        newfile.write('negative' + '\n')
newile.close()
print('generated summary table hairpins')

# make a summary file with results of the MK test for mature
newfile = open('MKtestmiRNAmatureNoSingleton.txt', 'w')
newfile.write('\t'.join(['Crm_miRNA', 'Pmirna', 'P4fold', 'Dmirna', 'D4fold', 'MK_P', 'Selection', 'MK_P_Corr', 'Selection_Corr']) + '\n')
for mirna in MKmature:
    newfile.write('\t'.join([mirna, str(MKmature[mirna][0]), str(MKmature[mirna][1]), str(MKmature[mirna][2]), str(MKmature[mirna][3]), str(MKmature[mirna][4])]) + '\t')
    if mirna in MatureNeutral:
        newfile.write('neutral' + '\t')
    elif mirna in MaturePositive:
        newfile.write('adaptive' + '\t')
    elif mirna in MatureNegative:
        newfile.write('negative' + '\t')
    newfile.write(str(PCorrMature[mirna]) + '\t')
    if mirna in MatureNeutralCorr:
        newfile.write('neutral' + '\n')
    elif mirna in MaturePositiveCorr:
        newfile.write('adaptive' + '\n')
    elif mirna in MatureNegativeCorr:
        newfile.write('negative' + '\n')
newile.close()
print('generated summary table matures')


# group mirnas based on level of conservation
Restricted = [mirna for mirna in famCons if famCons[mirna] == 'CrmCla']
Novel = [mirna for mirna in famCons if famCons[mirna] == 'Crm']
Conserved = [mirna for mirna in famCons if famCons[mirna] == 'Caeno']
print('restricted', len(Restricted))
print('novel', len(Novel))
print('conserved', len(Conserved))

# perform MK test based on entire conservation group
RestrictedCounts, NovelCounts, ConservedCounts = {'restricted': [0, 0, 0, 0]}, {'novel': [0, 0, 0, 0]}, {'conserved': [0, 0, 0, 0]} 
for mirna in MatureCounts:
    if mirna in Restricted:
        for i in range(4):
            RestrictedCounts['restricted'][i] += MatureCounts[mirna][i]
    elif mirna in Novel:
        for i in range(4):
            NovelCounts['novel'] += MatureCounts[mirna][i]
    elif mirna in Conserved:
        for i in range(4):
            ConservedCounts['conserved'] += MatureCounts[mirna][i]
print('sorted mirnas based on conservation')

MKRestricted = MK_test(RestrictedCounts, 'fisher')
MKNovel = MK_test(NovelCounts, 'fisher')
MKConserved = MK_test(ConservedCounts, 'fisher')

print(MKRestricted)
print(MKNovel)
print(MKConserved)


# make a dict with mirna for each conservation group
# use this dict to sample mirna for calulting distribution of alpha
# compute alpha
# plot alpha for each group?
# make summary file with each group




# make summary file

# see Lyu et al for organizing table




# site type site number D P PDAF.5% D/PDAF.5% MK test pvaluea ab (% of adaptive fixations)



# plot alpha for the different conservation group













