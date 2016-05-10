# -*- coding: utf-8 -*-
"""
Created on Fri Jun  5 12:26:07 2015

@author: Richard
"""
from accessories import *
import os


# open file file for reading
infile = open('cds_codon_table.txt', 'r')

# open file for writing
newfile = open('CDS_SNP_DIVERG.txt', 'w')

# get header
header = infile.readline().rstrip().split()
# insert column with KSR+PX reference allele
header.insert(13, 'KSR+PX_REF')
# add column with KSR+PX alternative allele
header.append('KSR+PX_ALT')

# write header to newfile
for item in header[:-1]:
    newfile.write(item + '\t')
newfile.write(header[-1] + '\n')

# loop over file
for line in infile:
    line = line.rstrip()
    if line != '':
        line = line.split()
        # compute the number of reference allele in KSR + PX
        ref = int(line[10]) + int(line[11])
        # compute the number of alternative alleles in KSR + PX
        alt = int(line[13]) + int(line[14])
        # insert ref in line
        line.insert(13, str(ref))
        # add alt to line
        line.append(str(alt))
        # write line to new file
        for item in line[:-1]:
            newfile.write(item + '\t')
        newfile.write(line[-1] + '\n')

# close files
infile.close()
newfile.close()
print('table updated with KSR+PX strains')


# make a list of alignment files
files = os.listdir('./pairs/Aligned_pairs/')
# make a list of alignment files
alignments = [filename for filename in files if filename[-8:] == '_aln.tfa']

# make a dict with crm and cla sequences
orthologs = {}
for filename in alignments:
    gene_name = filename[:filename.index('_CLA')]
    sequences = {}
    orthologs[gene_name] = {}
    sequences = convert_fasta('./pairs/Aligned_pairs/' + filename)
    for seq in sequences:
        if 'CRE' in seq:
            orthologs[gene_name]['crm'] = sequences[seq]
        else:
            orthologs[gene_name]['cla'] = sequences[seq]
print('dict with orthologs done')

# make a dict with the remanei CDS sequences from the snp file
rem_CDS = {}

# open file for reading
infile = open('CDS_SNP_DIVERG.txt', 'r')
# skip header
infile.readline()
# loop over file, populate dict with the crem CDS sequence {gene: CDS}
for line in infile:
    line = line.rstrip()
    if line != '':
        line = line.split()
        # get transcript name
        gene = line[2]
        # build CDS sequence
        if gene in rem_CDS:
            rem_CDS[gene] += line[5]
        else:
            rem_CDS[gene] = line[5]
# close file after reading
infile.close()
print('dict with remanei CDS sequences done')


# stop codons were removed in the alignments but are present in the SNP file
# add fake nucleotides at the end of remanei and latens CDS if stop codon is missing
# to bring the sequences to same length
for gene in orthologs:
    # check that gene has snps and is in the rem_CDS dict
    if gene in rem_CDS:
        # check if the remanei sequence in table has a stop codon
        if rem_CDS[gene][-3:] in {'TAG', 'TGA', 'TAA'}:
            # check if the rem sequences without stop codons differ
            if orthologs[gene]['crm'].replace('-', '') != rem_CDS[gene][:-3]:
                # if sequences differ raise ValueError
                raise ValueError('CDS sequences with stop are different')
            # check if the rem sequences have same length
            if len(orthologs[gene]['crm'].replace('-', '')) != len(rem_CDS[gene]):
                # if different length, add fake stop codon
                orthologs[gene]['crm'] += 'XXX'
                orthologs[gene]['cla'] += 'XXX'
        elif rem_CDS[gene][-3:] not in {'TAG', 'TGA', 'TAA'}:
            # check if rem sequences diff
            if orthologs[gene]['crm'].replace('-', '') != rem_CDS[gene]:
                raise ValueError('CDS sequences without STOP are different')
                # no need to do anything to the sequence if stop codon is not there
print('done adding fake stop codons')

# create a dict with remanei genes as key and a dict of codon_position: list of codons as value
# {gene: {i: [crem_codon, cla_codon]}}
codons = {}
for gene in orthologs:
    # get the remanei CDS sequence
    crm_cds = orthologs[gene]['crm']
    # get the latens CDS sequence
    cla_cds = orthologs[gene]['cla']
    # initiate gap variable
    gaps = 0
    # initialise j
    j = 0
    # initiate innder dict
    codons[gene] = {}
    # loop over the remanei CDS sequence
    for i in range(0, len(crm_cds), 3):
        if '-' in crm_cds[i:i+3]:
            gaps += crm_cds[i:i+3].count('-')
        else:
            j = i - gaps
            codons[gene][j] = [crm_cds[i:i+3], cla_cds[i:i+3]]
print('dict with codon positions done')

# open file for reading
infile = open('CDS_SNP_DIVERG.txt', 'r')
# create string to store the content of newfile
content = ''
# get header
header = infile.readline().rstrip().split()
# add columns to header
header.append('CLA_codon')
header.append('CLA_ref')
header.append('CLA_base_index')
header.append('CLA_codon_index')
header.append('CLA_position_in_codon')
header = '\t'.join(header)
# add header to content
content += (header + '\n')

# create a variable to track the gene ID
gene_ID = ''

# loop over file
for line in infile:
    # make a list of line
    line = line.rstrip()
    if line != '':
        line = line.split()
        # get gene name
        gene = line[2]        
        
        # check if gene has an ortholog with aligned sequences
        if gene not in orthologs:
            # no ortholog and ancestral state, add '???' for codon and '?' for ref
            line.append('???')
            line.append('?')
            # make a string of line
            line = '\t'.join(line)
            # add line to content
            content += (line + '\n')
        elif gene in orthologs:
            # ortholog exists, get the ancestral state
                       
            # make a list of codon positions, sort list
            positions = [pos for pos in codons[gene]]
            positions.sort()

            # check if we are looping in the same gene or if we start with a new gene
            if gene != gene_ID:
                # initialise counters
                i = 0
                j = 0
                k = 1
                # start a new gene
                # add first CLA codon
                line.append(codons[gene][j][1])
                # get the latens sequence gaped positions excluded
                cla_seq = ''.join([codons[gene][pos][1] for pos in positions])
                # add first nucleotide
                line.append(cla_seq[i])
                
                line.append(str(i))
                line.append(str(j))
                line.append(str(k))                
                
                # make string of line
                newline = '\t'.join(line)
                # add line to content
                content += (newline + '\n')
                # update i
                i += 1
                # update gene_ID
                gene_ID = gene
            elif gene == gene_ID:
                # still the same gene
                if k != 3:
                    line.append(codons[gene][j][1])
                    k += 1
                elif k ==3:
                    k = 1
                    j += 3
                    line.append(codons[gene][j][1])
                
                line.append(cla_seq[i])
                line.append(str(i))
                line.append(str(j))
                line.append(str(k))
                newline = '\t'.join(line)
                content += (newline + '\n')
                i += 1
                    
infile.close()
print('done adding ancestral state')


newfile = open('CDS_SNP_DIVERG.txt', 'w')
newfile.write(content)
newfile.close()
print('snp file updated')



