# -*- coding: utf-8 -*-
"""
Created on Mon May 30 20:23:16 2016

@author: Richard
"""

# use this script to generate a summary 




from protein_divergence import *
from chemoreceptors import *
import os




# get the list of codeml output file in each parition
TM_out = [filename for filename in os.listdir('./Partitions/Membrane/') if '.out' in filename]
Out_out = [filename for filename in os.listdir('./Partitions/Outside/') if '.out' in filename]
In_out = [filename for filename in os.listdir('./Partitions/Inside/') if '.out' in filename]
Extra_out = [filename for filename in os.listdir('./Partitions/Extra_membrane/') if '.out' in filename]
print('generated file lists')


# create dictionnaries to store the divergence of each partition {gene1 :[dN, dS, omega], gene2 : [dN, dS, omega]}
TM_diverg, Out_diverg, In_diverg, Extra_diverg = {}, {}, {}, {}

# loop over codeml output files
for filename in TM_out:
    # parse codeml outputfile, return a dict {gene :[ dN, dS. omega]}
    gene_diverg = parse_codeml_output('./Partitions/Membrane/' + filename)
    # populate divergence dict
    for gene in gene_diverg:
        TM_diverg[gene] = list(gene_diverg[gene])
print('got divergence in Transmembrane')

for filename in Out_out:
    # parse codeml outputfile, return a dict {gene :[ dN, dS. omega]}
    gene_diverg = parse_codeml_output('./Partitions/Outside/' + filename)
    # populate divergence dict
    for gene in gene_diverg:
        Out_diverg[gene] = list(gene_diverg[gene])
print('got divergence in Outside domain')

for filename in In_out:
    # parse codeml outputfile, return a dict {gene :[ dN, dS. omega]}
    gene_diverg = parse_codeml_output('./Partitions/Inside/' + filename)
    # populate divergence dict
    for gene in gene_diverg:
        In_diverg[gene] = list(gene_diverg[gene])
print('got divergence in Inside domain')

for filename in Extra_out:
    # parse codeml outputfile, return a dict {gene :[ dN, dS. omega]}
    gene_diverg = parse_codeml_output('./Partitions/Extra_membrane/' + filename)
    # populate divergence dict
    for gene in gene_diverg:
        Extra_diverg[gene] = list(gene_diverg[gene])
print('got divergence in Extramembrane')

# open file to store divergence in the different partitions
# do not consider partitions of genes that do not have TM and extra-TM domains
newfile = open('Chemo_divergence_partitions.txt', 'w')
# write header to file
header = '\t'.join(['Gene', 'TM_dN', 'TM_dS', 'TM_dN/dS', 'ExtraTM_dN',
                    'ExtraTM_dS', 'ExtraTM_dN/dS', 'Outside_dN',
                    'Outside_dS', 'Outside_dN/dS', 'Inside_dN',
                    'Inside_dS', 'Inside_dN/dS'])
newfile.write(header + '\n')

# loop over genes in TM_diverg
for gene in TM_diverg:
    line = [gene, str(TM_diverg[gene][0]), str(TM_diverg[gene][1]), str(TM_diverg[gene][2])]
    # check if gene has extra-cellular domain
    if gene in Extra_diverg:
        # gene has extra-transmembrane domain, add divergence to line
        for i in range(len(Extra_diverg[gene])):
            line.append(str(Extra_diverg[gene][i]))
    else:
        # gene does not have extra-membrane domain, add 'NA' to divergence
       line.extend(['NA', 'NA', 'NA'])
    # check if gene has outside domain
    if gene in Out_diverg:
        # gene has extra-cellular domain
        for i in range(len(Out_diverg[gene])):
            line.append(str(Out_diverg[gene][i]))
    else:
        # gene does not have extra-cellular domain
        line.extend(['NA', 'NA', 'NA'])        
    # check if gene has intra-cellular domain
    if gene in In_diverg:
        # gene has intra-cellular domain
        for i in range(len(In_diverg[gene])):
            line.append(str(In_diverg[gene][i]))
    else:
        # gene does not have intra-cellular domain, add 'NA' to line
        line.extend(['NA', 'NA', 'NA'])
    # write line to file
    line = '\t'.join(line)
    newfile.write(line + '\n')
    
# close file after writing
newfile.close()


