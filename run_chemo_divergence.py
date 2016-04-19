# -*- coding: utf-8 -*-
"""
Created on Mon Jun 22 15:36:59 2015

@author: Richard
"""

from protein_divergence import *
from chemoreceptors import *
import os


# get the alignment files in differnt partitions
TM_files = [filename for filename in os.listdir('./Partitions/Membrane/') if 'TM' in filename]
Out_files = [filename for filename in os.listdir('./Partitions/Outside/') if 'outside' in filename]
In_files = [filename for filename in os.listdir('./Partitions/Inside/') if 'inside' in filename]
Extra_files =[filename for filename in os.listdir('./Partitions/Extra_membrane/') if 'ExtraTM' in filename]

# generate codeml control files for each alignment files in each partition folder
for alignment_file in TM_files:
    generate_codeml_control_file(alignment_file, '../CREM_CLA_protein_divergence/codeml.ctl', './Partitions/Membrane/')
for alignment_file in Out_files:
    generate_codeml_control_file(alignment_file, '../CREM_CLA_protein_divergence/codeml.ctl', './Partitions/Outside/')
for alignment_file in In_files:
    generate_codeml_control_file(alignment_file, '../CREM_CLA_protein_divergence/codeml.ctl', './Partitions/Inside/')
for alignment_file in Extra_files:
    generate_codeml_control_file(alignment_file, '../CREM_CLA_protein_divergence/codeml.ctl', './Partitions/Extra_membrane/')

# run codeml on each partition
# get the the list of control files
TM_ctl = [filename for filename in os.listdir('./Partitions/Membrane/') if 'ctl' in filename]
Out_ctl = [filename for filename in os.listdir('./Partitions/Outside/') if 'ctl' in filename]
In_ctl = [filename for filename in os.listdir('./Partitions/Inside/') if 'ctl' in filename]
Extra_ctl = [filename for filename in os.listdir('./Partitions/Extra_membrane/') if 'ctl' in filename]

# copy tree files to directories
os.system('cp CREMCLATREE.tre.txt ./Partitions/Membrane/')
os.system('cp CREMCLATREE.tre.txt ./Partitions/Outside/')
os.system('cp CREMCLATREE.tre.txt ./Partitions/Inside/')
os.system('cp CREMCLATREE.tre.txt ./Partitions/Extra_membrane/')

# change directory to run codeml
os.chdir('./Partitions/Membrane/')
# loop over codml control files, run codeml
for filename in TM_ctl:
    os.system('codeml ' + filename)
# remove tree file
os.system('rm CREMCLATREE.tre.txt')

# change directory
os.chdir('../Outside')
# loop over control files, run codeml
for filename in Out_ctl:
    os.system('codeml ' + filename)
# remove tree file
os.system('rm CREMCLATREE.tre.txt')

# change directory
os.chdir('../Inside')
# loop over control files, run codeml
for filename in In_ctl:
    os.system('codeml ' + filename)
# remove tree file
os.system('rm CREMCLATREE.tre.txt')
    
# change directory
os.chdir('../Extra_membrane/')
# loop over control files, run codeml
for filename in Extra_ctl:
    os.system('codeml ' + filename)
# remove tree file
os.system('rm CREMCLATREE.tre.txt')


# change directory
os.chdir('../../')


# get the list of codeml output file in each parition
TM_out = [filename for filename in os.listdir('./Partitions/Membrane/') if '.out' in filename]
Out_out = [filename for filename in os.listdir('./Partitions/Outside/') if '.out' in filename]
In_out = [filename for filename in os.listdir('./Partitions/Inside/') if '.out' in filename]
Extra_out = [filename for filename in os.listdir('./Partitions/Extra_membrane/') if '.out' in filename]


# create dictionnaries to store the divergence of each partition {gene1 :[dN, dS, omega], gene2 : [dN, dS, omega]}
TM_diverg, Out_diverg, In_diverg, Extra_diverg = {}, {}, {}, {}

# loop over codeml output files
for filename in TM_out:
    # parse codeml outputfile, return a dict {gene :[ dN, dS. omega]}
    gene_diverg = parse_codeml_output('./Partitions/Membrane/' + filename)
    # populate divergence dict
    for gene in gene_diverg:
        TM_diverg[gene] = list(gene_diverg[gene])

for filename in Out_out:
    # parse codeml outputfile, return a dict {gene :[ dN, dS. omega]}
    gene_diverg = parse_codeml_output('./Partitions/Outside/' + filename)
    # populate divergence dict
    for gene in gene_diverg:
        Out_diverg[gene] = list(gene_diverg[gene])

for filename in In_out:
    # parse codeml outputfile, return a dict {gene :[ dN, dS. omega]}
    gene_diverg = parse_codeml_output('./Partitions/Inside/' + filename)
    # populate divergence dict
    for gene in gene_diverg:
        In_diverg[gene] = list(gene_diverg[gene])

for filename in Extra_out:
    # parse codeml outputfile, return a dict {gene :[ dN, dS. omega]}
    gene_diverg = parse_codeml_output('./Partitions/Extra_membrane/' + filename)
    # populate divergence dict
    for gene in gene_diverg:
        Extra_diverg[gene] = list(gene_diverg[gene])

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
       line.append('NA')
       line.append('NA')
       line.append('NA')
    # check if gene has outside domain
    if gene in Out_diverg:
        # gene has extra-cellular domain
        for i in range(len(Out_diverg[gene])):
            line.append(str(Out_diverg[gene][i]))
    else:
        # gene does not have extra-cellular domain
        line.append('NA')
        line.append('NA')
        line.append('NA')
    # check if gene has intra-cellular domain
    if gene in In_diverg:
        # gene has intra-cellular domain
        for i in range(len(In_diverg[gene])):
            line.append(str(In_diverg[gene][i]))
    else:
        # gene does not have intra-cellular domain, add 'NA' to line
        line.append('NA')
        line.append('NA')
        line.append('NA')
    # write line to file
    line = '\t'.join(line)
    newfile.write(line + '\n')
    
# close file after writing
newfile.close()


