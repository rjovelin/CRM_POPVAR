# -*- coding: utf-8 -*-
"""
Created on Wed May 11 11:29:46 2016

@author: RJovelin
"""

# use this script to filter weird divergence values in the remanei-latens codeml divergence file

# set threshold for dN, dS and omega, remove genes instead of removing values
# genes are removed if any divergence estimate is greater than following thresholds
dN = 2
dS = 2
omega = 7

# open file for reading
infile = open('CRM_CLA_protein_divergence.txt')
header = infile.readline().rstrip().split('\t')
# create a dict to store file cintent {gene : [data]}
genes = {}
for line in infile:
    line = line.rstrip()
    if line != '':
        line = line.split('\t')
        genes[line[0]] = line
infile.close()

# create a set of genes to delete
to_delete = set()
# loop over genes, delete genes if divergence is greater than threshold 
# delete genes for which omega is not defined
for i in genes:
    if genes[i][-1] == 'NA':
        to_delete.add(i)
    else:
        # compare divergence to trheshold
        if float(genes[i][4]) > dN:
            to_delete.add(i)
        if float(genes[i][5]) > dS:
            to_delete.add(i)
        if float(genes[i][6]) > omega:
            to_delete.add(i)

print('remove {0} genes'.format(len(to_delete)))
# delete genes
for i in to_delete:
    del genes[i]

# save data to new file
newfile = open('CRM_CLA_ProtDiverg_FILTERED.txt', 'w')
newfile.write('\t'.join(header) + '\n')
for i in genes:
    newfile.write('\t'.join(genes[i]) + '\n')
newfile.close()
