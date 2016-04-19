# -*- coding: utf-8 -*-
"""
Created on Tue Jun 23 15:58:57 2015

@author: Richard
"""

from chemoreceptors import *
from accessories import *
import os



# assign genes to chemoreceptor families
chemo = chemo_families('./PX356_protein_seq.tsv')

# remove ambiguous genes belonging to multiple families
chemo = remove_ambiguous_chemoreceptors(chemo)

# create a set of valid transcripts
transcripts = get_valid_transcripts('../CREM_CLA_protein_divergence/unique_transcripts.txt')

# get the CDS sequences
CDS = convert_fasta('../CREM_CLA_protein_divergence/noamb_PX356_all_CDS.fasta')

# create new directory
os.mkdir('Chemo_families')
# change directory
os.chdir('./Chemo_families/')

# create a directory named pairs
os.mkdir('./pairs/')

# generate fasta files
for family in chemo:
    # create new alignment file
    newfile = open('./pairs/' + family + '.tfa', 'w')
    # write all CDS for gene in this family
    for gene in chemo[family]:
        # check that gene is in transcripts
        if gene in transcripts:
            # check that gene doesn't have stop codon in the reference
            if '*' not in cds_translate(CDS[gene]):
                # write CDS to file in fasta
                newfile.write('>' + gene + '\n')
                newfile.write(CDS[gene] + '\n')
    # close file
    newfile.close()

# copy Tcoffee perl script to current directory
os.system('cp ../../CREM_CLA_protein_divergence/runTcoffee.pl ./')
    

# run t-coffee on chemo families
os.system('perl runTcoffee.pl')
print('done aligning sequences')


#
## convert alignment files to PAML format
## get a list with the file names
#files = os.listdir('./')
#for filename in files:
#    # exclude files for which genes have internal stop codons
#    if filename not in internal_stops:
#        alignment_file_to_PAML_format(filename, './')
#print('done generating codeml input files')
#
#
#
#





























# remove ambiguous genes


# generate fasta files with same-family genes


# align sequences


# make phylogenetic tree


# define groups based on tree topologies



# generate codeml ctl files


# run the M7. M8 and M8a tests


# perform LRT tests


# generate summary table and generate figures


# distribution of dn/ds values at sites for site in membrane and site outside of membrane domains