# -*- coding: utf-8 -*-
"""
Created on Fri May  8 17:01:03 2015

@author: Richard
"""

#!/usr/bin/env python3




from protein_divergence import *
from accessories import *
import os




# clean the reciprocal BLAST file to generate a file with 1:1 orthologs
clean_orthologs_file('mutual_blast_besthit.tab', '356_10172014.gff3', '534_10172014.gff3', 'crem_cla_orthologs.txt')

# create a directory named pairs
os.mkdir('./pairs/')

# save orthologous sequences into separate files
save_orthologous_seq_pairs_to_file('crem_cla_orthologs.txt', 'noamb_PX356_all_CDS.fasta', 'noamb_PX534_all_CDS.fasta', './pairs/')
print('done generating fasta files')

# run t-coffee
os.system('perl runTcoffee.pl')
print('done aligning sequences')


# make a set of filenames for which sequences have internal stop codons
internal_stops = set()
# get a list of files
files = os.listdir('./')
# loop over files, convert fasta to dict, translate CDS sequences, check for internal stop codons
for filename in files:
    # check if file is t-coffee alignment file
    if filename[-8:] == '_aln.tfa':
        # convert file to dict
        sequences = convert_fasta(filename)
        # translate each CDS sequence
        for gene in sequences:
            protein = cds_translate(sequences[gene])
            # check for internal stop codons
            if '*' in protein:
                # if yes, add the filename to set of files with internal stops
                internal_stops.add(filename)
print('number of ortholog pairs with internal stop codons: ', len(internal_stops))

# convert alignment files to PAML format
# get a list with the file names
files = os.listdir('./')
for filename in files:
    # exclude files for which genes have internal stop codons
    if filename not in internal_stops:
        alignment_file_to_PAML_format(filename, './')
print('done generating codeml input files')

# generate codeml control files
# get the list of alignment files
files = os.listdir('./')
for filename in files:
    if filename[-4:] == '.txt' and 'CRE' in filename and 'tre' not in filename:
        generate_codeml_control_file(filename, 'codeml.ctl', './')
print('done generating codeml control files')

# run codeml
# get the the list of control files
files = os.listdir('./')
for filename in files:
    if '.ctl' in filename and 'CRE' in filename:
        os.system('codeml ' + filename)
print('done running codeml')

# save divergent results to file
save_divergence_results_to_file('356_10172014.gff3', '534_10172014.gff3', 'CRM_CLA_protein_divergence.txt', './')
print('done computing sequence divergence')



