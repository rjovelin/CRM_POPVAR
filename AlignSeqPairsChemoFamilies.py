# -*- coding: utf-8 -*-
"""
Created on Thu Jun 25 21:34:28 2015

@author: Richard
"""


from chemoreceptors import *
from manipulate_sequences import *
from tcoffee_alignment import *
import os
from divergence import *
import numpy as np


# use this script to align pairs of protein sequences within each chemoreceptor family 


# assign genes to chemoreceptor families
chemo = chemo_families('../Genome_Files/PX356_protein_seq.tsv')
print('got chemo genes')

# remove ambiguous genes belonging to multiple families
chemo = remove_ambiguous_chemoreceptors(chemo)
print('removed genes assigned to multiple families')

# create a set of valid transcripts
transcripts = get_valid_transcripts('../Genome_Files/unique_transcripts.txt')
print('got valid transcripts')

# get the CDS sequences
CDS = convert_fasta('../Genome_Files/noamb_PX356_all_CDS.fasta')
print('got CDS sequences')

# create new directory
os.mkdir('Pairwise_Chemos')

# remove genes that are not in transcripts
chemo = clean_chemo_families(chemo, '../Genome_Files/unique_transcripts.txt')
print('got valid chemo genes')

# loop over families
for family in chemo:
    # create direcory with family name
    os.mkdir('./Pairwise_Chemos/' + family + '_family' + '/')
    # create a list of genes for the given family
    genes = [i for i in chemo[family]]
    # loop over genes in list, create files with pairs of protein sequences
    for i in range(len(genes)-1):
        for j in range(i+1, len(genes)):
            # create new file
            newfile = open('./Pairwise_Chemos/' + family + '_family/' + genes[i] + '_o_' + genes[j] + '.fasta', 'w')
            newfile.write('>' + genes[i] + '\n')
            newfile.write(cds_translate(CDS[genes[i]]) + '\n')
            newfile.write('>' + genes[j] + '\n')
            newfile.write(cds_translate(CDS[genes[j]]) + '\n')
            newfile.close()
print('generated sequence pairs fasta file')


# move to Pairwise folder
os.chdir('./Pairwise_Chemos/')

# run t-coffee on all files in each folder
# make a list of folders
folders = [folder for folder in os.listdir() if '_family' in folder]
# sort folders
folders.sort()
print('made list of folders')


# loop over directories, align pairs of sequences
for folder in folders:
    # move to folder so that tcoffee output files are saved in folder
    os.chdir(folder)
    # create a list of filenames
    files = [filename for filename in os.listdir() if 'fasta' in filename]
    print(folder, len(files))
    # run tcoffee
    for filename in files:
        os.system('t_coffee ' + filename)
    print('done aligning sequences for', folder)
    # move back to parent directory
    os.chdir('../')


# loop over directories, convert t-coffee format to fasta
for folder in folders:
    # create a list of t-coffee output files
    alignments = [filename for filename in os.listdir('./Pairwise_Chemos/' + folder) if '.aln' in filename]
    print(folder, len(alignments))
    # convert t-coffee format to text files
    for filename in alignments:
        convert_tcoffee_prot_to_fasta('./Pairwise_Chemos/' + folder + '/' + filename)
    print('fasta convertion done')

print('done aligning sequences')

