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
    # loop over genes in list, create files with pairs of protei sequences
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


# change directory
os.chdir('./Pairwise_Chemos/')

# run t-coffee on all files in each folder
# make a list of folders
folders = [folder for folder in os.listdir() if '_family' in folder]
# sort folders
folders.sort()

# loop over directories
for folder in folders:
    # change directory
    os.chdir(folder)
    # create a list of filenames
    files = [filename for filename in os.listdir() if 'fasta' in filename]
    print(folder, len(files))
    # run tcoffee
    for filename in files:
        os.system('t_coffee ' + filename)
    # create a list of t-coffee output files
    alignments = [filename for filename in os.listdir() if '.aln' in filename]
    # convert t-coffee format to text files
    for filename in alignments:
        convert_tcoffee_prot_to_fasta(filename)
print('done aligning sequences')


    
#    # create a list of filenames
#    files = [filename for filename in os.listdir() if '.txt' in filename]
#    # create a list to store the distances
#    distances = []
#    # loop over files, convert to fasta
#    for filename in files:
#        proteins = convert_fasta(filename)
#        # get the sequences
#        genes = [i for i in proteins]
#        seq1 = proteins[genes[0]]
#        seq2 = proteins[genes[1]]
#        # get p-distances
#        p_distance = pairwise_distance(seq1, seq2, 'protein') 
#        if p_distance != 'NA':
#            # multiply by 100 to get %
#            p_distance *= 100
#            distances.append(p_distance)
#    # create a histogram
#    hist = np.histogram(distances, range(0, 101, 10))
#    # get the folder name
#    folder_name = folder[:folder.index('_')] + '_hist'
#    # create new file in parent directory
#    newfile = open('../' + folder_name + '.txt', 'w')
#    for i in range(len(hist[0])):
#        newfile.write(str(hist[1][i]) + ':' + str(hist[1][i] + 10) + '\t' + str(hist[0][i]) + '\n')
#    # close file
#    newfile.close()
#    # go back to parent directory
#    os.chdir('../')
#    
#    
