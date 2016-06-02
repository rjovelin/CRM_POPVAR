# -*- coding: utf-8 -*-
"""
Created on Wed Jun  1 19:35:16 2016

@author: RJovelin
"""

from chemoreceptors import *
from manipulate_sequences import *
from tcoffee_alignment import *
import os
from divergence import *
import numpy as np


# use this script to generate a histogram with protein distance between pairs
# of chemoreceptor within each family


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
    
    # create a list of filenames
    files = [filename for filename in os.listdir() if '.txt' in filename]
    # create a list to store the distances
    distances = []
    # loop over files, convert to fasta
    for filename in files:
        proteins = convert_fasta(filename)
        # get the sequences
        genes = [i for i in proteins]
        seq1 = proteins[genes[0]]
        seq2 = proteins[genes[1]]
        # get p-distances
        p_distance = pairwise_distance(seq1, seq2, 'protein') 
        if p_distance != 'NA':
            # multiply by 100 to get %
            p_distance *= 100
            distances.append(p_distance)
    # create a histogram
    hist = np.histogram(distances, range(0, 101, 10))
    # get the folder name
    folder_name = folder[:folder.index('_')] + '_hist'
    # create new file in parent directory
    newfile = open('../' + folder_name + '.txt', 'w')
    for i in range(len(hist[0])):
        newfile.write(str(hist[1][i]) + ':' + str(hist[1][i] + 10) + '\t' + str(hist[0][i]) + '\n')
    # close file
    newfile.close()
    # go back to parent directory
    os.chdir('../')
    
    
