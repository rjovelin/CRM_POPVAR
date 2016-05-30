# -*- coding: utf-8 -*-
"""
Created on Sun Jun 21 11:17:26 2015

@author: Richard
"""

# use this script to predict membrane domains in chemoreceptors with Phoebius


from manipulate_sequences import *
from chemoreceptors import *
import os


# get the set of chemoreceptors from the iprscan outputfile
chemo = get_chemoreceptors('../Genome_Files/PX356_protein_seq.tsv') 
print('got chemo genes')

# create a set of valid transcripts
transcripts = get_valid_transcripts('../Genome_Files/unique_transcripts.txt')
print('got valid genes')

# create a set of valid chemoreceptors
GPCRs = set(gene for gene in chemo if gene in transcripts)
print('got valid GPCR genes')

# create a dict with the remanei CDS
CDS = convert_fasta('../Genome_Files/noamb_PX356_all_CDS.fasta')
print('got CDS sequences')

# create a directory to store the individual fasta files
os.mkdir('Chemo_genes')

# loop over genes in CDS
for gene in CDS:
    # check if gene if chemo and valid transcript
    if gene in GPCRs:
        # check that protein sequence doesn't have stop codon
        protein = cds_translate(CDS[gene])
        if '*' not in protein:
            # open file for writing
            newfile = open('Chemo_genes/' + gene + '_prot.fasta', 'w')
            # save protein to file
            newfile.write('>' + gene + '\n')
            newfile.write(protein + '\n')
            newfile.close()
print('generated protein fasta files')


# create list of fasta files
files = os.listdir('./Chemo_genes/')
print('got list of protein sequences')

# loop over files
for filename in files:
    #check that filename is fasta
    if '_prot.fasta' in filename:
        # get gene name
        gene = filename[:filename.index('_prot')]
        print(filename, gene, sep = '\t')
        # build command (phobius -png proba.png -gp proba.png -plp proba.txt input_file)
        os.system('phobius' + ' -png ' + './Chemo_genes/' + gene + '.png' + ' -gp ' + './Chemo_genes/' + gene + '.png' + ' -plp ' + './Chemo_genes/' + gene + '_proba.txt' + ' ./Chemo_genes/' + filename) 
print('done predicting membrane domains with phoebius')        
                
                    
                             
