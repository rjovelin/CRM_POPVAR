# -*- coding: utf-8 -*-
"""
Created on Thu May 28 10:28:18 2015

@author: Richard
"""

#!/usr/bin/env python3

from Manipulate_Sequences import *
import os

# use this function to clean up the file of 1:1 orthologs from transcripts that belong to a same gene 
def clean_cel_crem_orthologs_file(best_blast_hits, crem_gff, cel_gff, outputfile):
    '''
    (file, file, file) -> file
    Clean up the file with best blat hits from transcripts that belong to a same
    gene. Write to the outputfile the orthologous pairs between remanei and elegans
    by keeping a single transcript per gene in each species
    '''
    
#    # create dicts with {gene: [transcript1, transcript2]} in remanei and elegans
#    remanei = parent_gene(crem_gff)
#    elegans = celegans_gene_to_transcripts(cel_gff)
#    
    # create sets of genes for which transcripts are already used in ortholog pairs
    rem_already_used = set()
    cel_already_used = set()
    
    # create dicts of transcript : gene pairs {transcript : gene}
    remanei_ts = transcript_to_gene(crem_gff)
    elegans_ts = celegans_transcript_to_gene(cel_gff)
    
    # create a dict to store the ortholog pairs {crem_TS : cel_TS}
    orthologs = {}    
    
    # open best blast hits file
    infile = open(best_blast_hits, 'r')
    # skip header
    infile.readline()
    
    # loop over file
    for line in infile:
        line = line.rstrip()
        line = line.split()
        cel_TS = line[0]
        crem_TS = line[1]
        # check that the remanei gene has not been used
        if remanei_ts[crem_TS] not in rem_already_used:
            # check that the latens gene has not been used
            if elegans_ts[cel_TS] not in cel_already_used:
                # populate ortholog dict
                orthologs[crem_TS] = cel_TS
                # add corresponding genes to sets
                rem_already_used.add(remanei_ts[crem_TS])
                cel_already_used.add(elegans_ts[cel_TS])
    
    
    # open file for writing
    newfile = open(outputfile, 'w')
    # write header
    newfile.write('# 1:1 orthologs between remanei and elegans obtained by reciprocal BLAST hits\n')
    newfile.write('# genes are represented by a single transcript\n')
    newfile.write('# ambiguous sites in the reference genome are fixed\n')
    newfile.write('Cremanei' + '\t' + 'Celegans' + '\n')
    
    # write content to file
    for gene in orthologs:
        newfile.write(gene + '\t' + orthologs[gene] + '\n')
    
    # close files
    infile.close()
    newfile.close()


# use this function to generate the TargetScan input sequence file
def make_cremanei_celegans_targetscan_input_sequences(directory, outputfile):
    '''
    (str, file) -> file
    Generate the input sequence file for TargetScan using the aligned UTR sequences
    located in directory. Use only Cremanei and Celegans orthologs      
    '''

    # open outputfile for writing
    newfile = open(outputfile, 'w')

    # grab the aligned remanei and latens sequences and write to newfile
    files = os.listdir(directory)
    for filename in files:
        # check that the file contains the crem-cel orthologs
        if 'CRE' in filename and '.txt' in filename:
            # convert fasta file to a dictionnary for seq_name : sequence
            fasta_seq = convert_fasta(directory + filename)
            # get the name of the remanei transcript
            for seq_name in fasta_seq:
                if 'CRE_PX356' in seq_name:
                    transcript_name = seq_name
            # use the cremanei transcript name for the elegans ortholog
            # but add species ID to distinguish species
            # convert T to U
            for seq_name in fasta_seq:
                if 'CRE_PX356' in seq_name:
                    newfile.write(transcript_name + '\t' + '31234' + '\t' + fasta_seq[seq_name].upper().replace('T', 'U') + '\n')
                else:
                    newfile.write(transcript_name + '\t' + '6239' + '\t' + fasta_seq[seq_name].upper().replace('T', 'U') + '\n')
    
    
    newfile.close()


