# -*- coding: utf-8 -*-
"""
Created on Thu May 28 10:28:18 2015

@author: Richard
"""

#!/usr/bin/env python3

from accessories import *
import os

# Note that the elegans GFF file is different from the remanei and latens GFF files


# use this function to get the coordinates of the C. elegans 3' UTRs
def celegans_annotated_three_prime_coordinates(celegans_gff):
    '''
    (file) -> dict
    Returns a dictionnary with the 3' UTR coordinates of each celegans transcript
    Use a single UTR if multiple UTRs are annotated
    '''

    # create a dictionnary to stote the coordinates {transcript_name : [chromo, start, end, orienation]}
    annotated_three_prime = {}

    # open file for reading
    gff = open(celegans_gff, 'r')
    for line in gff:
        line = line.rstrip()
        if line != '':
            line = line.split()
            if line[1] == 'WormBase':
                if line[2] == 'three_prime_UTR':
                    transcript = line[8][line[8].index('Transcript:')+11 :]
                    chromo = line[0]
                    start = int(line[3])
                    end = int(line[4])
                    orientation = line[6]
                    annotated_three_prime[transcript] = [chromo, start, end, orientation]
    gff.close()
    return annotated_three_prime


# use this function to create a dict of transcript : gene pairs
def celegans_transcript_to_gene(celegans_gff):
    '''
    (file) -> dict
    Returns a dictionnary with celegans transcript : gene pairs from the gff annotation file
    '''
    #create a dictionnary of transcript : gene pairs
    transcripts_genes = {}

    # open file for reading
    gff = open(celegans_gff, 'r')
    for line in gff:
        line = line.rstrip()
        if line != '':
            line = line.split()
            if line[1] == 'WormBase':
                if line[2] == 'mRNA':
                    transcript = line[8][line[8].index('Transcript:')+11 : line[8].index(';')]
                    gene = line[8][line[8].index('Parent=Gene:')+12 : line[8].index(';', line[8].index('Parent'))]
                    transcripts_genes[transcript] = gene
    gff.close()
    return transcripts_genes


# use this function to create a dict if gene : list of rsnacripts airs
def celegans_gene_to_transcripts(celegans_gff):
    '''
    (file) -> dict
    Returns a dictionnary with celegans gene as key and a list of transcripts as value
    '''

    # get the dictionnary of transcripts : gene names pairs
    transcripts_genes = celegans_transcript_to_gene(celegans_gff)

    # create a reverse dictionnary of gene : [transcripts] pairs
    genes = {}
    for transcript in transcripts_genes:
        gene_name = transcripts_genes[transcript]
        if gene_name in genes:
            genes[gene_name].append(transcript)
        else:
            genes[gene_name] = [transcript]

    return genes


# use this function to get the UTR sequences of each transcript
def celegans_UTR_sequences(celegans_gff, assembly):
    '''
    (file, file, int) -> dict
    Returns a dictionnary with the transcript name as key and the UTR sequence
    from the assembly as value. Ignore transcripts lacking a UTR sequence
    (e.g. at the ends of chromosome or when adjacent to another gene)
    '''

    # convert assembly to fasta format
    genome = convert_fasta(assembly)
    
    # create a dictionnary with UTR coordinates
    three_prime = celegans_annotated_three_prime_coordinates(celegans_gff)
    
    # create a dict {transcript: sequence}
    UTR = {}
    # loop over transcripts
    for transcript in three_prime:
        # ignore transcripts that do not have a UTR (start = end)
        if not three_prime[transcript][1] == three_prime[transcript][2]:
            # get chromo
            chromo = three_prime[transcript][0]
            # get orientation
            orientation = three_prime[transcript][-1]
            # convert start to 0-based index (start  = start -1)
            start = three_prime[transcript][1] - 1
            # no need to convert end because end non-inclusive
            end = three_prime[transcript][2]
            # slice the squence
            UTR_seq = genome[chromo][start:end]
            # check orientation
            if orientation == '+':
                # populate dict
                UTR[transcript] = UTR_seq.upper()
            elif orientation == '-':
                # take reverse complement
                UTR_seq = reverse_complement(UTR_seq)
                # populate dict
                UTR[transcript] = UTR_seq.upper()
                
    return UTR


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


