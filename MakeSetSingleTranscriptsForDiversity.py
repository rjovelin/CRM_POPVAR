# -*- coding: utf-8 -*-
"""
Created on Thu Jun 11 15:06:16 2015

@author: Richard
"""


from manipulate_sequences import *
from genomic_coordinates import *
from miRNA_target import *

# use this script to make a set of unique transcript per gene for analyses of 
# polymorphism 


# create a dict with remanei CDS sequences
CDS = convert_fasta('noamb_PX356_all_CDS.fasta')

# get the set of genes that have orthologs
orthologs = set()
infile = open('crem_cla_orthologs.txt', 'r')
# skip header
line = infile.readline()
while '#' in line:
    line = infile.readline()
# loop over file, add remanei transcript to set of orthologs
for line in infile:
    if line.rstrip() != '':
        line = line.rstrip().split()
        gene = line[0]
        orthologs.add(gene)
# close file
infile.close()

# make a dict of {transcript: gene} pairs
transcripts = transcript_to_gene('356_10172014.gff3')

# make a dict of {gene : [list of transcripts]} pairs
genes = gene_to_transcripts('356_10172014.gff3')

# make a set of genes already used
already_picked = set()

# make a set of final transcripts to use in analyses of diversity
final_transcripts = set()

# add in priority the transcripts with 1:1 orthologs
# loop over genes in orthologs
for gene in orthologs:
    # add gene to final set of genes
    if transcripts[gene] not in already_picked:
        final_transcripts.add(gene)
        # add the parent gene to set of genes already picked
        already_picked.add(transcripts[gene])
    
# loop over genes in dict of gene : [transcript list]
for gene in genes:
    # if gene not already picked
    if gene not in already_picked:
        # if only 1 transcript, add the transcript to final set
        if len(genes[gene]) == 1:
            # only 1 transcript, add to final set
            final_transcripts.add(genes[gene][0])
            # add parent gene to set of already picked genes
            already_picked.add(gene)
        else:
            # multiple transcripts, take the transcript with longer CDS
            CDS_length = [(len(CDS[isoform]), isoform) for isoform in genes[gene]]
            # sort list according to CDS length
            CDS_length.sort()
            # pick the transcript with the longest CDS
            longest = CDS_length[-1][1]
            # add transcript with longest CDS to final set
            final_transcripts.add(longest)
            # add parent gene to already_picked
            already_picked.add(gene)
            
# open file for writing
newfile = open('unique_transcripts.txt', 'w')
for gene in final_transcripts:
    newfile.write(gene + '\n')
newfile.close()

            
        
            
    
    
    