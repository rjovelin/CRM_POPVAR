# -*- coding: utf-8 -*-
"""
Created on Fri May  6 11:43:14 2016

@author: RJovelin
"""



def transcript_to_gene(caeno_gff):
    '''
    (file) -> dict
    Returns a dictionnary with transcript : gene pairs from the gff annotation file
    '''
    #create a dictionnary of transcript : gene pairs
    transcripts_genes = {}

    # open file for reading
    gff = open(caeno_gff, 'r')
    for line in gff:
        line = line.rstrip()
        if line != '':
            line = line.split()
            if len(line) > 8:
                if line[2] == 'mRNA':
                    transcript = line[8][line[8].index('ID=')+3 : line[8].index(';')]
                    gene = line[8][line[8].index('Parent')+7 : line[8].index(';', line[8].index('Parent'))]
                    transcripts_genes[transcript] = gene
    gff.close()
    return transcripts_genes


def gene_to_transcripts(caeno_gff):
    '''
    (file) -> dict
    Returns a dictionnary with gene as key and a list of transcripts as value
    '''

    # get the dictionnary of transcripts : gene names pairs
    transcripts_genes = transcript_to_gene(caeno_gff)

    # create a reverse dictionnary of gene : [transcripts] pairs
    genes = {}
    for transcript in transcripts_genes:
        gene_name = transcripts_genes[transcript]
        if gene_name in genes:
            genes[gene_name].append(transcript)
        else:
            genes[gene_name] = [transcript]

    return genes


# use this function to generate a set of genes with indels to exclude from analysis
def get_genes_with_indels(indel_transcripts):
    '''
    (file) -> set
    Take a file with transcrtipts that have an indel in their CDS and
    a return a set of such transcripts
    '''
    
    # create set of genes
    indels = set()
    
    # open file for reading
    infile = open(indel_transcripts, 'r')
    # loop over file
    for line in infile:
        line = line.rstrip()
        if line != '':
            indels.add(line)
            
    # close file
    infile.close()
    
    return indels
    

# use this function to generate a set of transcripts to use for analyses of diversity
def get_valid_transcripts(unique_transcripts):
    '''
    (file) -> set
    Take the file with transcript names and return a set of transcripts
    that each have a unique parent gene to use in analyses of diversity
    '''
    
    # open file for reading
    infile = open(unique_transcripts, 'r')
    # make a set of transcript names
    transcripts = set()
    for line in infile:
        line = line.rstrip()
        transcripts.add(line)
    
    # close file
    infile.close()
    
    return transcripts
    


