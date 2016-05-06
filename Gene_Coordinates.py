# -*- coding: utf-8 -*-
"""
Created on Fri May  6 11:43:14 2016

@author: RJovelin
"""


# use this function to map transcript name to gene names
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


# use this function to map all transcript names to gene names
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
    

# use this function to map the remanei chromosomes to elegans chromosomes
def remanei_elegans_chromosomes(chromosome_file):
    '''
    (file) -> dict
    Take the chromosome_file and return a dictionnary of remanei sacffolds/
    linkage groups as key and the corresponding elegans chromosome as value    
    '''
    
    # create dict {rem_LG: elegans_LG}
    LG = {}
    
    # open file for reading
    infile = open(chromosome_file, 'r')
    # skip header
    infile.readline()
    # loop over file, populate dict
    for line in infile:
        line = line.rstrip()
        if line != '':
            line = line.split(',')
            rem_chromo = line[0]
            elegans_chromo = line[1]
            LG[rem_chromo] = elegans_chromo
    # close file after reading
    infile.close()
    
    return LG

    
# use this function to partition a collection of genes into lists of X-linked and autosomal genes 
def X_autosomal_genes(genes, caeno_gff, chromosome_file):
    '''
    (collection, file, file) -> (list, list)
    Take a collection of genes (set or list), the GFF file, and a file with
    remanei and elegans corresponding chromosomes and return a tuple with 
    a list of X-linked genes and a list of autosomal genes    
    '''
    
    # create a dict of remanei : elegans chromosomes
    chromosomes = remanei_elegans_chromosomes(chromosome_file)    
    
    # get the genes coordinates {gene1 : [chromo, start, end, sense]}
    genes_coordinates = get_genes_coordinates(caeno_gff)
    
    # create lists of genes
    X_genes, autosomal_genes = [], []
    
    # loop over genes
    for gene in genes:
        # get chromo
        chromo = genes_coordinates[gene][0]
        # check that chromo has a elegans equivalent
        if chromo in chromosomes:
            # get the elegans chromosome
            if chromosomes[chromo] == 'X':
                # chromo is X chromosome
                X_genes.append(gene)
            else:
                # gene is autosomal
                autosomal_genes.append(gene)
                
    return X_genes, autosomal_genes