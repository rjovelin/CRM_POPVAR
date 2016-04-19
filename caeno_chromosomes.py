# -*- coding: utf-8 -*-
"""
Created on Tue Jun 16 12:55:59 2015

@author: Richard
"""

from miRNA_target import *

# use this function to map the remanei chromosomes to elegans chromosomes
def remanei_elegans_chromosomes(chromosome_file):
    '''
    (file) -> dict
    Take the chromosome_file and return a dicttionnary of remanei sacffolds/
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