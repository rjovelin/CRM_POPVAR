# -*- coding: utf-8 -*-
"""
Created on Wed Jun 24 16:26:02 2015

@author: Richard
"""


from scipy import stats
from manipulate_sequences import *




# use this file to map GO identifiers to description
def GO_to_description(GO_file):
    '''
    (file) -> dict
    Take the file of Gene Ontologies and return a dictionnary with GO identifier
    as key and a description of the identifier
    '''
    
    # create dict of {GO_id : description}
    GO = {}
    
    # open file for reading
    infile = open(GO_file, 'r')
    # loop over file
    for line in infile:
        # check if line contains GO identifier
        if line.startswith('id: GO:'):
            # parse line and populate dict with GO identifier
            line = line.rstrip().split()
            GO_id = line[-1]
            GO[GO_id] = ''
        elif line.startswith('name:'):
            # parse line to get GO_id description, add value to dict
            name = line.rstrip()[line.index(': ')+2:]
            GO[GO_id] = name
            
    # close file after reading
    infile.close()
    
    return GO
            
        
# use this function to count the number of genes 
def GO_from_iprscan(iprscan_file):
    '''
    (file) -> dict
    Take the iprscan output file and return a dictionnary with Gene Ontology
    identifier as key and a set of associated genes    
    '''
    
    # create dict
    GO = {}
    
    # open file for reading
    infile = open(iprscan_file, 'r')
    # loop over file
    for line in infile:
        # check if gene has GO
        if 'GO:' in line:
            line = line.rstrip().split('\t')
            # check if gene has multiple GO
            if line[-1].count('GO:') > 1:
                # need to get all GO ids
                ontology = line[-1].split('|')
            else:
                # only 1 GO id
                ontology = [line[-1]]
                
            # loop over GO ids and populate dict
            for GO_id in ontology:
                # check if GO_id is key 
                if GO_id in GO:
                    # GO_id is key, add gene to set
                    GO[GO_id].add(line[0])
                else:
                    # initialize set
                    GO[GO_id] = set()
                    # add gene
                    GO[GO_id].add(line[0])
    
    # close file after reading
    infile.close()
    
    return GO
    
# use this function to remove non valid transcripts from dict of GO terms : associated genes
def GO_to_genes(iprscan_file, unique_transcripts):
    '''
    (file, file) -> dict
    Take the output file from Interproscan and the file with set of valid
    transcripts and return a dictionnary with Gene Ontology identifier as 
    key and set of genes as value, that include only the set of valid genes
    (ie, all transcripts in the set map to a unique parent gene)
    '''
    
    # get the set of valid genes
    transcripts = get_valid_transcripts(unique_transcripts)
    
    # get the dict of GO terms : gene set pairs
    GO = GO_from_iprscan(iprscan_file)
    
    # loop over GO terms
    for GO_term in GO:
        # create a list of genes to delete
        to_remove = []
        # loop over genes ssociated with GO_term
        for gene in GO[GO_term]:
            # if gene not valid, add to list to remove
            if gene not in transcripts:
                to_remove.append(gene)
        # check if there are genes to remove
        if len(to_remove) != '':
            # loop over genes to remove and delete them from set
            for gene in to_remove:
                GO[GO_term].remove(gene)
    
    # delete GO_terms that have no associated genes
    to_remove = []
    for GO_term in GO:
        if len(GO[GO_term]) == 0:
            to_remove.append(GO_term)
    if len(to_remove) != 0:
        for GO_term in to_remove:
            del GO[GO_term]
       
    return GO


# use this function to test the over-representation of Gene Ontology terms 
def test_overepresentation(Gene_list, GO_genes_dict):
    '''
    (list, dict) -> dict
    Take a list of gene of interest and the dictionnary of GO terms : associated genes
    (with only valid transcripts), and return a dictionnary of GO_term and probability
    that the term is over-represented in the gene list
    Precondition: the gene list contains only valid transcripts
    (ie. transcripts mapping to unique parent genes)
    '''
    
    # create a dict {GO_term: P_val}
    GO_pval = {}    
    
    # compute the total number genes with GO terms M
    total_genes = set()
    for GO_term in GO_genes_dict:
        for gene in GO_genes_dict[GO_term]:
            total_genes.add(gene)
    M = len(total_genes)
    
    # create a set of GO terms in the sample gene list
    sample_GO = set()
    
    # compute the total number of genes with GO terms in the sample gene list N
    sample_genes = set()
        
    # loop over GO terms in dict
    for GO_term in GO_genes_dict:
        # loop of genes
        for gene in GO_genes_dict[GO_term]:
            if gene in Gene_list:
                # add GO_term associated to genes in sample list
                sample_GO.add(GO_term)
                # add genes that have GO terms in sample list
                sample_genes.add(gene)
    N = len(sample_genes)
            
     
    # loop over GO_term in set of GO terms in gene list
    for GO_term in sample_GO:
        # compute the total number of genes associated with GO_term n
        n = len(GO_genes_dict[GO_term])
        # compute the number of genes associated with GO_term in sample k
        k = len(set(gene for gene in GO_genes_dict[GO_term] if gene in Gene_list))
        # compute the P-value of over-representation of that given GO term 
        p_val = stats.hypergeom.cdf(k, M, n, N)
        # populate dict
        GO_pval[GO_term] = p_val
        
    return GO_pval
    
    
    
            
    