# -*- coding: utf-8 -*-
"""
Created on Mon Jun 15 13:37:54 2015

@author: Richard
"""

import numpy as np




# use this function to het the expression level in remanei males and females
def expression_male_female(expression_male_female_file):
    '''
    (file) -> dict
    Returns a dictionnary for each gene in the expression file with
    a tuple including the average expression in fermale, in male, and p-value
    of the mean difference between female and male
    '''

    # create dictionnary {gene : (female, male, p-val)}
    expression = {}

    # open file for reading
    infile = open(expression_male_female_file, 'r')
    # skip header
    infile.readline()
    # read file and store relevant info in dict
    for line in infile:
        line = line.rstrip()
        if line != '':
            line = line.split()
            # modify gene so that it matches in the GFF file
            gene = 'CRE_PX356_' + line[0]
            # get expression level in males and females 
            female = float(line[4])
            male = float(line[8])
            # get the p-value comparing expression difference between males and females
            p_val = float(line[9])
            expression[gene] = (female, male, p_val)

    infile.close()
    return expression


# use this function to assign gene names to remanei gene IDs
def crem_wormbase_gene_ID(gene_ID_file):
    '''
    (file) -> dict
    Returns a dictionnary with the wormbase gene ID as key and
    the remanei gene name as value
    '''

    # create dictionnary {WBGene : CRE0000}
    identifiers = {}

    # open file for reading
    ID = open(gene_ID_file, 'r')

    # go through the file
    for line in ID:
        line = line.rstrip()
        if line != '':
            line = line.split(',')
            # use only live genes
            if line[-1] == 'Live':
                if line[1].startswith('WBG') and line[-2].startswith('CRE'):
                    identifiers[line[1]] = line[-2]

    ID.close()
    return identifiers


# use this function to get average expression level across development
def expression_developmental_stages(expression_file, gene_ID_file):
    '''
    (file, file) -> dict
    Returns a dictionnary with cremanei gene name and average expression
    level measured across embryonic development
    '''

    # get the cremanei gene names for each remanei Wormbase ID
    identifiers = crem_wormbase_gene_ID(gene_ID_file)

    # make a dict with average expression for each WBGene ID
    expression_WBG = {}
    
    # open expression file for reading
    infile = open(expression_file, 'r')
    infile.readline()

    # go through the file
    for line in infile:
        line = line.rstrip()
        if line != '':
            # convert str to list
            line = line.split()
            # extract WBGene ID from list
            gene = line.pop(0)
            # convert strings in list to floats
            line = list(map(lambda x: float(x), line))
            # compute mean expression
            mean = np.mean(line)
            expression_WBG[gene] = mean

    # create dict {CRE_gene_ID: mean eexpression}
    expression = {}
    for gene in identifiers:
        # check if gene has expression
        if gene in expression_WBG:
            # get the remanei gene name
            gene_name = identifiers[gene]
            # modify remanei gene name to match with gene names in the GFF file
            gene_name = 'CRE_PX356_' + gene_name
            expression[gene_name] = expression_WBG[gene]
    # close file after reading
    infile.close()
    return expression


            
            
