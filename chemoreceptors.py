# -*- coding: utf-8 -*-
"""
Created on Sat Jun 20 19:26:00 2015

@author: Richard
"""

from manipulate_sequences import *
from mK_test import *
from divergence import *


# use this function to get the gene IDs of all predicted GPCRs
def get_chemoreceptors(iprscan_file):
    '''
    (file) -> set
    Take the iprscan outputfile and return a set of genes that are predicted
    to be G-protein-coupled receptors
    '''
    
    # create a set of genes
    chemo = set()
    
    # open file for reading
    infile = open(iprscan_file, 'r')
    
    # loop over file
    for line in infile:
        if 'GPCR' in line:
            chemo.add(line.rstrip().split()[0])
        if 'chemoreceptor' in line:
            chemo.add(line.rstrip().split()[0])
    
    # close file
    infile.close()
    
    return chemo


# use this function to get a set of valid chemoreceptor genes
def clean_chemoreceptor_genes(iprscan_file, unique_transcripts):
    '''
    (file, file) -> set
    Take the iprscan output file, the file with valid transcripts and return
    a set of G-protein coupled receptors with a single transcript per gene
    '''
    
    # get the set of valid transcripts
    transcripts = get_valid_transcripts(unique_transcripts)
    # get the set of chemoreceptors
    chemo = get_chemoreceptors(iprscan_file)
    # create a set of valid chemoreceptors
    GPCR = set(gene for gene in chemo if gene in transcripts)
    
    return GPCR

        
# use this function to create a dict with aligned codons
def get_aligned_codons(CDS_alignment_file, directory):
    '''
    (file) -> dict
    Take a codon-based alignment file between remanei and latens and the
    directory in which file is located and return a dictionnary with codon
    index corresponding to the remanei codon and a list of cremanei and clatens
    aligned codons (ie. ignore gaps in the remanei sequence)
    '''
    
    # convert CDS file to dict
    orthologs = convert_fasta(directory + CDS_alignment_file)
    
    # create a dict {i: [crem_codon, cla_codon]}}
    # use the codon position in the remanei CDS sequence, ignoring gaps

    codons = {}

    # get the remanei and latens CDS sequences
    for gene in orthologs:
        if 'CRE' in gene:
            crm_gene = gene
            crm_cds = orthologs[crm_gene]
        else:
            cla_gene = gene
            cla_cds = orthologs[cla_gene]
        
    # initiate gap variable
    gaps = 0
    # initialise j
    j = 0
    # loop over the remanei CDS sequence
    for i in range(0, len(crm_cds), 3):
        if '-' in crm_cds[i:i+3]:
            gaps += crm_cds[i:i+3].count('-')
        else:
            j = i - gaps
            codons[j] = [crm_cds[i:i+3], cla_cds[i:i+3]]
    
    return codons



# use this fucntion to parse the phobious outfile with posterior probabilities
def parse_phobius_output(proba_file, directory):
    '''
    (file) -> dict
    Take the phobious output file with probabilities that amino acid are 
    inside the cell, outside or transmembrane and the directory where the file
    is located and return a dict with the amino acid index and list with
    amino acid and probabilities
    '''
    
    # create dict {codon-index: [aa, inside, outside, membrane, signal]}
    proba = {}
    
    # open file for reading
    infile = open(directory + proba_file, 'r')
    # skip 2 first lines
    infile.readline()
    infile.readline()
    # loop over file
    for line in infile:
        line = line.rstrip()
        if line != '':
            line = line.split()
            # get index 0-based of corresponding codons
            i = (int(line[0]) - 1) * 3
            # get amino acid
            aa = line[1]
            # get probabilities
            inside = float(line[2])
            outside = float(line[3])
            membrane = float(line[4])
            signal = float(line[5])
            # populate dict
            proba[i] = [aa, inside, outside, membrane, signal]
            
    # close file
    infile.close()
    
    return proba
    
    
# use this function to parse the codeml output file
def parse_codeml_output(codeml_outputfile):
    '''
    (file) -> dict
    Extract the dN and dS values from the codml output files, and return a 
    dictionnary with transcript name as key and a list with dN, dS, dN/dS as value
    '''
    
    # create dict {gene : [dN, dS, dN/dS]}
    divergence = {}
    
    # iniate crm_gene
    crm_gene = ''
    
    # open file for reading
    infile = open(codeml_outputfile, 'r')
    # loop over file
    for line in infile:
        # get remanei gene name
        if line.startswith('CRE'):
            # line starts with remanei gene name
            # get gene name only if not already found
            if crm_gene == '':
                line = line.rstrip().split()
                crm_gene = line[0]
        elif 'tree length for dN' in line:
            line = line.split()
            dN = float(line[-1])
        elif 'tree length for dS' in line:
            line = line.split()
            dS = float(line[-1])
    if dS != 0:
        omega = dN / dS
    elif dS == 0:
        omega = 'NA'
    
    divergence[crm_gene] = [dN, dS, omega]
    
    # close file
    infile.close()
    
    return divergence
            


# use this function to get the chemoreceptor genes in each chemoreceptor family
def chemo_families(iprscan_file):
    '''
    (file) -> dict
    Take the iprscan outputfile and return a dictionnary with chemoreceptor gene
    family as key and a set of chemoreceptor genes as value
    '''
    
    # create dictionnary
    chemo = {}

    # open file for reading
    infile = open(iprscan_file, 'r')
    
    # loop over file
    for line in infile:
        # check if gene is chemoreceptor and that method is Pfam (only Pfam returns families)
        if 'chemoreceptor' in line and 'Pfam' in line:
            line = line.rstrip().split('\t')            
            # Sre, Srg and Srab family needs to be parsed differently
            if 'Sre' in line[5]:
                family = line[5].split()[2]
            elif 'Srg' in line[5]:
                family = line[5].split()[0]
            elif '7TM GPCR' in line[5] and 'ab' in line[5]:
                family = 'Srab'
            else:
                family = line[5].split()[-1]
                
            # check that family is valid chemoreceptor family (some genes are not assigned to families)
            if family in {'Srbc', 'Sra', 'Srab', 'Sru', 'Sri', 'Srx', 'Srsx', 'Str', 'Srg',
                          'Srz', 'Srh', 'Srd', 'Srt', 'Srb', 'Srj', 'Srv', 'Sre', 'Srw'}:
                # check if family is key in dict
                if family in chemo:
                    # add gene to set
                    chemo[family].add(line[0])
                else:
                    # initiate key: gene set pair
                    chemo[family] = set()
                    chemo[family].add(line[0])
    
    # close file
    infile.close()
    
    return chemo
    
 
# use this function to find ambiguous chemoreceptors assigned to multiple families
def find_ambiguous_chemoreceptors(chemo):
    '''
    (dict) -> set
    Take the dictionnary of chemoreceptor family : set of genes pairs
    and return a set of ambiguous chemoreceptors that are assigned to multiple
    families
    '''

    # create a set of genes that are in more than 1 family
    overlap = set()

    # create a list of families
    families = [fam for fam in chemo]
    
    # compare each family, find overlapping genes and add them to set
    for i in range(len(families)-1):
        for j in range(i+1, len(families)):
            # check that genes are present in both families
            if len(chemo[families[i]].intersection(chemo[families[j]])) != 0:
                # add all overlapping genes to set
                for gene in chemo[families[i]].intersection(chemo[families[j]]):
                    overlap.add(gene)
    
    return overlap

# use this function to discard ambiguous chemoreceptor genes
def remove_ambiguous_chemoreceptors(chemo):
    '''
    (dict) -> dict
    Take the dictionnary of chamoreceptor family : gene set pairs and return
    a modified dictionnary in which ambiguous cjemoreceptors belonging to multiple
    families have been removed
    '''
    
    # find ambiguous chemoreceptors    
    ambiguous = find_ambiguous_chemoreceptors(chemo)
    
    # loop over chemoreceptor families
    for family in chemo:
        # loop over ambiguous gene
        for gene in ambiguous:
            # if gene in family, remove gene from family
            if gene in chemo[family]:
                chemo[family].remove(gene)
    # return modified dict
    return chemo
                
# use this function to get rif of non-valid transcripts in the chemo_fam dict
def clean_chemo_families(chemo, unique_transcripts):
    '''
    (dict, file) -> dict
    Take the dict of chemoreceptor families : gene set pairs and the file
    with unique transcripts and return a modified dictionnary that contains
    only the valid chemoreceptor genes (ie. a a single transcript mapping to a
    single parent gene)
    '''
    
    # get the valid transcripts
    transcripts = get_valid_transcripts(unique_transcripts)
    
    # loop over families in chemo
    for family in chemo:
        # create a list of genes to remove
        to_remove = []
        # loop over genes in given family
        for gene in chemo[family]:
            # check if gene in transcripts
            if gene not in transcripts:
                to_remove.append(gene)
        # check if genes need to be removed
        if len(to_remove) != '':
            for gene in to_remove:
                chemo[family].remove(gene)
    
    # remove families with no genes
    to_remove = []
    for family in chemo:
        if len(chemo[family]) == 0:
            to_remove.append(family)
    if len(to_remove) != 0:
        del chemo[family]
        
    return chemo
            
    
# use this function to reverse the dict of chemo family : gene pairs    
def chem_genes_to_family(chemo):
    '''
    (dict) -> dict
    Take the dictionnary of chemoreceptor family : gene set pairs and reverse it
    to return a dictionnary of gene : family pair
    '''
    
    # create new dict
    genes = {}
    
    # loop over families
    for family in chemo:
        # loop over genes and populate new dict
        for gene in chemo[family]:
            genes[gene] = family
            
    return genes
    
    
# use this function to compare nucleotide divergence between transmembrane domains and extra-mebrane domains    
def parse_chemo_divergence_table(chemo_divergence_file, rate, threshold):
    '''
    (file, str, num) -> (list, list)
    Take the divergence table of chemoreceptor domains, a threshold to eliminate
    high values and rate ('dN', 'dS' or 'omega') to specifiy the type of
    nucleotide divergence, and return a tuple with lists of nucleotide divergence
    for transmebrane domains and for extra-membrane domains
    '''
        
    # create lists to store nucleotide divergence
    TM, Extra_TM = [], []    
        
    # open file for reading
    infile = open(chemo_divergence_file, 'r')
    # skip header
    infile.readline()
    
    # loop over file
    for line in infile:
        if line.rstrip() != '':
            line = line.rstrip().split()
            gene = line[0]
            # check which nucleotide rate to use
            if rate == 'dN':
                # use non-synonymous nucleotide divergence
                TM_diverg = float(line[1])
                ExTM_diverg = float(line[4])
            elif rate == 'dS':
                # use synonymous nucleotide divergence
                TM_diverg = float(line[2])
                ExTM_diverg = float(line[5])
            elif rate == 'omega':
                # use ratio of dN/dS
                # check that omega is defined
                if line[3] != 'NA' and line[6] != 'NA':
                    # omega is defined
                    TM_diverg = float(line[3])
                    ExTM_diverg = float(line[6])
            # check threshold, do not consider genes with dN and dS > threshold
            if float(line[1]) <= threshold and float(line[4]) <= threshold and float(line[2]) <= threshold and float(line[5]) <= threshold: 
                # store divergence in lists
                if TM_diverg <= threshold and ExTM_diverg <= threshold:
                    TM.append(TM_diverg)
                    Extra_TM.append(ExTM_diverg)
                
    # close file
    infile.close()
    
    return TM, Extra_TM


# use this function to parse the file of divergence between remanei and latens
def parse_divergence_file(divergence_file):
    '''
    (file) -> dict
    Take the file of nucleotide divergence between C. remanei and C. latens and
    return a dict with remanei gene as key and list of dN, dS and dN/dS as value
    '''
    
    # create dict
    divergence = {}
    
    # open file for reading
    infile = open(divergence_file, 'r')
    # skip header
    infile.readline()
    # loop over file
    for line in infile:
        if line.rstrip() != '':
            line = line.rstrip().split()
            # get gene
            gene = line[0]
            # get dN, dS and dN/dS
            dN = float(line[4])
            dS = float(line[5])
            # check if omega is defined
            if line[6] == 'NA':
                # omega is not defined
                omega = line[6]
            else:
                # omega is defined
                omega = float(line[6])
            # populate dict
            divergence[gene] = [dN, dS, omega]
    
    # close file
    infile.close()
    
    return divergence
    
        
# use this function to compare nucleotide divergence among chemoreceptor families
def chemoreceptor_family_divergence(divergence_file, iprscan_file, rate):
    '''
    (file, file, str) -> dict
    Take the summary file of nucleotide divergence between remanei and latens,
    the iprscan file and the nucleotide rate (dN, dS or omega) and return a
    dictionnary with chemoreceptor family as key and a list with dN, dS or
    omega as value
    '''

    # get the dict of chemo families
    families = chemo_families(iprscan_file)
    # remove ambiguous genes
    families = remove_ambiguous_chemoreceptors(families)
    # reverse dictionnary {gene : family}
    chemo_genes = chem_genes_to_family(families)
    # parse divergence file
    divergence = parse_divergence_file(divergence_file)
    
    # create dictionnary 
    fam_divergence = {}
    
    # loop over genes in divergence
    for gene in divergence:
        # check if gene is chemoreceptor
        if gene in chemo_genes:
            # check which rate to use
            if rate == 'dN':
                nucleotide_diverg = divergence[gene][0]
            elif rate == 'dS':
                nucleotide_diverg = divergence[gene][1]
            elif rate == 'omega':
                nucleotide_diverg = divergence[gene][-1]
            # check if divergence is defined
            if nucleotide_diverg != 'NA':
                # populate dict
                fam = chemo_genes[gene]
                # check if family is key
                if fam in fam_divergence:
                    # family is key, add divergence to list
                    fam_divergence[fam].append(nucleotide_diverg)
                else:
                    # family is not key, initiate list
                    fam_divergence[fam] = [nucleotide_diverg]
                    
    return fam_divergence
    









  