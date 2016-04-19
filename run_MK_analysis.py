# -*- coding: utf-8 -*-
"""
Created on Wed Jul  1 16:51:27 2015

@author: Richard
"""

from accessories import *
from MK_test import *
from crm_expression import *
from multiple_testing import *
from Gene_ontologies import *
import numpy as np

# create a set of valid transcripts
transcripts = get_valid_transcripts('unique_transcripts.txt')

# create a set of genes with orthologs
orthos = set()

# open divergence file
infile = open('CRM_CLA_prot_diverg_filtered.txt', 'r') 
# skiper header
infile.readline()
# loop over file
for line in infile:
    if line.rstrip() != '':
        line = line.rstrip().split()
        orthos.add(line[0])
# close file
infile.close()


# get expression level for all genes
expression = expression_developmental_stages('WBPaper00041190.cre.mr.csv', '../Premature_stops/c_remanei.PRJNA53967.WS248.geneIDs.txt')


# make a list of expression level
level = [expression[gene] for gene in expression]
# compute quartiles
Q1 = np.percentile(level, 25)
Q2 = np.median(level)
Q3 = np.percentile(level, 75)

# create dicts genes: expression pairs for different expression groups
# expression groups are based on the distribution of expression level, 
# using quartiles as limits
low, midlow, midhigh, high = {}, {}, {}, {}
for gene in expression:
    if expression[gene] < Q1:
        low[gene] = expression[gene]
    elif expression[gene] >= Q1 and expression[gene] < Q2:
        midlow[gene] = expression[gene]
    elif expression[gene] >= Q2 and expression[gene] < Q3:
        midhigh[gene] = expression[gene]
    elif expression[gene] >= Q3:
        high[gene] = expression[gene]
        
# open new file to record results of MK analysis
newfile = open('summary_MK_analysis.txt', 'w')

# record expression
newfile.write('Expression groups:\n')
newfile.write('-' * 19 + '\n')
newfile.write('Expression_group\tExpression_range\tN_genes\n')
newfile.write('\t'.join(['Low', '0-{0}'.format(Q1), str(len(low))]) + '\n')
newfile.write('\t'.join(['Midlow', '{0}-{1}'.format(Q1, Q2), str(len(midlow))]) + '\n')
newfile.write('\t'.join(['Midhigh', '{0}-{1}'.format(Q2, Q3), str(len(midhigh))]) + '\n')
newfile.write('\t'.join(['High', '>{0}'.format(Q3)]) + '\n')
newfile.write('\n')
        
# count divergence and polymorphisms in Ontario strains
# eliminate singleton polymorphisms
# dict in form {gene : [PN, PS, DN, DS]}        
MK = count_polym_diverg('CDS_SNP_DIVERG.txt', 'KSR_PX', 'count', 0, 1)

# remove genes that are not in orthos
to_remove = [gene for gene in MK if gene not in orthos]
for gene in to_remove:
    del MK[gene]


# perform MK test using Fisher exact test
mk = MK_test(MK, 'fisher')

# make a list of genes with significant MK test
significant_MK = [gene for gene in mk if mk[gene][-1] < 0.05]


# count genes with significant MK in expression groups
low_significant = [gene for gene in significant_MK if gene in low]
midlow_significant = [gene for gene in significant_MK if gene in midlow]
midhigh_significant  = [gene for gene in significant_MK if gene in midhigh]
high_significant = [gene for gene in significant_MK if gene in high]

# determine genes with significant MK that are under positive or negative selection
negative, positive = [], []
for gene in MK:
    # check if gene has significant departure from neutrality
    if gene in significant_MK:
        # check if PS is defined
        if MK[gene][1] != 0:
            # compute ratio PN / PS
            polymorphism = MK[gene][0] / MK[gene][1]
        elif MK[gene][1] == 0:
            # add small number to PS
            polymorphism = MK[gene][0] / (MK[gene][1] + 0.1)
        # check if DS is defined
        if MK[gene][3] != 0:
            # compute DN / DS
            divergence = MK[gene][2] / MK[gene][3]
        elif MK[gene][3] == 0:
            # add small number to DS
            divergence = MK[gene][2] / (MK[gene][3] + 0.1)
        # check if divergence is lower r higher than polymorphism
        if divergence > polymorphism:
            # positive selection
            positive.append(gene)
        elif divergence < polymorphism:
            negative.append(gene)
            
# determine the number of genes under positive and negative selection for expression groups
low_positive = [gene for gene in positive if gene in low]
midlow_positive = [gene for gene in positive if gene in midlow]
midhigh_positive = [gene for gene in positive if gene in midhigh]
high_positive = [gene for gene in positive if gene in high]            
   
low_negative = [gene for gene in negative if gene in low]
midlow_negative = [gene for gene in negative if gene in midlow]
midhigh_negative = [gene for gene in negative if gene in midhigh]
high_negative = [gene for gene in high if gene in high]


# apply a Bejamini-Hochberg correction for multiple testing
mk_p = {}
for gene in mk:
    mk_p[gene] = mk[gene][-1]
pvals = [(p, gene) for gene, p in mk_p.items()]
corrected = Benjamini_Hochberg_correction(pvals)

# determine the number of genes in expression groups for significant genes after correction
# use a 10% FDR
FDR = 0.1
low_positive_corrected = [gene for gene in corrected if gene in low_positive and corrected[gene] < FDR]
midlow_positive_corrected = [gene for gene in corrected if gene in midlow_positive and corrected[gene] < FDR]
midhigh_positive_corrected = [gene for gene in corrected if gene in midhigh_positive and corrected[gene] < FDR]
high_positive_corrected = [gene for gene in corrected if gene in high_positive and corrected[gene] < FDR]

low_negative_corrected = [gene for gene in corrected if gene in low_negative and corrected[gene] < FDR]
midlow_negative_corrected = [gene for gene in corrected if gene in midlow_negative and corrected[gene] < FDR]
midhigh_neggative_corrected = [gene for gene in corrected if gene in midhigh_negative and corrected[gene] < FDR]
high_negative_corrected = [gene for gene in corrected if gene in high_negative and corrected[gene] < FDR]

# get a list of genes with significant MK test after correction
significant_corrected = [gene for gene in corrected if corrected[gene] < FDR]
# get a list of genes with positive selection after correction
positive_corrected = [gene for gene in corrected if gene in positive and corrected[gene] < FDR]
# get a list of genes with negative selection after correction
negative_corrected = [gene for gene in corrected if gene in negative and corrected[gene] < FDR]







# create a dict of GO terms : genes
GO = GO_to_genes('../Chemoreceptors/PX356_protein_seq.tsv', 'unique_transcripts.txt')
# test over=representation among genes with positive selection
GO_positive_pval = test_overepresentation(positive_corrected, GO)
# test over-representation among genes with purifying selection
GO_negative_pval = test_overepresentation(negative_corrected, GO)


positive_over = [gene for gene in GO_positive_pval if GO_positive_pval[gene] < FDR]
negative_over = [gene for gene in GO_negative_pval if GO_negative_pval[gene] < FDR]

print(len(positive_over))
print(len(negative_over))



#>>> significant = {}
#>>> for gene in corrected:
#...     
#KeyboardInterrupt
#>>> corrected[genes[0]]
#0.32820927617136486
#>>> for gene in corrected:
#...     if corrected[gene] < 0.05:
#...         significant.add(gene)
#... 
#Traceback (most recent call last):
#  File "<stdin>", line 3, in <module>
#AttributeError: 'dict' object has no attribute 'add'
#>>> significant = set()
#>>> for gene in corrected:
#...     if corrected[gene] < 0.05:
#...         significant.add(gene)
#... 
#>>> len(significant)
#279
#>>> len(MK)
#15379
#>>> len(mk)
#15379
#>>> from accessories import *
#>>> transcripts = get_valid_transcripts('unique_transcripts.txt')
#>>> positive = [gene for gene in significant if gene in transcripts]
#>>> len(positive)
#279
#>>> 



#print('significant - low - midlow - midhigh - high:')
#print(len(low_significant), len(midlow_significant), len(midhigh_significant), len(high_significant))
#
#print('positive  - low - midlow - midhigh - high:')
#print(len(low_positive), len(midlow_positive), len(midhigh_positive), len(high_positive))
#   
#
#print('significant_corrected', len(significant_corrected))
#print('low - midlow - midhigh - positive high after correction')
#print(len(low_positive_corrected), len(midlow_positive_corrected), len(midhigh_positive_corrected), len(high_positive_corrected))
#print('low - midlow - midhigh - negative high after correction')
#print(len(low_negative_corrected), len(midlow_negative_corrected), len(midhigh_neggative_corrected), len(high_negative_corrected))
#










# applying correction for multiple testing
# make a list of tuples (p-value, gene)






# what are the genes with significant MK test, negative selection, positive selection



# how many genes in expression groups for significant MK genes, 
# for genes negative selection, for genes with positive selection


# same questions for significant genes after correction




# close file after writing
newfile.close()



