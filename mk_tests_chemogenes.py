# -*- coding: utf-8 -*-
"""
Created on Sat Jun  4 19:09:25 2016

@author: Richard
"""


# use this script to 


# use Agg backend on server without X server
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib import rc
rc('mathtext', default='regular')
# import modules
import numpy as np
from scipy import stats
import math
import os
import sys
# import custom modules
from chemoreceptors import *
from manipulate_sequences import *
from mk_test import *
from crm_expression import *
from multiple_testing import *
from gene_ontologies import *


# create a set of valid transcripts
transcripts = get_valid_transcripts('../Genome_Files/unique_transcripts.txt')
print('got list of valid transcripts')

# get the chemoreceptor genes in each chemoreceptor family
Families = chemo_families('../Genome_Files/PX356_protein_seq.tsv')
print('assigned GPCRs to gene families')

# remove non valid genes
for family in Families:
    to_remove = []
    for gene in Families[family]:
        if gene not in transcripts:
            to_remove.append(gene)
    for gene in to_remove:
        Families[family].remove(gene)
print('removed non genes in chemo families')

# get the set of chemoreceptors from the iprscan outputfile
chemo = get_chemoreceptors('../Genome_Files//PX356_protein_seq.tsv') 
print('parsed chemo genes')

# create a set of valid chemoreceptors
GPCRs = set(gene for gene in chemo if gene in transcripts)
print('GPCRS', len(GPCRs))
print('made set of valid chemo genes')

NonGPCRs = set(gene for gene in transcripts if gene not in chemo)
print('NonGPCRs', len(NonGPCRs))
print('made set of valid non-chemo genes')


# count divergence and polymorphisms in Ontario strains
# eliminate singleton polymorphisms
# dict in form {gene : [PN, PS, DN, DS]}        
MK = count_polym_diverg('../Genome_Files/CDS_SNP_DIVERG.txt', 'KSR_PX', 'count', 0, 1)

# perform MK test using Fisher exact test
mk = MK_test(MK, 'fisher')



# make a list of genes with significant MK test
significant_MK = [gene for gene in mk if mk[gene][-1] < 0.05]


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
            
print('significant', len(significant_MK))
print('positive', len(positive))
print('negative', len(negative))





# apply a Bejamini-Hochberg correction for multiple testing
mk_p = {}
for gene in mk:
    mk_p[gene] = mk[gene][-1]
pvals = [(p, gene) for gene, p in mk_p.items()]
corrected = Benjamini_Hochberg_correction(pvals)

## determine the number of genes in expression groups for significant genes after correction
## use a 10% FDR
#FDR = 0.1
#low_positive_corrected = [gene for gene in corrected if gene in low_positive and corrected[gene] < FDR]
#midlow_positive_corrected = [gene for gene in corrected if gene in midlow_positive and corrected[gene] < FDR]
#midhigh_positive_corrected = [gene for gene in corrected if gene in midhigh_positive and corrected[gene] < FDR]
#high_positive_corrected = [gene for gene in corrected if gene in high_positive and corrected[gene] < FDR]
#
#low_negative_corrected = [gene for gene in corrected if gene in low_negative and corrected[gene] < FDR]
#midlow_negative_corrected = [gene for gene in corrected if gene in midlow_negative and corrected[gene] < FDR]
#midhigh_neggative_corrected = [gene for gene in corrected if gene in midhigh_negative and corrected[gene] < FDR]
#high_negative_corrected = [gene for gene in corrected if gene in high_negative and corrected[gene] < FDR]
#
## get a list of genes with significant MK test after correction
#significant_corrected = [gene for gene in corrected if corrected[gene] < FDR]
## get a list of genes with positive selection after correction
#positive_corrected = [gene for gene in corrected if gene in positive and corrected[gene] < FDR]
## get a list of genes with negative selection after correction
#negative_corrected = [gene for gene in corrected if gene in negative and corrected[gene] < FDR]
#
#
#
#
#
#
#
## create a dict of GO terms : genes
#GO = GO_to_genes('../Chemoreceptors/PX356_protein_seq.tsv', 'unique_transcripts.txt')
## test over=representation among genes with positive selection
#GO_positive_pval = test_overepresentation(positive_corrected, GO)
## test over-representation among genes with purifying selection
#GO_negative_pval = test_overepresentation(negative_corrected, GO)
#
#
#positive_over = [gene for gene in GO_positive_pval if GO_positive_pval[gene] < FDR]
#negative_over = [gene for gene in GO_negative_pval if GO_negative_pval[gene] < FDR]
#
#print(len(positive_over))
#print(len(negative_over))



