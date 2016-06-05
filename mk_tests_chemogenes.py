# -*- coding: utf-8 -*-
"""
Created on Sat Jun  4 19:09:25 2016

@author: Richard
"""


# use this script to generate table with number of genes departing from neutrality with MK test
# and alpha, the proportion of amino acid fixed by adaptive evolution estimated by the method of 
# Smith and Eyre-Walker Science 2002


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

# create a set of non-chemo gene
NonGPCRs = set(gene for gene in transcripts if gene not in chemo)
print('NonGPCRs', len(NonGPCRs))
print('made set of valid non-chemo genes')

# count divergence and polymorphisms in Ontario strains
# eliminate singleton polymorphisms
# dict in form {gene : [PN, PS, DN, DS]}        
MK = count_polym_diverg('../Genome_Files/CDS_SNP_DIVERG.txt', 'KSR_PX', 'count', 0, 1)

# remove genes that are not valid transcripts
to_remove = [gene for gene in MK if gene not in transcripts]
for gene in to_remove:
    del MK[gene]
print('removed {0} genes from MK analysis'.format(len(to_remove)))

# perform MK test using Fisher exact test
mk = MK_test(MK, 'fisher')

# make a list of genes with significant MK test
significant = [gene for gene in mk if mk[gene][-1] < 0.05]
nonsignificant = [gene for gene in mk if mk[gene][-1] >= 0.05]

# determine genes with significant MK that are under positive or negative selection
negative, positive = NaturalSelection(MK, significant)

print('significant', len(significant))
print('positive', len(positive))
print('negative', len(negative))
assert len(mk) == len(positive) + len(negative) + len(nonsignificant)

# create lists of GPCR and non-GPCR genes
PositiveChemo = [gene for gene in positive if gene in GPCRs]
NegativeChemo = [gene for gene in negative if gene in GPCRs]
NonsignificantChemo = [gene for gene in nonsignificant if gene in GPCRs]

PositiveNC = [gene for gene in positive if gene in NonGPCRs]
NegativeNC = [gene for gene in negative if gene in NonGPCRs]
NonsignificantNC = [gene for gene in nonsignificant if gene in NonGPCRs]

print('chemo +', len(PositiveChemo))
print('chemo -', len(NegativeChemo))
print('chemo NS', len(NonsignificantChemo))
print('NC +', len(PositiveNC))
print('NC -', len(NegativeNC))
print('NC NS', len(NonsignificantNC))
assert len(mk) == len(PositiveChemo) + len(NegativeChemo) + len(NonsignificantChemo) + \
                  len(PositiveNC) + len(NegativeNC) + len(NonsignificantNC)

# apply a Bejamini-Hochberg correction for multiple testing {gene: corrected_Pval}
mk_p = {}
for gene in mk:
    mk_p[gene] = mk[gene][-1]
pvals = [(p, gene) for gene, p in mk_p.items()]
corrected = Benjamini_Hochberg_correction(pvals)

# count the number of genes with significant MK after BJ correction with 10% FDR 
FDR = 0.1
signifcorr = [gene for gene in corrected if corrected[gene] < FDR]
nonsignifcorr = [gene for gene in corrected if corrected[gene] >= FDR]
print('FDR', FDR, len(signifcorr))
print('FDR', FDR, len(nonsignifcorr))

# determine genes with significant MK that are under positive or negative selection
negativecorr, positivecorr = NaturalSelection(MK, signifcorr)

# create lists of GPCR and non-GPCR genes
PositiveChemoCorr = [gene for gene in positivecorr if gene in GPCRs]
NegativeChemoCorr = [gene for gene in negativecorr if gene in GPCRs]
NonsignificantChemoCorr = [gene for gene in nonsignifcorr if gene in GPCRs]

PositiveNCCorr = [gene for gene in positivecorr if gene in NonGPCRs]
NegativeNCCorr = [gene for gene in negativecorr if gene in NonGPCRs]
NonsignificantNCCorr = [gene for gene in nonsignifcorr if gene in NonGPCRs]

print('chemo +', len(PositiveChemoCorr))
print('chemo -', len(NegativeChemoCorr))
print('chemo NS', len(NonsignificantChemoCorr))
print('NC +', len(PositiveNCCorr))
print('NC -', len(NegativeNCCorr))
print('NC NS', len(NonsignificantNCCorr))
assert len(mk) == len(PositiveChemoCorr) + len(NegativeChemoCorr) + len(NonsignificantChemoCorr) + \
                  len(PositiveNCCorr) + len(NegativeNCCorr) + len(NonsignificantNCCorr)


# compute alpha according to the Smith-EyreWalker 2002 method

# generate dicts with polymorphism amd divergence counts for chemo and nonchemo genes
ChemoPolymDivCounts, NCPolymDivCounts = {}, {}
for gene in MK:
    if gene in GPCRs:
        ChemoPolymDivCounts[gene] = list(MK[gene])
    elif gene in NonGPCRs:
        NCPolymDivCounts[gene] = list(MK[gene])

# compute average alpha for chemo and non-chemo genes

for i in range(10):
    ChemoAlpha = ComputeAlphaSEW2002(ChemoPolymDivCounts, i)
    NCAlpha = ComputeAlphaSEW2002(NCPolymDivCounts, i)
    print('alpha GPCR', i, ChemoAlpha)
    print('alpha NC', i, NCAlpha)


# Bootstrap genes to compute 






