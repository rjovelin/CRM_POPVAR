# -*- coding: utf-8 -*-
"""
Created on Sun Jun  5 15:24:32 2016

@author: Richard
"""


# use this script to plot the Pn/Ps ratio and Dn/Ds ratio for chemo and non chemo genes
# and the distribution of these 2 ratios per gene


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
import random
# import custom modules
from chemoreceptors import *
from manipulate_sequences import *
from mk_test import *


# create a set of valid transcripts
transcripts = get_valid_transcripts('../Genome_Files/unique_transcripts.txt')
print('got list of valid transcripts')

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

# generate dicts with polymorphism amd divergence counts for chemo and nonchemo genes
ChemoPolymDivCounts, NCPolymDivCounts = {}, {}
for gene in MK:
    if gene in GPCRs:
        ChemoPolymDivCounts[gene] = list(MK[gene])
    elif gene in NonGPCRs:
        NCPolymDivCounts[gene] = list(MK[gene])


# create lists of Pn/Ps for chemo and non chemo genes
chemoPolymRatio = sum([ChemoPolymDivCounts[gene][0] for gene in ChemoPolymDivCounts]) / sum([ChemoPolymDivCounts[gene][1] for gene in ChemoPolymDivCounts])
NCPolymRatio = sum([NCPolymDivCounts[gene][0] for gene in NCPolymDivCounts]) / sum([NCPolymDivCounts[gene][1] for gene in NCPolymDivCounts])
# create lists of Dn/Ds for chemo and non chemo genes
chemoDivRatio = sum([ChemoPolymDivCounts[gene][2] for gene in ChemoPolymDivCounts]) / sum([ChemoPolymDivCounts[gene][3] for gene in ChemoPolymDivCounts])
NCDivRatio = sum([NCPolymDivCounts[gene][2] for gene in NCPolymDivCounts]) / sum([NCPolymDivCounts[gene][3] for gene in NCPolymDivCounts])

# create lists of ratio for individial genes if Ps and Ds are > 0





