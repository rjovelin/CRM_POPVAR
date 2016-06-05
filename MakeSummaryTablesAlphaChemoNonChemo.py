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
import random
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

# record alpha removing genes with PS = 0
ChemoAlpha = ComputeAlphaSEW2002(ChemoPolymDivCounts, 1)
NCAlpha = ComputeAlphaSEW2002(NCPolymDivCounts, 1)
# Bootstrap genes to generate distrobution of alpha
ChemoBootstrap = BootstrapAlphaSEW2002(ChemoPolymDivCounts, 1, 1000, 500)
print('boostraped alpha for chemo genes')
NCBootstrap = BootstrapAlphaSEW2002(NCPolymDivCounts, 1, 1000, 5000)
print('bootstraped alpha for non-chemo genes')

# compte 95% confidence interval
SEMChemoBootstrap = np.std(ChemoBootstrap) / math.sqrt(len(ChemoBootstrap))
ChemoLCI = np.mean(ChemoBootstrap) - (1.96 * SEMChemoBootstrap)
ChemoHCI = np.mean(ChemoBootstrap) + (1.96 * SEMChemoBootstrap)
SEMNCBootstrap = np.std(NCBootstrap) / math.sqrt(len(NCBootstrap))
NCLCI = np.mean(NCBootstrap) - (1.96 * SEMNCBootstrap)
NCHCI = np.mean(NCBootstrap) + (1.96 * SEMNCBootstrap)

print('chemo', ComputeAlphaSEW2002(ChemoPolymDivCounts, 1), '({0}, {1})'.format(ChemoLCI, ChemoHCI))
print('non-chemo', ComputeAlphaSEW2002(NCPolymDivCounts, 1), '({0}, {1})'.format(NCLCI, NCHCI))


# make summary table
newfile = open('TableSummaryChemoNonChemoAlpha.txt', 'w')
header = ['Genes', 'Adaptive', 'Negative', 'Neutral', 'AdaptiveCorrected', 'NegativeCorrected', 'NeutralCorrected', 'alpha (95% CI)']
newfile.write((len('\t'.join(header)) + 1) * '-' + '\n')
newfile.write('\t'.join(header) + '\n')
newfile.write((len('\t'.join(header)) + 1) * '-' + '\n') 
line1 = ['GPCRs', len(PositiveChemo), len(NegativeChemo), len(NonsignificantChemo),
         len(PositiveChemoCorr), len(NegativeChemoCorr), len(NonsignificantChemoCorr),
         '{0} ({1}, {2})'.format(round(ChemoAlpha, 3), round(ChemoLCI, 3), round(ChemoHCI, 3))]
line1 = list(map(lambda x: str(x), line1))
newfile.write('\t'.join(line1) + '\n')

line2 = ['Non-chemo', len(PositiveNC), len(NegativeNC), len(NonsignificantNC),
         len(PositiveNCCorr), len(NegativeNCCorr), len(NonsignificantNCCorr),
         '{0} ({1}, {2})'.format(round(NCAlpha, 3), round(NCLCI, 3), round(NCHCI, 3))]
line2 = list(map(lambda x: str(x), line2))
newfile.write('\t'.join(line2) + '\n')         
newfile.write((len('\t'.join(header)) + 1) * '-' + '\n')
newfile.write('\n' * 3)

# perform chi square of independance for genes under selection for chemo and non-chemo genes
header = ['Comparison', 'Chi2', 'P']
comps = ['+ vs NS', '- vs NS', '+ vs -', '+ vs NS corr', '- vs NS corr', '+ vs - corr']
data = [[[len(PositiveChemo), len(NonsignificantChemo)], [len(PositiveNC), len(NonsignificantNC)]],
        [[len(NegativeChemo), len(NonsignificantChemo)], [len(NegativeNC), len(NonsignificantNC)]],
        [[len(PositiveChemo), len(NegativeChemo)], [len(PositiveNC), len(NegativeNC)]],
        [[len(PositiveChemoCorr), len(NonsignificantChemoCorr)], [len(PositiveNCCorr), len(NonsignificantNCCorr)]],
        [[len(NegativeChemoCorr), len(NonsignificantChemoCorr)], [len(NegativeNCCorr), len(NonsignificantNCCorr)]],
        [[len(PositiveChemoCorr), len(NegativeChemoCorr)], [len(PositiveNCCorr), len(NegativeNCCorr)]]
        ]

newfile.write((len('\t'.join(header)) + 1) * '-' + '\n')
newfile.write('\t'.join(header) + '\n')
newfile.write((len('\t'.join(header)) + 1) * '-' + '\n')
for i in range(len(comps)):
    chi2 = stats.chi2_contingency(data[i])[0]
    p = stats.chi2_contingency(data[i])[1]
    newfile.write('\t'.join([comps[i], str(round(chi2, 4)), str(round(p, 4))]) + '\n')
newfile.write((len('\t'.join(header)) + 1) * '-' + '\n') 

# close file after writing
newfile.close()

# compute average alpha for chemo and non-chemo genes for different syn polym cutoffs
newfile = open('TableSynPolymCutoffAlphaChemoNC.txt', 'w')
header = ['x', 'alpha GPCRs', 'alpha Non-chemo']
newfile.write((len('\t'.join(header)) + 1) * '-' + '\n')
newfile.write('\t'.join(header) + '\n')
newfile.write((len('\t'.join(header)) + 1) * '-' + '\n')

for i in range(10):
    ChemoAlpha = ComputeAlphaSEW2002(ChemoPolymDivCounts, i)
    NCAlpha = ComputeAlphaSEW2002(NCPolymDivCounts, i)
    newfile.write('\t'.join([str(i), str(round(ChemoAlpha, 4)), str(round(NCAlpha, 4))]) + '\n')
newfile.write((len('\t'.join(header)) + 1) * '-' + '\n')
newfile.close()
