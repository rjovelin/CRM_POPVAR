# -*- coding: utf-8 -*-
"""
Created on Sun Jun  5 23:31:07 2016

@author: Richard
"""

# use this script to 



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
from premature_stops import *


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


# make a dictionary with relative positions of PTC SNPs for chemo and non-chemo genes
ChemoPTC = RelativePositionPTC('../Genome_Files/CDS_SNP_DIVERG.txt',
                               GPCRs, '../Genome_Files/noamb_PX356_all_CDS.fasta')
NCPTC = RelativePositionPTC('../Genome_Files/CDS_SNP_DIVERG.txt',
                               NonGPCRs, '../Genome_Files/noamb_PX356_all_CDS.fasta')




#########################################

# express relative position in % of the CDS length
for i in range(len(ChemoPTC)):
    ChemoPTC[i] = ChemoPTC[i] * 100
for i in range(len(NCPTC)):
    NCPTC[i] = NCUPTC[i] * 100

# make a histogram
ChemoHist = np.histogram(ChemoUpstream, range(0, 110, 10))
NCHist = np.histogram(NCUpstream, range(0, 110, 10))

# transform the gene counts to proportions
ChemoFreq = [i / sum(ChemoHist[0]) for i in ChemoHist[0]]
NCFreq = [i / sum(NCHist[0]) for i in NCHist[0]]

print(len(NCFreq))
print(NCFreq)

# test uniformity of the distribution of the 5' most upstream PTC
chemotest = stats.chisquare(ChemoHist[0])
NCtest = stats.chisquare(NCHist[0])
print('chemo Pval uniformity', chemotest[1])
print('non-chemo Pval uniformity', NCtest[1])

# create figure
fig = plt.figure(1, figsize = (4, 2))
# add a plot to figure (1 row, 1 column, 1 plot)
ax = fig.add_subplot(1, 1, 1)  

# set width of bar
width = 0.1

# plot positions chemo
graph1 = ax.bar([i / 10 for i in range(10)], ChemoFreq, width, color = '#de2d26', edgecolor = '#de2d26', linewidth = 1, alpha = 0.7)
# plot positions non-chemo
graph2 = ax.bar([i / 10 for i in range(10)], NCFreq, width, color = '#3182bd', edgecolor = '#3182bd', linewidth = 1, alpha = 0.7)

ax.set_ylabel('Proportion of genes with a PTC', size = 10, ha = 'center', fontname = 'Arial')

# determine tick position on x axis
xpos =  [i / 10 for i in range(10)] + [1]
xtext = [0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
xtext = list(map(lambda x : str(x), xtext))
# set up tick positions and labels
plt.xticks(xpos, xtext, fontsize = 10, fontname = 'Arial')

# set x axis label
ax.set_xlabel('Decile of CDS length', size = 10, ha = 'center', fontname = 'Arial')

# do not show lines around figure, keep bottow line  
ax.spines["top"].set_visible(False)    
ax.spines["bottom"].set_visible(True)    
ax.spines["right"].set_visible(False)    
ax.spines["left"].set_visible(False)      
# offset the spines
for spine in ax.spines.values():
  spine.set_position(('outward', 5))
  
# add a light grey horizontal grid to the plot, semi-transparent, 
ax.yaxis.grid(True, linestyle='--', which='major', color='lightgrey', alpha=0.5, linewidth = 0.5)  
# hide these grids behind plot objects
ax.set_axisbelow(True)

# do not show ticks on 1st graph
ax.tick_params(
    axis='x',       # changes apply to the x-axis and y-axis (other option : x, y)
    which='both',      # both major and minor ticks are affected
    bottom='on',      # ticks along the bottom edge are off
    top='off',         # ticks along the top edge are off
    right = 'off',
    left = 'off',          
    labelbottom='on', # labels along the bottom edge are off 
    colors = 'black',
    labelsize = 10,
    direction = 'out') # ticks are outside the frame when bottom = 'on

# do not show ticks
ax.tick_params(
    axis='y',       # changes apply to the x-axis and y-axis (other option : x, y)
    which='both',      # both major and minor ticks are affected
    bottom='off',      # ticks along the bottom edge are off
    top='off',         # ticks along the top edge are off
    right = 'off',
    left = 'off',          
    labelbottom='off', # labels along the bottom edge are off 
    colors = 'black',
    labelsize = 10,
    direction = 'out') # ticks are outside the frame when bottom = 'on

for label in ax.get_yticklabels():
    label.set_fontname('Arial')

# create legend
chemoptc = mpatches.Patch(facecolor = '#de2d26' , edgecolor = 'grey', linewidth = 1, label= 'GPCR', alpha = 0.7)
NCptc = mpatches.Patch(facecolor = '#3182bd', edgecolor = 'grey', linewidth = 1, label = 'NC', alpha = 0.7)
plt.legend(handles=[chemoptc, NCptc], loc = 1, fontsize = 8, frameon = False)

# add margin on the x-axis
plt.margins(0.05)

fig.savefig('testfile.pdf', bbox_inches = 'tight')





#########################################


# express relative position in % of the CDS length
for i in range(len(upstream_PTC)):
    upstream_PTC[i] = upstream_PTC[i] * 100
# make a histogram
upstream_hist = np.histogram(upstream_PTC, range(0, 101, 10))
# open file to write the CDS truncations
newfile = open('CDS_truncation.txt', 'w')
newfile.write('PTC_relative_position' + '\t' + 'PTC_counts' + '\n')
for i in range(len(upstream_hist[0])):
    newfile.write(str(upstream_hist[1][i]) + ':' + str(upstream_hist[1][i] + 10) + '\t' + str(upstream_hist[0][i]) + '\n')
newfile.close()



# test uniformity of the distribution of the 5' most upstream PTC
test_uniform_upstream = stats.chisquare(upstream_hist[0])
# write results of test to summary file
summary_file.write('\n')
summary_file.write('test of uniformity of the distribution of 5\' most upstream PTC\n') 
summary_file.write('-' * 63 + '\n')
summary_file.write('chi2' + '\t' + 'p-val' + '\n')
summary_file.write(str(test_uniform_upstream[0]) + '\t' + str(test_uniform_upstream[1]) + '\n')





# compare expression level between PTC and non-PTC genes
# use all PTC genes, including genes with indels

# get the set of genes with PTC
PTC_genes = count_genes_with_indels_premature_stops('../CREM_CLA_protein_divergence/CDS_SNP_DIVERG.txt', '../CREM_CLA_protein_divergence/unique_transcripts.txt')
# get the set of valid transcripts
transcripts = get_valid_transcripts('../CREM_CLA_protein_divergence/unique_transcripts.txt')
# create a set of non-PTC genes
non_PTC_genes = set(gene for gene in transcripts if gene not in PTC_genes)
# create a dict of {gene: expression leve;}
crm_expression = expression_developmental_stages('../CREM_CLA_protein_divergence/WBPaper00041190.cre.mr.csv', 'c_remanei.PRJNA53967.WS248.geneIDs.txt')
# create lists with expression level for PTC and non-PTC genes
PTC_expression = [crm_expression[gene] for gene in PTC_genes if gene in crm_expression]
non_PTC_expression = [crm_expression[gene] for gene in non_PTC_genes if gene in crm_expression]
# test differences in gene expression between PTC and non-PTC genes
wilcoxon, p_val = stats.ranksums(PTC_expression, non_PTC_expression)
# write mean expression, SEM and differences to file

summary_file.write('\n')
summary_file.write('Expression level (includes PTC genes with indels):\n')
summary_file.write('-' * 51 + '\n')
summary_file.write('\t'.join(['genes', 'N', 'mean_expression', 'SEM']) + '\n')
# write the sample size, mean expression, SEM
summary_file.write('\t'.join(['PTC_genes', str(len(PTC_expression)), str(np.mean(PTC_expression)), str(np.std(PTC_expression) / math.sqrt(len(PTC_expression)))]) + '\n')
summary_file.write('\t'.join(['non_PTC_genes', str(len(non_PTC_expression)), str(np.mean(non_PTC_expression)), str(np.std(non_PTC_expression) / math.sqrt(len(non_PTC_expression)))]) + '\n')
# write results of statistical test
summary_file.write('Wilcoxon_sum_rank_test:' + '\t' + str(wilcoxon) + '\t' + str(p_val) + '\n')





# use this funtion to get the relative position of all PTC SNPs, including genes with alleles
def RelativePositionPTC(snp_file, GeneList, CDS_fasta):




