# -*- coding: utf-8 -*-
"""
Created on Sun Jun  5 23:29:57 2016

@author: Richard
"""

# use this script to make a multiple plot figure comparing PTC between chemo and non-chemo genes 

# visit this page for plotting multiple subplots in a singlr figure
#http://matplotlib.org/users/gridspec.html


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
from crm_expression import *


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


# compare the proportion of PTC genes among chemo and non-chemo genes

# count the number of genes PTC SNPs, including genes with indels
PTCGenes = count_genes_with_indels_premature_stops('../Genome_Files/CDS_SNP_DIVERG.txt',
                                                      '../Genome_Files/unique_transcripts.txt')
                                                      
# sort PTC genes into chemo and non-chemo
ChemoPTC = [gene for gene in PTCGenes if gene in GPCRs]
print('chemo PTC', len(ChemoPTC))
NCPTC = [gene for gene in PTCGenes if gene in NonGPCRs]
print('non-chemo PTC', len(NCPTC))
# sort nonPTC genes into chemo and non-chemo genes
ChemoNonPTC = [gene for gene in GPCRs if gene not in PTCGenes]
print('chemo non PTC', len(ChemoNonPTC))
NCNonPTC = [gene for gene in NonGPCRs if gene not in PTCGenes]
print('non-chemo non-PTC', len(NCNonPTC))
assert len(transcripts) == len(ChemoPTC) + len(NCPTC) + len(ChemoNonPTC) + len(NCNonPTC), 'counts do not add up'

# test for difference in PTC genes among GPCRs and nonGPCR genes
Pchi2PTC = stats.chi2_contingency([[len(ChemoPTC), len(ChemoNonPTC)], [len(NCPTC), len(NCNonPTC)]])[1]
print('chi2', Pchi2PTC)

# get proportions of genes
ChemoPTCProp = len(ChemoPTC) / len(GPCRs)
ChemoNonPTCProp = len(ChemoNonPTC) / len(GPCRs)
NCPTCProp = len(NCPTC) / len(NonGPCRs)
NCNonPTCProp = len(NCNonPTC) / len(NonGPCRs)


# compare expression level between PTC and non-PTC genes
# use all PTC genes, including genes with indels

# create a dict of {gene: expression level}
Expression = expression_developmental_stages('../Genome_Files/WBPaper00041190.cre.mr.csv',
                                             '../Genome_Files/c_remanei.PRJNA53967.WS248.geneIDs.txt')

# create lists with expression for PTC and non-PTC chemo and non-chemo genes
ChemoPTCExp = [Expression[gene] for gene in Expression if gene in GPCRs and gene in PTCGenes]
ChemoNonPTCExp = [Expression[gene] for gene in Expression if gene in GPCRs and gene not in PTCGenes]
NCPTCExp = [Expression[gene] for gene in Expression if gene in NonGPCRs and gene in PTCGenes] 
NCNonPTCExp = [Expression[gene] for gene in Expression if gene in NonGPCRs and gene not in PTCGenes]

# test differences in gene expression between chemo PTC and non-PTC genes
PchemoPtc = stats.ranksums(ChemoPTCExp, ChemoNonPTCExp)[1]
PNCPtc = stats.ranksums(NCPTCExp, NCNonPTCExp)[1]
print('chemo expression PTC', np.mean(ChemoPTCExp), np.mean(ChemoNonPTCExp), PchemoPtc)
print('nonchemo expression PTC', np.mean(NCPTCExp), np.mean(NCNonPTCExp), PNCPtc)


# compute the distribution of position of the 5' miost upstream PTC

# get the relative position of the 5' most PTC allele in chemo and non-chemo genes
ChemoUpstream = position_first_PTC('../Genome_Files/CDS_SNP_DIVERG.txt', GPCRs,
                                   '../Genome_Files/transcripts_indels_CDS.txt', '../Genome_Files/noamb_PX356_all_CDS.fasta')
print('grabbed positions of 5\' PTC alleles for chemo genes')

NCUpstream = position_first_PTC('../Genome_Files/CDS_SNP_DIVERG.txt', NonGPCRs,
                                '../Genome_Files/transcripts_indels_CDS.txt', '../Genome_Files/noamb_PX356_all_CDS.fasta')
print('grabbed positions of 5\' PTC alleles for non-chemo genes')

# express relative position in % of the CDS length
for i in range(len(ChemoUpstream)):
    ChemoUpstream[i] = ChemoUpstream[i] * 100
for i in range(len(NCUpstream)):
    NCUpstream[i] = NCUpstream[i] * 100

# make a histogram
ChemoHist = np.histogram(ChemoUpstream, range(0, 110, 10))
NCHist = np.histogram(NCUpstream, range(0, 110, 10))

# transform the gene counts to proportions
ChemoFreq = [i / sum(ChemoHist[0]) for i in ChemoHist[0]]
NCFreq = [i / sum(NCHist[0]) for i in NCHist[0]]


# test uniformity of the distribution of the 5' most upstream PTC
chemotest = stats.chisquare(ChemoHist[0])
NCtest = stats.chisquare(NCHist[0])
print('chemo Pval uniformity', chemotest[1])
print('non-chemo Pval uniformity', NCtest[1])

# compare the CDF of PTC position along chemo and non-chemo genes

# make a dictionary with relative positions of PTC SNPs for chemo and non-chemo genes
ChemoPTCPos = RelativePositionPTC('../Genome_Files/CDS_SNP_DIVERG.txt',
                               GPCRs, '../Genome_Files/noamb_PX356_all_CDS.fasta')
print('recorded all PTC positions for chemo genes')

NCPTCPos = RelativePositionPTC('../Genome_Files/CDS_SNP_DIVERG.txt',
                               NonGPCRs, '../Genome_Files/noamb_PX356_all_CDS.fasta')
print('recorded all PTC postions for non-chemo genes')

# express relative position in % of the CDS length
for i in ChemoPTCPos:
    ChemoPTCPos[i] = ChemoPTCPos[i] * 100
for i in NCPTCPos:
    NCPTCPos[i] = NCPTCPos[i] * 100

# make a histogram
ChemoPosHist = np.histogram(ChemoUpstream, range(0, 105, 5))
NCPosHist = np.histogram(NCUpstream, range(0, 105, 5))

# transform the gene counts to proportions
ChemoPosFreq = [i / sum(ChemoPosHist[0]) for i in ChemoPosHist[0]]
NCPosFreq = [i / sum(NCPosHist[0]) for i in NCPosHist[0]]
# sort values
ChemoPosFreq.sort()
NCPosFreq.sort()
print('values are sorted')
print(len(NCPosFreq))
print(NCPosFreq)


# create figure
fig = plt.figure(1, figsize = (6, 6))


# create a function to format the subplots
def FormatAx(AxNum, XTicksPos, XTicklabels, Xscale, Data, figure, width, colorscheme,
             XLabel, YLabel, isXLabel = True, Legend = False, FirstData = '', AnnotateP = False):
    
    # use subplot2grid to specify position of subplots (# row, # column)
    if AxNum == 1:
        ax = plt.subplot2grid((3,2), (0,0))
        # plot proportions of genes
        ax.bar(Xscale, Data, width, color = colorscheme, edgecolor = 'black', linewidth = 1)
    elif AxNum == 2:
        ax = plt.subplot2grid((3,2), (0,0))
        # plot proportion of genes on top of the other proportions 
        ax.bar(Xscale, Data, width, color = colorscheme, edgecolor = 'black', linewidth = 1, bottom = FirstData)
    elif AxNum == 3:
        ax = plt.subplot2grid((3,2), (0,1))
        # plot bar graphs comparing expression Means = Data[0], SEM = Data[1] 
        ax.bar(Xscale, Data[0], width, color = colorscheme, yerr = Data[1], edgecolor = 'black', linewidth = 1,
               error_kw=dict(elinewidth=1, ecolor='black', markeredgewidth = 1))
    elif AxNum == 4:
        # plot the % truncation 
        ax = plt.subplot2grid((3,2), (1,0))
        ax.bar(Xscale, Data, width, color = colorscheme, edgecolor = 'black', linewidth = 1)
    elif AxNum == 5:
        # plot the % truncation for chemo genes
        ax = plt.subplot2grid((3,2), (2,0))
        ax.bar(Xscale, Data, width, color = colorscheme, edgecolor = 'black', linewidth = 1)
    elif AxNum >= 6:
        ax = plt.subplot2grid((3,2), (1,1), rowspan = 2)
        ax.step(Data, np.linspace(0, 1, len(Data), endpoint = False), linewidth = 1.2, color = colorscheme, alpha = 0.7)
    
        
    # set Y label
    ax.set_ylabel(YLabel, size = 10, ha = 'center', fontname = 'Arial')
    # set x ticks
    plt.xticks(XTicksPos, XTicklabels, fontsize = 10, fontname = 'Arial')
    # set x axis label
    if isXLabel == True:
        ax.set_xlabel(XLabel, size = 10, ha = 'center', fontname = 'Arial')

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

#    if Legend == True:
#        # create legend
#        chemoptc = mpatches.Patch(facecolor = '#de2d26' , edgecolor = 'grey', linewidth = 1, label= 'GPCR', alpha = 0.7)
#        NCptc = mpatches.Patch(facecolor = '#3182bd', edgecolor = 'grey', linewidth = 1, label = 'NC', alpha = 0.7)
#        plt.legend(handles=[chemoptc, NCptc], loc = 1, fontsize = 8, frameon = False)

    # add margin on the x-axis
    #plt.margins(0.05)

    return ax





# plot proportions of PTC
PTCProp = [ChemoPTCProp, NCPTCProp]
ax1 = FormatAx(1, [0.1, 0.3], ['GPCRs', 'NC'], [0, 0.2], [ChemoPTCProp, NCPTCProp], fig, 0.2, ['#de2d26', '#3182bd'],
               '', 'Proportion of genes', isXLabel = False, Legend = False, FirstData = '', AnnotateP = False)
print('plotted graph 1')

ax2 = FormatAx(2, [0.1, 0.3], ['GPCRs', 'NC'], [0, 0.2], [ChemoNonPTCProp, NCNonPTCProp], fig, 0.2, ['#fee0d2', '#deebf7'],
               '', 'Proportion of genes', isXLabel = False, Legend = False, FirstData = PTCProp, AnnotateP = False)
print('plotted graph 2')

# plot expression level
# create a list of expression data
expdata = [[np.mean(ChemoPTCExp), np.mean(ChemoNonPTCExp), np.mean(NCPTCExp), np.mean(NCNonPTCExp)],
          [np.std(ChemoPTCExp) / math.sqrt(len(ChemoPTCExp)), np.std(ChemoNonPTCExp) / math.sqrt(len(ChemoNonPTCExp)),
           np.std(NCPTCExp) / math.sqrt(len(NCPTCExp)), np.std(NCNonPTCExp) / math.sqrt(len(NCNonPTCExp))]]
ax3 = FormatAx(3, [0.2, 0.7], ['GPCRs', 'NC'], [0, 0.2, 0.5, 0.7], expdata, fig, 0.2, ['#de2d26', '#fee0d2', '#3182bd', '#deebf7'],
               '', 'Expression level', isXLabel = False, Legend = False, FirstData = '', AnnotateP = False)

print('plotted graph 3')

# plot truncation for chemo genes
ax4 = FormatAx(4, [i / 10 for i in range(10)] + [1], ['0', '10', '20', '30', '40', '50', '60', '70', '80', '90', '100'], [i / 10 for i in range(10)], ChemoFreq, fig, 0.1, '#de2d26',
               'Decile of CDS length', 'Proportion of genes with a PTC', isXLabel = True, Legend = False, FirstData = '', AnnotateP = False)

print('plotted graph 4')

ax5 = FormatAx(5, [i / 10 for i in range(10)] + [1], ['0', '10', '20', '30', '40', '50', '60', '70', '80', '90', '100'], [i / 10 for i in range(10)], NCFreq, fig, 0.1, '#3182bd',
               'Decile of CDS length', 'Proportion of genes with a PTC', isXLabel = True, Legend = False, FirstData = '', AnnotateP = False)

print('plotted graph 5')

# plot the CDS of PTC allele counts
ax6 = FormatAx(6, [i / 10 for i in range(11)], ['0', '10', '20', '30', '40', '50', '60', '70', '80', '90', '100'], [i / 10 for i in range(20)], ChemoPosFreq, fig, 0.1, '#de2d26',
               'CDS length', 'Proportions of PTC alleles', isXLabel = True, Legend = False, FirstData = '', AnnotateP = False)

print('plotted graph 6')


ax7 = FormatAx(7, [i / 10 for i in range(11)], ['0', '10', '20', '30', '40', '50', '60', '70', '80', '90', '100'], [i / 10 for i in range(20)], NCPosFreq, fig, 0.1, '#3182bd',
               'CDS length', 'Proportions of PTC alleles', isXLabel = True, Legend = False, FirstData = '', AnnotateP = False)

print('plotted graph 7')


############ code below works
#
## set width of bar
#width = 0.1
#
## plot positions chemo
#graph1 = ax.bar([i / 10 for i in range(10)], ChemoFreq, width, color = '#de2d26', edgecolor = '#de2d26', linewidth = 1, alpha = 0.7)
## plot positions non-chemo
#graph2 = ax.bar([i / 10 for i in range(10)], NCFreq, width, color = '#3182bd', edgecolor = '#3182bd', linewidth = 1, alpha = 0.7)
#
#ax.set_ylabel('Proportion of genes with a PTC', size = 10, ha = 'center', fontname = 'Arial')
#
## determine tick position on x axis
#xpos =  [i / 10 for i in range(10)] + [1]
#xtext = [0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
#xtext = list(map(lambda x : str(x), xtext))
## set up tick positions and labels
#plt.xticks(xpos, xtext, fontsize = 10, fontname = 'Arial')
#
## set x axis label
#ax.set_xlabel('Decile of CDS length', size = 10, ha = 'center', fontname = 'Arial')
#
## do not show lines around figure, keep bottow line  
#ax.spines["top"].set_visible(False)    
#ax.spines["bottom"].set_visible(True)    
#ax.spines["right"].set_visible(False)    
#ax.spines["left"].set_visible(False)      
## offset the spines
#for spine in ax.spines.values():
#  spine.set_position(('outward', 5))
#  
## add a light grey horizontal grid to the plot, semi-transparent, 
#ax.yaxis.grid(True, linestyle='--', which='major', color='lightgrey', alpha=0.5, linewidth = 0.5)  
## hide these grids behind plot objects
#ax.set_axisbelow(True)
#
## do not show ticks on 1st graph
#ax.tick_params(
#    axis='x',       # changes apply to the x-axis and y-axis (other option : x, y)
#    which='both',      # both major and minor ticks are affected
#    bottom='on',      # ticks along the bottom edge are off
#    top='off',         # ticks along the top edge are off
#    right = 'off',
#    left = 'off',          
#    labelbottom='on', # labels along the bottom edge are off 
#    colors = 'black',
#    labelsize = 10,
#    direction = 'out') # ticks are outside the frame when bottom = 'on
#
## do not show ticks
#ax.tick_params(
#    axis='y',       # changes apply to the x-axis and y-axis (other option : x, y)
#    which='both',      # both major and minor ticks are affected
#    bottom='off',      # ticks along the bottom edge are off
#    top='off',         # ticks along the top edge are off
#    right = 'off',
#    left = 'off',          
#    labelbottom='off', # labels along the bottom edge are off 
#    colors = 'black',
#    labelsize = 10,
#    direction = 'out') # ticks are outside the frame when bottom = 'on
#
#for label in ax.get_yticklabels():
#    label.set_fontname('Arial')
#
## create legend
#chemoptc = mpatches.Patch(facecolor = '#de2d26' , edgecolor = 'grey', linewidth = 1, label= 'GPCR', alpha = 0.7)
#NCptc = mpatches.Patch(facecolor = '#3182bd', edgecolor = 'grey', linewidth = 1, label = 'NC', alpha = 0.7)
#plt.legend(handles=[chemoptc, NCptc], loc = 1, fontsize = 8, frameon = False)
#
## add margin on the x-axis
#plt.margins(0.05)

# make sure subplots do not overlap
plt.tight_layout()


fig.savefig('testfile.pdf', bbox_inches = 'tight')



