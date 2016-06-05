# -*- coding: utf-8 -*-
"""
Created on Sun Jun  5 14:32:34 2016

@author: Richard
"""


# use this script to plot the cumulative distribution function MAF for
# REP, SYN and PTC alleles for chemo and non-chemo genes

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
import sys
import json
# import custom modules
from manipulate_sequences import *
from miRNA_target import *
from genomic_coordinates import *
from sites_with_coverage import *
from divergence import *
from premature_stops import *
from randomize_SNPs import *
from get_coding_sequences import *
from mk_test import * 
from chemoreceptors import *


#######

# this is a modified function from MAF_SNP to compute MAF for a list of genes of interest    
def GenesMAF(snp_file, GeneList, indel_transcripts, site_type, threshold):
    '''
    (file, list, file, str, int) -> list
    Take the file of snps, a list of genes of interest, the file with
    the transcripts having indels, the type of mutation to record (REP, PTC
    or SYN), and a threshold of the minimum sample size to accept and return
    a list with minor allele frequencies in the KSR and PX strains combined
    Precondition: codons with more than 1 SNP are not considered
    '''
    
    # create a list to store the frequencies
    MAF = []    
        
    # get the set of transcripts with indels in their CDS
    indels = get_genes_with_indels(indel_transcripts)
    
    # create a set of valid nucleotides
    valid_base = {'A', 'T', 'C', 'G'}    
        
    # open file for reading
    infile = open(snp_file, 'r')
    # skip header
    infile.readline()
    # loop over file
    for line in infile:
        line = line.rstrip()
        if line != '':
            line = line.split()
            # check that gene is in GeneList
            gene = line[2]
            if gene in GeneList:
                # check only valid snps
                if line[6] == 'snp' and line[5] != line[7] and line[5] in valid_base and line[7] in valid_base:
                    # get ref codon and alternative codons
                    ref_codon, alt_codon = line[3], line[8]
                    # do not consider Ns in codons
                    if 'N' not in ref_codon and 'N' not in alt_codon:
                        # get allele counts
                        ref_count, alt_count = int(line[13]), int(line[17])
                        # verify that SNP is in KSR+PX
                        if alt_count != 0 and ref_count != 0:
                            # check that sample size is greater than threshold
                            if ref_count + alt_count >= threshold:
                                # check that ref and alt codons have only 1 difference
                                if diff_codon(ref_codon, alt_codon) <= 1:
                                    # don't count reference stop codons
                                    if genetic_code[ref_codon] != '*':
                                        # check site_type
                                        if site_type == 'PTC':
                                            # check that gene is not in indels
                                            if gene not in indels:
                                                # count premature stops
                                                # check that alternative allele is stop codon
                                                if genetic_code[alt_codon] == '*':
                                                    # check which allele is minor allele
                                                    if alt_count <= ref_count:
                                                        # alternative alelle is minor
                                                        # add minor frequency to list
                                                        MAF.append(alt_count / (alt_count + ref_count))
                                                    elif alt_count > ref_count:
                                                        # reference allele is minor
                                                        # add minor frequency to list
                                                        MAF.append(ref_count / (alt_count + ref_count))
                                        elif site_type == 'REP':
                                            # count replacement site
                                            # verify that change is non-synonymous but not PTC
                                            if genetic_code[ref_codon] != genetic_code[alt_codon] and genetic_code[alt_codon] != '*':
                                                # site is nonsynonymous
                                                # check which allele is minor 
                                                if alt_count <= ref_count:
                                                    # alternative allele is minor
                                                    # add minor frequency to list
                                                    MAF.append(alt_count / (alt_count + ref_count))
                                                elif alt_count > ref_count:
                                                    # reference codon is minor
                                                    # add minor frequency to list
                                                    MAF.append( ref_count / (alt_count + ref_count))
                                        elif site_type == 'SYN':
                                            # count synonymous sites
                                            # verify that site is synonymous
                                            if genetic_code[ref_codon] == genetic_code[alt_codon]:
                                                # site is synonymous
                                                # check which allele is minor
                                                if alt_count <= ref_count:
                                                    # alternative allele is minor
                                                    # add minor frequency to list
                                                    MAF.append(alt_count / (ref_count + alt_count))
                                                elif alt_count > ref_count:
                                                    # reference allele is minor
                                                    # add minor frequency to list
                                                    MAF.append(ref_count / (ref_count + alt_count))
                                       
    # close file after reading
    infile.close()
    
    return MAF

#######    


# create a set of valid transcripts
transcripts = get_valid_transcripts('../Genome_Files/unique_transcripts.txt')
print('got list of valid transcripts')

# get the set of chemoreceptors from the iprscan outputfile
chemo = get_chemoreceptors('../Genome_Files/PX356_protein_seq.tsv') 
print('parsed chemo genes')

# create a set of valid chemoreceptors
GPCRs = set(gene for gene in chemo if gene in transcripts)
print('GPCRS', len(GPCRs))
print('made set of valid chemo genes')

# create a set of non-chemo gene
NonGPCRs = set(gene for gene in transcripts if gene not in chemo)
print('NonGPCRs', len(NonGPCRs))
print('made set of valid non-chemo genes')
assert len(transcripts) == len(GPCRs) + len(NonGPCRs)

# make a list with MAF of replacement SNPs, excluding sites with sample size < 10
ChemoMAFREP = GenesMAF('../Genome_Files/CDS_SNP_DIVERG.txt', GPCRs, 
                  '../Genome_Files/transcripts_indels_CDS.txt', 'REP', 10)
print('REP Chemo', len(ChemoMAFREP))
print('got MAF Rep for chemo genes')

# make a list with MAF of synonymous sites, excluding sites with sample size < 10
ChemoMAFSYN = GenesMAF('../Genome_Files/CDS_SNP_DIVERG.txt', GPCRs, 
                  '../Genome_Files/transcripts_indels_CDS.txt', 'SYN', 10)
print('SYN Chemo', len(ChemoMAFSYN))
print('got MAF Syn for chemo genes')

# make a list with MAF of PTC sites, excluding sites with sample size < 10
ChemoMAFPTC = GenesMAF('../Genome_Files/CDS_SNP_DIVERG.txt', GPCRs, 
                  '../Genome_Files/transcripts_indels_CDS.txt', 'PTC', 10)
print('PTC Chemo', len(ChemoMAFPTC))
print('got MAF PTC for chemo genes')


# make a list with MAF of replacement SNPs, excluding sites with sample size < 10
NCMAFREP = GenesMAF('../Genome_Files/CDS_SNP_DIVERG.txt', NonGPCRs, 
                  '../Genome_Files/transcripts_indels_CDS.txt', 'REP', 10)
print('REP NC', len(NCMAFREP))
print('got MAF Rep for non-chemo genes')

# make a list with MAF of synonymous sites, excluding sites with sample size < 10
NCMAFSYN = GenesMAF('../Genome_Files/CDS_SNP_DIVERG.txt', NonGPCRs, 
                  '../Genome_Files/transcripts_indels_CDS.txt', 'SYN', 10)
print('SYN NC', len(NCMAFSYN))
print('got MAF Syn for non-chemo genes')

# make a list with MAF of PTC sites, excluding sites with sample size < 10
NCMAFPTC = GenesMAF('../Genome_Files/CDS_SNP_DIVERG.txt', NonGPCRs, 
                  '../Genome_Files/transcripts_indels_CDS.txt', 'PTC', 10)
print('PTC NC', len(NCMAFPTC))
print('got MAF PTC for non-chemo genes')


# sort MAF values
ChemoMAFREP.sort()
ChemoMAFSYN.sort()
ChemoMAFPTC.sort()
NCMAFREP.sort()
NCMAFSYN.sort()
NCMAFPTC.sort()
print('values are sorted')

# create figure
fig = plt.figure(1, figsize = (3.5, 2.5))
# add a plot to figure (1 row, 1 column, 1 plot)
ax = fig.add_subplot(1, 1, 1)  

# plot MAF synonymous sites
graph1 = ax.step(ChemoMAFSYN, np.linspace(0, 1, len(ChemoMAFSYN), endpoint=False), linewidth = 1.2, color = '#fb6a4a', alpha = 0.7)
# plot MAF replacement sites
graph2 = ax.step(ChemoMAFREP, np.linspace(0, 1, len(ChemoMAFREP), endpoint=False), linewidth = 1.2, color = '#cb181d', alpha = 0.7)
# plot MAF PTC
graph3 = ax.step(ChemoMAFPTC, np.linspace(0, 1, len(ChemoMAFPTC), endpoint=False), linewidth = 1.2, color = '#fcae91', alpha = 0.7)
# plot MAF syn NC
graph4 = ax.step(NCMAFSYN, np.linspace(0, 1, len(NCMAFSYN), endpoint=False), linewidth = 1.2, color = '#6baed6', alpha = 0.7)
# plot MAF rep NC
graph5 = ax.step(NCMAFREP, np.linspace(0, 1, len(NCMAFREP), endpoint=False), linewidth = 1.2, color = '#2171b5', alpha = 0.7)
# plot MAF ptc NC
graph6 = ax.step(NCMAFPTC, np.linspace(0, 1, len(NCMAFPTC), endpoint=False), linewidth = 1.2, color = '#bdd7e7', alpha = 0.7)
print('plotted CDF')


# add label for the Y axis
ax.set_ylabel('Proportion of SNPs', size = 10, ha = 'center', fontname = 'Arial')
# set x axis label
ax.set_xlabel('Minor Allele Frequency', size = 10, ha = 'center', fontname = 'Arial')

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

# add lines
lns = graph1+graph2+graph3+graph4+graph5+graph6
# get labels
labs = ['GPCR Syn.', 'GPCR Rep.', 'GPCR PTC', 'NC Syn.', 'NC Rep.', 'NC PTC']
# plot legend
ax.legend(lns, labs, loc=2, fontsize = 8, frameon = False)

fig.savefig('testfile.pdf', bbox_inches = 'tight')