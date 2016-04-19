# -*- coding: utf-8 -*-
"""
Created on Sat Sep 12 19:59:58 2015

@author: Richard
"""

from caeno_chromosomes import *
from divergence import *
from accessories import *
from chemoreceptors import *
import math
import numpy as np
from scipy import stats


# compute theta at synonymous sites
theta_syn = compute_theta_diversity('CDS_SNP_DIVERG.txt', 'unique_transcripts.txt', 'SYN', 10)
print('computed theta at synonymous sites')
# compute theta at nonsynonymous sites
theta_rep = compute_theta_diversity('CDS_SNP_DIVERG.txt', 'unique_transcripts.txt', 'REP', 10)
print('computed theta at replacement sites')
# sort genes according to their chromosomal origin
# make list of genes in the theta dicts                                  
X_genes_syn, autosomal_genes_syn = X_autosomal_genes(theta_syn, '356_10172014.gff3', 'scaffold_chromosome.csv')
X_genes_rep, autosomal_genes_rep = X_autosomal_genes(theta_rep, '356_10172014.gff3', 'scaffold_chromosome.csv')
# make sets of genes from lists
X_genes_syn = set(X_genes_syn)
X_genes_rep = set(X_genes_rep)
autosomal_genes_syn = set(autosomal_genes_syn)
autosomal_genes_rep = set(autosomal_genes_rep)
# make sets of sex-linked and autosomal genes by combining sets
X_genes = X_genes_syn.intersection(X_genes_rep)
autosomal_genes = autosomal_genes_syn.intersection(autosomal_genes_rep)
print('sorted genes according to linkage')
print('X-linked', len(X_genes))
print('autosomal ', len(autosomal_genes))


# generate dicts with gene : theta pairs with genes sorted according to linkage
theta_syn_X, theta_syn_auto, theta_rep_X, theta_rep_auto = {}, {}, {}, {}
# sort syn values
for gene in theta_syn:
    if gene in X_genes:
        theta_syn_X[gene] = theta_syn[gene]
    elif gene in autosomal_genes:
        theta_syn_auto[gene] = theta_syn[gene]

# sort rep values
for gene in theta_rep:
    if gene in X_genes:
        theta_rep_X[gene] = theta_rep[gene]
    elif gene in autosomal_genes:
        theta_rep_auto[gene] = theta_rep[gene]

# generate dict with rep/ syn ratio
theta_omega = {}
for gene in theta_rep:
    if gene in theta_syn and theta_syn[gene] != 0:
        theta_omega[gene] = theta_rep[gene] / theta_syn[gene]

# sort omega values
theta_omega_X, theta_omega_auto = {}, {}
for gene in theta_omega:
    if gene in X_genes:
        theta_omega_X[gene] = theta_omega[gene]
    elif gene in autosomal_genes:
        theta_omega_auto[gene] = theta_omega[gene]

# save theta values to files
# generate files with genes : theta
newfile = open('theta_X_rep.txt', 'w')
for gene in theta_rep_X:
    newfile.write(gene + '\t' + str(theta_rep_X[gene]) + '\n')
newfile.close()

newfile = open('theta_auto_rep.txt', 'w')
for gene in theta_rep_auto:
    newfile.write(gene + '\t' + str(theta_rep_auto[gene]) + '\n')
newfile.close()

newfile = open('theta_X_syn.txt', 'w')
for gene in theta_syn_X:
    newfile.write(gene + '\t' + str(theta_syn_X[gene]) + '\n')
newfile.close()

newfile = open('theta_auto_syn.txt', 'w')
for gene in theta_syn_auto:
    newfile.write(gene + '\t' + str(theta_syn_auto[gene]) + '\n')
newfile.close()

newfile = open('theta_X_omega.txt', 'w')
for gene in theta_omega_X:
    newfile.write(gene + '\t' + str(theta_omega_X[gene]) + '\n')
newfile.close()

newfile = open('theta_auto_omega.txt', 'w')
for gene in theta_omega_auto:
    newfile.write(gene + '\t' + str(theta_omega_auto[gene]) + '\n')
newfile.close()

# make list of theta values
SYN_X = [theta_syn_X[gene] for gene in theta_syn_X]
REP_X = [theta_rep_X[gene] for gene in theta_rep_X]
SYN_auto = [theta_syn_auto[gene] for gene in theta_syn_auto]
REP_auto = [theta_rep_auto[gene] for gene in theta_rep_auto]
OMEGA_X = [theta_omega_X[gene] for gene in theta_omega_X]
OMEGA_auto = [theta_omega_auto[gene] for gene in theta_omega_auto]

# open summary file
summary_file = open('summary_X-auto_comparison.txt', 'w')

# compare diversity between X-linked genes and autosomal genes
summary_file.write('comparison of nucleotide diversity between X-linked and autosomal genes\n')
summary_file.write('-' * 72 + '\n')

# create a list of theta
theta_sites = [SYN_X, SYN_auto, REP_X, REP_auto, OMEGA_X, OMEGA_auto]
# create a list of corrsponding site type
site_types = ['SYN_X', 'SYN_auto', 'REP_X', 'REP_auto', 'OMEGA_X', 'OMEGA_auto']

# write header
summary_file.write('sites' + '\t' + 'N' + '\t' + 'mean_theta' + '\t' + 'SEM' + '\t' + 'median_theta' + '\n')
for i in range(len(theta_sites)):
    summary_file.write('\t'.join([site_types[i], str(len(theta_sites[i])), str(np.mean(theta_sites[i])), str(np.std(theta_sites[i]) / math.sqrt(len(theta_sites[i]))), str(np.median(theta_sites[i]))]) + '\n')
    
summary_file.write('\n')

summary_file.write('Wilcoxon rank sum test of mean differences\n')
summary_file.write('-' * 43 + '\n')
summary_file.write('sites' + '\t' + 'wilcoxon' + '\t' + 'P' + '\n')

# loop over list, compute pairwise difference
for i in range(0, len(site_types), 2):
    wilcoxon, p = stats.ranksums(theta_sites[i], theta_sites[i+1])
    summary_file.write('\t'.join([site_types[i] + '_vs_' + site_types[i+1], str(wilcoxon), str(p)]) + '\n')


# compute divergence for X- linked and autos dN, dS and omega
# parse divergence file
divergence = parse_divergence_file('CRM_CLA_prot_diverg_filtered.txt')
# create dicts of divergence values for each gene
dN_X, dN_auto, dS_X, dS_auto, omega_X, omega_auto = {}, {}, {}, {}, {}, {}
# populate dict
for gene in divergence:
    if gene in X_genes:
        dN_X[gene] = divergence[gene][0]
        dS_X[gene] = divergence[gene][1]
        # filter omega values, keep values < 5
        if divergence[gene][2] != 'NA' and divergence[gene][2] < 5:
            omega_X[gene] = divergence[gene][2]
    elif gene in autosomal_genes:
        dN_auto[gene] = divergence[gene][0]
        dS_auto[gene] = divergence[gene][1]
        # filter omega values, keep values < 5
        if divergence[gene][2] != 'NA' and divergence[gene][2] < 5:
            omega_auto[gene] = divergence[gene][2]             
            
# write divergence values to file
newfile = open('dN_X.txt', 'w')
for gene in dN_X:
    newfile.write(gene + '\t' + str(dN_X[gene]) + '\n')
newfile.close()

newfile = open('dN_auto.txt', 'w')
for gene in dN_auto:
    newfile.write(gene + '\t' + str(dN_auto[gene]) + '\n')
newfile.close()

newfile = open('dS_X.txt', 'w')
for gene in dS_X:
    newfile.write(gene + '\t' + str(dS_X[gene]) + '\n')
newfile.close()

newfile = open('dS_auto.txt', 'w')
for gene in dS_auto:
    newfile.write(gene + '\t' + str(dS_auto[gene]) + '\n')
newfile.close()

newfile = open('dN_dS_X.txt', 'w')
for gene in omega_X:
    newfile.write(gene + '\t' + str(omega_X[gene]) + '\n')
newfile.close()

newfile = open('dN_dS_auto.txt', 'w')
for gene in omega_auto:
    newfile.write(gene + '\t' + str(omega_auto[gene]) + '\n')
newfile.close()


summary_file.write('\n')

# compare divergence between X-linked genes and autosomal genes
summary_file.write('comparison of interspecies divergence between X-linked and autosomal genes\n')
summary_file.write('-' * 75 + '\n')

# create list of divergence values
DN_X = [dN_X[gene] for gene in dN_X]
DN_auto = [dN_auto[gene] for gene in dN_auto]
DS_X = [dS_X[gene] for gene in dS_X]
DS_auto = [dS_auto[gene] for gene in dS_auto]
OMEGA_X = [omega_X[gene] for gene in omega_X]
OMEGA_auto = [omega_auto[gene] for gene in omega_auto]

# create a list of theta
div_sites = [DS_X, DS_auto, DN_X, DN_auto, OMEGA_X, OMEGA_auto]
# create a list of corrsponding site type
site_types = ['dS_X', 'dS_auto', 'dN_X', 'dN_auto', 'OMEGA_X', 'OMEGA_auto']

# write header
summary_file.write('sites' + '\t' + 'N' + '\t' + 'mean_divergence' + '\t' + 'SEM' + '\t' + 'median_divergence' + '\n')
for i in range(len(div_sites)):
    summary_file.write('\t'.join([site_types[i], str(len(div_sites[i])), str(np.mean(div_sites[i])), str(np.std(div_sites[i]) / math.sqrt(len(div_sites[i]))), str(np.median(div_sites[i]))]) + '\n')
    
summary_file.write('\n')

summary_file.write('Wilcoxon rank sum test of mean differences\n')
summary_file.write('-' * 43 + '\n')
summary_file.write('sites' + '\t' + 'wilcoxon' + '\t' + 'P' + '\n')

# loop over list, compute pairwise difference
for i in range(0, len(site_types), 2):
    wilcoxon, p = stats.ranksums(div_sites[i], div_sites[i+1])
    summary_file.write('\t'.join([site_types[i] + '_vs_' + site_types[i+1], str(wilcoxon), str(p)]) + '\n')

summary_file.close()