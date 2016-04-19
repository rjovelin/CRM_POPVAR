# -*- coding: utf-8 -*-
"""
Created on Mon Jun 29 19:27:18 2015

@author: Richard
"""



from chemoreceptors import *
from accessories import *
from diversity_chemo_partitions import *
import os
import numpy as np
from scipy import stats
import math


# get the set of chemoreceptors from the iprscan outputfile
chemo = get_chemoreceptors('./PX356_protein_seq.tsv') 

# create a set of valid transcripts
transcripts = get_valid_transcripts('../CREM_CLA_protein_divergence/unique_transcripts.txt')

# create a set of valid chemoreceptors
GPCRs = set(gene for gene in chemo if gene in transcripts)

# create a dict proba_domain {gene: {codon_index: [list of probabilities]}}
proba_domain = {}

# create list of phobius output files
files = [filename for filename in os.listdir('./Chemo_genes/') if '_proba.txt' in filename]

# loop over genes in GPCRs
for gene in GPCRs:
    # check has predictions
    for filename in files:
        if gene == filename[:filename.index('_proba')]:
            # initialize dict
            proba_domain[gene] = {}
            # get the dict of probabilities
            probabilities = parse_phobius_output(gene + '_proba.txt', './Chemo_genes/')
            # populate dict
            for key, val in probabilities.items():
                proba_domain[gene][key] = val

# count SNPs, degenerate sites for partitions for nonsynonymous sites
TM_rep, EX_rep = count_SNPs_degenerate_sites_chemo_paritions('../CREM_CLA_protein_divergence/CDS_SNP_DIVERG.txt', proba_domain, 'REP', 0.95)
# compute theta at replacement sites for partitions
# accept a minimum of 15 sites
TM_theta_rep, EX_theta_rep = compute_theta_chemo_partitions(TM_rep, EX_rep, 'REP', 15)

# count SNPs, degenerate sites for partitions for synonymous sites
TM_syn, EX_syn = count_SNPs_degenerate_sites_chemo_paritions('../CREM_CLA_protein_divergence/CDS_SNP_DIVERG.txt', proba_domain, 'SYN', 0.95)
# compute theta at synonymous sites for partitions
# accept a minum of 15 sites
TM_theta_syn, EX_theta_syn = compute_theta_chemo_partitions(TM_syn, EX_syn, 'SYN', 15)

# count SNPs for the CDS in each partition
TM_cod, EX_cod = count_SNPs_degenerate_sites_chemo_paritions('../CREM_CLA_protein_divergence/CDS_SNP_DIVERG.txt', proba_domain, 'coding', 0.95)
# compute theta for CDS, accept a minimum of 15 sites
TM_theta_cod, EX_theta_cod = compute_theta_chemo_partitions(TM_cod, EX_cod, 'coding', 15)



# create lists for theta in different partitions keeping the same gene between lists
theta_rep_TM = []
theta_rep_EX = []
for gene in TM_theta_rep:
    theta_rep_TM.append(TM_theta_rep[gene])
    theta_rep_EX.append(EX_theta_rep[gene])
    
theta_syn_TM = []
theta_syn_EX = []
for gene in TM_theta_syn:
    theta_syn_TM.append(TM_theta_syn[gene])
    theta_syn_EX.append(EX_theta_syn[gene])
    
theta_cod_TM = []
theta_cod_EX = []
for gene in TM_theta_cod:
    theta_cod_TM.append(TM_theta_cod[gene])
    theta_cod_EX.append(EX_theta_cod[gene])
    

# open file, dump theta for replacement sites
newfile = open('chemo_theta_rep_partitions.txt', 'w')
for theta in theta_rep_TM:
    newfile.write(str(theta) + '\t' + 'TM' + '\n')
for theta in theta_rep_EX:
    newfile.write(str(theta) + '\t' + 'EX' + '\n')
newfile.close()

# open file, dump theta for synonymous sites
newfile = open('chemo_theta_syn_partitions.txt', 'w')
for theta in theta_syn_TM:
    newfile.write(str(theta) + '\t' + 'TM' + '\n')
for theta in theta_syn_EX:
    newfile.write(str(theta) + '\t' + 'EX' + '\n')
newfile.close()

# open file, dump theta for CDS
newfile = open('chemo_theta_coding_partitions.txt', 'w')
for theta in theta_cod_TM:
    newfile.write(str(theta) + '\t' + 'TM' + '\n')
for theta in theta_cod_EX:
    newfile.write(str(theta) + '\t' + 'EX' + '\n')
newfile.close()


# create dicts {gene: theta_rep / theta_syn} for each partition
TM_omega , EX_omega = {}, {}
# loop over genes in TM_theta_rep
for gene in TM_theta_rep:
    # check if gene TM_theta_syn
    if gene in TM_theta_syn:
        # check that theta syn different than 0
        if TM_theta_syn[gene] != 0:
            TM_omega[gene] = TM_theta_rep[gene] / TM_theta_syn[gene]
# loop over genes in EX_theta_rep
for gene in EX_theta_rep:
    # check if gene in EX_theta_syn
    if gene in EX_theta_syn:
        # check if theta syn is defined
        if EX_theta_syn[gene] != 0:
            EX_omega[gene] = EX_theta_rep[gene] / EX_theta_syn[gene]

# create lists of theta ratio for each partition, keeping the same gene order
theta_omega_TM , theta_omega_EX = [], []
for gene in TM_omega:
    # check that gene in EX_omega
    if gene in EX_omega:
        theta_omega_TM.append(TM_omega[gene])
        theta_omega_EX.append(EX_omega[gene])

# open file, dump theta ratios
newfile = open('chemo_theta_omega_partitions.txt', 'w')
for theta in theta_omega_TM:
    newfile.write(str(theta) + '\t' + 'TM' + '\n')
for theta in theta_omega_EX:
    newfile.write(str(theta) + '\t' + 'EX' + '\n')
newfile.close()


# open file to write summary of analysis
newfile = open('summary_chemo_partitions_diversity.txt', 'w')

newfile.write('comparison theta nonsynonymous TM & EX:\n')
newfile.write('-' * 40 + '\n')

newfile.write('\t'.join(['\t', 'N', 'mean', 'standard_error']) + '\n')
newfile.write('\t'.join(['theta_nonsyn_TM', str(len(theta_rep_TM)), str(np.mean(theta_rep_TM)), str(np.std(theta_rep_TM) / math.sqrt(len(theta_rep_TM)))]) + '\n')
newfile.write('\t'.join(['theta_nonsyn_EX', str(len(theta_rep_EX)), str(np.mean(theta_rep_EX)), str(np.std(theta_rep_EX) / math.sqrt(len(theta_rep_EX)))]) + '\n')
# test difference using paired test
newfile.write('Wilcoxon paired test: ' + '\t' + str(stats.wilcoxon(theta_rep_TM, theta_rep_EX)[1]) + '\n')

newfile.write('\n')


newfile.write('comparison theta synonymous TM & EX:\n')
newfile.write('-' * 37 + '\n')

newfile.write('\t'.join(['\t', 'N', 'mean', 'standard_error']) + '\n')
newfile.write('\t'.join(['theta_nonsyn_TM', str(len(theta_syn_TM)), str(np.mean(theta_syn_TM)), str(np.std(theta_syn_TM) / math.sqrt(len(theta_syn_TM)))]) + '\n')
newfile.write('\t'.join(['theta_nonsyn_EX', str(len(theta_syn_EX)), str(np.mean(theta_syn_EX)), str(np.std(theta_syn_EX) / math.sqrt(len(theta_syn_EX)))]) + '\n')
# test difference using paired test
newfile.write('Wilcoxon paired test: ' + '\t' + str(stats.wilcoxon(theta_syn_TM, theta_syn_EX)[1]) + '\n')

newfile.write('\n')

newfile.write('comparison theta nonsynonymous / theta syn TM & EX:\n')
newfile.write('-' * 45 + '\n')

newfile.write('\t'.join(['\t', 'N', 'mean', 'standard_error']) + '\n')
newfile.write('\t'.join(['theta_nonsyn_TM', str(len(theta_omega_TM)), str(np.mean(theta_omega_TM)), str(np.std(theta_omega_TM) / math.sqrt(len(theta_omega_TM)))]) + '\n')
newfile.write('\t'.join(['theta_nonsyn_EX', str(len(theta_omega_EX)), str(np.mean(theta_omega_EX)), str(np.std(theta_omega_EX) / math.sqrt(len(theta_omega_EX)))]) + '\n')
# test difference using paired test
newfile.write('Wilcoxon paired test: ' + '\t' + str(stats.wilcoxon(theta_omega_TM, theta_omega_EX)[1]) + '\n')

newfile.write('\n')

newfile.write('comparison theta coding sites TM & EX:\n')
newfile.write('-' * 39 + '\n')

newfile.write('\t'.join(['\t', 'N', 'mean', 'standard_error']) + '\n')
newfile.write('\t'.join(['theta_nonsyn_TM', str(len(theta_cod_TM)), str(np.mean(theta_cod_TM)), str(np.std(theta_cod_TM) / math.sqrt(len(theta_cod_TM)))]) + '\n')
newfile.write('\t'.join(['theta_nonsyn_EX', str(len(theta_cod_EX)), str(np.mean(theta_cod_EX)), str(np.std(theta_cod_EX) / math.sqrt(len(theta_cod_EX)))]) + '\n')
newfile.write('Wilcoxon paired test: ' + '\t' + str(stats.wilcoxon(theta_cod_TM, theta_cod_EX)[1]) + '\n')


newfile.close()