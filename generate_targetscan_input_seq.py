# -*- coding: utf-8 -*-
"""
Created on Thu May 28 14:01:09 2015

@author: Richard
"""


import os
from miRNA_target import *
from genomic_coordinates import *
from protein_divergence import *
from TcoffeeAlignment import *


# compute threshold based on the distribution of elegans UTR length
UTR_length = celegans_three_prime_UTR_length('c_elegans.PRJNA13758.WS248.annotations.gff3')
threshold = get_percentile(UTR_length, 99)


# generate the Targetscan input file for all the remanei transcripts

# create a dict of transcript : gene  pairs
crm_transcripts = transcript_to_gene('../CREM_CLA_protein_divergence/356_10172014.gff3')
# create a dict of remanei UTR sequences
crm_UTR = fetch_UTR_sequences('../CREM_CLA_protein_divergence/356_10172014.gff3', '../CREM_CLA_protein_divergence/noamb_356_v1_4.txt', threshold)
# create targetscan input sequence file
make_cremanei_targetscan_input_sequences(crm_UTR, crm_transcripts, '../CREM_CLA_protein_divergence/crem_cla_orthologs.txt', 'Crm_UTR_seq_targetscan.txt')
print('remanei targetscan input done')



# generate the Targetscan input file for the remanei and latens orthologs

# create a folder to store the ortholog pairs between remanei and latens
os.mkdir('./Crm_Cla_UTR_sequences/')

# create a dict of latens UTR sequences
cla_UTR = fetch_UTR_sequences('../CREM_CLA_protein_divergence/534_10172014.gff3', '../CREM_CLA_protein_divergence/noamb_534_v1.txt', threshold)

# create a dict of orthologous pairs {crem_TS: cla_TS}
crm_cla_orthos = orthologous_pairs('../CREM_CLA_protein_divergence/crem_cla_orthologs.txt')

# save UTRs of orthologous pairs to separate files
for gene in crm_cla_orthos:
    # check if remanei gene and its latens ortholog both have UTR
    if gene in crm_UTR and crm_cla_orthos[gene] in cla_UTR:
        # get the UTR sequences
        crm_seq = crm_UTR[gene]
        cla_seq = cla_UTR[crm_cla_orthos[gene]]
        # open file for writing
        newfile = open('./Crm_Cla_UTR_sequences/' + gene + '_UTR' + '.fas', 'w')
        # write remanei sequence header
        newfile.write('>' + gene + '\n')
        # write remanei sequence
        newfile.write(crm_seq + '\n')
        # write latens sequence header
        newfile.write('>' + crm_cla_orthos[gene] + '\n')
        # write latens sequence
        newfile.write(cla_seq + '\n')
        # close file after writing
        newfile.close()
    
# align remanei and latens UTR sequences
run_tcoffee_noncoding('./Crm_Cla_UTR_sequences/')
# t-coffee outputs are saved in the current directory so move outputfiles
os.system('mv CRE_PX356*UTR* ./Crm_Cla_UTR_sequences/')
# convert the t-coffee output to fasta format and save in text files
generate_fasta_from_tcoffee('./Crm_Cla_UTR_sequences/')

# generate TargetScan input sequence file with remanei and latens aligned sequences
make_cremanei_clatens_targetscan_input_sequences('./Crm_Cla_UTR_sequences/', 'Crm_Cla_UTR_seq_targetscan.txt')
print('remanei - latens targetscan input done')



# generate the Targetscan input file for the remanei and elegans orthologs
# use only transcripts that are in the remanei targetscan input file

# create a folder to store the ortholog pairs between remanei and elegans
os.mkdir('./Crm_Celegans_UTR_sequences/')

# create a dict of elegans UTR sequences
elegans_UTR = celegans_UTR_sequences('c_elegans.PRJNA13758.WS248.annotations.gff3', 'Celegans_WS248_genome.txt')

# create a dict of orthologous pairs {crem_TS: cel_TS}
crm_cel_orthos = orthologous_pairs('crm_cel_orthologs.txt')

# create a set of remanei transcripts that are used to predict target sites
valid_transcripts = set()
# open remanei targetscn input file
infile = open('Crm_UTR_seq_targetscan.txt', 'r')
# loop over file
for line in infile:
    line = line.rstrip()
    if line != '':
        line = line.split()
        valid_transcripts.add(line[0])
infile.close()

# save UTRs of orthologous pairs to separate files
# go over remanei transcripts that have an elegans ortholog
for gene in crm_cel_orthos:
    # check that gene is valid
    if gene in valid_transcripts:
        # check that gene and its elegans ortholog have both UTR sequences
        if gene in crm_UTR and crm_cel_orthos[gene] in elegans_UTR:
            # get the UTR sequences
            crm_seq = crm_UTR[gene]
            cel_seq = elegans_UTR[crm_cel_orthos[gene]]
            # open file for writing
            newfile = open('./Crm_Celegans_UTR_sequences/' + gene + '_UTR' + '.fas', 'w')
            # wwrite remanei header
            newfile.write('>' + gene + '\n')
            # write remanei sequence
            newfile.write(crm_UTR[gene] + '\n')
            # write elegans header
            newfile.write('>' + crm_cel_orthos[gene] + '\n')
            # write elegans sequence
            newfile.write(elegans_UTR[crm_cel_orthos[gene]] + '\n')
            # close file after writing
            newfile.close()
            

# align remanei and latens UTR sequences
run_tcoffee_noncoding('./Crm_Celegans_UTR_sequences/')
# t-coffee outputs are saved in the current directory so move outputfiles
os.system('mv CRE_PX356*UTR* ./Crm_Celegans_UTR_sequences/')
# convert the t-coffee output to fasta format and save in text files
generate_fasta_from_tcoffee('./Crm_Celegans_UTR_sequences/')

# generate TargetScan input sequence file with remanei and latens aligned sequences
make_cremanei_celegans_targetscan_input_sequences('./Crm_Celegans_UTR_sequences/', 'Crm_Cel_UTR_seq_targetscan.txt')
print('remanei - elegans targetscan input sequences done')





