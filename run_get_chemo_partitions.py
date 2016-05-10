# -*- coding: utf-8 -*-
"""
Created on Sun Jun 21 11:17:26 2015

@author: Richard
"""

from accessories import *
from chemoreceptors import *
import os


# get the set of chemoreceptors from the iprscan outputfile
chemo = get_chemoreceptors('./PX356_protein_seq.tsv') 

# create a set of valid transcripts
transcripts = get_valid_transcripts('../CREM_CLA_protein_divergence/unique_transcripts.txt')

# create a set of valid chemoreceptors
GPCRs = set(gene for gene in chemo if gene in transcripts)

# create a dict with the remanei CDS
CDS = convert_fasta('../CREM_CLA_protein_divergence/noamb_PX356_all_CDS.fasta')

# create a directory to store the individual fasta files
os.mkdir('Chemo_genes')

# loop over genes in CDS
for gene in CDS:
    # check if gene if chemo and valid transcript
    if gene in GPCRs:
        # check that protein sequence doesn't have stop codon
        protein = cds_translate(CDS[gene])
        if '*' not in protein:
            # open file for writing
            newfile = open('Chemo_genes/' + gene + '_prot.fasta', 'w')
            # save protein to file
            newfile.write('>' + gene + '\n')
            newfile.write(protein + '\n')
            newfile.close()

# create list of fasta files
files = os.listdir('./Chemo_genes/')

# loop over files
for filename in files:
    #check that filename is fasta
    if '_prot.fasta' in filename:
        # get gene name
        gene = filename[:filename.index('_prot')]
        # build command (phobius -png proba.png -gp proba.png -plp proba.txt input_file)
        os.system('phobius' + ' -png ' + './Chemo_genes/' + gene + '.png' + ' -gp ' + './Chemo_genes/' + gene + '.png' + ' -plp ' + './Chemo_genes/' + gene + '_proba.txt' + ' ./Chemo_genes/' + filename) 
        
# create directories to store the aligned partitions
os.mkdir('Partitions')
os.mkdir('./Partitions/Membrane/')
os.mkdir('./Partitions/Extra_membrane/')
os.mkdir('./Partitions/Inside/')
os.mkdir('./Partitions/Outside/')

# make a list of files in alignment directory
ali_files = os.listdir('../CREM_CLA_protein_divergence/pairs/Aligned_pairs/')

# make a list of alignment files
alignments = [filename for filename in ali_files if filename[-8:] == '_aln.tfa']


# loop over genes in GPCRs
for gene in GPCRs:
    # loop over alignment
    for filename in alignments:
        # check that gene in filename
        if gene == filename[:filename.index('_CLA')]:
            # get the aligned codons 
            codons = get_aligned_codons(filename, '../CREM_CLA_protein_divergence/pairs/Aligned_pairs/')
            # get the dict of probabilities
            probabilities = parse_phobius_output(gene + '_proba.txt', './Chemo_genes/')
            # get a list of codon index
            codon_index = [i for i in codons]
            # sort list
            codon_index.sort()
            # get the list of amino acid index
            aa_index = [i for i in probabilities]
            # sort list
            aa_index.sort()
            # check that the list of index are the same
            if aa_index != codon_index:
                print(gene, codon_index, aa_index)
                raise ValueError('Codon and AA index lists are different')
            # create sequences to store the different partitions
            crm_TM, crm_intra, crm_extra, crm_not_TM, cla_TM, cla_intra, cla_extra, cla_not_TM = '', '', '', '', '', '', '', ''
            # loop over the aa_indices
            for i in aa_index:
                # check that sequences in each dict is the same
                if probabilities[i][0] != cds_translate(codons[i][0]):
                    raise ValueError('Protein sequences in ortholog and probability dicts are different')
                # check probabilities and build sequences
                if probabilities[i][1] >= 0.95:
                    # build intra and not membrane sequences
                    crm_intra += codons[i][0]
                    cla_intra += codons[i][1]
                    crm_not_TM += codons[i][0]
                    cla_not_TM += codons[i][1]
                elif probabilities[i][2] >= 0.95:
                    # build outside and not membrane sequences
                    crm_extra += codons[i][0]
                    cla_extra += codons[i][1]
                    crm_not_TM += codons[i][0]
                    cla_not_TM += codons[i][1]
                elif probabilities[i][3] >= 0.95:
                    # build membrane sequences
                    crm_TM += codons[i][0]
                    cla_TM += codons[i][1]
                elif probabilities[i][4] >= 0.95:
                    # build not_membrane sequences
                    crm_not_TM += codons[i][0]
                    cla_not_TM += codons[i][1]
            # get cla_gene name
            cla_gene = filename[filename.index('CLA'):filename.index('_aln')]
            # check that remanei sequence is not empty and that latens sequence has at least 5 codons
            if len(crm_TM) != 0 and len(crm_not_TM) != 0 and len(cla_TM.replace('-', '')) >= 5 and len(cla_not_TM.replace('-', '')) >= 5:
                # gene has both membrane and extra-membrane residues
                # open file for writing
                newfile = open('./Partitions/Membrane/' + gene + '_TM.txt', 'w')
                # write alignment file in codeml input format
                newfile.write('2' + ' ' + str(len(crm_TM)) + '\n')
                newfile.write('>' + gene + '\n')
                newfile.write(crm_TM + '\n')
                newfile.write('>' + cla_gene + '\n')
                newfile.write(cla_TM + '\n')
                newfile.close()
                # open file for writing
                newfile = open('./Partitions/Extra_membrane/' + gene + '_ExtraTM.txt', 'w')
                # write alignment file in codeml input format
                newfile.write('2' + ' ' + str(len(crm_not_TM)) + '\n')
                newfile.write('>' + gene + '\n')
                newfile.write(crm_not_TM + '\n')
                newfile.write('>' + cla_gene + '\n')
                newfile.write(cla_not_TM + '\n')
                newfile.close()
            if len(crm_intra) != 0 and len(cla_intra.replace('-', '')) >= 5:
                # gene has intra-cellular domain
                # open file for writing
                newfile = open('./Partitions/Inside/' + gene + '_inside.txt', 'w')
                # write alignment file in codeml input format
                newfile.write('2' + ' ' + str(len(crm_intra)) + '\n')
                newfile.write('>' + gene + '\n')
                newfile.write(crm_intra + '\n')
                newfile.write('>' + cla_gene + '\n')
                newfile.write(cla_intra + '\n')
                newfile.close()
            if len(crm_extra) != 0 and len(cla_extra.replace('-', '')) >= 5:
                # gene has extra-cellular domain
                # open file for writing 
                newfile = open('./Partitions/Outside/' + gene + '_outside.txt', 'w')
                # write alignment file in codeml input format
                newfile.write('2' + ' ' + str(len(crm_extra)) + '\n')
                newfile.write('>' + gene + '\n')
                newfile.write(crm_extra + '\n')
                newfile.write('>' + cla_gene + '\n')
                newfile.write(cla_extra + '\n')
                newfile.close()
                
            
                
                    
                             
