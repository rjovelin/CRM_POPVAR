# -*- coding: utf-8 -*-
"""
Created on Mon Aug 17 14:07:49 2015

@author: RJovelin
"""

from manipulate_sequences import *
from sites_with_coverage import *
from miRNA_target import *
from piRNAs import *
from genomic_coordinates import *
import random
import os

# use this function to randomly sample SNPs with replacement
def SNP_MAF_randomization(SNP_positions, N_SNPs, N_replicates):
    '''
    (dict, int, int) -> dict
    Take a dict with chromo as key of inner dictionaries with SNP position as key
    and list of allele counts as value, the number of SNPs to sample, and the
    number of replicates, and return a dictionary of replicate number : list of  
    minor allele frequencies pairs
    Precondition: the SNP_positions dict has positions of polymorphic sites only
    '''
    
    # make a dictionary with int as key and chromo as value
    # {int : chromo}
    LG = {}
    i = 0
    for chromo in SNP_positions:
        LG[i] = chromo
        i += 1
    
    # make a dictionary with chromo as key of inner dictionaries with int : SNP position pairs 
    # {chromo : int : SNP positions}
    SNPs = {}
    for chromo in SNP_positions:
        # initiate integer key
        i = 0
        # create inner dict
        SNPs[chromo] = {}
        # loop over chromo, populate dict
        for position in SNP_positions[chromo]:
            SNPs[chromo][i] = position
            i += 1
    
    # create a dict to store the MAF distribution for each repliacate
    MAF = {}
    
    # get N_ SNPs for each replicate
    j = N_replicates
    while j != 0:
        # initiate list
        MAF[j] = []
        # sample N_sites at random, compute MAF and store in list
        k = N_SNPs
        while k != 0:
            # pick a random chromosome
            m = random.randint(0, len(LG)-1)
            chromo = LG[m]
            # pick a random SNP on that chromo
            n = random.randint(0, len(SNPs[chromo]) -1)
            site = SNPs[chromo][n]
            snp = SNP_positions[chromo][site]
            ref_count = snp[2]
            alt_count = snp[3]
            if alt_count <= ref_count:
                # alternative allele is minor
                # add MAF to list
                MAF[j].append(alt_count  / (ref_count + alt_count))
            elif alt_count > ref_count:
                # reference allele is minor
                # add MAF to list
                MAF[j].append(ref_count / (ref_count + alt_count))
            # update SNP counter
            k -= 1
        # update replicate number
        j -= 1
        
    return MAF
                    
            

# use this function to get indices of SNPs flanking pirnas or mirnas
def get_small_rna_flanking_SNPs(chromo_sites, rna_coord, flanking_size):
    '''
    (dict, dict, int) -> dict
    Take the dictionary with allele counts at each position with coverage (SNPs
    and non-SNPs), a dict with RNAs (miRNAs or piRNAs) coordinates on each chromo,
    the number of flanking sites to consider and return a dictionary with chromo as key of inner dictionaries of indices of SNPs only 
    for sites flanking small RNAs and exluding SNPs in small RNAs
    Precondition: All position are 0-based indices
    '''
    
    # get the position of small rnas on each chromo
    # {chromo: {set of position}}
    rna_pos = get_small_rna_sites(rna_coord)
    
    # create a dict to store the indices of flanking sites
    rna_flanking = {}

    # loop over chromo in rna_coord
    for chromo in rna_coord:
        # check if chromo in chromo_sites
        if chromo in chromo_sites:
            # loop over rna on chromo
            for i in range(len(rna_coord[chromo])):
                # get start, end position
                start = rna_coord[chromo][i][0]
                end = rna_coord[chromo][i][1]
                # get start and end position of flanking sites
                flank_start = start - flanking_size
                flank_end = end + flanking_size
                # create a list of flanking positions
                positions = [j for j in range(flank_start, start)]
                for j in range(end, flank_end):
                    positions.append(j)
                # check if chromo in rna flanking
                if chromo not in rna_flanking:
                    rna_flanking[chromo] = {}
                # add positions : allele counts to rna flanking
                for j in positions:
                    # check if site in is chromo_sites
                    if j in chromo_sites[chromo]:
                        # add only polymorphic sites
                        if chromo_sites[chromo][j][2] != 0 and chromo_sites[chromo][j][3] != 0:
                            # do not add sites falling in neighbioring small rnas
                            if j not in rna_pos[chromo]:
                                rna_flanking[chromo][j] = list(chromo_sites[chromo][j])
                                            
    # remove chromo with empty positions
    to_remove = []
    for chromo in rna_flanking:
        if len(rna_flanking[chromo]) == 0:
            to_remove.append(chromo)
    if len(to_remove) != 0:
        for chromo in to_remove:
            del rna_flanking[chromo]
            
    return rna_flanking
    
    
# use this function to get indices of SNPs in UTR
def get_UTR_SNPs(chromo_sites, caeno_gff, elegans_gff, genome_fasta, quantile, unique_transcripts, mirna_target_coord_file):
    '''
    (dict, file, file, file, int, file, file)
    Take the dictionary with allele counts at each position with coverage, 
    the remanei gff file, the elegans gff file, the remanei genome fasta file,
    the quantile of the distribution of elegans annotated UTRs, the file with 
    valid transcripts (ie. a single transcript per gene), the file with miRNA
    target coordinates and return a dictionary with chromo as key and indices
    of SNPs in UTR with their allele counts excluding miRNA target sites
    Precondition: All position are 0-based indices
    '''
    
    # compute threshold based on the distribution of elegans UTR length
    UTR_length = celegans_three_prime_UTR_length(elegans_gff)
    threshold = get_percentile(UTR_length, quantile)
    # get UTR coord {TS1 : [chromo, start, end, orientation]}
    three_prime = get_three_prime_UTR_positions(caeno_gff, genome_fasta, threshold)    
    # get a set of valid transcripts
    valid_transcripts = get_valid_transcripts(unique_transcripts)
    
    # remove UTR of non-valid transscripts
    to_remove = [gene for gene in three_prime if gene not in valid_transcripts]
    if len(to_remove) != 0:
        for gene in to_remove:
            del three_prime[gene]
    
    # create a dict UTR_pos
    UTR_pos = {}
    # loop over genes in three_prime
    for gene in three_prime:
        # get chromo
        chromo = three_prime[gene][0]
        # convert to 0-based
        start = three_prime[gene][1] -1
        end = three_prime[gene][2]
        # check if chromo in UTR_pos
        if chromo in UTR_pos:
            for j in range(start, end):
                UTR_pos[chromo].add(j)
        else:
            UTR_pos[chromo] = set()
            for j in range(start, end):
                UTR_pos[chromo].add(j)
    
    # get the positions of the miRNA targets {chromo: [[start, end, orientation]]}
    targets = get_miRNA_target_loci(mirna_target_coord_file, unique_transcripts, 'all')
    
    # create a dict with the positions of the targets on each chromo
    # {chromo: {set of positions}}
    targets_pos = {}
    for chromo in targets:
        if chromo not in targets_pos:
            targets_pos[chromo] = set()
        for i in range(len(targets[chromo])):
            start = targets[chromo][i][0]
            end = targets[chromo][i][1]
            for j in range(start, end):
                targets_pos[chromo].add(j)
    
    # create a dict with positions of SNPs in UTR
    UTR_snps = {}    
    
    # loop over chromo in UTR_pos
    for chromo in UTR_pos:
        # check if chromo in chromo_sites
        if chromo in chromo_sites:
            # check if chromo in UTR_snps
            if chromo not in UTR_snps:
                # initate inner dict
                UTR_snps[chromo] = {}
            # loop over positions on chromo
            for i in UTR_pos[chromo]:
                # check if site has coverage
                if i in chromo_sites[chromo]:
                    # do not record sites falling at miRNA targets
                    if chromo not in targets_pos:
                        # record all polymorphic sites
                        # check if site is polymorphic
                        if chromo_sites[chromo][i][2] != 0 and chromo_sites[chromo][i][3] != 0:
                            # site is polymorphic
                            UTR_snps[chromo][i] = list(chromo_sites[chromo][i])
                    elif chromo in targets_pos:
                        # do not record sites in targets
                        if i not in targets_pos[chromo]:
                            # check if site is polymorphic
                            if chromo_sites[chromo][i][2] != 0 and chromo_sites[chromo][i][3] != 0:
                                # site is polymorphic
                                UTR_snps[chromo][i] = list(chromo_sites[chromo][i])
    
    # remove chromo if no SNPs on chromo
    to_remove = [chromo for chromo in UTR_snps if len(UTR_snps[chromo]) == 0]
    if len(to_remove) != 0:
        for chromo in to_remove:
            del UTR_snps[chromo]
    
    return UTR_snps


# use this function to count the number of SNPs in each MAF bin
def SNP_proportions_MAF_bin(MAF_list):
    '''
    (list) -> dict
    Take a list of MAF frequencies and return a dictionnary with SNP proportions
    in each MAF bin in %
    '''
    
    # create a dict {MAF_lower bound : SNP_count}
    MAF = {}
    # initate empty lists
    for i in range(0, 5, 1):
        MAF[i/10] = 0
        
    # loop over MAF values
    for i in MAF_list:
        # check MAF value, update counter
        if i < 0.1:
            MAF[0.0] += 1
        elif 0.1 <= i < 0.2:
            MAF[0.1] += 1
        elif 0.2 <= i < 0.3:
            MAF[0.2] += 1
        elif 0.3 <= i < 0.4:
            MAF[0.3] += 1
        elif i >= 0.4:
            MAF[0.4] += 1
            
    # get total SNP count
    total = len(MAF_list)
    
    # loop over MAF bins
    for i in MAF:
        MAF[i] = (MAF[i] / total) * 100 
    
    return MAF    
    
    
# use this function to parition MAF vallues of each replicate into MAF bins
def get_MAF_distribution_from_replicates(MAF_replicates):
    '''
    (dict) -> dict
    Take the dictionnary with the MAF distribution of each replicate from the 
    random sampling and return a dictionary with MAF lower bound as key and 
    a list of SNP proportions (in %) in each MAF bin pooled from all replicates
    '''
    
    # create a dictionary {MAF_lower_limnit : [list of MAF values]}
    MAF = {}
    # initate empty lists
    for i in range(0, 5, 1):
        MAF[i/10] = []
       
    # loop over replicate in MAF replciates
    for i in MAF_replicates:
        # get the SNP proportions in each MAF bin
        # {MAF_lower bound : SNP_proportion}
        maf = SNP_proportions_MAF_bin(MAF_replicates[i])
        # loop over maf_lower bound and add proportion to list
        for j in maf:
            MAF[j].append(maf[j])
    
    return MAF
    
    
    
    
# use this function to get the DAF of SNPs in UTR    
def get_DAF_UTR_non_target(chromo_sites, crm_cla_target_sites_file, UTR_alignments_folder, genome_fasta):
    '''
    (dict, file, str, file) -> dict
    Take the dictionary with allele counts at sites with coverage and minimum
    sample size, the file with mirna target sites coordinates (and UTR coordinates
    in genome) for genes that have orthologs between remanei and latens,
    the folder where the remanei-latens UTR alignments are located, and fasta file
    of the remanei genome and return a dictionary with chromo as key and inner dictionaries
    with SNP positions : DAF pairs in UTR (including target sites)
    Preconditions: all positions are 0-based
    '''
    
    # convert genome fasta file to dict
    genome = convert_fasta(genome_fasta)
    
    # create a list of UTR alignment files
    files = [i for i in os.listdir(UTR_alignments_folder) if 'CRE' in i and '.txt' in i]
    
    # create a dict to store the UTR positions and associated DAF on each chromo
    # {chromo: {positions: DAF}}
    DAF = {}
    
    # open file for reading
    infile = open(crm_cla_target_sites_file, 'r')
    # skip header
    infile.readline()
    # loop over file
    for line in infile:
        line = line.rstrip()
        if line != '':
            line = line.split()
            # get transcript
            gene = line[0]
            # find the UTR alignment file
            for filename in files:
                if gene == filename[:filename.index('_UTR')]:
                    utr_file = filename
                    break
            # convert the alignment to dict
            UTR_ali = convert_fasta(UTR_alignments_folder + utr_file)
            # find the remanei and latens utr
            for seqname in UTR_ali:
                if 'CRE' in seqname:
                    crm_utr = seqname
                elif 'CLA' in seqname:
                    cla_utr = seqname
            # get UTR coord on chromo 0-based
            UTR_start = int(line[9]) - 1
            UTR_end = int(line[10])
            # get orientation
            orientation = line[15]
            # get chromo
            chromo = line[7]
            # get the positions of UTR on chromo
            positions = [i for i in range(UTR_start, UTR_end)]
            # reverse order if orientation is -
            if orientation == '-':
                positions.sort()
                positions.reverse()
            # get the remanei and latens aligned UTR sequences, without gaped positions in remanei
            # remanei UTR sequence from alignment without gaps == UTR sequence grabed from chromo
            crmutrseq = ''
            clautrseq = ''
            for i in range(len(UTR_ali[crm_utr])):
                if UTR_ali[crm_utr][i] != '-':
                    crmutrseq += UTR_ali[crm_utr][i]
                    clautrseq += UTR_ali[cla_utr][i]
            # check UTR lengths
            assert len(crmutrseq) == len(clautrseq), 'length of remanei UTR gap-free is different from latens aligned UTR'
            assert len(crmutrseq) == len(positions), 'length of remanei UTR gap-free is different from UTR length from chromo'
            # check remanei UTR sequences
            checkutrseq = ''
            if orientation == '+':
                for i in range(len(positions)):
                    checkutrseq += genome[chromo][positions[i]]
            elif orientation == '-':
                # need to take the complement 
                for i in range(len(positions)):
                    checkutrseq += seq_complement(genome[chromo][positions[i]])
            assert checkutrseq == crmutrseq, 'remanei UTR gap-free and UTR from chromo are different'
                      
            # check if chromo in chromo sites
            if chromo in chromo_sites:
                # loop oer positions
                for i in range(len(positions)):
                    # check if site has coverage or minimum sample size
                    if positions[i] in chromo_sites[chromo]:
                        # get ref and alt counts
                        ref_count = chromo_sites[chromo][positions[i]][2]
                        alt_count = chromo_sites[chromo][positions[i]][3]
                        # check orientation to get ref and alt alleles
                        if orientation == '+':
                            ref = chromo_sites[chromo][positions[i]][0]
                            alt = chromo_sites[chromo][positions[i]][1]
                        elif orientation == '-':
                            # take the complement of ref and alt
                            ref = seq_complement(chromo_sites[chromo][positions[i]][0])
                            alt = seq_complement(chromo_sites[chromo][positions[i]][1])
                        # check that SNP can be polarozed (ancestral == ref or alt)
                        if clautrseq[i] == ref or clautrseq[i] == alt:
                            # SNP can be polarized
                            # check if site is polymorphic
                            if ref_count != 0 and alt_count != 0:
                                # site is polymorphic, find derived and ancestral allele
                                # double check that ref and alt are different
                                assert ref != alt, 'ref and alt counts are different but alleles are the same'
                                if clautrseq[i] == ref:
                                    # reference allele is ancestral, alt is derived
                                    # compute DAF
                                    freq = alt_count / (ref_count + alt_count)
                                elif clautrseq[i] == alt:
                                    # alternative allele is ancestral, ref is derived
                                    # compute DAF
                                    freq = ref_count / (ref_count + alt_count)
                                # add site position : freq to dict                                
                                # check if chromo is key in DAF dict
                                if chromo in DAF:
                                    DAF[chromo][positions[i]] = freq
                                else:
                                    # initiate inner dict
                                    DAF[chromo] = {}
                                    DAF[chromo][positions[i]] = freq
                                    
    # close file after reading
    infile.close()
    
    return DAF
    
    
   
# use this function to randomly sample SNPs with replacement
def SNP_DAF_randomization(SNP_positions, N_SNPs, N_replicates):
    '''
    (dict, int, int) -> dict
    Take a dict with chromo as key of inner dictionaries with SNP position as key
    and DAF as value, the number of SNPs to sample, and the
    number of replicates, and return a dictionary of replicate number : list of  
    derived allele frequencies pairs
    Precondition: the SNP_positions dict has positions of polymorphic sites only
    '''
    
    # make a dictionary with int as key and chromo as value
    # {int : chromo}
    LG = {}
    i = 0
    for chromo in SNP_positions:
        LG[i] = chromo
        i += 1
    
    # make a dictionary with chromo as key of inner dictionaries with int : SNP position pairs 
    # {chromo : int : SNP positions}
    SNPs = {}
    for chromo in SNP_positions:
        # initiate integer key
        i = 0
        # create inner dict
        SNPs[chromo] = {}
        # loop over chromo, populate dict
        for position in SNP_positions[chromo]:
            SNPs[chromo][i] = position
            i += 1
    
    # create a dict to store the DAF distribution for each repliacate
    DAF = {}
    
    # get N_ SNPs for each replicate
    j = N_replicates
    while j != 0:
        # initiate list
        DAF[j] = []
        # sample N_sites at random, add DAF to list
        k = N_SNPs
        while k != 0:
            # pick a random chromosome
            m = random.randint(0, len(LG)-1)
            chromo = LG[m]
            # pick a random SNP on that chromo
            n = random.randint(0, len(SNPs[chromo]) -1)
            site = SNPs[chromo][n]
            # get daf at that site
            daf = SNP_positions[chromo][site]
            # add daf to list
            DAF[j].append(daf)
            # update SNP counter
            k -= 1
        # update replicate number
        j -= 1
        
    return DAF
    
    
    
# use this function to parition DAF vallues of each replicate into DAF bins
def get_DAF_distribution_from_replicates(DAF_replicates):
    '''
    (dict) -> dict
    Take the dictionnary with the DAF distribution of each replicate from the 
    random sampling and return a dictionary with DAF lower bound as key and 
    a list of SNP proportions (in %) in each DAF bin pooled from all replicates
    '''
    
    # create a dictionary {DAF_lower_limnit : [list of DAF values]}
    DAF = {}
    # initate empty lists
    for i in range(0, 10, 1):
        DAF[i/10] = []
       
    # loop over replicate in DAF replciates
    for i in DAF_replicates:
        # get the SNP proportions in each DAF bin
        # {DAF_lower bound : SNP_proportion}
        daf = SNP_proportions_DAF_bin(DAF_replicates[i])
        # loop over maf_lower bound and add proportion to list
        for j in daf:
            DAF[j].append(daf[j])
    
    return DAF
    
    
# use this function to count the number of SNPs in each DAF bin
def SNP_proportions_DAF_bin(DAF_list):
    '''
    (list) -> dict
    Take a list of DAF frequencies and return a dictionnary with SNP proportions
    in each DAF bin in %
    '''
    
    # create a dict {DAF_lower bound : SNP_count}
    DAF = {}
    # initate empty lists
    for i in range(0, 10, 1):
        DAF[i/10] = 0
        
    # loop over DAF values
    for i in DAF_list:
        # check DAF value, update counter
        if i < 0.1:
            DAF[0.0] += 1
        elif 0.1 <= i < 0.2:
            DAF[0.1] += 1
        elif 0.2 <= i < 0.3:
            DAF[0.2] += 1
        elif 0.3 <= i < 0.4:
            DAF[0.3] += 1
        elif 0.4 <= i < 0.5:
            DAF[0.4] += 1
        elif 0.5 <= i < 0.6:
            DAF[0.5] += 1
        elif 0.6 <= i < 0.7:
            DAF[0.6] += 1
        elif 0.7 <= i < 0.8:
            DAF[0.7] += 1
        elif 0.8 <= i < 0.9:
            DAF[0.8] += 1
        elif i >= 0.9:
            DAF[0.9] += 1
            
    # get total SNP count
    total = len(DAF_list)
    
    # loop over DAF bins
    for i in DAF:
        DAF[i] = (DAF[i] / total) * 100 
    
    return DAF    