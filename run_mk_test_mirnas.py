# -*- coding: utf-8 -*-
"""
Created on Wed Jun 29 16:07:21 2016

@author: RJovelin
"""

# use this script to compute the MK test for miRNAs

import os
from manipulate_sequences import *
from divergence import *
from miRNA_target import *
from sites_with_coverage import *


# convert genome to dict
CrmGenome = convert_fasta('../Genome_Files/noamb_356_v1_4.txt')
print('converted genome to dict')

# create a set of valid transcripts (1 transcript mapped to 1 gene)
valid_transcripts = get_valid_transcripts('../Genome_Files/unique_transcripts.txt')
print('made a set of valid transcripts')

# get miRNA coordinates {chromo: [[start, end, orientation]]}
mirna_coord = get_mirna_loci('CRM_miRNAsCoordinatesFinal.txt')
print('got miRNA coordinates')

# get mature coordinates {chromo: [[start, end orientation]]}
mature_coord = get_mirna_loci('CRM_MatureCoordinatesFinal.txt')
print('got mature miR coordinates')

# create a dict with remanei mirna hairpin name and remanei-latens aligned hairpins
hairpins = {}
# get the hairpin alignment files
hairpinfiles = [i for i in os.listdir('Crm_Cla_miRNA_orthos') if 'hairpin' in i and '.txt' in i]
for filename in hairpinfiles:
    # extract crm mir name from file
    crmname = filename[:filename.index('_')]
    # convert file to dict
    fastafile = convert_fasta('Crm_Cla_miRNA_orthos' + '/' + filename)
    # initialize dict
    hairpins[crmname] = {}
    for name in fastafile:
        hairpins[crmname][name] = fastafile[name]
print('got aligned hairpins')
print('aligned hairpins', len(hairpins))


# create a dict with remanei mirna mature name and remanei-latens aligned mature
matures = {}
# get the mature alignment files
maturefiles = [i for i in os.listdir('Crm_Cla_miRNA_orthos') if 'mature' in i and '.txt' in i]
for filename in maturefiles:
    # extract crm mature name from file
    crmname = filename[:filename.index('_mature')]
    # convert file to dict
    fastafile = convert_fasta('Crm_Cla_miRNA_orthos' + '/' + filename)
    # initialize dict
    matures[crmname] = {}
    for name in fastafile:
        matures[crmname][name] = fastafile[name]
print('got aligned matures')
print('aligned matures', len(matures))


# get miRNA coordinates {chromo: [[start, end, orientation]]}
mirna_coord = get_mirna_loci('CRM_miRNAsCoordinatesFinal.txt')
print('got miRNA coordinates')


# get mature coordinates {chromo: [[start, end orientation]]}
mature_coord = get_mirna_loci('CRM_MatureCoordinatesFinal.txt')
print('got mature miR coordinates')


# make a dict with family level conservation for all miRNAs {name : conservarion}
famCons = {}
infile = open('CRM_miRNAsCoordinatesFinal.txt')
infile.readline()
for line in infile:
    line = line.rstrip()
    if line != '':
        line = line.split('\t')
        name, conservation = line[0], line[-1]
        famCons[name] = conservation
infile.close()
print('got mirna conservation level')

# create a dict with coordinates of mature sequences
miR_coord = {}
infile = open('CRM_MatureCoordinatesFinal.txt')
infile.readline()
for line in infile:
    line = line.rstrip()
    if line != '':
        line = line.split('\t')
        name, chromo, start, end, orientation = line[0], line[1], int(line[2]) -1, int(line[3]), line[4]
        if orientation == '+':
            assert CrmGenome[chromo][start : end] == line[-1], 'mirna does not match mature on +'
        elif orientation == '-':
            assert reverse_complement(CrmGenome[chromo][start : end]) == line[-1], 'mirna does not match mature on -'
        miR_coord[name] = [chromo, start, end, orientation]
infile.close()
print('got mature coordinates')


# create a dict with coordinates of hairpin sequences
hairpin_coord = {}
infile = open('CRM_miRNAsCoordinatesFinal.txt')
infile.readline()
for line in infile:
    line = line.rstrip()
    if line != '':
        line = line.split('\t')
        name, chromo, start, end, orientation = line[0], line[1], int(line[2]) - 1, int(line[3]), line[4]
        # remove arm from name
        name = name[: name.index('_')]
        if orientation == '+':
            assert CrmGenome[chromo][start : end] == line[5], 'mirna does not match hairpin on +'
        elif orientation == '-':
            assert reverse_complement(CrmGenome[chromo][start : end]) == line[5], 'mirna does not match hairpin on -'
        hairpin_coord[name] = [chromo, start, end, orientation]
infile.close()
print('got hairpin coordinates')


# get the allele counts for all sites with coverage
# {chromo: {site : [ref_allele, alt_allele, count_ref, count alt]}]}
chromo_sites = get_non_coding_snps('../SNP_files/', 10)
print('got allele counts at all sites')

# check that mirna names in aligned sequences dict are in coordinates dict
for mirna in hairpins:
    assert mirna in hairpin_coord, 'mirna with aligned sequences does not have coordinates'
for mirna in matures:
    assert mirna in miR_coord, 'miR with aligned sequences does not have coordinates'        


# create a dict to record the number of fixed diffs and polmorphisms for each mirna hairpin and mature
# {mirna: [D, P]}
hairpin_diffs, mature_diffs = {}, {}

# loop over aligned hairpins
for mirna in hairpins:
    # get chromo, start, end and orientation
    chromo, start, end, orientation = hairpin_coord[mirna][0], hairpin_coord[mirna][1], hairpin_coord[mirna][2], hairpin_coord[mirna][3] 
    # extract sequence from genome
    sequence = CrmGenome[chromo][start: end]
    if orientation == '-':
        sequence = reverse_complement(sequence)
    # get mirna positions on chromo
    positions = [i for i in range(start, end)]    
    # get the position in decreasing order if orientation is -    
    if orientation == '-':
        positions.reverse()
    # get the remanei and latens mirna name, and correspsonding sequences
    for name in hairpins[mirna]:
        if name.startswith('crm'):
            crmmirna = name
            assert crmirna == mirna, 'mirna names do not match'
            crmseq = hairpins[mirna][crmmirna]
        elif name.startswith('cla'):
            clamirna = name
            claseq = hairpins[mirna][clamirna]
    # set up a gap counter
    gaps = 0
    # loop over crm sequence
    # keep track of index to look for SNP data when gaps are present
    for i in range(len(crmseq)):
        # get ancestral allele and remanei allele
        ancestral, crmallele = claseq[i], crmseq[i]
        if crmallele != '-':
            # get the index to look in positions
            j = i - gaps
            # check that nucleotide in crmseq correspond to sequence from genome
            assert sequence[j] == crmallele, 'nucleotides do not match between extracted sequence and aligned sequence'
            # check that index in positions is correct
            # check that ref allele in SNP dictionary is correct
            if orientation == '+':
                assert CrmGenome[chromo][positions[j]] == crmallele, 'no match with nucleotide extracted with list index on +'
                if positions[j] in chromo_sites[chromo]:
                    assert chromo_sites[chromo][positions[j]][0] == crmallele, 'no match with ref allele in SNP dict in +'
                else:
                    print(j, positions[j], orientation)
            elif orientation == '-':
                assert seq_complement(CrmGenome[chromo][positions[j]]) == crmallele, 'no match with nucleotide extracted with list index on -'
                if positions[j] in chromo_sites[chromo]:
                    assert seq_complement(chromo_sites[chromo][positions[j]][0]) == crmallele, 'no match with ref allele in SNP dict in -'
                else:
                    print(j, positions[j], orientation)
        else:
            gaps += 1
        # determine if site a fized difference or a polymorphism
        # check that ancestral allele is valid base
        if ancestral in 'ATCG' and crmallele in 'ATCG':
            # check that site has coverage            
            if positions[j] in chromo_sites[chromo]:
                # fixed diff if ref_count != 0 and alt_count = 0 and ref != ancestral 
                # fixed diff is ref_count = 0 and alt_count != 0 and alt != ancestral
                # polymorphism if (ref_count != 0 and alt_count != 0) and (ref = amcestral or alt = ancestral)
                # get ref and alt counts
                ref_count, alt_count = chromo_sites[chromo][positions[j]][2], chromo_sites[chromo][positions[j]][3]
                # get reference and alternative alleles
                ref, alt = chromo_sites[chromo][positions[j]][0], chromo_sites[chromo][positions[j]][1]
                
                if ref_count != 0 and alt_count == 0 and ref != ancestral:
                    # fixed difference, populate dict
                    if mirna in hairpin_diffs:
                        hairpin_diffs[mirna][0] += 1
                    else:
                        hairpin_diffs[mirna] = [0, 0]
                elif ref_count == 0 and alt_count != 0 and alt != ancestral:
                    # fixed difference, populate dict
                    if mirna in hairpin_diffs:
                        hairpin_diffs[mirna][0] += 1
                    else:
                        hairpin_diffs[mirna] = [0, 0]
                elif (ref_count != 0 and alt_count != 0) and (ref == ancestral or alt == ancestral):
                    # polymorphism, populate dict
                    if mirna in hairpin_diffs:
                        hairpin_diffs[mirna][1] += 1
                    else:
                        hairpin_diffs[mirna] = [0, 0]
            





# add option to consider singletons or not


# add code to consider only positions with sample size > 10


# make it a function so it can be use for mature mir as well


# consider only 4-fold degenerate sites in entire genome for neutral control


# group mirnas based on level of conservation
# perform MK test based on each mirna individually
# perform MK test based on entire conservation group by pooling numbers
# see Lyu et al for organizing table
# compute alpha for all hairpins, all matures, and each conservation group



# site type site number D P PDAF.5% D/PDAF.5% MK test pvaluea ab (% of adaptive fixations)


# need a script to count P and D for non-coding sites



# extract from Lyu et al Plos genetics 2014
# New MicroRNAs in Drosophila—Birth, Death and Cycles of Adaptive Evolution

'''
We used the McDonald-Kreitman test (MK test) [32] framework
to detect positive selection in miRNAs from each age group
based on the polymorphisms within D. melanogaster and the
divergence between D. melanogaster and D. simulans. Precursor or
mature sequences of each miRNA group were combined and
treated as the functional category, while the 4-fold degenerate sites
in the whole genome were used as the neutral control. The
divergence is calculated by counting the number of changed
nucleotide sites between D. melanogaster (dm3) and D. simulans
(droSim1) based on the UCSC whole genome alignment.
Polymorphism data was retrieved from Drosophila Population
Genomics Project (DPGP, http://www.dpgp.org/, release 1.0).

SNPs that were detected on more than thirty individuals and
exhibited a derived allele frequency (DAF) > 5% were used for the
MK test.
The proportion of adaptively fixed mutations (a) was estimated
as previously described [75]. To estimate the evolutionary fate of
each miRNA, we first screened for adaptive miRNAs among the
238 candidates by using each miRNA’s precursor together with
the 50 bp flanking sequences on both sides as the functional sites.
The p-values of multiple MK tests were adjusted by the BenjaminiHochberg
method [76] and the adaptive significance of each
candidate is re-validated by using the precursor alone in the MK
test. We then identified the conservative miRNAs by comparing
the number of substitutions in the miRNA precursors (KmiR) with
the number of substitutions in the synonymous sites (KS) between
D.melanogaster and D.simulans. miRNAs with KmiR/KS < 0.5 were
considered to be conservatively evolving. Kimura’s 2-parameter
model [72] and the Nei-Gojobori model [77] were used to
calculate KmiR and KS, respectively. Finally, excluding the
adaptive and conservative miRNAs, the remaining were considered
to be in transition between adaptive to conservative/death
'''





# use this function to count the number of replacement and synonymous changes
def count_polym_diverg(snp_file, strains, rare_sites, cutoff, raw_count):
    '''
    (file, str, str, float, int) -> dict 
    Take the file with the SNPs in coding sequences, a given focal group of strains
    and return a dictionnary with gene as key and a list of containing PN, PS,
    DN, DS as value. If rare_sites is freq then polymorphic sites with a
    frequency < cutoff are ignored, and if rare_sites is count then sites that
    are < than raw_count in the sample are ignored
    '''
    
    # PN: number of replacement polymorphisms
    # PS: number of synonymous polymorphisms
    # DN: number of fixed replacements
    # DS: number of fixed synonymous changes    
    
    # create a set of stop codons
    stop_codons = {'TAG', 'TGA', 'TAA'} 
    
    # initialize dict
    SNPs = {}
    # initialize gene : empty list pairs
    # open file for reading
    infile = open(snp_file, 'r')
    # skip header
    infile.readline()
    # loop over file, get the gene name as key and initialize list
    # list contains PN, PS, DN, DS
    for line in infile:
        if line.rstrip() != '':
            line = line.rstrip().split()
            gene = line[2]
            if gene not in SNPs:
                # {gene : [PN, PS, DN, DS]}
                SNPs[gene] = [0, 0, 0, 0]
    # close file after reading
    infile.close()
    
    # open file for reading
    infile = open(snp_file, 'r')
    # skip header
    infile.readline()
    # loop over file
    for line in infile:
        line = line.rstrip()
        if line != '':
            line = line.split()
            # get dict key
            gene = line[2]
            # record only sites for which ancestral state is defined
            if line[19] in {'A', 'C', 'T', 'G'}:
                # record only sites labeled snp or no_snp and valid snp type
                if line[6] in {'snp', 'no_snp'} and line[9] in {'NA', 'COD', 'SYN', 'REP'}:
                    # get reference and alternative codons and latens codon
                    ref_codon, alt_codon, cla_codon = line[3], line[8], line[18]
                    # get reference, alternative and latens alleles
                    ref, alt, cla_base = line[5], line[7], line[19]
                    # get SNP type
                    snp = line[9]
                    # do not consider stop codons
                    if ref_codon not in stop_codons and alt_codon not in stop_codons and cla_codon not in stop_codons:
                        # do not consider codons with 'Ns'
                        if 'N' not in ref_codon and 'N' not in alt_codon and 'N' not in cla_codon:
                            # check focal group, consider only KSR strains or KRS and PX combined
                            if strains == 'KSR':
                                # get reference allele count
                                ref_count = int(line[10])
                                # get alternative allele count
                                alt_count = int(line[14])
                            elif strains == 'KSR_PX':
                                # get reference allele count
                                ref_count = int(line[13])
                                # get alternative allele count
                                alt_count = int(line[17])
                            # Consider sites that have a sample size >= 10
                            if ref_count + alt_count >= 10:
                                # do not consider codons with more than 1 substitutions
                                if diff_codon(cla_codon, alt_codon) <= 1:
                                    # determine if site is polymorphic or fixed
                                    # fixed diff if alternative allele fixed and different from latens 
                                    if ref_count == 0 and alt_count != 0 and cla_base != alt:
                                        # fixed difference between latens and remanei
                                        # check if change is synonymous or nonsynonymous
                                        if genetic_code[cla_codon] != genetic_code[alt_codon]:
                                            # fixed nonsynonymous change
                                            SNPs[gene][2] += 1
                                        elif genetic_code[cla_codon] == genetic_code[alt_codon]:
                                            # fixed synonymous change
                                            SNPs[gene][3] += 1
                                    elif ref_count != 0 and alt_count == 0 and cla_base != ref:
                                        # fixed difference between latens and remanei
                                        # do not consider codons with more than 1 substitutions
                                        # check if change is synonymous or nonsynonymous
                                        if genetic_code[cla_codon] != genetic_code[ref_codon]:
                                            # fixed nonsynonymous change
                                            SNPs[gene][2] += 1
                                        elif genetic_code[cla_codon] == genetic_code[ref_codon]:
                                            # fixed synonymous change
                                            SNPs[gene][3] += 1
                                    elif ref_count != 0 and alt_count != 0 and (ref == cla_base or alt == cla_base):
                                        # site is polymorphic
                                        # set up boolean to identify polymorphic site after filtering based on MAF or raw count
                                        PolymorphicSite = False
                                        # check if cutoff or raw_count applies
                                        if rare_sites == 'freq':
                                            # use frequency cutoff
                                            # compare the snp MAF frequency to cutoff
                                            if ref_count >= alt_count:
                                                freq = alt_count / (ref_count + alt_count)
                                            elif ref_count < alt_count:
                                                freq = ref_count / (ref_count + alt_count)
                                            if freq >= cutoff:
                                                # record polymorphic site
                                                PolymorphicSite = True
                                        elif rare_sites == 'count':
                                            # use allele count to filter sites
                                            # check that allele with lowest count > raw_count threshold
                                            if ref_count >= alt_count and alt_count > raw_count:
                                                # alt_count is greater than minimum required threshold
                                                # record polymorphic site
                                                PolymorphicSite = True
                                            elif ref_count < alt_count and ref_count > raw_count:
                                                # ref_count is greater than minimum required thereshold
                                                # record polymorphic site
                                                PolymorphicSite = True
                                        # determine the type of mutation if polymorphic site is to be recorded
                                        if PolymorphicSite == True:
                                            if snp == 'SYN' and genetic_code[ref_codon] == genetic_code[alt_codon]:
                                                # mutation is synonymous
                                                SNPs[gene][1] += 1
                                            elif snp == 'REP' and genetic_code[ref_codon] != genetic_code[alt_codon]:
                                                # mutation is nonsynonymous
                                                SNPs[gene][0] += 1
                                                    
    # close file after readling
    infile.close()
    
    # remove genes with no polymorhisms or no divergence
    # because MK test cannot be computed for these genes
    
    to_remove = []
    for gene in SNPs:
        if SNPs[gene][0] == 0 and SNPs[gene][1] == 0:
            to_remove.append(gene)
        elif SNPs[gene][2] == 0 and SNPs[gene][3] == 0:
            to_remove.append(gene)
        
    for gene in to_remove:
      del SNPs[gene]  
    
    # return dict
    return SNPs