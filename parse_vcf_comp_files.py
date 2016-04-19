# -*- coding: utf-8 -*-
"""
Created on Tue Jul 14 13:32:33 2015

@author: RJovelin
"""

import os



# use this function to check that all chromos are recorded
def check_chromo(directory, gff_file):
    '''
    (str, file) -> (set, set)
    Take the directory containing the vcf files and the remanei GFF annotation file
    and return a tuple with sets of chromos extracted respectively from the GFF
    file and from the vcf_comp files    
    '''
    
    # make a list of files in directory
    files = [filename for filename in os.listdir(directory) if 'ksrpx' in filename and 'filtered.comp' in filename]
    # extract the chromo from filenames
    chromo_comp = set()
    for filename in files:
        chromo = filename[filename.index('.')+1: filename.index('.filtered')]
        chromo_comp.add(chromo)
        
    # make a set of chromos from the GFF file
    chromo_gff = set()
    # open file for reading
    infile = open(gff_file, 'r')
    # loop over the file
    for line in infile:
        # chromo sequence are at the end of the file, skip fasta sequences
        if not line.startswith('>'):
            if 'scaffold' in line or 'linkage_group' in line:
                line = line.rstrip().split()
                chromo = line[0]
                chromo_gff.add(chromo)
    
    return chromo_gff, chromo_comp
        
# use this function to make a dict from the vcf comp file 
def from_vcf_comp_to_dict(vcf_comp_file, repeat_positions, remove_repeats):
    '''
    (file, dict, bool) -> (str, dict)
    Take the filtered vcf file with SNP info for KSR and PX strains,
    a dictionary with chromo: set of indices for all repeats on that chromo, 
    a boolean (True or False) specifying if SNPs in repeat positions should be
    removed or not, return a tuple with the chromosome and a dictionnary with
    positions as key and a list of allele counts as value
    Precondition: the file contains SNP info for KSR and PX strains only
    '''
    
    # open file for reading
    infile = open(vcf_comp_file, 'r')
    
    # create dict to store the SNP counts {pos: [ref_PX, ref_KSR, alt_PX, alt_KSR]}
    snps = {}
    
    # create a set of valid nucleotides
    valid_bases = {'A', 'C', 'G', 'T', 'a', 't', 'c', 'g'}
    
    # loop over file
    for line in infile:
        if line.rstrip() != '':
            line = line.rstrip().split()
            # get chromosome
            chromo = line[0]
            # get position 0-based index
            position = int(line[1]) - 1
            # parse strains
            strains = line[2].split('-')
            # make a set of alleles
            alleles = set()
            # get the reference allele
            ref_allele = strains[0].split(';')[1].upper()
            # get the alternative allele
            for i in range(len(strains)):
                alt_allele = strains[i].split(';')[2].upper()
                if alt_allele != '.':
                    break
            # set up counter
            px_ref = 0
            px_alt = 0
            ksr_ref = 0
            ksr_alt = 0
            
            # check that reference allele is valid nucleotide
            if ref_allele in valid_bases:
                # loop over strains
                for i in range(len(strains)):
                    # parse each string
                    individual = strains[i].split(';')
                    if 'PX' in individual[0]:
                        # strain is from PX, check if allele is alternative or reference
                        if individual[2] == '.':
                            # reference allele, update counter
                            px_ref += 1
                            # add reference allele to set of alleles
                            alleles.add(individual[1].upper())
                        elif individual[2] != '.':
                            # alternative allele, update counter
                            px_alt += 1
                            # add alternative allele to set of alleles
                            alleles.add(individual[2].upper())
                    else:
                        # strain is from KSR, check if allele is alternative or reference
                        if individual[2] == '.':
                            # reference allele, update counter
                            ksr_ref += 1
                            # add reference allele to set of alleles
                            alleles.add(individual[1].upper())
                        elif individual[2] != '.':
                            # alternative allele, update counter
                            ksr_alt += 1
                            # add alternative allele to set of alleles
                            alleles.add(individual[2].upper())
                # check if number of alleles is at most 2 and only contains valid nucleotides (no N)
                if len(alleles) <= 2 and alleles.issubset(valid_bases): 
                    # check if reference allele is different than alternative allele
                    if alt_allele == '.':
                        alt_allele = ref_allele
                    # populate dict {pos: [ref_allele, alt_allele, ref_PX, ref_KSR, alt_PX, alt_KSR]}
                    snps[position] = [ref_allele, alt_allele, px_ref, px_alt, ksr_ref, ksr_alt]
    
    # check if remove_repeats == True   
    if remove_repeats == True:
        # check that SNPs are on chromo and that repeats are on chromo
        if len(snps) != 0 and chromo in repeat_positions:
            # remove SNPs falling in repeats
            for i in repeat_positions[chromo]:
                # check if positions in snp dict
                if i in snps:
                    # remove position
                    del snps[i]
    
    # check that is not empty (some files are empty (eg 'ksrpx.scaffold_3050.filtered.comp))
    if len(snps) != 0:
        print(chromo, len(snps))
        return chromo, snps
    else:
        print('', 0)
        return('', 0)
                
                    
# use this function to save the SNP info for each vcf_comp_file into a separate outputfile
def generate_snp_file(chromo, snps, outputfile):
    '''
    (str, dict, file) -> file
    Take a given chromo and a dictionnary with SNP counts in PX and KSR strains 
    for each position on chromo and save the SNP information to outputfile    
    '''
    
    # open outputfile for writing
    newfile = open(outputfile, 'w')
    
    header = '\t'.join(['chromo', 'pos', 'ref', 'alt', 'flag_snp', 'px_ref', 'ksr_ref', 'px+ksr_ref', 'px_alt', 'ksr_alt', 'px+ksr_alt'])
    # write header to file
    newfile.write(header + '\n')    
    
    # make a list of positions in dict
    positions = [j for j in snps]
    # sort list
    positions.sort()
    print('positions sorted')
    print(chromo, len(positions))
    
    # loop over the sorted positions
    for i in positions:
        # get ref and alt alleles
        ref_allele = snps[i][0]
        alt_allele = snps[i][1]
        # get 1-based position
        pos = i + 1
        px_ref = snps[i][2]
        px_alt = snps[i][3]
        ksr_ref = snps[i][4]
        ksr_alt = snps[i][5]
        pxksr_ref = px_ref + ksr_ref
        pxksr_alt = px_alt + ksr_alt
        if ref_allele == alt_allele:
            SNP = 'no_snp'
        elif ref_allele != alt_allele:
            SNP = 'snp'
                
        # write info to file
        line = '\t'.join([chromo, str(pos), ref_allele, alt_allele, SNP, str(px_ref), str(ksr_ref), str(pxksr_ref), str(px_alt), str(ksr_alt), str(pxksr_alt)])
        newfile.write(line + '\n')
    
    # close file after reading    
    newfile.close()
    
     
# use this function to make a dict from the vcf comp file that include PB, KSR and PX strains 
def parse_vcf_comp_PB_ON(vcf_comp_file, repeat_positions, remove_repeats):
    '''
    (file, dict, bool) -> (str, dict)
    Take the filtered vcf file with SNP info for KSR, PX strains aned PB strains,
    a dictionary with chromo: set of indices for all repeats on that chromo, 
    a boolean (True or False) specifying if SNPs in repeat positions should be
    removed or not, return a tuple with the chromosome and a dictionnary with
    positions as key and a list of allele counts as value
    '''
    
    # open file for reading
    infile = open(vcf_comp_file, 'r')
    
    # create dict to store the SNP counts {pos: [ref, alt, ref_PX, alt_PX, ref_KSR, alt_KSR, ref_PB, alt_PB]}  
    snps = {}
    
    # create a set of valid nucleotides
    valid_bases = {'A', 'C', 'G', 'T', 'a', 't', 'c', 'g'}
    
    # loop over file
    for line in infile:
        if line.rstrip() != '':
            line = line.rstrip().split()
            # get chromosome
            chromo = line[0]
            # get position 0-based index
            position = int(line[1]) - 1
            # parse strains
            strains = line[2].split('-')
            # make a set of alleles
            alleles = set()
            # get the reference allele
            ref_allele = strains[0].split(';')[1].upper()
            # get the alternative allele
            for i in range(len(strains)):
                alt_allele = strains[i].split(';')[2].upper()
                if alt_allele != '.':
                    break
            # set up counter
            px_ref = 0
            px_alt = 0
            ksr_ref = 0
            ksr_alt = 0
            pb_ref = 0
            pb_alt = 0
            
            # check that reference allele is valid nucleotide
            if ref_allele in valid_bases:
                # loop over strains
                for i in range(len(strains)):
                    # parse each string
                    individual = strains[i].split(';')
                    if 'PX' in individual[0]:
                        # strain is from PX, check if allele is alternative or reference
                        if individual[2] == '.':
                            # reference allele, update counter
                            px_ref += 1
                            # add reference allele to set of alleles
                            alleles.add(individual[1].upper())
                        elif individual[2] != '.':
                            # alternative allele, update counter
                            px_alt += 1
                            # add alternative allele to set of alleles
                            alleles.add(individual[2].upper())
                    elif 'PB' in individual[0]:
                        # strain is from PB, check if allele is alternative or reference
                        if individual[2] == '.':
                            # reference allele, update counter
                            pb_ref += 1
                            # add reference allele to set of alleles
                            alleles.add(individual[1].upper())
                        elif individual[2] != '.':
                            # alternative allele, update counter
                            pb_alt += 1
                            # add alternative allele to set of alleles
                            alleles.add(individual[2].upper())
                    else:
                        # strain is from KSR, check if allele is alternative or reference
                        if individual[2] == '.':
                            # reference allele, update counter
                            ksr_ref += 1
                            # add reference allele to set of alleles
                            alleles.add(individual[1].upper())
                        elif individual[2] != '.':
                            # alternative allele, update counter
                            ksr_alt += 1
                            # add alternative allele to set of alleles
                            alleles.add(individual[2].upper())
                # check if number of alleles is at most 2 and only contains valid nucleotides (no N)
                if len(alleles) <= 2 and alleles.issubset(valid_bases): 
                    # check if reference allele is different than alternative allele
                    if alt_allele == '.':
                        alt_allele = ref_allele
                    # populate dict {pos: [ref, alt, ref_PX, alt_PX, ref_KSR, alt_KSR, ref_PB, alt_PB]} 
                    snps[position] = [ref_allele, alt_allele, px_ref, px_alt, ksr_ref, ksr_alt, pb_ref, pb_alt]
    
    # check if remove_repeats == True   
    if remove_repeats == True:
        # check that SNPs are on chromo and that repeats are on chromo
        if len(snps) != 0 and chromo in repeat_positions:
            # remove SNPs falling in repeats
            for i in repeat_positions[chromo]:
                # check if positions in snp dict
                if i in snps:
                    # remove position
                    del snps[i]
    
    # check that is not empty (some files are empty (eg 'ksrpx.scaffold_3050.filtered.comp))
    if len(snps) != 0:
        print(chromo, len(snps))
        return chromo, snps
    else:
        print('', 0)
        return('', 0)
                
                    
# use this function to save the SNP info for each vcf_comp_file into a separate outputfile
def generate_PB_ON_snp_file(chromo, snps, outputfile):
    '''
    (str, dict, file) -> file
    Take a given chromo and a dictionnary with SNP counts in PB, PX and KSR strains 
    for each position on chromo and save the SNP information to outputfile    
    '''
    
    # open outputfile for writing
    newfile = open(outputfile, 'w')
    
    header = '\t'.join(['chromo', 'pos', 'ref', 'alt', 'flag_snp', 'px_ref', 'ksr_ref', 'px+ksr_ref', 'pb_ref', 'px_alt', 'ksr_alt', 'px+ksr_alt', 'pb_alt'])
    # write header to file
    newfile.write(header + '\n')    
    
    # make a list of positions in dict
    positions = [j for j in snps]
    # sort list
    positions.sort()
    print('positions sorted')
    print(chromo, len(positions))
    
    # loop over the sorted positions
    for i in positions:
        # get ref and alt alleles
        ref_allele = snps[i][0]
        alt_allele = snps[i][1]
        # get 1-based position
        pos = i + 1
        px_ref = snps[i][2]
        px_alt = snps[i][3]
        ksr_ref = snps[i][4]
        ksr_alt = snps[i][5]
        pb_ref = snps[i][6]
        pb_alt = snps[i][7]
        pxksr_ref = px_ref + ksr_ref
        pxksr_alt = px_alt + ksr_alt
        if ref_allele == alt_allele:
            SNP = 'no_snp'
        elif ref_allele != alt_allele:
            SNP = 'snp'
                
        # write info to file
        line = '\t'.join([chromo, str(pos), ref_allele, alt_allele, SNP, str(px_ref), str(ksr_ref), str(pxksr_ref), str(pb_ref), str(px_alt), str(ksr_alt), str(pxksr_alt), str(pb_alt)])
        newfile.write(line + '\n')
    
    # close file after reading    
    newfile.close()

 
  
    
    