# -*- coding: utf-8 -*-
"""
Created on Thu May 28 18:34:40 2015

@author: Richard
"""



from manipulate_sequences import *
from miRNA_target import *
from genomic_coordinates import *



def check_seq_position():
    '''
    
    '''
    
    gaps = 0

    for i in range(len(b)):
        if b[i] != '-':
            j = i - gaps
        else:
            gaps += 1
        print(i, j, a[j], b[i], end = '\t')
        if b[i] == '-':
            print(b[i])
        else:
            print(a[j] == b[i])









# use this function to verify that sites in species 1 and 2 are conserved
def is_site_conserved(UTR_sp1, UTR_sp2, start, end):
    '''
    (str, int, int, str, int, int) -> bool
    Take the aligned UTR sequences from species 1 and 2, and the coordinates
    in 0-based index of sites in the alignment and return True is sites are
    conserved and False otherwise
    '''
    
    site_sp1 = UTR_sp1[start: end].upper()
    site_sp2 = UTR_sp2[start: end].upper()
    
    return site_sp1 == site_sp2    
    
    

# use this function to verify that site is 8mer-1A
def is_site_8mer_1A(site_type, start, end, UTR_seq):
    '''
    (str, str, int, int) -> bool
    Take the site type, the coordinates of the site in the the UTR sequence
    (0-index based), and return True is site type is 8mer-1a, otherwise
    return False
    '''
    
    # 8mer-1A: exact match at 2-8 (the seed + position 8) followed by an 'A'
    
    # check that site type is 8mer-1a    
    if site_type == '8mer-1a':
        # check that site has correct length and that last position is A
        if (end - start) == 8 and UTR_seq.upper()[end -1] == 'A':
            return True
        else:
            return False
        
# use this function to verify that site is 8mer-1U    
def is_site_8mer_1U(site_type, start, end, UTR_seq):
    '''
    (str, str, int, int) -> bool
    Take the site type, the coordinates of the site in the the UTR sequence
    (0-index based), and return True is site type is 8mer-1u, otherwise
    return False
    '''
        
    # 8mer-1U: exact match at 2-8 (the seed + position 8) followed by an 'U'  
        
    # check if site is 8mer-1u
    if site_type == '8mer-1u':
        # check that site has correct length and that last position is U
        if (end - start) == 8 and UTR_seq.upper()[end - 1] == 'U' :
            return True
        else:
            return False
    
# use this function to verify that site is 7mer-m8
def is_site_7mer_m8(site_type, start, end, UTR_seq):
    '''
    (str, str, int, int) -> bool
    Take the site type, the coordinates of the site in the the UTR sequence
    (0-index based), and return True is site type is 7mer-m8, otherwise
    return False
    '''
    
    # 7mer-m8: exact match at 2-8 (the seed + position 8)
    
    # check that site_type is 7mer-m8
    if site_type == '7mer-m8':
        # check that site has correct length 
        if end - start == 7:
            return True
        else:
            return False


# use this function to verify that site is 7mer-1A
def is_site_7mer_1A(site_type, start, end, UTR_seq):
    '''
    (str, str, int, int) -> bool
    Take the site type, the coordinates of the site in the the UTR sequence
    (0-index based), and return True is site type is 7mer-1a, otherwise
    return False
    '''
    
    # 7mer-1A: exact match at 2-7 (the seed) followed by an 'A'
    
    # check that site is 7mir-1a
    if site_type == '7mer-1a':
        # check that site has correct length and that last position is 'A'
        if end - start == 7 and UTR_seq[end - 1] == 'A':
            return True
        else:
            return False


# use this function to verify that site is 6mer
def is_site_6mer(site_type, start, end, UTR_seq):
    '''
    (str, str, int, int) -> bool
    Take the site type, the coordinates of the site in the the UTR sequence
    (0-index based), and return True is site type is 6mer, otherwise
    return False
    '''
    
    # 6mer: exact match at 2-7 of the mature miRNA (the seed)
    
    # check that site is 6mer
    if site_type == '6mer':
        # check that site has correct length
        if end - start  == 6:
            return True
        else:
            False


# use this function to verify that site is 6mer-1a
def is_site_6mer_1A(site_type, start, end, UTR_seq):
    '''
    (str, str, int, int) -> bool
    Take the site type, the coordinates of the site in the the UTR sequence
    (0-index based), and return True is site type is 6mer_1a, otherwise
    return False
    '''

    # 6mer-1A: exact match at 2-6 followed by an 'A'
    
    # check that site is 6mer-1a
    if site_type == '6mer-1a':
        # check that site has correct length and last position is A
        if end - start == 6 and UTR_seq[end -1] == 'A':
            return True
        else:
            return False


# use this function to verify that the site matches the description
def check_site_type(site_type, start, end, UTR_seq):
    '''
    (str, int, int, str) -> bool
    Given a type of miRNA site, and its coordinates in UTR_seq in 0-index based,
    verify that the site in the site conforms to expected length and type
    '''
    
    if site_type == '8mer-1a':
        return is_site_8mer_1A(site_type, start, end, UTR_seq)
    elif site_type == '8mer-1u':
        return is_site_8mer_1U(site_type, start, end, UTR_seq)    
    elif site_type == '7mer-m8':
        return is_site_7mer_m8(site_type, start, end, UTR_seq)
    elif site_type == '7mer-1a':
        return is_site_7mer_1A(site_type, start, end, UTR_seq)
    elif site_type == '6mer':
        return is_site_6mer(site_type, start, end, UTR_seq)
    elif site_type == '6mer-1a':
        return is_site_6mer_1A(site_type, start, end, UTR_seq)   
    
    
# use this function to generate dicts f transcript: aligned UTRs
def get_aligned_UTR_from_TargetScan_input(UTR_seq_targetscan, other_species):
    '''
    (file, str) -> (dict, dict)
    Take the Targetscan input file with aligned UTR sequences between remanei
    and other_species elegans or latens and return a tuple with dictionnaries
    of remanei UTR sequences and other_species UTR sequences. Keep gaps in all
    sequences.    
    '''
    
    # get the species_ID
    if other_species == 'elegans':
        species2 = '6239'
    elif other_species == 'latens':
        species2 = '1503980'
    
    # open file for reading
    infile = open(UTR_seq_targetscan, 'r')
    
    # create dict to store the remanei and other_species UTR sequences
    # transcripts have the remanei transcript name in both dicts
    remanei_UTR = {}
    species_UTR = {}
    
    # read file and populate dict
    for line in infile:
        line = line.rstrip()
        if line != '':
            line = line.split()
            # check species_ID
            if line[1] == species2:
                # populate other_species dict
                species_UTR[line[0]] = line[-1]
            elif line[1] == '31234':
                remanei_UTR[line[0]] = line[-1]
    
    # close file after reading
    infile.close()
    
    return remanei_UTR, species_UTR




# use this function to get the indices of target sites on chromosome
# from the site positions in the targetscan output
def convert_site_position_UTR_to_site_position_chromo(site_start, site_end, UTR_start, UTR_end, UTR_sense, length_chromo):
    '''
    (int, int, int, int, str, int) -> (int, int)
    Take the coordinates of a target site on UTR predicted by targetscan,
    the coordinates of the UTR and its orientation on chromosome, the chromosome
    length and return the coordinates of the site relative to chromosome    
    Precondition: all coordinates are 0-index based    
    '''

    # site_start is starting index of site on UTR
    # site_end is ending index of site on UTR    
    # UTR_start is starting index of UTR in chromo
    # UTR_end is ending index of UTR in chromo
    # UTR_sense is orientation of UTR on chromo
    # chromo_length is the length of chromo
    
    
    # compute site length and UTR length
    site_length = site_end - site_start
    UTR_length = UTR_end - UTR_start    
    
    # check  UTR orientation
    if UTR_sense == '+':
        # site coordinates on chromosome can be directly deduced
        # start = start_index_site + start_index_UTR
        start = site_start + UTR_start
        # end = start + length_site, also end = end_index_site + start_index_UTR
        end = start + site_length
        return start, end
    elif UTR_sense == '-':
        # sites are predicted based on reverse_complement of UTR
        # need to adjust coordinates of site on chromo
        # get the indices of site in UTR sequence in - sense
        site_start_UTR = UTR_length - site_end
        # site indices on chromo can be deduced from site indices on UTR(-)
        start = site_start_UTR + UTR_start
        end = start + site_length
        return start, end


  
        
# use this function to parse the targetscan output files with remanei and latens sites
def parse_crm_cla_sites(UTR_seq_targetscan, predicted_targets, seeds):
    '''
    (file, file, dict) -> dict
    Take the targetscan input file with aligned UTR sequences, the targetscan
    output file with predicted targets and a dict with seeds: mirna pairs
    and return a dict with site as key and list of seeds targeting that
    site as value
    '''
    
    # targetscan calls sites conserved when gaps are destroy a site 
    # these sites are considered not conserved here
        
    # create dicts for each species {transcript: gaped UTR seq}
    # transcript names are the remanei name in each species
    remanei_UTR, latens_UTR = get_aligned_UTR_from_TargetScan_input(UTR_seq_targetscan, 'latens')
        
    # create a dict of {mirna name : seed_seq}
    mirna_seed = {}
    for seed in seeds:
        for mirna in seeds[seed]:
            mirna_seed[mirna] = seed
    
    # create a dict to store the target sites
    targets = {}

    # open file for reading
    infile = open(predicted_targets, 'r')
    # skip header
    infile.readline()
    # go through file
    for line in infile:
        line = line.rstrip()
        if line != '':
            line = line.split()
            # record only the remanei sites
            if line[2] == '31234':
                # get the site type
                site_type = line[8]
                # get 0-ibdex based MSA positions
                MSA_start = int(line[3]) - 1
                MSA_end = int(line[4])
                # get transcript name           
                transcript = line[0]
                # get UTR_seq
                crm_UTR_seq = remanei_UTR[transcript]
                cla_UTR_seq = latens_UTR[transcript]
                # verify site type has correct length (ie. no gaps)
                if check_site_type(site_type, MSA_start, MSA_end, crm_UTR_seq) == True:
                    # get coordinates in UTR in 0-index based
                    start = int(line[5]) - 1
                    end = int(line[6])
                    # get the seed sequence
                    seed_seq = mirna_seed[line[1]]
                    # check if site is conserved or not
                    if is_site_conserved(crm_UTR_seq, cla_UTR_seq, MSA_start, MSA_end) == True:
                        # if site conserved, add string to show conservation in remanei and latens
                        site_info = [transcript, site_type, str(start), str(end), str(MSA_start), str(MSA_end), 'crm:cla']
                    elif is_site_conserved(crm_UTR_seq, cla_UTR_seq, MSA_start, MSA_end) == False:
                        # if site not conserved, use string 'crm' to show site present only in remanei
                        site_info = [transcript, site_type, str(start), str(end), str(MSA_start), str(MSA_end), 'crm']
                    # some sites may be shared by different mirna families (eg: 7mer-1a, 6mer)
                    # use site_info as key and list of seed_seq as value
                    site_info = '|'.join(site_info)
                    # populate dict
                    if site_info in targets:
                        targets[site_info].append(seed_seq)
                    else:
                        targets[site_info] = [seed_seq]
                        
    # close file after reading
    infile.close()    
    
    return targets
    

# use this function to generate the summary table of remanei and latens sites
# use this table to compute divergence at mirna target sites between species
def make_table_cremanei_clatens_sites(caeno_gff, assembly, threshold, seeds, targets, outputfile):
    '''
    (file, file, int, dict, dict, file) -> file
    Take the remanei genome and its GGF annotation files, a threshold used to
    grab the UTR sequences, a dict with seeds :mirna pairs, a dict with
    site info: seeds pairs and generate a summary table of target sites
    '''
    
    # targets is a dict with site info : seed pairs
    # {'transcript|site_type|start|end|MSA_start|MSA_end|conservation': [seed1, seed2]}
    # seeds is a dict of seeds : mirna pairs 
       
    # get the coordinates of the UTR/downstream seq
    # {TS1 : [chromo, start, end, orientation]}
    UTR_coord = get_three_prime_UTR_positions(caeno_gff, assembly, threshold)
    
    # convert genome fasta file to dict
    # convert the fasta assembly to a dictionnary with chromo / scaffold as key and sequence as value
    genome = convert_fasta(assembly)
    
    # get the remanei transcript coordinates
    # {gene1 : [chromo, start, end, orientation]}
    transcripts_coord = get_genes_coordinates(caeno_gff)
    
    # create a dict of remanei transcripts with annotated UTRs
    # {transcript_name : [chromo, start, end, orienation]}
    annotated_UTR = grab_annotated_three_prime_coordinates(caeno_gff)
    
    # open file for writing
    newfile = open(outputfile, 'w')
       
    # site positions in the UTR and in the MSA correspond to UTR in the + sense
    # even if UTR orientation is - on chromo     
    header = '\t'.join(['transcript', 'seed', 'N_mirnas', 'type_of_site',
                        'conservation_level', 'site_MSA_start(+)', 'site_MSA_end(+)',
                        'chromo', 'downstream/UTR', 'UTR_start', 'UTR_end',
                        'site_UTR_start(+)', 'site_UTR_end(+)', 'site_chromo_start',
                        'site_chromo_end', 'orientation'])
                        
    # write header to file
    newfile.write(header + '\n')                        
        
    # loop over sites in dict
    for site in targets:
        # get site features by parsing the site string
        # positions are 0-based index
        site_info = site.split('|')
        transcript = site_info[0]
        site_type = site_info[1]
        start = int(site_info[2])
        end = int(site_info[3])
        MSA_start = int(site_info[4])
        MSA_end = int(site_info[5])
        conservation = site_info[6]
        # compute the number of miRNA regulator for this site
        N_regulators = 0
        for seed_seq in targets[site]:
            N_regulators += len(seeds[seed_seq])
        # get all seeds targeting the given site
        mir_fams = ':'.join(targets[site])
        # get chromosome
        chromo = transcripts_coord[transcript][0]
        # get orientation
        orientation = transcripts_coord[transcript][-1]
        # is UTR annotated UTR or downstream sequence
        if transcript in annotated_UTR:
            UTR_status = 'UTR'
        else:
            UTR_status = 'downstream'
        # get the site positions in 1-based index  
        # get positions of site in the aligned sequence
        site_MSA_start  = MSA_start + 1
        site_MSA_end = MSA_end
        # get positions of site in the UTR 1-based index for file
        site_UTR_start = start + 1
        site_UTR_end = end
        # get UTR coordinates in 0-based index
        UTR_start = UTR_coord[transcript][1] - 1
        UTR_end = UTR_coord[transcript][2]
        # get the positions of site on chromo
        # get chromo length
        length_chromo = len(genome[chromo])
        site_chromo_start, site_chromo_end = convert_site_position_UTR_to_site_position_chromo(start, end, UTR_start, UTR_end, orientation, length_chromo)
        # convert site_chromo_start and UTR_start to 1-based index
        site_chromo_start += 1
        UTR_start += 1
        
        # write content to outputfile
        content = '\t'.join([transcript, mir_fams, str(N_regulators), site_type, conservation, str(site_MSA_start),
        str(site_MSA_end), chromo, UTR_status, str(UTR_start), str(UTR_end), str(site_UTR_start),
        str(site_UTR_end), str(site_chromo_start), str(site_chromo_end), orientation])
        
        newfile.write(content + '\n')
        
    
    # close file after writing
    newfile.close()
        
        

    
# use this function to create dicts from the targetscan output
def parse_targetscan_output(UTR_seq_targetscan, predicted_targets, seeds, other_species):
    '''
    (file, file, dict, str) -> dict
    Take the targetscan input sequence file, the targetscan output, a dict of
    seeds :mirna pairs, and a species name ('remanei', 'elegans' or 'latens')
    and return a dict with site as key and a list of seeds targeting that site
    as value. Only records conserved sites if species is not remanei
    '''
    
    # targetscan calls sites conserved when gaps are destroy a site 
    # these sites are considered not conserved here
    
    # check if other_species if different than remanei
    # create dicts for each species {transcript: gaped UTR seq}
    # transcript names are the remanei name in each species
    if other_species != 'remanei':
        # create dicts with remanei and latens or elegans UTR aligned sequences
        remanei_UTR, species_UTR = get_aligned_UTR_from_TargetScan_input(UTR_seq_targetscan, other_species)
    elif other_species == 'remanei':
        # create a dict with the remanei sequences
        remanei_UTR = {}        
        # open file for reading
        infile = open(UTR_seq_targetscan, 'r')
        # read file and populate dict
        for line in infile:
            line = line.rstrip()
            if line != '':
                line = line.split()
                remanei_UTR[line[0]] = line[-1]
        # close file after reading
        infile.close()
    
#    # create a dict of seed sequence : list of mirnas sharing the seed
#    seeds = seed_mirnas(mature_fasta)
       
    # create a dict of {mirna name : seed_seq}
    mirna_seed = {}
    for seed in seeds:
        for mirna in seeds[seed]:
            mirna_seed[mirna] = seed
    
    
    # create a dict to store the target sites
    targets = {}

    # open file for reading
    infile = open(predicted_targets, 'r')
    # skip header
    infile.readline()
    # go through file
    for line in infile:
        line = line.rstrip()
        if line != '':
            line = line.split()
            # record only the remanei sites
            if line[2] == '31234':
                # get the site type
                site_type = line[8]
                # get 0-ibdex based MSA positions
                MSA_start = int(line[3]) - 1
                MSA_end = int(line[4])
                # get transcript name           
                transcript = line[0]
                # get UTR_seq
                crm_UTR_seq = remanei_UTR[transcript]
                if other_species != 'remanei':
                    # get the UTR sequence of the other species
                    species_UTR_seq = species_UTR[transcript]
                # verify site type has correct length (ie. no gaps)
                if check_site_type(site_type, MSA_start, MSA_end, crm_UTR_seq) == True:
                    # get coordinates in UTR in 0-index based
                    start = int(line[5]) - 1
                    end = int(line[6])
                    # get the seed sequence
                    seed_seq = mirna_seed[line[1]]
                    # check if site is conserved or not
                    # record site: site is transcript:start:end:site_type
                    if other_species != 'remanei':
                        # only record conserved sites between remanei and other_species
                        if is_site_conserved(crm_UTR_seq, species_UTR_seq, MSA_start, MSA_end) == True:
                            # if site conserved, add string to show conservation in remanei and other_species
                            site_info = [transcript, str(start), str(end), site_type]
                            # some sites may be shared by different mirna families (eg: 7mer-1a, 6mer)
                            # use site_info as key and list of seed_seq as value
                            site_info = '|'.join(site_info)
                            # populate dict
                            if site_info in targets:
                                targets[site_info].append(seed_seq)
                            else:
                                targets[site_info] = [seed_seq]
                    elif other_species == 'remanei':
                        # record all sites
                        site_info = [transcript, str(start), str(end), site_type]
                        # some sites may be shared by different mirna families (eg: 7mer-1a, 6mer)
                        # use site_info as key and list of seed_seq as value
                        site_info = '|'.join(site_info)
                        # populate dict
                        if site_info in targets:
                            targets[site_info].append(seed_seq)
                        else:
                            targets[site_info] = [seed_seq]
                        
    # close file after reading
    infile.close()    
    
    return targets


# use this function to make a summary table with remanei target sites
# to use to analyze SNPs in targets
def summary_table_cremanei_target_sites(crm_sites, crm_cla_conserved_sites, crm_cel_conserved_sites, assembly, caeno_gff, seeds, threshold, outputfile):
    '''
    (dict, dict, dict, file, file, dict, int, file)
    Take the parsed target prediction files, the genome and its annotation files, 
    a dict with seeds: mirna pairs, a threshold of UTR length, and save information
    about the remanei miRNA target sites, adding conservation level by comparing
    the target sites obtained with the remanei transcripts only, remanei and latens
    orthologs and remanei and elegans orthologs, and by adding the coordinates of
    the site on chromosome    
    '''
     

    # targetscan calls sites conserved when gaps are destroy a site 
    # these sites are considered not conserved here    
    
    # get the coordinates of the UTR/downstream seq
    # {TS1 : [chromo, start, end, orientation]}
    UTR_coord = get_three_prime_UTR_positions(caeno_gff, assembly, threshold)
    
    # convert genome fasta file to dict
    # convert the fasta assembly to a dictionnary with chromo / scaffold as key and sequence as value
    genome = convert_fasta(assembly)
         
    # get the remanei transcript coordinates
    # {gene1 : [chromo, start, end, orientation]}
    transcripts_coord = get_genes_coordinates(caeno_gff)
    
    # create a dict of remanei transcripts with annotated UTRs
    # {transcript_name : [chromo, start, end, orienation]}
    annotated_UTR = grab_annotated_three_prime_coordinates(caeno_gff)
    
    # open file for writing
    newfile = open(outputfile, 'w')
    
    # site positions in the UTR and in the MSA correspond to UTR in the + sense
    # even if UTR orientation is - on chromo     
        
    header = '\t'.join(['transcript', 'seed', 'N_mirnas', 'type_of_site',
                        'conservation_level', 'chromo', 'downstream/UTR',
                        'UTR_start', 'UTR_end', 'site_UTR_start(+)',
                        'site_UTR_end(+)', 'site_chromo_start',
                        'site_chromo_end', 'orientation'])
                        
    # write header to file
    newfile.write(header + '\n')                        
        
    # site_info = [transcript, str(start), str(end), site_type]
        
    # loop over sites in dict of remanei target sites
    for site in crm_sites:
        # get site features by parsing the site string
        site_info = site.split('|')
        # get positions (0-based index)
        transcript = site_info[0]
        start = int(site_info[1])
        end = int(site_info[2])
        site_type = site_info[-1]        
        
        # get the conservation level
        if site in crm_cla_conserved_sites and site in crm_cel_conserved_sites:
            conservation = 'crm:cla:cel'
        elif site in crm_cla_conserved_sites and site not in crm_cel_conserved_sites:
            conservation = 'crm:cla'
        elif site not in crm_cla_conserved_sites and site in crm_cel_conserved_sites:
            conservation = 'crm:cel'
        elif site not in crm_cla_conserved_sites and site not in crm_cel_conserved_sites:
            conservation = 'crm'
        
        # compute the number of miRNA regulator for this site
        N_regulators = 0
        for seed_seq in crm_sites[site]:
            N_regulators += len(seeds[seed_seq])
        # get all seeds targeting the given site
        mir_fams = ':'.join(crm_sites[site])
        # get chromosome
        chromo = transcripts_coord[transcript][0]
        # get orientation
        orientation = transcripts_coord[transcript][-1]
        # is UTR annotated UTR or downstream sequence
        if transcript in annotated_UTR:
            UTR_status = 'UTR'
        else:
            UTR_status = 'downstream'
        # get positions of site in the UTR, 1-based index for file
        site_UTR_start = start + 1
        site_UTR_end = end
        # get UTR coordinates in 0-based index
        UTR_start = UTR_coord[transcript][1] - 1
        UTR_end = UTR_coord[transcript][2]
        # get the positions of site on chromo
        # get chromo length
        length_chromo = len(genome[chromo])
        site_chromo_start, site_chromo_end = convert_site_position_UTR_to_site_position_chromo(start, end, UTR_start, UTR_end, orientation, length_chromo)
        # convert site_chromo_start and UTR_start to 1-based index
        site_chromo_start += 1
        UTR_start += 1
        
        # write content to outputfile
        content = '\t'.join([transcript, mir_fams, str(N_regulators), site_type, conservation, 
                             chromo, UTR_status, str(UTR_start), str(UTR_end), str(site_UTR_start),
                             str(site_UTR_end), str(site_chromo_start), str(site_chromo_end), orientation])
        
        newfile.write(content + '\n')
        
    
    # close file after writing
    newfile.close()
        


# use this function to get the coordinates of the latens-specific miRNA target sites
def grab_latens_specific_sites(UTR_seq_targetscan, predicted_targets, mature_fasta):
    '''
    (file, file, file) -> dict
    Take the targetscan input file with aligned UTR sequences, the targetscan
    output file with predicted targets, the fasta file of remanei mature mirnas
    and return a dict with site as key and list of seeds for targeting that site as value
    '''
    
    # targetscan calls sites conserved when gaps are destroy a site 
    # these sites are considered not conserved here
        
    # create dicts for each species {transcript: gaped UTR seq}
    # transcript names are the remanei name in each species
    remanei_UTR, latens_UTR = get_aligned_UTR_from_TargetScan_input(UTR_seq_targetscan, 'latens')
        
    # create a dict of seed sequence : list of mirnas sharing the seed
    seeds = seed_mirnas(mature_fasta)
       
    # create a dict of {mirna name : seed_seq}
    mirna_seed = {}
    for seed in seeds:
        for mirna in seeds[seed]:
            mirna_seed[mirna] = seed
    
    # create a dict to store the target sites
    targets = {}

    # open file for reading
    infile = open(predicted_targets, 'r')
    # skip header
    infile.readline()
    # go through file
    for line in infile:
        line = line.rstrip()
        if line != '':
            line = line.split()
            # record only the latens sites
            if line[2] == '1503980':
                # get the site type
                site_type = line[8]
                # get 0-ibdex based MSA positions
                MSA_start = int(line[3]) - 1
                MSA_end = int(line[4])
                # get transcript name           
                transcript = line[0]
                # get UTR_seq
                crm_UTR_seq = remanei_UTR[transcript]
                cla_UTR_seq = latens_UTR[transcript]
                # verify site type has correct length (ie. no gaps)
                if check_site_type(site_type, MSA_start, MSA_end, crm_UTR_seq) == True:
                    # get coordinates in UTR in 0-index based
                    start = int(line[5]) - 1
                    end = int(line[6])
                    # get the seed sequence
                    seed_seq = mirna_seed[line[1]]
                    # record only latens-specific miRNA targets
                    if is_site_conserved(crm_UTR_seq, cla_UTR_seq, MSA_start, MSA_end) == False:
                        # if site not conserved, use string 'cla' to show site present only in latens
                        site_info = [transcript, site_type, str(start), str(end), str(MSA_start), str(MSA_end), 'cla']
                        # some sites may be shared by different mirna families (eg: 7mer-1a, 6mer)
                        # use site_info as key and list of seed_seq as value
                        site_info = '|'.join(site_info)
                        # populate dict
                        if site_info in targets:
                            targets[site_info].append(seed_seq)
                        else:
                            targets[site_info] = [seed_seq]
                        
    # close file after reading
    infile.close()    
    
    return targets
    

# use this function to generate a summary table of latens specific sites
# use this table to find the non-remanei targets with 1 substition 
def make_table_clatens_specific_sites(cla_specific_targets, mature_fasta, caeno_gff, outputfile):
    '''
    (dict, file, file, file) -> file
    Take file with coordinates of the latens-specific miRNA targets, 
    the fasta file of mature remanei miRNAs, the remanei GFF file, and    
    generate a summary table of coordinates if latens-specific
    miRNA target sites
    '''
    
    # get the remanei transcript coordinates
    # {gene1 : [chromo, start, end, orientation]}
    transcripts_coord = get_genes_coordinates(caeno_gff)
    
    # create a dict of seed sequence : list of mirnas sharing the seed
    seeds = seed_mirnas(mature_fasta)    
    # open file for writing
    newfile = open(outputfile, 'w')
    # site positions in the UTR and in the MSA correspond to UTR in the + sense
    # even if UTR orientation is - on chromo     
    header = '\t'.join(['transcript', 'seed', 'N_mirnas', 'type_of_site',
                        'conservation_level', 'site_MSA_start(+)', 'site_MSA_end(+)',
                        'site_UTR_start(+)', 'site_UTR_end(+)', 'transcript_orientation'])
    # write header to file
    newfile.write(header + '\n')                        
    
    # loop over sites in dict
    for site in cla_specific_targets:
        # get site features by parsing the site string
        # positions are 0-based index
        site_info = site.split('|')
        transcript = site_info[0]
        site_type = site_info[1]
        start = int(site_info[2])
        end = int(site_info[3])
        MSA_start = int(site_info[4])
        MSA_end = int(site_info[5])
        conservation = site_info[6]
        # compute the number of miRNA regulator for this site
        N_regulators = 0
        for seed_seq in cla_specific_targets[site]:
            N_regulators += len(seeds[seed_seq])
        # get all seeds targeting the given site
        mir_fams = ':'.join(cla_specific_targets[site])
        # get the site positions in 1-based index  
        # get positions of site in the aligned sequence
        site_MSA_start  = MSA_start + 1
        site_MSA_end = MSA_end
        # get positions of site in the UTR 1-based index for file
        site_UTR_start = start + 1
        site_UTR_end = end
                
        # write content to outputfile
        content = '\t'.join([transcript, mir_fams, str(N_regulators), site_type, conservation, str(site_MSA_start),
        str(site_MSA_end), str(site_UTR_start), str(site_UTR_end), transcripts_coord[transcript][-1]])
        
        newfile.write(content + '\n')
        
    # close file after writing
    newfile.close()


