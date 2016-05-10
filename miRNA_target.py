# -*- coding: utf-8 -*-
"""
Created on Tue May 26 11:09:39 2015

@author: Richard
"""



   
# use this function to generate the targetscan input file for the Cremanei transcripts only
def make_cremanei_targetscan_input_sequences(Crm_UTR_sequences, Crm_transcripts_genes, orthologs_file, outputfile):
    '''
    (dict, dict, file, file) -> file
    With the dict of remanei transcripts' UTR sequences, generate the inputfile
    for TargetScan. Using the ortholog file and the dict of transcripts : genes pairs
    use a single transcript per gene. Use the latens ortholog if present or
    use any other transcript if no ortholog exists.
    '''
    
    # Crem_transcripts_genes is a dict of transcripts : genes pairs
    # Crm_UTR_sequences is a dict of transcript: sequence pairs    
    # ortholog file contains the 1: 1 orthologs between remanei and latens
    
    # make a set of transcripts with 1:1 orthologs in C. latens
    orthologs = set()
    # open file for reading
    orthofile = open(orthologs_file, 'r')
    # header contains multiple lines
    # read first line
    line = orthofile.readline()
    # keep reading the entire header
    while not line.startswith('Cremanei'):
        line = orthofile.readline()
    # go through file, get the name of the remanei transcript
    for line in orthofile:
        line = line.rstrip()
        if line != '':
            line = line.split()
            orthologs.add(line[0])
    # close file after reading
    orthofile.close()

    # make a set of genes already picked
    picked_genes = set()
    
    # open outfile for writing (transcript_ID species_ID UTR_seq)
    newfile = open(outputfile, 'w')
    
    # go through transcripts in UTR dict
    # record all transcripts with orthologs first
    for transcript in Crm_UTR_sequences:
        if transcript in orthologs:
            # only 1 transcript per gene is recorded in ortholog file, 
            # but reinforce 1 transcript per gene
            if Crm_transcripts_genes[transcript] not in picked_genes:
                # gene has not been picked, can record transcript + UTR
                # convert T to U
                newfile.write(transcript + '\t' + '31234' + '\t' + Crm_UTR_sequences[transcript].replace('T', 'U') + '\n')
                # add the gene to set of already picked genes to enure 1 transcript per gene
                picked_genes.add(Crm_transcripts_genes[transcript])
    
    # after recording the transcripts with orthologs, record other transcripts
    for transcript in Crm_UTR_sequences:
        # keeping 1 transcript per gene
        if Crm_transcripts_genes[transcript] not in picked_genes:
            # record transcript, convert T to U
            newfile.write(transcript + '\t' + '31234' + '\t' + Crm_UTR_sequences[transcript].replace('T', 'U') + '\n')
            # add the gene to set of already picked genes to enure 1 transcript per gene
            picked_genes.add(Crm_transcripts_genes[transcript])
    
    
    newfile.close()  

   
# use this function to generate the targetscan input file with aligned Cremanei-Clatens sequences  
def make_cremanei_clatens_targetscan_input_sequences(directory, outputfile):
    '''
    (str, file) -> file
    Generate the input sequence file for TargetScan using the aligned UTR sequences
    located in directory. Use only Cremanei and Clatens orthologs      
    '''

    # open outputfile for writing
    newfile = open(outputfile, 'w')

    # grab the aligned remanei and latens sequences and write to newfile
    files = os.listdir(directory)
    for filename in files:
        # check that the file contains the crem-cla orthologs
        if 'CRE' in filename and '.txt' in filename:
            # convert fasta file to a dictionnary for seq_name : sequence
            fasta_seq = convert_fasta(directory + filename)
            # get the name of the remanei transcript
            for seq_name in fasta_seq:
                if 'CRE_PX356' in seq_name:
                    transcript_name = seq_name
            # use the cremanei transcript name for the latens ortholog
            # but add species ID to distinguish species
            # convert T to U
            for seq_name in fasta_seq:
                if 'CRE_PX356' in seq_name:
                    newfile.write(transcript_name + '\t' + '31234' + '\t' + fasta_seq[seq_name].upper().replace('T', 'U') + '\n')
                elif 'CLA' in seq_name:
                    newfile.write(transcript_name + '\t' + '1503980' + '\t' + fasta_seq[seq_name].upper().replace('T', 'U') + '\n')
    
    
    newfile.close()
                
         
# use this function to generate the TargetScan input sequence file with remanei-elegans aligned sequences
def make_cremanei_celegans_targetscan_input_sequences(directory, outputfile):
    '''
    (str, file) -> file
    Generate the input sequence file for TargetScan using the aligned UTR sequences
    located in directory. Use only Cremanei and Celegans orthologs      
    '''

    # open outputfile for writing
    newfile = open(outputfile, 'w')

    # grab the aligned remanei and latens sequences and write to newfile
    files = os.listdir(directory)
    for filename in files:
        # check that the file contains the crem-cel orthologs
        if 'CRE' in filename and '.txt' in filename:
            # convert fasta file to a dictionnary for seq_name : sequence
            fasta_seq = convert_fasta(directory + filename)
            # get the name of the remanei transcript
            for seq_name in fasta_seq:
                if 'CRE_PX356' in seq_name:
                    transcript_name = seq_name
            # use the cremanei transcript name for the elegans ortholog
            # but add species ID to distinguish species
            # convert T to U
            for seq_name in fasta_seq:
                if 'CRE_PX356' in seq_name:
                    newfile.write(transcript_name + '\t' + '31234' + '\t' + fasta_seq[seq_name].upper().replace('T', 'U') + '\n')
                else:
                    newfile.write(transcript_name + '\t' + '6239' + '\t' + fasta_seq[seq_name].upper().replace('T', 'U') + '\n')
    
    
    newfile.close()         
         
# use this function to get the mirna coordinates on each chromo
def get_mirna_loci(mirna_coord_file):
    '''
    (file) -> dict
    Take the file with mirna genomic coordinates (hairpin or mature)
    and return a dictionary with chromo as key and a list of list coordinates
    for each mirna on that chromo
    '''
    
    # create a dict {chromo: [[start, end, orientation]]}
    mirnas_coord = {}
    # openfile for reading
    infile = open(mirna_coord_file, 'r')
    # skip header
    infile.readline()
    # loop over file
    for line in infile:
        line = line.rstrip()
        if line != '':
            line = line.split()
            chromo = line[1]
            # get positions in 0-based indices
            start = int(line[2]) -1
            end = int(line[3])
            orientation = line[4]
            # check if chromo in dict
            if chromo in mirnas_coord:
                mirnas_coord[chromo].append([start, end, orientation])
            else:
                mirnas_coord[chromo] = [[start, end, orientation]]
    infile.close()

    return mirnas_coord
    

# use this function to grab the coordinates of the miRNA targets 
def get_miRNA_target_coord(mirna_target_coord_file):
    '''
    (file) -> dict
    Take the file of miRNA target coordinates and return a dict with gene as key
    and a list with miRNA coordinates for all miRNAs targetting that gene
    '''
    
    # create a dict to store target site information
    # {gene: [[chromo, start, end, orientation, seed, N_mirnas, site_type, conservation, utr]]}
    targets = {}
        
    # open file for reading
    infile = open(mirna_target_coord_file, 'r')
    # skip header
    infile.readline()
    # loop over file
    for line in infile:
        line = line.rstrip()
        if line != '':
            line = line.split()
            # get transcript
            transcript = line[0]
            # get miRNA seed
            seed = line[1]
            # get the number of mirnas targeting that site
            N_mirnas = int(line[2])
            # get the type of site
            site_type = line[3]
            # get conservation level
            conservation = line[4]
            # get chromo
            chromo = line[5]
            # get UTR sequence type
            utr = line[6]
            # get start and end positions on chromo
            start = int(line[11]) -1
            end = int(line[12])
            # get orientation
            orientation = line[-1]
            # check if transcript is key
            if transcript in targets:
                targets[transcript].append([chromo, start, end, orientation, seed, N_mirnas, site_type, conservation, utr])
            else:
                targets[transcript] = [[chromo, start, end, orientation, seed, N_mirnas, site_type, conservation, utr]]
                
    # close file after reading
    infile.close()
    
    return targets
    
    

# use this function to get the coordinates of each target site on each chromo
def get_miRNA_target_loci(mirna_target_coord_file, unique_transcripts, conservation_level):
    '''
    (file, file, str) -> dict
    Take the file with miRNA target coordinates, the file with valid transcripts
    (a single transcript per gene), a string with the conservation level 
    (all, crm, crm-cla, crm-cla-cel) and return a dict with chromo as key and
    a list of target site coordinates for each target on that chromo
    '''
       
    # get the set of valid transcripts
    valid_transcripts = get_valid_transcripts(unique_transcripts)
        
    # create a dict {chromo: [[start, end, orientation]]}
    targets = {}

    # open file for reading
    infile = open(mirna_target_coord_file, 'r')
    # skip header
    infile.readline()
    # loop over the file
    for line in infile:
        line = line.rstrip()
        if line != '':
            line = line.split()
            # get transcript
            gene = line[0]
            # check that gene is valid
            if gene in valid_transcripts:
                # get conservation
                conservation = line[4]
                # get chromo
                chromo = line[5]
                # get orientation
                orientation = line[-1]
                # get site positions on chromo, 0-based
                start = int(line[11]) -1
                end = int(line[12])
                # check conservation level
                if conservation_level == 'all':
                    # record all site, regardeless of conservation
                    # populate dict, check if chromo is key
                    if chromo in targets:
                        targets[chromo].append([start, end, orientation])
                    else:
                        targets[chromo] = [[start, end, orientation]]
                elif conservation_level == 'crm':
                    # record only remanei-specific sites
                    if conservation == 'crm':
                        # populate dict, check if chromo is key
                        if chromo in targets:
                            targets[chromo].append([start, end, orientation])
                        else:
                            targets[chromo] = [[start, end, orientation]]
                elif conservation_level == 'crm-cla':
                    # record remanei-latens conserved target sites
                    if conservation == 'crm:cla':
                        # poopulate dict, check if chromo is key
                        if chromo in targets:
                            targets[chromo].append([start, end, orientation])
                        else:
                            targets[chromo] = [[start, end, orientation]]
                elif conservation_level == 'crm-cla-cel':
                    # record remanei-latens-elegans conserved target sites
                    if conservation == 'crm:cla:cel':
                        # populate dict, check if chromo is key
                        if chromo in targets:
                            targets[chromo].append([start, end, orientation])
                        else:
                            targets[chromo] = [[start, end, orientation]]
                        
    # close file after reading
    infile.close()
    
    return targets



# use this function to get the positions all miRNA target sites on each chromo
def get_all_miRNA_target_sites(mirna_target_coord_file):
    '''
    (file) -> dict
    Take the file with miRNA target coordinates and return a duictionary with 
    chromo as key and a set with all the positions falling in miRNA targets
    Preconditions: positions of miRNA target sites are 0-based
    '''
       
    # create a dict {chromo: {set of 0-based mirna target positions}}
    targets = {}

    # open file for reading
    infile = open(mirna_target_coord_file, 'r')
    # skip header
    infile.readline()
    # loop over the file
    for line in infile:
        line = line.rstrip()
        if line != '':
            line = line.split()
            # get chromo
            chromo = line[5]
            # get site positions on chromo, 0-based
            start = int(line[11]) -1
            end = int(line[12])
            # check if chromo is key in targets
            if chromo not in targets:
                # initiate set
                targets[chromo] = set()
            # add all the target site positions to set
            for i in range(start, end):
                targets[chromo].add(i)
            
    # close file after reading
    infile.close()
    
    return targets



# use this function to slice the 7bp seed from the mature sequence 
def grab_seed(mature_seq):
    '''
    (str) -> str
    Slice the miRNA mature sequence to return the 7bp seed motif
    '''
    seed = mature_seq[1:8]
    return seed


# use this function to generate a dict of seeds and list pf mirnas pairs
def seed_mirnas(mature_fasta):
    '''
    (file) -> dict
    Take a fasta file of miRNA mature sequences and return a dictionnary with
    seed as key and a list of mirnas from the same family (sharing the same seed)
    as value
    '''
    
    # convert fasta file to 
    mirnas = convert_fasta(mature_fasta)
    # create a dict of seed : [mir1, mir2]
    seeds = {}
    # loop over mature sequences
    for name in mirnas:
        # get seed sequence
        seed_seq = grab_seed(mirnas[name])
        # populate dict
        if seed_seq in seeds:
            seeds[seed_seq].append(name)
        else:
            seeds[seed_seq] = [name]
    return seeds


# use this function to generate a set of seeds for a given species, or group of mature sequences
def get_all_seeds_in_species(mature_fasta):
    '''
    (file) -> set
    Take a fasta file of miRNA mature sequences and return a set of miRNA 
    seed sequences present in a given species
    '''
    
    # create a dict of seed sequence and list of mirna pairs
    seeds = seed_mirnas(mature_fasta)
    # create a set of seeds present in species
    seed_species = {i for i in seeds}
    return seed_species


# use this function to check if a miRNA family is conserved in a group of miRNAs or species
def is_miRNA_family_conserved(seed_seq, seeds_species):
    '''
    (str, set) -> bool
    Return True if the seed sequence is conserved and present in the set of 
    seeds from species, return False otherwise
    '''
    
    if seed_seq in seeds_species:
        return True
    else:
        return False

    

# use this function to get the seed sequences of all miRNAs
def get_seed_sequences(mirna_fam_file):
    '''
    (file) -> list
    Take the miRNA seed family file and return a list with all seeds with U 
    replaced by T    
    '''
    
    # create list
    seeds = []
    
    # open file for reading
    infile = open(mirna_fam_file, 'r')
    # loop over file
    for line in infile:
        line = line.rstrip()
        if line != '':
            line = line.split()
            # replace U by T and add seed to list
            seeds.append(line[1].replace('U', 'T'))
    
    # close file
    infile.close()
    return seeds
    


# use this function to check if a potential miRNA site is matching a seed sequence
def is_matching_seed(site, site_type, seed):
    '''
    (str, str, str)
    Take a DNA sequence, potentially a miRNA target site of given type, and a
    known miRNA seed and return True if the site match the seed or False if it doesn not
    Precondition: site contains only valid nucleotides
    '''
    # convert to upper case
    seed = seed.upper()
    site = site.upper()    
    # convert seed from RNA to DNA
    if 'U' in seed:
        seed = seed.replace('U', 'T')
        
    # check the site_type
    if site_type == '8mer-1a' or site_type == '8mer-1u':
        if reverse_complement(site)[1:] == seed:
            return True
        else:
            return False
    elif site_type == '7mer-m8':
        if reverse_complement(site) == seed:
            return True
        else:
            return False
    elif site_type == '7mer-1a':
        if reverse_complement(site)[1:] == seed[:1]:
            return True
        else:
            return False


       
# use this function to get the DAF of the miRNA targets
def get_DAF_miRNA_targets(chromo_sites, genome_fasta, crm_cla_target_sites_file, UTR_alignments_folder, conservation_scores, conservation_level, unique_transcripts):
    '''
    (dict, file, file, str, dict, str, file) -> list
    Take the dictionary with allele counts for sites with coverage
    (and minimum sample size), the fasta file of the genome sequence,
    the file with remanei target coordinates for sites of genes that have an
    ortholog in latens, the folder with remanei-latens UTR alignments,
    the dictionary with site conservation, a string specifying the level of
    conservation to consider (all, crm, crm-cla, crm-cla-cel), the file with valid transcripts
    (ie. 1 transcript per gene) and return a list of DAF for SNPs in target sites
    Precondition: consider only target sites with coverage at all positions
    '''
    
    # convert genome sequence to dict
    genome = convert_fasta(genome_fasta)     
        
    # create a list of UTR alignment files
    files = [i for i in os.listdir(UTR_alignments_folder) if 'CRE' in i and '.txt' in i]
    
    # get the set of valid transcripts
    valid_transcripts = get_valid_transcripts(unique_transcripts)    
    
    # create a list to store the DAF
    DAF = []
    
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
            # check if valid transcript
            if gene in valid_transcripts:
                # get seed
                seed = line[1]
                # get site_type
                site_type = line[3]
                # get chromo
                chromo = line[7]
                # get coordinates 0-based on chromo
                chromo_start = int(line[13]) -1
                chromo_end = int(line[14])
                orientation = line[15]
                # get coordinates in the multi-sequence UTR alignment
                # the positions correspond to + orientation
                msa_start = int(line[5]) -1
                msa_end = int(line[6])
                # build site 'transcript;seed;site_type;site_start_chromo;site_end_chromo;orientation'
                site = gene + ';' + seed + ';' + site_type + ';' + str(chromo_start) + ';' + str(chromo_end) + ';' + orientation
                # get the target site on chromo
                target_chromo = genome[chromo][chromo_start:chromo_end]
                # check target site orientation on chromo
                if orientation == '-':
                    # take reverse complement of target chromo
                    target_chromo = reverse_complement(target_chromo)
                # get the target site in the aligned remanei UTR
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
                # get the target site in remanei UTR
                crm_target_utr = UTR_ali[crm_utr][msa_start:msa_end]
                if crm_target_utr != target_chromo:
                    print(target_chromo)
                    print(crm_target_utr)
                    print(chromo)
                    print(gene)
                    print(site_type)
                    print(utr_file)
                    print(chromo_start, chromo_end)
                    print(orientation)
                    print(msa_start, msa_end)
                    
                # verify that target site i chromo is the same as target site on UTR
                assert crm_target_utr == target_chromo, 'target sites in UTR and chromo do not match'
                # get the target site in latens UTR
                cla_target = UTR_ali[cla_utr][msa_start: msa_end]
                # get the allele counts for the target site positions on chromo
                # check that chromo is key in chromo_sites
                if chromo in chromo_sites:
                    positions = [i for i in range(chromo_start, chromo_end) if i in chromo_sites[chromo]]
                    # do not consider target sites that do not have coverage (or minimum sample size) on all sites
                    if len(positions) == len(target_chromo):
                        # check orientation
                        if orientation == '-':
                            # reverse sort the positions to align chromo positions to the cla target
                            positions.sort()
                            positions.reverse()
                        # check that target site obtained from chromo positons correspond to crem utr target
                        check_target = ''
                        for i in positions:
                            # check orientation
                            if orientation == '+':
                                # get the reference allele in {[chromo: {site : [ref_allele, alt_allele, count_ref, count alt]}]
                                check_target += chromo_sites[chromo][i][0]
                            elif orientation == '-':
                                # take the complement of the reference allele
                                check_target += seq_complement(chromo_sites[chromo][i][0])
                        assert check_target == crm_target_utr, 'remanei target sites from chromo indices and UTR indices do not match'
                        # loop over positions
                        for i in range(len(positions)):
                            # get allele counts
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
                            # check that SNP can be polarized (ancestral allele  == ref or alt)
                            if cla_target[i] == ref or cla_target[i] == alt:
                                # SNP can be polarized
                                # check if site is polymorphic                           
                                if ref_count != 0 and alt_count != 0:
                                    # site is polymorphic, find derived and ancestral allele
                                    # check that ref and alt alleles are different
                                    assert ref != alt, 'ref and alt counts are different but alleles are the same'
                                    if cla_target[i] == ref:
                                        # reference allele is ancestral, alt is derived
                                        # compute DAF
                                        freq = alt_count / (ref_count + alt_count)
                                    elif cla_target[i] == alt:
                                        # reference allele is derived, alt is ancestral
                                        freq = ref_count / (ref_count + alt_count)
                                    # check conservation level
                                    if conservation_level == 'all':
                                        # add freq to list
                                        DAF.append(freq)
                                    elif conservation_level == 'crm':
                                        # record DAF for remanei-specific targets
                                        # check conservation level
                                        if conservation_scores[site] == 'crm':
                                            # use conservation in the crm-target-coord
                                            # not the crm-cla-target-coord because
                                            # conservation with cel is not labeled in this file
                                            # add freq to list
                                            DAF.append(freq)
                                    elif conservation_level == 'crm-cla':
                                        # record DAF for remanei-latens conserved targets
                                        # check conservation level
                                        if conservation_scores[site] == 'crm:cla':
                                            # use conservation in the crm-target-coord
                                            # not the crm-cla-target-coord because
                                            # conservation with cel is not labeled in this file
                                            # add freq to list
                                            DAF.append(freq)
                                    elif conservation_level == 'crm-cla-cel':
                                        # record DAF for site conserved in remanei, latens and elegans
                                        # need to find conservation from dict
                                        if conservation_scores[site] == 'crm:cla:cel':
                                            # use conservation in the crm-target-coord
                                            # not the crm-cla-target-coord because
                                            # conservation with cel is not labeled in this file
                                            # add freq to list
                                            DAF.append(freq)
                                        
    # close file after reading
    infile.close()
    
    return DAF

   

# use this function to get the conservation level of remanei target sites
def get_conservation_score(crm_target_sites_file, unique_transcripts):
    '''
    (file, file) -> dict
    Take the file with all remanei mirna target sites and return a dictionary 
    with site as key and conservation level as value.
    Consider sites only in valid transcripts (ie. 1 transcript per gene)
    '''
    
    # create a set of valid transcripts
    valid_transcripts = get_valid_transcripts(unique_transcripts)
    
    # create a dict to store the conservation level for each site
    # site is a string defined as 'transcript;seed;site_type;site_start_chromo;site_end_chromo;orientation'
    targets = {}
        
    # open file for reading
    infile = open(crm_target_sites_file, 'r')
    # skip header
    infile.readline()
    # loop over file
    for line in infile:
        line = line.rstrip()
        if line != '':
            line = line.split()
            # get transcript
            gene = line[0]
            # check if transcript is valid
            if gene in valid_transcripts:
                # get seed
                seed = line[1]
                # get site_type
                site_type = line[3]
                # get 0-based coordinates on chromo
                start = int(line[11]) -1
                end = int(line[12])
                orientation = line[13]
                # get conservation
                conservation = line[4]
                # build site
                site = gene + ';' + seed + ';' + site_type + ';' + str(start) + ';' + str(end) + ';' + orientation
                targets[site] = conservation
            
    # close file
    infile.close()
    
    return targets
    
        

# ise this function to get the name of the UTR alignment file of a gene if instrest     
def get_UTR_ali_file(gene, files):
    '''
    (str, list) -> str
    Take a gene of interest and a list of UTR alignments files and return the 
    file name corresponding to the alignment of the gene of interest
    Precondition: the gene name is in the file name
    '''
    # loop over file names in files
    for filename in files:
        # check if the gene is in the file name
        if gene == filename[:filename.index('_UTR')]:
            utr_file = filename
            break
    return utr_file    
    

# use this function to find the remanei and latens sequences in a dict with 2 sequences
def get_crm_cla_seqnames(UTR_ali):
    '''
    (dict) -> (str, str)
    Take a dictionary with sequence name : sequence pairs for a single remanei
    sequence and a single latens sequence and return a tuple with the remanei
    and latens squence names
    '''
    # find the remanei and latens utr
    for seqname in UTR_ali:
        if 'CRE' in seqname:
            crm_utr = seqname
        elif 'CLA' in seqname:
            cla_utr = seqname
    return crm_utr, cla_utr

# use this fucntion to get the sequence of the target sites using sites with allele counts
def get_target_seq_from_sites_with_allele_counts(chromo_sites, chromo, positions, orientation):
    '''
    (dict, str, list, str) -> str
    Take the dictionary of site positions with allele counts, the chromosome,
    a list of positions corresponding to the target sites on chromo and the orientation
    of the target on chromo and return the target site sequence
    Preconditions: Positions are 0-based, the target sites has coverage
    and minimum sample size for all its sites, the indices in positions list
    are already sorted so that target is orientated 5'-3'
    '''
       
    # initate str
    target_seq = ''
    # loop over sorted positions
    for i in positions:
        # check orientation
        if orientation == '+':
            # get the reference allele in {[chromo: {site : [ref_allele, alt_allele, count_ref, count alt]}]
            target_seq += chromo_sites[chromo][i][0]
        elif orientation == '-':
            # take the complement of the reference allele
            target_seq += seq_complement(chromo_sites[chromo][i][0])
            
    return target_seq


# use this function to check if a target site has a single SNP
def get_SNP_count(chromo_sites, chromo, positions):
    '''
    (dict, str, list) -> int
    Take the dictionary with site positions : allele counts, the chromo where
    the sequence of interest is located and a sorted list with positions of the
    sequence of interest on chromo and return the number of SNPs in the sequence
    Preconditions: positions are 0-based, the sequence has coverage and minumum sample size
    for all its sites and the list of its positions is already ordered such that
    the sequence is orientated 5'-3' (positions are in decreasing order is the sequence orientation is '-')    
    '''
    
    # set snp count
    SNPs = 0    
    # loop over positions
    for i in range(len(positions)):
        # get the ref and alt allele counts
        ref_count = chromo_sites[chromo][positions[i]][2]
        alt_count = chromo_sites[chromo][positions[i]][3]
        # check if site is polymorphic
        if ref_count != 0 and alt_count != 0:
            SNPs += 1
    return SNPs


    
 



# use this function to compute the derived allele frequency at a single snp
def compute_snp_daf(ref, alt, anc, ref_count, alt_count):
    '''
    (str, str, str, int, int) -> float
    Take the reference and alternative alleles at a snp position, the ancestral
    allele at that position and the counts of reference and alternative alleles
    and return the derived allele frequency
    Precondition: the SNP can be polarized: ancestral allele is the same as
    the reference of the alternative allele
    '''
    
    if anc == ref:
        # ref is ancestral, alt is derived
        freq = alt_count / (ref_count + alt_count)
    elif anc == alt:
        # alt is ancestral, ref is derived
        freq = ref_count / (ref_count + alt_count)
        
    return freq
    
    

# use this function to get the find the targets with a single high DAF SNP
def find_miRNA_targets_with_high_DAF(chromo_sites, genome_fasta, crm_cla_target_sites_file, UTR_alignments_folder, unique_transcripts, daf_lower, daf_upper):
    '''
    (dict, file, file, str, str, float, float) -> dict
    Take the dictionary with allele counts for sites with coverage
    (and minimum sample size), the fasta file of the genome sequence,
    the file with remanei target coordinates for sites of genes that have an
    ortholog in latens, the folder with remanei-latens UTR alignments,
    the file with valid transcripts (ie. 1 transcript per gene), 2 float numbers
    specifying the range of DAF to consider and return a dictionary with target
    gene as key and a list of list of features of all the target sites of that gene
    Precondition: consider only target sites with coverage at all positions
    '''
    
    valid_bases = {'A', 'T', 'C', 'G'}    
    
    # convert genome sequence to dict
    genome = convert_fasta(genome_fasta)     
        
    # create a list of UTR alignment files
    files = [i for i in os.listdir(UTR_alignments_folder) if 'CRE' in i and '.txt' in i]
    
    # get the set of valid transcripts
    valid_transcripts = get_valid_transcripts(unique_transcripts)    
    
    # create a dict to store the target
    # {gene : [[ref, alt, anc, daf, seed, crm_target, alt_target, cla_target, site_type, start, end, orientation]]}
    targets = {}
    
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
            # check if valid transcript
            if gene in valid_transcripts:
                # get seed
                seed = line[1]
                # get site_type
                site_type = line[3]
                # get chromo
                chromo = line[7]
                # get coordinates 0-based on chromo
                chromo_start = int(line[13]) -1
                chromo_end = int(line[14])
                orientation = line[15]
                # get coordinates in the multi-sequence UTR alignment
                # the positions correspond to + orientation
                msa_start = int(line[5]) -1
                msa_end = int(line[6])
                # get the target site on chromo
                target_chromo = genome[chromo][chromo_start:chromo_end]
                # check target site orientation on chromo
                if orientation == '-':
                    # take reverse complement of target chromo
                    target_chromo = reverse_complement(target_chromo)
                # get the target site in the aligned remanei UTR
                # find the UTR alignment file
                utr_file = get_UTR_ali_file(gene, files)                
                # convert the alignment to dict
                UTR_ali = convert_fasta(UTR_alignments_folder + utr_file)
                # find the remanei and latens utr
                crm_utr, cla_utr = get_crm_cla_seqnames(UTR_ali)
                # get the target site in remanei UTR
                crm_target_utr = UTR_ali[crm_utr][msa_start:msa_end]
                # verify that target site i chromo is the same as target site on UTR
                assert crm_target_utr == target_chromo, 'target sites in UTR and chromo do not match'
                # get the target site in latens UTR
                cla_target = UTR_ali[cla_utr][msa_start: msa_end]
                # check that both remanei and latens targets have only valid nucleotides
                if set(crm_target_utr).issubset(valid_bases) and set(cla_target).issubset(valid_bases):
                    # check that chromo is key in chromo_sites
                    if chromo in chromo_sites:
                        # get the positions of thet remanei target site on chromo
                        positions = [i for i in range(chromo_start, chromo_end) if i in chromo_sites[chromo]]
                        # do not consider target sites that do not have coverage (or minimum sample size) on all sites
                        if len(positions) == len(target_chromo):
                            # check orientation
                            if orientation == '-':
                                # reverse sort the positions to align chromo positions to the cla target
                                positions.sort()
                                positions.reverse()
                            # check that target site obtained from chromo positons correspond to crem utr target
                            check_target = get_target_seq_from_sites_with_allele_counts(chromo_sites, chromo, positions, orientation)
                            assert check_target == crm_target_utr, 'remanei target sites from chromo indices and UTR indices do not match'
                            # get SNP count in target
                            snp_count = get_SNP_count(chromo_sites, chromo, positions)
                            # keep only targets with a single snp
                            if snp_count == 1:
                                # set up boolean
                                polarized = False
                                # find position of the snp [list of snp positions]
                                snp_pos = find_SNP_positions(chromo_sites, chromo, positions)
                                # check that there is a single snp
                                assert len(snp_pos) == 1, 'snp count is not 1'
                                # get the position on chromo with the snp
                                snp_chromo_pos = positions[snp_pos[0]]
                                # get reference  and alternative alleles
                                # check orientation to get ref and alt alleles
                                if orientation == '+':
                                    ref = chromo_sites[chromo][snp_chromo_pos][0]
                                    alt = chromo_sites[chromo][snp_chromo_pos][1]
                                elif orientation == '-':
                                    # take the complement of ref and alt
                                    ref = seq_complement(chromo_sites[chromo][snp_chromo_pos][0])
                                    alt = seq_complement(chromo_sites[chromo][snp_chromo_pos][1])                                
                                # get allele counts
                                ref_count = chromo_sites[chromo][snp_chromo_pos][2]
                                alt_count = chromo_sites[chromo][snp_chromo_pos][3]
                                # find the number of differences between cla and crm target
                                target_diff = match_diff(crm_target_utr, cla_target)
                                # accept only 1 or 0 mismatches
                                if target_diff <= 1:
                                    # check that SNP can be poloraized
                                    if target_diff == 0:
                                        # check if ancestral allele is reference or alternative
                                        if cla_target[snp_pos[0]] == ref or cla_target[snp_pos[0]] == alt:
                                            polarized = True
                                    elif target_diff == 1:
                                        # check if the mismatch is at the same position of the snp
                                        subs_pos = find_mismatch_positions(crm_target_utr, cla_target)
                                        # check that there is a single mismatch
                                        assert len(subs_pos) == 1, 'more than 1 mismatch between crm and cla targets'
                                        # check that mismatch is at the SNP position
                                        if subs_pos == snp_pos:
                                            # mismatch and snp at same site
                                            # check if ancestral allele is reference or alternative
                                            if cla_target[snp_pos[0]] == ref or cla_target[snp_pos[0]] == alt:
                                                polarized = True
                                    # if snp can be polarized, compute daf
                                    if polarized == True:
                                        # get alternative target
                                        if snp_pos[0] == 0:
                                            alt_target = alt + crm_target_utr[1:]
                                        else:
                                            alt_target = crm_target_utr[:snp_pos[0]] + alt + crm_target_utr[snp_pos[0]+1:]
                                        daf = compute_snp_daf(ref, alt, cla_target[snp_pos[0]], ref_count, alt_count)
                                        # check if daf is within daf_lower and daf_upper
                                        if daf_lower <= daf <= daf_upper:
                                            # keep target site and populate dict
                                            # {gene : [ref, alt, anc, daf, seed, crm_target, cla_target, site_type, start, end, orientation]}
                                            # check if gene is key in dict
                                            if gene in targets:
                                                targets[gene].append([ref, alt, cla_target[snp_pos[0]], daf, seed, crm_target_utr, alt_target, cla_target, site_type, chromo_start, chromo_end, orientation])
                                            else:
                                                targets[gene] = [[ref, alt, cla_target[snp_pos[0]], daf, seed, crm_target_utr, alt_target, cla_target, site_type, chromo_start, chromo_end, orientation]]
    
                                       
    # close file after reading
    infile.close()
    
    return targets


    
# use this function to get the find the targets with a single high DAF SNP
def count_targets_with_single_snp(chromo_sites, genome_fasta, crm_cla_target_sites_file, unique_transcripts):
    '''
    (dict, file, file, file) -> int
    Take the dictionary with allele counts for sites with coverage
    (and minimum sample size), the fasta file of the genome sequence,
    the file with remanei target coordinates for sites of genes that have an
    ortholog in latens, the file with valid transcripts (ie. 1 transcript per gene),
    and return the number of miRNA target sites with a single SNP among sites
    in UTR of genes that have an ortholog between remanei and latens
    Precondition: consider only target sites with coverage at all positions
    '''
    
    # convert genome sequence to dict
    genome = convert_fasta(genome_fasta)     
        
    # get the set of valid transcripts
    valid_transcripts = get_valid_transcripts(unique_transcripts)    
    
    # create a dict to store the target
    # {gene : [[ref, alt, anc, daf, seed, crm_target, alt_target, cla_target, site_type, start, end, orientation]]}
    targets = 0
    
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
            # check if valid transcript
            if gene in valid_transcripts:
                # get chromo
                chromo = line[7]
                # get coordinates 0-based on chromo
                chromo_start = int(line[13]) -1
                chromo_end = int(line[14])
                orientation = line[15]
                # get the target site on chromo
                target_chromo = genome[chromo][chromo_start:chromo_end]
                # check target site orientation on chromo
                if orientation == '-':
                    # take reverse complement of target chromo
                    target_chromo = reverse_complement(target_chromo)
                if chromo in chromo_sites:
                    # get the positions of thet remanei target site on chromo
                    positions = [i for i in range(chromo_start, chromo_end) if i in chromo_sites[chromo]]
                    # do not consider target sites that do not have coverage (or minimum sample size) on all sites
                    if len(positions) == len(target_chromo):
                        # get SNP count in target
                        snp_count = get_SNP_count(chromo_sites, chromo, positions)
                        # keep only targets with a single snp
                        if snp_count == 1:
                            targets += 1
                                
    infile.close()
    return targets 
    




        
        

    
   
# use this function to the find the cla targets with a single high DAF SNP in non-target remanei
def find_non_remanei_targets_in_latens_with_high_DAF(chromo_sites, mature_fasta, genome_fasta, caeno_gff, threshold, cla_specific_target_coord_file, UTR_alignments_folder, unique_transcripts, daf_lower, daf_upper):
    '''
    (dict, file, file, file, int, file, str, str, float, float) -> dict
    Take the dictionary with allele counts for sites with coverage
    (and minimum sample size), the file with remanei miRNA seed pairs, 
    the fasta file of the genome sequence, the remanei GFF file, the threshold
    of UTR length obtained from the distribution of elegans UTRs,
    the file with latens specific miRNA target coordinates for genes that have an
    orthologs, the folder with remanei-latens UTR alignments,
    the file with valid transcripts (ie. 1 transcript per gene), 2 float numbers
    specifying the range of DAF to consider and return a dictionary with target
    gene as key and a list of list of features of all the target sites of that gene
    Precondition: consider only target sites with coverage at all positions
    '''
    
    # get the coordinates of the UTR in remanei
    # {TS1 : [chromo, start, end, orientation]}
    # coordinates are 1-based in this dict
    UTR_coord = get_three_prime_UTR_positions(caeno_gff, genome_fasta, threshold)
        
    # create a dict of seed sequence : list of mirnas sharing the seed
    seed_mature_pairs = seed_mirnas(mature_fasta) 
    # make a list of remanei seeds
    seeds = [i for i in seed_mature_pairs]
    
    valid_bases = {'A', 'T', 'C', 'G'}    
    
    # convert genome sequence to dict
    genome = convert_fasta(genome_fasta)     
        
    # create a list of UTR alignment files
    files = [i for i in os.listdir(UTR_alignments_folder) if 'CRE' in i and '.txt' in i]
    
    # get the set of valid transcripts
    valid_transcripts = get_valid_transcripts(unique_transcripts)    
    
    # create a dict to store the target
    # {gene : [[ref, alt, anc, daf, seed, crm_target, alt_target, cla_target, site_type, start, end, orientation]]}
    targets = {}
    
    # open file for reading
    infile = open(cla_specific_target_coord_file, 'r')
    # skip header
    infile.readline()
    # loop over file
    for line in infile:
        line = line.rstrip()
        if line != '':
            line = line.split()
            # get transcript
            gene = line[0]
            # check if valid transcript
            if gene in valid_transcripts:
                # get seed
                seed = line[1]
                # get site_type
                site_type = line[3]
                # get coordinates in the multi-sequence UTR alignment
                # the positions correspond to + orientation
                msa_start = int(line[5]) -1
                msa_end = int(line[6])
                # get chromo
                chromo = UTR_coord[gene][0]
                # get the orientation of the remanei UTR in genome
                orientation = line[-1]
                # find the UTR alignment file
                utr_file = get_UTR_ali_file(gene, files)
                # convert the alignment to dict
                UTR_ali = convert_fasta(UTR_alignments_folder + utr_file)
                # find the remanei and latens utr
                crm_utr, cla_utr = get_crm_cla_seqnames(UTR_ali)
                # get the latens target
                cla_target = UTR_ali[cla_utr][msa_start:msa_end]
                # get the remanei seq aligned to the latens target
                crm_target_utr = UTR_ali[crm_utr][msa_start:msa_end]
                # check that both remanei and latens targets have only valid nucleotides
                if set(crm_target_utr).issubset(valid_bases) and set(cla_target).issubset(valid_bases):
                    # check that the remanei seq is not a miRNA target
                    # loop of remanei target, check if potential crm site is target
                    for seed_seq in seeds:
                        if is_matching_seed(crm_target_utr, site_type, seed_seq) == True:
                            valid_crm_target = True
                            break
                        else:
                            valid_crm_target = False
                    # do not consider latens targets that are remanei targets
                    if valid_crm_target == False:
                        # get the coordinates of the remanei target from UTR aligned to UTR
                        start_non_aligned_utr = convert_MSA_coord_to_UTR_coord(msa_start, UTR_ali[crm_utr])
                        end_non_aligned_utr = start_non_aligned_utr + len(crm_target_utr)
                        # get the coordinates of the remanei target in genome
                        chromo_start, chromo_end = convert_site_position_UTR_to_site_position_chromo(start_non_aligned_utr, end_non_aligned_utr,
                                                                                                     UTR_coord[gene][1]-1, UTR_coord[gene][2], orientation, len(genome[chromo]))
                        # check that crm target from aligned UTR and from chromo are the same
                        crm_target_chromo = genome[chromo][chromo_start: chromo_end]
                        if orientation == '-':
                            crm_target_chromo = reverse_complement(crm_target_chromo)
                        assert crm_target_utr == crm_target_chromo, 'crm targets from aligned UTR and from chromo don\'t match'
                        # get the number of SNPs in crm
                        if chromo in chromo_sites:
                            # get the positions of crm target on chromo 
                            positions = [i for i in range(chromo_start, chromo_end) if i in chromo_sites[chromo]]
                            # do not consider target sites that do not have coverage or minimum sample size for all sites
                            if len(positions) == len(crm_target_utr):
                                # check orientation
                                if orientation == '-':
                                    # reverse sort the positions to align chromo positions to the cla target
                                    positions.sort()
                                    positions.reverse()
                                # check that target site obtained from chromo positons correspond to crem utr target
                                check_target = get_target_seq_from_sites_with_allele_counts(chromo_sites, chromo, positions, orientation)
                                assert check_target == crm_target_utr, 'remanei target sites from chromo indices and UTR indices do not match'
                                # get SNP count in target 
                                snp_count = get_SNP_count(chromo_sites, chromo, positions)
                                # keep only targets with a single snp
                                if snp_count == 1:
                                    # set up boolean
                                    polarized = False
                                    # find position of the snp [list of snp positions]
                                    snp_pos = find_SNP_positions(chromo_sites, chromo, positions)
                                    # check that there is a single snp
                                    assert len(snp_pos) == 1, 'snp count is not 1'
                                    # get the position of the SNP on chromo
                                    snp_chromo_pos = positions[snp_pos[0]]
                                    # get reference  and alternative alleles
                                    # check orientation to get ref and alt alleles
                                    if orientation == '+':
                                        ref = chromo_sites[chromo][snp_chromo_pos][0]
                                        alt = chromo_sites[chromo][snp_chromo_pos][1]
                                    elif orientation == '-':
                                        # take the complement of ref and alt
                                        ref = seq_complement(chromo_sites[chromo][snp_chromo_pos][0])
                                        alt = seq_complement(chromo_sites[chromo][snp_chromo_pos][1])                                
                                    # get allele counts
                                    ref_count = chromo_sites[chromo][snp_chromo_pos][2]
                                    alt_count = chromo_sites[chromo][snp_chromo_pos][3]
                                    # find the number of differences between cla and crm target
                                    target_diff = match_diff(crm_target_utr, cla_target)
                                    # accept only 1 or 0 mismatches
                                    if target_diff <= 1:
                                        # check that SNP can be poloraized
                                        if target_diff == 0:
                                            # check if ancestral allele is reference or alternative
                                            if cla_target[snp_pos[0]] == ref or cla_target[snp_pos[0]] == alt:
                                                polarized = True
                                        elif target_diff == 1:
                                            # check if the mismatch is at the same position of the snp
                                            subs_pos = find_mismatch_positions(crm_target_utr, cla_target)
                                            # check that there is a single mismatch
                                            assert len(subs_pos) == 1, 'more than 1 mismatch between crm and cla targets'
                                            # check that mismatch is at the SNP position
                                            if subs_pos == snp_pos:
                                                # mismatch and snp at same site
                                                # check if ancestral allele is reference or alternative
                                                if cla_target[snp_pos[0]] == ref or cla_target[snp_pos[0]] == alt:
                                                    polarized = True
                                        # if snp can be polarized, compute daf
                                        if polarized == True:
                                            # get alternative target
                                            if snp_pos[0] == 0:
                                                alt_target = alt + crm_target_utr[1:]
                                            else:
                                                alt_target = crm_target_utr[:snp_pos[0]] + alt + crm_target_utr[snp_pos[0]+1:]
                                            daf = compute_snp_daf(ref, alt, cla_target[snp_pos[0]], ref_count, alt_count)
                                            # check if daf is within daf_lower and daf_upper
                                            if daf_lower <= daf <= daf_upper:
                                                # keep target site and populate dict
                                                # {gene : [ref, alt, anc, daf, seed, crm_target, cla_target, site_type, start, end, orientation]}
                                                # check if gene is key in dict
                                                if gene in targets:
                                                    targets[gene].append([ref, alt, cla_target[snp_pos[0]], daf, seed, crm_target_utr, alt_target, cla_target, site_type, chromo_start, chromo_end, orientation])
                                                else:
                                                    targets[gene] = [[ref, alt, cla_target[snp_pos[0]], daf, seed, crm_target_utr, alt_target, cla_target, site_type, chromo_start, chromo_end, orientation]]
    
    # close file after reading
    infile.close()
    
    return targets
    
    
    
    
# use this function to transform site info for a 8mer site into a 7mer-m8 site    
def from_8mer_to_heptamer(site_info):
    '''
    (list) -> list
    Take the list of features 8mer a target site with high DAF SNP and return
    a modified list with features corresponding to a 7mer-m8 target site if
    ref and alt target have a SNP or return an empty list if ref target and alt
    target are the same (ie. the SNP was located in the 'A' of the 8mer-1a)
    Precondition: target site is 8mer
    '''
    
    # check that site is 8mer
    assert '8mer' in site_info[8], 'target site is not 8mer'
    
    # convert ref, alt and ancestral target sites from 8mer to 7mer-m8
    for i in range(5,8):
        site_info[i] = site_info[i][:-1]
    # convert target site type
    site_info[8] = '7mer-m8'        
    # adjust coordinates
    # check orientation
    if site_info[-1] == -1:
        site_info[9] = site_info[9] -1
    elif site_info[-1] == '+':
        site_info[10] = site_info[10] -1
    # check if there is a SNP between alt and ref targets
    if site_info[5] == site_info[6]:
        # no SNP after removing the A in the 8mer-1a
        return []
    else:
        # SNP is within the 7mer-m8
        return site_info
    
       
# use this function to combine targetscan 7mer-m8 with heptamers from transformed 8mer-1a
def merge_8mer_with_7merm8(target_high_DAF_SNPs):
    '''
    (dict) -> dict
    Take the dictionary with gene: list of high DAF SNPs and return a modified 
    dictionary, taking into account only 7mer-m8 target sites and heptamer
    seed matching sites transformed from 8mer sites to 7mer-m8 sites
    '''
    Gstr = lambda x:  str(x)
    
    # remove any 7mer-1a site    
    # loop over genes
    for gene in target_high_DAF_SNPs:
        # create a list of sites to remove
        to_remove = []
        # loop over sites for the given target gene
        for site in target_high_DAF_SNPs[gene]:
            if site[8] == '7mer-1a':
                to_remove.append(site)
        # loop over sites to remove
        for site in to_remove:
            target_high_DAF_SNPs[gene].remove(site)
    
    # remove gene without target sites
    to_remove = []
    for gene in target_high_DAF_SNPs:
        if len(target_high_DAF_SNPs[gene]) == 0:
            to_remove.append(gene)
    for gene in to_remove:
        del target_high_DAF_SNPs[gene]
    
    # loop over genes
    for gene in target_high_DAF_SNPs:
        # loop over site
        for i in range(len(target_high_DAF_SNPs[gene])):
            # check if site is 8mer
            if '8mer' in target_high_DAF_SNPs[gene][i][8]:
                # convert site features from 8mer to 7mer-m8
                target_high_DAF_SNPs[gene][i] = from_8mer_to_heptamer(target_high_DAF_SNPs[gene][i])
    # remove empty lists resulting from the conversion of 8mer to 7mer-m8 when SNPs were removed
    for gene in target_high_DAF_SNPs:
        to_remove = []
        for site in target_high_DAF_SNPs[gene]:
            if site  == []:
                to_remove.append(site)
        for site in to_remove:
            target_high_DAF_SNPs[gene].remove(site)
            
    # remove doublons resultins from 8mer to 7mer-m8 conversion
    for gene in target_high_DAF_SNPs:
        # make a string with each item in list
        for i in range(len(target_high_DAF_SNPs[gene])):
            target_high_DAF_SNPs[gene][i] = ';'.join(list(map(Gstr, target_high_DAF_SNPs[gene][i])))
        # convert list to set
        target_high_DAF_SNPs[gene] = set(target_high_DAF_SNPs[gene])
        # convert back to list
        target_high_DAF_SNPs[gene] = list(target_high_DAF_SNPs[gene])
        # split str
        for i in range(len(target_high_DAF_SNPs[gene])):
            target_high_DAF_SNPs[gene][i] = target_high_DAF_SNPs[gene][i].split(';')

    # convert back str to nums
    for gene in target_high_DAF_SNPs:
        # loop over sites
        for i in range(len(target_high_DAF_SNPs[gene])):
            target_high_DAF_SNPs[gene][i][3] = float(target_high_DAF_SNPs[gene][i][3])
            target_high_DAF_SNPs[gene][i][9] = int(target_high_DAF_SNPs[gene][i][9])
            target_high_DAF_SNPs[gene][i][10] = int(target_high_DAF_SNPs[gene][i][10])

    return target_high_DAF_SNPs
        
        
    
# use this function to infer target site loss and gain among sites with high DAF SNPs    
def infer_target_gain_loss(crm_target_DAF, cla_target_DAF, mature_seeds):
    '''
    (dict, dict, list) -> dict, dict, dict
    Take the dictionary of remanei target sites with high DAF, and the dictionary
    of latens target sites with high DAF, the list of remanei seed sequences
    and return 3 dictionaries with target gene as key: one dictionary when SNPs
    lead to site loss, one dictionary when SNPs lead to site gain from a
    non-target site and one dictionary when SNP to site gain from a target site
    Precondition: all sites are 7mer-m8 heptamer seed match
    '''
    
    
    # create dicts
    gain, loss, site_to_site = {}, {}, {}
    
    # loop over genes in crm dict
    for gene in crm_target_DAF:
        # loop over target sites for the given gene
        for site in crm_target_DAF[gene]:
            # get ref, alt, ancestral alleles
            ref, alt, ancestral = site[0], site[1], site[2]
            crm_target_ref, crm_target_alt, cla_target = site[5], site[6], site[7]
            # get site type
            site_type = site[8]
            assert site_type == '7mer-m8', 'site type is different than 7mer-m8'
            # find derved allele
            if ref == ancestral:
                # alt is derived, check if target_alt match a seed or not
                for seed_seq in mature_seeds:
                    new_site = is_matching_seed(crm_target_alt, '7mer-m8', seed_seq)
                    if new_site == True:
                        break
                if new_site == True:
                    # SNP created a new target site, populate dict
                    if gene in site_to_site:
                        site_to_site[gene].append(site)
                    else:
                        site_to_site[gene] = [site]
                elif new_site == False:
                    # SNP destroyed target site, populate dict
                    if gene in loss:
                        loss[gene].append(site)
                    else:
                        loss[gene] = [site]
            elif alt == ancestral:
                # ref is derived, site gain
                # check if target_alt match a seed or not
                for seed_seq in mature_seeds:
                    new_site = is_matching_seed(crm_target_alt, '7mer-m8', seed_seq)
                    if new_site == True:
                        break
                if new_site == True:
                    # SNP create a new site from a site, populate dict
                    if gene in site_to_site:
                        site_to_site[gene].append(site)
                    else:
                        site_to_site[gene] = [site]
                elif new_site == False:
                    # SNP created a site de novo, populate dict
                    if gene in gain:
                        gain[gene].append(site)
                    else:
                        gain[gene] = [site]
                        
    # loop over genes in cla_target
    for gene in cla_target_DAF:
        # loop over the target sites for the given target gene
        for site in cla_target_DAF[gene]:
            # get ref, alt, ancestral alleles
            ref, alt, ancestral = site[0], site[1], site[2]
            crm_target_ref, crm_target_alt, cla_target = site[5], site[6], site[7]
            # get site type
            site_type = site[8]
            assert site_type == '7mer-m8', 'site type is different than 7mer-m8'
            # find derived allele
            if ref == ancestral:
                # alt is derived, check if alt is site or not
                for seed_seq in mature_seeds:
                    new_site = is_matching_seed(crm_target_alt, '7mer-m8', seed_seq)
                    if new_site == True:
                        break
                if new_site == True:
                    # SNP created a site from a site, populate dict
                    if gene in site_to_site:
                        site_to_site[gene].append(site)
                    else:
                        site_to_site[gene] = [site]
                elif new_site == False:
                    # SNP destroyed a site, populate dict
                    if gene in loss:
                        loss[gene].append(site)
                    else:
                        loss[gene] = [site]
            elif alt == ancestral:
                # alt is target, ref is derived, check if ref is target
                for seed_seq in mature_seeds:
                    new_site = is_matching_seed(crm_target_ref, '7mer-m8', seed_seq)
                    if new_site == True:
                        break
                if new_site == True:
                    # SNP created a site from site, populate dict
                    if gene in site_to_site:
                        site_to_site[gene].append(site)
                    else:
                        site_to_site[gene] = [site]
                elif new_site == False:
                    # SNP destroyed a target site, populate dict
                    if gene in loss:
                        loss[gene].append(site)
                    else:
                        loss[gene] = [site]
                     
                
            
    return loss, gain, site_to_site            
        

    
if __name__ == '__main__':
    from genomic_coordinates import *
    from manipulate_sequences import *
    from divergence import *
    from parse_targetscan_output import *
    import os