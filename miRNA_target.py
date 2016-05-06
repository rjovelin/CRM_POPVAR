# -*- coding: utf-8 -*-
"""
Created on Tue May 26 11:09:39 2015

@author: Richard
"""

#!/usr/bin/env python3


from Cel_UTR import *
from accessories import *
from piRNAs import *
from miRfam_targetscan_input import *
from parse_targetscan_output import *
import os

   
# use this function to get the coorinates of the annotated 3'UTR
def grab_annotated_three_prime_coordinates(caeno_gff):
    '''
    (file) -> dict
    Returns a dictionnary with the 3' UTR coordinates of each transcript
    '''

    # create a dictionnary to stote the coordinates {transcript_name : [chromo, start, end, orienation]}
    annotated_three_prime = {}

    # open file for reading
    gff = open(caeno_gff, 'r')
    for line in gff:
        line = line.rstrip()
        if line != '':
            line = line.split()
            if len(line) > 8:
                if line[2] == 'three_prime_UTR':
                    transcript = line[8][line[8].index('ID=')+3 : line[8].index(':')]
                    chromo = line[0]
                    start = int(line[3])
                    end = int(line[4])
                    orientation = line[6]
                    annotated_three_prime[transcript] = [chromo, start, end, orientation]
    gff.close()
    return annotated_three_prime


# use this function to get the coodinates of all transcripts
def get_genes_coordinates(caeno_gff):
    '''
    (file) -> dict
    Returns a dictionnary with the transcripts as key and a list of coordinates [chromosome, start, end, orientarion] as value
    '''

    # create a dictionnary to store the gene coordinates {gene1 : [chromo, start, end, sense]}
    ts_positions = {}

    # open annotation file for reading
    gff = open(caeno_gff, 'r')
    for line in gff:
        line = line.rstrip()
        if line != '':
            line = line.split()
            if len(line) > 8:
                if line[2] == 'mRNA':
                    chromo = line[0]
                    start = int(line[3])
                    end = int(line[4])
                    sense = line[6]
                    transcript = line[8][line[8].index('ID=')+3 : line[8].index(';')]
                    ts_positions[transcript] = [chromo, start, end, sense]
    gff.close()
    return ts_positions


# use this function to create a dictionnary of transcript features that can
# be searched by chromosome and then by coordinates
def chromo_positions(caeno_gff):
    '''
    (file) -> dict
    Return a dictionnary of scaffold/chromosome as key and dictionnaries
    with a tuple of gene coordinates as key and a list of features as value    
    '''
    
    # get the transcripts coordinates
    ts_positions = get_genes_coordinates(caeno_gff)
    
    # create a dict
    # {chromo : {(start, end): [transcript, start, end, orientation]}}
    LG_pos = {}
    
    # warning: some transcripts have the same coordinates
    # need to grab the same transcripts when inverting the dict
    
    # warning: to make this function reusable, create a sorted list of genes
    # from ts_positions to loop over 
    proteins = [gene for gene in ts_positions]
    proteins.sort()
    # reverse so that genes with CRE_ID are looped over first
    proteins.reverse()
    
    # create dict of chromo as key and inner dicts of coordinates and list
    # of transcripts with identical positions on chromo
    # {chromo: {coordinates: [transcript1, transcript2]}}
    coord = {}
       
    for gene in proteins:
        # get positions  
        positions = (ts_positions[gene][1], ts_positions[gene][2])
        # get chromo
        chromo = ts_positions[gene][0]
        # populate dict
        if chromo in coord:
            if positions in coord[chromo]:
                coord[chromo][positions].append(gene)
            else:
                coord[chromo][positions] = [gene]
        else:
            coord[chromo] = {}
            coord[chromo][positions] = [gene]
    
    # create a list of transcripts with unique positions on chromosomes
    unique_genes = []
    # loop over chromo in coord
    for chromo in coord:
        # loop over positions in chromo
        for position in coord[chromo]:
            # get the first gene of the list
            unique_genes.append(coord[chromo][position][0])
    # sort list 
    unique_genes.sort()
    
    # loop over transcripts with unique positions
    for gene in unique_genes:
        # get features
        LG = ts_positions[gene][0]
        start = ts_positions[gene][1]
        end = ts_positions[gene][2]
        orientation = ts_positions[gene][-1]
        # populate dict
        if LG in LG_pos:
            # if key is main dict, add key: value pair for inner dict
            LG_pos[LG][(start, end)] = [gene, start, end, orientation]
        else:
            # if key not in main dict, initialize inner dict
            LG_pos[LG] = {}
            LG_pos[LG][(start, end)] = [gene, start, end, orientation]
            
    return LG_pos

# use this function to find the next downstream gene for a given transcript
def find_next_downstream_gene(transcripts_genes, LG_pos, chromosome, position, coordinates):
    '''
    (dict, dict, str, int, list) -> int
    For a given given gene at index position in the list coordinates,
    return the index of the next downstream gene in the list
    '''
    
    # transcript_genes is a dict of transcript : gene pairs
    # LG_pos is a dict {chromo : {(start, end): [transcript, start, end, orientation]}}
    
    # get the orientation of focal transcript
    orientation = LG_pos[chromosome][coordinates[position]][-1]
    # get the transcript name of the focal transcript
    transcript_name = LG_pos[chromosome][coordinates[position]][0]
    # get the parent gene of the focal transcript
    gene = transcripts_genes[transcript_name]
    
    # set up boolean
    next_gene_not_found = True
    
    # if orientation == +, loop over coordinates starting at position and find next gene
    if orientation == '+':
        i = position + 1
        
        # loop over coordinates starting at position until all transcripts
        # have been evaluated until next gene is found
        while i != len(coordinates) and next_gene_not_found:
            # check if transcripts are from the same gene
            transcript_ID = LG_pos[chromosome][coordinates[i]][0]
            gene_2 = transcripts_genes[transcript_ID]
            # if genes are different, return the index of the next gene
            if gene != gene_2:
                next_gene_not_found = False
                return i
            else:
                # keep looking
                i += 1
        
        # if all the positions in coordinates have been evaluated, then all
        # transcripts after focal transcript are from the same gene
        if i == len(coordinates) and next_gene_not_found == True:
            return None 
    
    # if orienttion  == -, loop in reverse order
    elif orientation == '-':
        i = position - 1
        # loop over coordinates starting at position until all transcripts
        # have been evaluated until next gene is found
        while i != -1 and next_gene_not_found == True:
            # check if transcripts are from the same gene
            transcript_ID = LG_pos[chromosome][coordinates[i]][0]
            gene_2 = transcripts_genes[transcript_ID]
            # if genes are different, return the index of the next gene
            if gene != gene_2:
                next_gene_not_found = False
                return i
            else:
                # keep looking
                i -= 1
        
        # if all the positions in coordinates have been evaluated, then all
        # transcripts before focal transcript are from the same gene
        if i == -1 and next_gene_not_found == True:
            return None



# use this function to get the start and end positions of downstream sequence
# for a single transcript with no downstream distinct gene
def get_UTR_coord_border_chromo(orientation, coordinates, k, threshold, genome, chromo):
    '''
    (str, list, int, int, dict, str) -> (int, int)
    Return a tuple with the start and end position of the downstream sequence
    of the gene at position k in coordinates with given orientation, using the
    threshold length of sequence to extract and a dict of chromo: sequence pairs
    Precondition: coordinates is sorted
    '''
    
    # orientation is sense of transcript on chromosome: '+' or '-'
    # coordinates is a list of tuple with (start, end) positions for all transcripts
    # threshold is the maximum UTR length to extract
    # genome is a dict with chromo as key and sequence as value
    
    
    if orientation == '+':
        # check if len(downstream) > len(chromo), take sequence up to end of chromo
        if coordinates[k][1] + threshold >= len(genome[chromo]):
            # check if transcripts end at the end of chromo
            if coordinates[k][1] == len(genome[chromo]):
                # start = end, no downstream sequence
                start = coordinates[k][1]
                end = coordinates[k][1]
            else:
                # start = end of transcript + 1, end = end of chromo
                start = coordinates[k][1] + 1
                end = len(genome[chromo])
        # check if len(downstream) < len(chromo), take sequence of length theshold                 
        elif coordinates[k][1] + threshold < len(genome[chromo]):
            # start = end of transcript + 1, end = end of transcript + threshold
            start = coordinates[k][1] + 1
            end = start + threshold -1
    # if '-', compare length of downstream sequence to start of transcript
    elif orientation == '-':
        # check if distance from start of chromo to start of transcript < threshold
        if coordinates[k][0] <= threshold:
            # start = start of chromosome
            start  = 1
            # check if transcript is at begining of chromosome
            if coordinates[k][0] -1 > 1:
                # downstream sequence exists, end = start of transcript - 1
                end = coordinates[k][0] - 1
            else:
                # no downstream sequence
                end = 1
        # check if distance from start of chromo to start of transcript > threshold                 
        elif coordinates[k][0] > threshold:
            # end = start of transcript - 1, start = end - threshold + 1
            end = coordinates[k][0] - 1
            start = end - threshold + 1
    
    return start, end 
    


# use this function to get the (start, end) positions of downstream sequence
# of a gene that has a neighboring gene
def get_UTR_coord_between_neighbors(orientation, coordinates, k, m, threshold):
    '''
    (str, list, int, int, int) -> (int, int)
    Return a tuple with the start and end position of the downstream sequence
    of the gene at position k in coordinates with given orientation, by comparing the 
    positions of gene k with the positions of its neighbor at position m in 
    coordinates, using the threshold length of sequence to extract
    Precondition: coordinates is sorted
    '''
    
    # orientation is sense of transcript on chromosome: '+' or '-'
    # coordinates is a list of tuple with (start, end) positions for all transcripts
    # threshold is the maximum UTR length to extract
    # genome is a dict with chromo as key and sequence as value
    # gene at position k is the focal gene, gene at position m is the next downstream gene
    
    # assume genes do not overlap
    overlap = False
    
    # check if genes overlap
    if k < m:
        # downstream is after in the list, k orientation is '+'
        assert orientation == '+'
        # check if gene k and gene m overlap
        if coordinates[k][1] > coordinates[m][0]:
            # genes overlap, start of m is within gene k
            overlap = True
    elif k > m:
        # downstream is before in the list, k orientation is '-'
        assert orientation == '-'
        # check if gene k and gene m overlap
        if coordinates[k][0] < coordinates[m][1]:
            # genes overlap, start of k is within gene m
            overlap = True
    
    # if gene overlap, no downstream sequence
    if overlap == True:
        start = coordinates[k][1]
        end = coordinates[k][1]
    # if gene do not overlap, compare the distance between genes
    elif overlap == False:
        # check if the downstream gene is before or after gene k in the list
        if k < m:
            # compare the distance between gene k and gene m to threshold
            distance = coordinates[m][0] - coordinates[k][1] - 1
            if distance <= threshold:
                # take entire distance
                # check that distance is longer than 1 nucleotide
                if distance > 1:
                    start = coordinates[k][1] + 1
                    end = coordinates[m][0] -1
                else:
                    # no downstream sequence, genes are adjacent
                    start = coordinates[k][1]
                    end = coordinates[k][1]
            elif distance > threshold:
                # take downstream sequence of length threshold
                start = coordinates[k][1] + 1
                end = start + threshold - 1
        # check if downstream is before of after gene k in the list        
        elif k > m:
            # compare the distance between gene k and gene m to threshold
            distance = coordinates[k][0] - coordinates[m][1] - 1
            if distance <= threshold:
                # take entire distance
                # check that distance is longer than 1 nucleotide
                if distance > 1:
                    start = coordinates[m][1] + 1
                    end = coordinates[k][0] - 1
                else:
                    # no downstream sequence, genes are adjacent
                    start = coordinates[k][0]
                    end = coordinates[k][0]
            elif distance > threshold:
                # take downstream sequence of length threshold
                end = coordinates[k][0] - 1
                start = end - threshold + 1
    
    return start, end
        
   
            
# use this function to get the coordinates of the UTR/downstream sequences
# of each transcript
def get_three_prime_UTR_positions(caeno_gff, assembly, threshold):
    '''
    (file) -> dict
    Returns a dictionnary with a list containing the positions of the annotated
    3' UTRs or the downstream sequences up to length threshold, its orientation
    and the parent gene for each transcript in the gff file
    Precondition: Coordinates are start = 0 and end = 0 for overlapping genes
    '''

    # convert the fasta assembly to a dictionnary with chromo / scaffold as key and sequence as value
    genome = convert_fasta(assembly)
        
    # create a dict for each transcript, searchable by chromosome and then
    # by coordinates {chromo: {(start, end): [transcript, start, end, orienation]}}
    LG_pos = chromo_positions(caeno_gff)

    # get the transcripts : gene names pairs
    transcripts_genes = transcript_to_gene(caeno_gff)
    
    # get the coordinates of the annotated 3' UTR {TS : [chromo, start, end, orienation]}
    annotated_three_prime = grab_annotated_three_prime_coordinates(caeno_gff)
            
    # create a dictionnary to store the UTR/downstream positions of each transcript
    # {TS1 : [chromo, start, end, orientation]}
    three_prime = {}

    # go through each chromosome that have transcripts
    for chromo in LG_pos:
        # make a list of (start, end) coordinates for each transcript on chromo
        coordinates = [i for i in LG_pos[chromo]]
        # sort coordinates
        coordinates.sort()
        # check if more than 1 gene on chromo
        if len(coordinates) == 1:
            # check orientation
            transcript = LG_pos[chromo][coordinates[0]][0]
            orientation = LG_pos[chromo][coordinates[0]][-1]
            # get the start, end positions
            start, end = get_UTR_coord_border_chromo(orientation, coordinates, 0, threshold, genome, chromo)
            # check if the transcript has annotated UTR
            if transcript in annotated_three_prime:
                # use the coordinates of the annotated 3' UTR
                three_prime[transcript] = list(annotated_three_prime[transcript])
            elif transcript not in annotated_three_prime:
                # use the coordinates of the downstream sequence
                three_prime[transcript] = [chromo, start, end, orientation]
                           
        # if more than 1 transcript on chromosome:
        elif len(coordinates) > 1:
            # go through the genes in coordinates
            # take the downstream sequence of each transcript up to the next gene
            # ie. downstream seq of one transcript may include coding sequence
            # of another transcript of the same gene
            for i in range(len(coordinates)):
                # get transcript name
                transcript = LG_pos[chromo][coordinates[i]][0]
                # get orientation
                orientation = LG_pos[chromo][coordinates[i]][-1]
                # find the index of the next distinct neighboring gene
                j = find_next_downstream_gene(transcripts_genes, LG_pos, chromo, i, coordinates)
                # check if transcript has a downstream gene
                if j == None:
                    # no downstream gene
                    # transcript is first or last on chromo
                    # or downstream transcripts belong to the same gene
                    # get start end positions
                    start, end = get_UTR_coord_border_chromo(orientation, coordinates, i, threshold, genome, chromo)
                elif str(j).isdigit():
                    # downstream gene exists, get start and end positions
                    start, end = get_UTR_coord_between_neighbors(orientation, coordinates, i, j, threshold)
                # check if transcript is annonated in three_prime
                if transcript in annotated_three_prime:
                    # use the coordinates of the annotated 3' UTR
                    three_prime[transcript] = list(annotated_three_prime[transcript])
                elif transcript not in annotated_three_prime:
                    # use the coordinates of the downstream sequence
                    three_prime[transcript] = [chromo, start, end, orientation]


    return three_prime                        
                        
                        
# use this function to get the UTR/dowsntream sequences of each transcript
def fetch_UTR_sequences(caeno_gff, assembly, threshold):
    '''
    (file, file, int) -> dict
    Returns a dictionnary with the transcript name as key and the UTR sequence
    from the assembly as value. Ignore transcripts lacking a UTR sequence
    (e.g. at the ends of chromosome or when adjacent to another gene)
    '''

    # convert assembly to fasta format
    genome = convert_fasta(assembly)
    
    # create a dictionnary with UTR coordinates
    three_prime = get_three_prime_UTR_positions(caeno_gff, assembly, threshold)
    
    # create a dict {transcript: sequence}
    UTR = {}
    # loop over transcripts
    for transcript in three_prime:
        # ignore transcripts that do not have a UTR (start = end)
        if not three_prime[transcript][1] == three_prime[transcript][2]:
            # get chromo
            chromo = three_prime[transcript][0]
            # get orientation
            orientation = three_prime[transcript][-1]
            # convert start to 0-based index (start  = start -1)
            start = three_prime[transcript][1] - 1
            # no need to convert end because end non-inclusive
            end = three_prime[transcript][2]
            # slice the squence
            UTR_seq = genome[chromo][start:end]
            # check orientation
            if orientation == '+':
                # populate dict
                UTR[transcript] = UTR_seq.upper()
            elif orientation == '-':
                # take reverse complement
                UTR_seq = reverse_complement(UTR_seq)
                # populate dict
                UTR[transcript] = UTR_seq.upper()
                
    return UTR
            

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
                
         
# use this function to generate the TargetScan input sequence file for cremanei-celegans analysis
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
    Take the file with mirna coordinates and return a dictionary with chromo
    as key and a list of list coordinates for each mirna on that chromo
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
    

#use this function to get the coordinates of mature miRNAs on each chromo
def get_mature_loci(mature_coord_file):
    '''
    (file) -> dict
    Take the file with genomic coordinates of mature miRNAs and return a 
    dictionary with chromo as key and a list of list coordinates for each mature 
    miR of that chrom
    '''
    
    # create a dict {chromo: [[start, end orientation]]}
    mature_coord = {}
    
    # open file for reading
    infile = open(mature_coord_file, 'r')
    # skip header
    infile.readline()
    infile.readline()
    # loop over file
    for line in infile:
        line = line.rstrip()
        if line != '':
            line = line.split()
            # get chromo
            chromo = line[2]
            # get start, end position in 0-based
            start = int(line[3]) - 1
            end = int(line[4])
            # get orientation
            orientation = line[-1]
            # check if chromo in dict
            if chromo in mature_coord:
                mature_coord[chromo].append([start, end, orientation])
            else:
                mature_coord[chromo] = [[start, end, orientation]]
    # close file after reading
    infile.close()
    
    return mature_coord
    

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

# use this function to find the indices where 2 sequences have a difference
def find_mismatch_positions(seq1, seq2):
    '''
    (str, str) -> list
    Take 2 sequences and return a list of indices where the 2 sequences differ
    Preconditions: the 2 sequences are aligned and/or have the same length
    '''
    
    # create a list to store the indices
    pos = []
    # loop over seq1, compare each position to seq2
    for i in range(len(seq1)):
        # compare seq1 and se2
        if seq1.upper()[i] != seq2.upper()[i]:
            pos.append(i)
            
    return pos
    
 
# use this function to find the indices of a SNP
def find_SNP_positions(chromo_sites, chromo, positions):
    '''
    (dict, str, list) -> list
    Take the dictionary with site positions : allele counts, the chromo where
    the sequence of interest is located and a sorted list with positions of the
    sequence of interest on chromo and return a list with the indices of the
    SNPs in the sequence
    Precondition: positions are 0-based, the sequence has coverage and
    minumum sample size for all its sites and the list of its positions is
    already ordered such that the sequence is orientated 5'-3'
    (positions are in decreasing order is the sequence orientation is '-')    
    '''
    
    # create a list to store the indisces of the snps
    SNPs = []    
    # loop over positions
    for i in range(len(positions)):
        # get the ref and alt allele counts
        ref_count = chromo_sites[chromo][positions[i]][2]
        alt_count = chromo_sites[chromo][positions[i]][3]
        # check if site is polymorphic
        if ref_count != 0 and alt_count != 0:
            # site is polymorphic, add i to list
            SNPs.append(i)
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
        
        
# use this function to get the start position of a target site in the non-aligned UTR 
# from its start position in the aligned UTR
def convert_MSA_coord_to_UTR_coord(start, aligned_UTR_seq):
    '''
    (int, str) -> int
    Take a the start position of a target site in the the aligned UTR sequence
    return its start position in the non-aligned UTR (ie. same sequence, gapped removed) 
    Precondition: positions are 0-based
    '''
    
    gaps = 0
    for i in range(len(aligned_UTR_seq[:start])):
        if aligned_UTR_seq[i] == '-':
            gaps += 1
    
    return start - gaps
            
    
   
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
        

    
