# -*- coding: utf-8 -*-
"""
Created on Fri May  6 11:43:14 2016

@author: RJovelin
"""

from Manipulate_Sequences import *
from mirna_targets import *


# use this function to map transcript name to gene names
def transcript_to_gene(caeno_gff):
    '''
    (file) -> dict
    Returns a dictionnary with transcript : gene pairs from the gff annotation file
    '''
    #create a dictionnary of transcript : gene pairs
    transcripts_genes = {}

    # open file for reading
    gff = open(caeno_gff, 'r')
    for line in gff:
        line = line.rstrip()
        if line != '':
            line = line.split()
            if len(line) > 8:
                if line[2] == 'mRNA':
                    transcript = line[8][line[8].index('ID=')+3 : line[8].index(';')]
                    gene = line[8][line[8].index('Parent')+7 : line[8].index(';', line[8].index('Parent'))]
                    transcripts_genes[transcript] = gene
    gff.close()
    return transcripts_genes


# use this function to map all transcript names to gene names
def gene_to_transcripts(caeno_gff):
    '''
    (file) -> dict
    Returns a dictionnary with gene as key and a list of transcripts as value
    '''

    # get the dictionnary of transcripts : gene names pairs
    transcripts_genes = transcript_to_gene(caeno_gff)

    # create a reverse dictionnary of gene : [transcripts] pairs
    genes = {}
    for transcript in transcripts_genes:
        gene_name = transcripts_genes[transcript]
        if gene_name in genes:
            genes[gene_name].append(transcript)
        else:
            genes[gene_name] = [transcript]

    return genes


# remanei and elegans GFF have different formats
# use this function to create a dict of transcript : gene pairs
def celegans_transcript_to_gene(celegans_gff):
    '''
    (file) -> dict
    Returns a dictionnary with celegans transcript : gene pairs from the gff annotation file
    '''
    #create a dictionnary of transcript : gene pairs
    transcripts_genes = {}

    # open file for reading
    gff = open(celegans_gff, 'r')
    for line in gff:
        line = line.rstrip()
        if line != '':
            line = line.split()
            if line[1] == 'WormBase':
                if line[2] == 'mRNA':
                    transcript = line[8][line[8].index('Transcript:')+11 : line[8].index(';')]
                    gene = line[8][line[8].index('Parent=Gene:')+12 : line[8].index(';', line[8].index('Parent'))]
                    transcripts_genes[transcript] = gene
    gff.close()
    return transcripts_genes


# remanei and elegans GFF have different formats
# use this function to create a dict if gene : list of rsnacripts airs
def celegans_gene_to_transcripts(celegans_gff):
    '''
    (file) -> dict
    Returns a dictionnary with celegans gene as key and a list of transcripts as value
    '''

    # get the dictionnary of transcripts : gene names pairs
    transcripts_genes = celegans_transcript_to_gene(celegans_gff)

    # create a reverse dictionnary of gene : [transcripts] pairs
    genes = {}
    for transcript in transcripts_genes:
        gene_name = transcripts_genes[transcript]
        if gene_name in genes:
            genes[gene_name].append(transcript)
        else:
            genes[gene_name] = [transcript]

    return genes



# use this function to generate a set of genes with indels to exclude from analysis
def get_genes_with_indels(indel_transcripts):
    '''
    (file) -> set
    Take a file with transcrtipts that have an indel in their CDS and
    a return a set of such transcripts
    '''
    
    # create set of genes
    indels = set()
    
    # open file for reading
    infile = open(indel_transcripts, 'r')
    # loop over file
    for line in infile:
        line = line.rstrip()
        if line != '':
            indels.add(line)
            
    # close file
    infile.close()
    
    return indels
    

# use this function to generate a set of transcripts to use for analyses of diversity
def get_valid_transcripts(unique_transcripts):
    '''
    (file) -> set
    Take the file with transcript names and return a set of transcripts
    that each have a unique parent gene to use in analyses of diversity
    '''
    
    # open file for reading
    infile = open(unique_transcripts, 'r')
    # make a set of transcript names
    transcripts = set()
    for line in infile:
        line = line.rstrip()
        transcripts.add(line)
    
    # close file
    infile.close()
    
    return transcripts
    

# use this function to map the remanei chromosomes to elegans chromosomes
def remanei_elegans_chromosomes(chromosome_file):
    '''
    (file) -> dict
    Take the chromosome_file and return a dictionnary of remanei sacffolds/
    linkage groups as key and the corresponding elegans chromosome as value    
    '''
    
    # create dict {rem_LG: elegans_LG}
    LG = {}
    
    # open file for reading
    infile = open(chromosome_file, 'r')
    # skip header
    infile.readline()
    # loop over file, populate dict
    for line in infile:
        line = line.rstrip()
        if line != '':
            line = line.split(',')
            rem_chromo = line[0]
            elegans_chromo = line[1]
            LG[rem_chromo] = elegans_chromo
    # close file after reading
    infile.close()
    
    return LG

    
# use this function to partition a collection of genes into lists of X-linked and autosomal genes 
def X_autosomal_genes(genes, caeno_gff, chromosome_file):
    '''
    (collection, file, file) -> (list, list)
    Take a collection of genes (set or list), the GFF file, and a file with
    remanei and elegans corresponding chromosomes and return a tuple with 
    a list of X-linked genes and a list of autosomal genes    
    '''
    
    # create a dict of remanei : elegans chromosomes
    chromosomes = remanei_elegans_chromosomes(chromosome_file)    
    
    # get the genes coordinates {gene1 : [chromo, start, end, sense]}
    genes_coordinates = get_genes_coordinates(caeno_gff)
    
    # create lists of genes
    X_genes, autosomal_genes = [], []
    
    # loop over genes
    for gene in genes:
        # get chromo
        chromo = genes_coordinates[gene][0]
        # check that chromo has a elegans equivalent
        if chromo in chromosomes:
            # get the elegans chromosome
            if chromosomes[chromo] == 'X':
                # chromo is X chromosome
                X_genes.append(gene)
            else:
                # gene is autosomal
                autosomal_genes.append(gene)
                
    return X_genes, autosomal_genes
    
    
# use this function to get the length of all annotated 3' UTRs in C. elegans    
def celegans_three_prime_UTR_length(celegans_gff):
    '''
    (file) -> list
    Returns a list with the length of the annotated 3' UTRs in C. elegans gff annotation file
    '''
    #create list to store UTR length
    UTR = []
    # opengff file for reading
    cel = open(celegans_gff, 'r')
    # go through the file, extract UTR length and store in list
    for line in cel:
        line = line.rstrip()
        if line != '':
            line = line.split('\t')
            if len(line) >= 5:
                if line[1] == 'WormBase':
                    if line[2] == 'three_prime_UTR':
                        UTR_length = (int(line[4]) - int(line[3])) + 1
                        UTR.append(UTR_length)
    cel.close()
    return UTR

# use this function to get the value corresponding to percentile
def get_percentile(L, percentile):
    '''
    (list, int) -> num
    Return the value corresponding to the percentile from the list of values L
    Precondition: percentile is not in %
    '''

    # order the list
    L.sort()
    # use the nearest rank method to find the percentile rank
    Q = percentile / 100 * len(L)
    if int(Q+1) >= len(L):
        Qposition = int(Q-1)
    else:
        Qposition = int(Q+1)
    return L[Qposition]


# use this function to get the coordinates of the C. elegans 3' UTRs
def celegans_annotated_three_prime_coordinates(celegans_gff):
    '''
    (file) -> dict
    Returns a dictionnary with the 3' UTR coordinates of each celegans transcript
    Use a single UTR if multiple UTRs are annotated
    '''

    # create a dictionnary to stote the coordinates {transcript_name : [chromo, start, end, orienation]}
    annotated_three_prime = {}

    # open file for reading
    gff = open(celegans_gff, 'r')
    for line in gff:
        line = line.rstrip()
        if line != '':
            line = line.split()
            if line[1] == 'WormBase':
                if line[2] == 'three_prime_UTR':
                    transcript = line[8][line[8].index('Transcript:')+11 :]
                    chromo = line[0]
                    # get positions 0-based
                    start = int(line[3]) - 1
                    end = int(line[4])
                    orientation = line[6]
                    annotated_three_prime[transcript] = [chromo, start, end, orientation]
    gff.close()
    return annotated_three_prime
 
   
# use this function to get the UTR sequences of each transcript
def celegans_UTR_sequences(celegans_gff, assembly):
    '''
    (file, file, int) -> dict
    Returns a dictionnary with the transcript name as key and the UTR sequence
    from the assembly as value. Ignore transcripts lacking a UTR sequence
    (e.g. at the ends of chromosome or when adjacent to another gene)
    '''

    # convert assembly to fasta format
    genome = convert_fasta(assembly)
    
    # create a dictionnary with UTR coordinates
    three_prime = celegans_annotated_three_prime_coordinates(celegans_gff)
    
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
            # get positions (already 0-based)
            start = three_prime[transcript][1] 
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
    

# use this function to get the coordinates of the annotated 3'UTR in remanei or latens
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


# use this function to get the coodinates of all transcripts in remanei or latens
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
    Coordinates are 0-based    
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
    Take a dict of transcript: gene pairs, a dict of transcript genomic features,
    and for a given given gene at index position in the list coordinates,
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



# use this function to get the coordinates of the UTR/downstream sequences of each transcript
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


   
# use this function to convert coordinates in a sequence to genomic coordinates
def convert_seq_coord_genome_coord(site_start, site_end, sequence_start, sequence_end, sequence_sense, length_chromo):
    '''
    (int, int, int, int, str, int) -> (int, int)
    Take the coordinates of a feature in a given sequence,
    the coordinates of that sequence in the genome, its orientation,
    the chromosome length where the sequence is located and return the
    coordinates of the feature relative to chromosome    
    Precondition: all coordinates are 0-index based and the sequence is continuous   
    '''

    # compute site length and sequence length
    site_length = site_end - site_start
    sequence_length = sequence_end - sequence_start    
    
    # check  sequence orientation
    if sequence_sense == '+':
        # site coordinates on chromosome can be directly deduced
        # start = start_index_site + start_index_sequence
        start = site_start + sequence_start
        # end = start + length_site, also end = end_index_site + start_index_sequence
        end = start + site_length
        return start, end
    elif sequence_sense == '-':
        # sites are predicted based on reverse_complement of sequence
        # need to adjust coordinates of site on chromo
        # get the indices of site in sequence in - sense
        site_start_sequence = sequence_length - site_end
        # site indices on chromo can be deduced from site indices on sequence(-)
        start = site_start_sequence + sequence_start
        end = start + site_length
        return start, end
        

# use this function to get positions not falling in site categories
def keep_intergenic_positions(genome, CDS_pos, UTR_pos, intron_pos):
    '''
    (dict, dict, dict, dict) -> dict
    Take the dictionary of genome sequence, the dict of CDS coordinates,
    predicted UTR coordinates, intron coordinates and return a dictionnary of 
    intergenic positions for each chromo after removing sites falling in the other sites
    Precondition: all dicts are in the form {chromo: {set of positions}}
    '''
    
    # make a dict with positions in intergenic regions
    # {chromo: {set of positions}}
    intergenic_pos = {}
    for chromo in genome:
        intergenic_pos[chromo] = set(i for i in range(len(genome[chromo])))
    # loop over CDS pos, remove positions falling in CDS
    for chromo in CDS_pos:
        # loop over position on that chromo
        for i in CDS_pos[chromo]:
            if i in intergenic_pos[chromo]:
                intergenic_pos[chromo].remove(i)
    # loop over UTR pos, remove positions falling in UTR
    for chromo in UTR_pos:
        # loop over positions on that chromo
        for i in UTR_pos[chromo]:
            if i in intergenic_pos[chromo]:
                intergenic_pos[chromo].remove(i)
    # loop over intron pos
    for chromo in intron_pos:
        # loop over positoons on that chromo
        for i in intron_pos[chromo]:
            if i in intergenic_pos[chromo]:
                intergenic_pos[chromo].remove(i)
            
    return intergenic_pos