# -*- coding: utf-8 -*-
"""
Created on Fri Jul 10 09:31:03 2015

@author: Richard
"""

from manipulate_sequences import *
from repeats_TEs import *
from get_coding_sequences import *
from genomic_coordinates import *





# use this function to verify the piRNA sequences in WS230
def verify_piRNA_sequences(pirnas_file, WS230_genome_file):
    '''
    (file, file ) -> dict
    Take the file with piRNA sequences and coordinates and the WS230 genome
    sequence from which coordinates were derived, and return a dict with
    remanei: sequence pairs for which the piRNA genomic sequence doesn't match
    the piRNA sequence extracted from the DNA sequence in pirnas_file
    '''
    
    # convert genome into dict
    genome = convert_fasta(WS230_genome_file)
    # convert pirnas_file into dict
    pirnas = convert_fasta(pirnas_file)
    
    # create dict to store the pirnas with mismatch sequences
    different_seq = {}
    
    # loop over pirna in pirnas dict
    for seqname in pirnas:
        # extract coordinates from seqname
        coordinates = seqname.split('|')
        # extract chromosome, orientation, start and end
        chromo = coordinates[0]
        orientation = coordinates[1]
        start = int(coordinates[2]) -1
        end = int(coordinates[-1])
        # extract genome pirnas sequence
        genomic_pirna = genome[chromo][start:end].upper()
        if orientation == '-':
            genomic_pirna = reverse_complement(genomic_pirna)
        # extra pirna sequence from DNA sequence
        pirna_seq = pirnas[seqname][70:].upper()
        # compare genomic pirna sequence and extracted pirna sequence
        if pirna_seq != genomic_pirna:
            # populate dict
            different_seq[seqname] = [pirna_seq, genomic_pirna]
    
    return different_seq
    
    
# use this function to map piRNAs from WS230 to PX356 genome 
def map_remanei_pirnas(genome_file, pirnas_file):
    '''
    (file, file) -> dict
    Take the PX356 genome assembly and the file with piRNA sequences and their 
    coordinates in the WS230 assembly and return a dict with number as key and 
    a list of coordinates where the pirna and its 70 bp sequences match
    '''
    
    # convert genome file to dict
    genome = convert_fasta(genome_file)
    # convert pirnas_file to dict
    pirnas = convert_fasta(pirnas_file)
    # create a dict to store the coordinates {i: [[chromo, sens, start, end], [chromo, sens, start, end]]}
    pirna_coordinates = {}
    # initialize key
    i = 0
    # get a list of numbers multiple of 100
    num = [k for k in range(0, len(pirnas), 100)]
    
    
    # create a dict to store the position already recorded to avoid 
    # recording the same position for identical piRNAs
    already_recorded = {}
    # initialise list values
    for chromo in genome:
        already_recorded[chromo] = []
    
    # loop over pirnas
    for seqname in pirnas:
        # initialize empty list
        pirna_coordinates[i] = []
        # get pirna+upstream seq
        pirna_up = pirnas[seqname].upper()
        # get pirna_seq
        pirna_seq = pirna_up[70:].upper()
        # loop over genome
        for chromo in genome:
            # count number of times pirna + upstream sequence is found in chromo
            N_found = genome[chromo].count(pirna_up)
            # if N_found != 0: get all the indices with a sequence match
            # set up index to look for a match
            j = -1
            while N_found != 0:
                # get the index of pirna + 70 bp upstream
                j = genome[chromo].index(pirna_up, j+1)
                # create a list of coordinates for the piRNA sequence only (without the upstream sequence)
                coord = [chromo, '+', j + 70, j + 70 + len(pirna_seq)]
                # check if coordinates is already recorded
                if coord not in already_recorded[chromo]:
                    # record piRNA coordinates in dict [chromo, orientation, start, end]
                    pirna_coordinates[i].append(coord)
                    # add coordinates to list of recorded coord
                    already_recorded[chromo].append(coord)
                N_found -= 1
            # count the number of times reverse complement of pirna + upstream is found in chromo
            pirna_up_rv = reverse_complement(pirna_up)
            N_found = genome[chromo].count(pirna_up_rv)
            # if rev compl of pirna is found, get all indices where sequence matches
            # sex up index to look for for a match
            j = -1
            while N_found != 0:
                j = genome[chromo].index(pirna_up_rv, j+1)
                # create a list of coordinates for the piRNA sequence only (without the upstream sequence)
                coord = [chromo, '-', j, j + len(pirna_seq)]
                # check if coord already recorded
                if coord not in already_recorded[chromo]:
                    # record coordinates in dict [chromo, orientation, start, end]
                    pirna_coordinates[i].append(coord)
                    # add coord to list of recorded coord
                    already_recorded[chromo].append(coord)
                N_found -= 1
        # update key counter
        i += 1
        if i in num:
            print(i)
        
    # remove pirnas that didn't match the PX356 genome
    to_remove = []
    for i in pirna_coordinates:
        if len(pirna_coordinates[i]) == 0:
            to_remove.append(i)
    for i in to_remove:
        del pirna_coordinates[i]
    
    print(len(to_remove), 'pirnas did not match')
    
    return pirna_coordinates
    
    
# use this function to generate a file with coordinates of piRNAs
def write_pirna_coordinates_to_file(pirna_coordinates, outputfile):
    '''
    (dict) -> file
    Take the dictionnary of piRNA coordinates and write the coordinates of each
    piRNA in file, including all genomic locations when a a single sequence maps
    in multiple locations
    '''
    
    # open file for writing
    newfile = open(outputfile, 'w')
    # write header
    header = '\t'.join(['piRNA_name', 'chromosome', 'orientation', 'start', 'end'])
    newfile.write(header + '\n')
    # create a counter variable to assign names to piRNAs
    i = 0
    # loop over piRNAs in dictionnary
    for pirna in pirna_coordinates:
        # loop over the coordinates
        for coord in pirna_coordinates[pirna]:
            # get the pirna name
            name = 'piRNA_' + str(i)
            # get chromo, orientation
            chromo = coord[0]
            orientation = coord[1]
            # convert start and end positions to 1-based indices
            start = coord[2] + 1
            end = coord[3]
            # write coordinates to file
            line = '\t'.join([name, chromo, orientation, str(start), str(end)])
            newfile.write(line + '\n')
            # update counter variable
            i += 1
    
    # close file after writing
    newfile.close()
            
        
# use this function to get the piRNA sequences from the PX356 genome 
def grab_piRNA_sequences(pirna_coord_file, genome_file):
    '''
    (file, file) -> dict
    Take the file with piRNA coordinates in the PX356 remanei genome sequence
    and the genome sequence in fasta format and return a dictionnary of piRNA name:
    sequence pairs
    Note: Some piRNAs have multiple locations (N pirna loci different than N pirna sequences)    
    '''
    
    # convert genome to dict
    genome = convert_fasta(genome_file)
    
    # create a dict to store pirnas {pirna_num: sequence}
    pirnas = {}
    
    # open file for reading
    infile = open(pirna_coord_file, 'r')
    # skip header
    infile.readline()
    # loop over files
    for line in infile:
        if line.rstrip() != '':
            line = line.rstrip().split()
            # get chromo
            chromo = line[1]
            # get orientation
            orientation = line[2]
            # get start and end , 0-index based
            start = int(line[3]) - 1
            end = int(line[4])
            pirna_seq = genome[chromo][start:end]
            if orientation == '-':
                pirna_seq = reverse_complement(pirna_seq)
            pirnas[line[0]] = pirna_seq
            
    # close file after reading
    infile.close()
    
    return pirnas


 
# use this function to count the number of pirnas on each chromo  
def count_pirna_per_chromo(pirna_coord_file):
    '''
    (file) -> dict
    Take the file with the piRNA coordinates in PX356 and return a dictionnary
    with chromo as key and the count of piRNAs on chromo as value
    '''
    
    # create dict {chromo: number of pirnas}
    pirna_count = {}
    
    # open file for reading
    infile = open(pirna_coord_file, 'r')
    # skip header
    infile.readline()
    
    # loop over file
    for line in infile:
        if line.rstrip() != '':
            line = line.rstrip().split()
            # get chromo
            chromo = line[1]
            # check if chromo is key in dict
            if chromo in pirna_count:
                # chromo already key, add 1
                pirna_count[chromo] += 1
            else:
                # chromo not key, add key to dict and set counter to 1
                pirna_count[chromo] = 1
    # close file after reading
    infile.close()
    
    return pirna_count
    



# use this function to get the starting position of each pirna
def chromo_pirna_start(pirna_coord_file):
    '''
    (file) -> dict
    Take the file with piRNA coordinates in PX356 and return a dictionnary
    with chromo as key and a list of start positions for each piRNA on that 
    chromo
    '''
    
    # create dict
    pirna_start = {}
    
    # open file for reading
    infile = open(pirna_coord_file, 'r')
    # skip header
    infile.readline()
    
    # loop over file
    for line in infile:
        if line.startswith('p'):
            line = line.rstrip().split()
            chromo = line[1]
            start = int(line[3]) - 1
            # check if chromo is key in dict
            if chromo in pirna_start:
                # add start position
                pirna_start[chromo].append(start)
            else:
                # initialize list
                pirna_start[chromo] = [start]
                
    # close file after reading
    infile.close()
    
    return pirna_start
    
    
# use this function to generate save the number of pirnas per window
def cluster_pirnas(pirna_start, genome, chromo, window_size, outputfile):
    '''
    (dict, dict, str, int, file) -> file
    Take a dictionnary with chromo: list of start position pairs, a dict
    with chromo : sequence pairs, the focal chromosome, the size of the window,
    and save the number of pirnas per window in outputfile
    '''
    
    # make list of size window that contains only 0s:
    # each value in the list is the count of position for the range [0 - window[ etc
    range_counts = [0] * (len(genome[chromo]) // window_size)
    
    # loop over starting positions
    for start in pirna_start[chromo]:
        # determine the index in the list range_count where the position should be added
        which_range = start // window_size
        if which_range == len(range_counts):
            which_range -= 1
        # count pirnas
        range_counts[which_range] += 1
        
    # open file for writing
    newfile = open(outputfile, 'w')
    # write header
    newfile.write('Range' + '\t' + 'Lower_point' + '\t' + 'Midpoint' + '\t' + 'Higher_point' + '\t' + 'Count' + '\t' + '\n')
        
    # loop over indices of list
    for i in range(len(range_counts)):
        newfile.write('[' + str(i * window_size) + '-' + str((i * window_size) + window_size -1) + ']' + '\t')
        newfile.write(str(i * window_size) + '\t')
        newfile.write(str(int(((i * window_size) + window_size) / 2)) + '\t')
        newfile.write(str((i * window_size) + window_size) + '\t')
        newfile.write(str(range_counts[i]) + '\n')
        
    newfile.close()







# use this function to get the coordinates of each pirna on each chromosomes
def get_pirna_loci(pirna_coord_file):
    '''
    (file) -> dict
    Take the file with pirna coordinates and return a dictionary with chromo as
    key a list of pirna coordinates as value
    '''
    
    # create a dict with chromo as key and a list of pirna coordinates as value
    # {chromo: [[start, end, orienation]]}
    pirna_loci = {}
    
    # open file for reading
    infile = open(pirna_coord_file, 'r')
    # skip header
    infile.readline()
    # loop over files
    for line in infile:
        if line.rstrip() != '':
            line = line.rstrip().split()
            # get chromo
            chromo = line[1]
            # get orientation
            orientation = line[2]
            # get start and end , 0-index based
            start = int(line[3]) - 1
            end = int(line[4])
            # check if chromo in dict
            if chromo in pirna_loci:
                pirna_loci[chromo].append([start, end, orientation])
            else:
                pirna_loci[chromo] = [[start, end, orientation]]
                
    # close file after reading
    infile.close()
    
    return pirna_loci









  

# use this function to find potential piRNA targets in TEs
def find_pirna_TE_targets(pirna_coord_file, repeatmasker_output, genome_fasta, mismatch):
    '''
    (file, file, file, int) -> dict
    Take the file with piRNA coordinates, the output file of Repeat Masker,
    the fasta file of the genome sequence, and the number of allowed mismatch 
    a piRNA and a target, and return a dict with repeat as key and a list of
    list of coordinates for all targets in that repeat 
    '''
    
    # convert fasta to dict
    genome = convert_fasta(genome_fasta)

    # get the pirna sequences
    pirnas = grab_piRNA_sequences(pirna_coord_file, genome_fasta)
    
    # make a dictionary of {repeat class : set of repeat names} 
    repfam = repeat_family(repeatmasker_output)
    
    # invert dict to create dict {repeat name : repeat class} pairs
    famrep = {}
    for fam in repfam:
        for repeat in repfam[fam]:
            famrep[repeat] = fam
    
    # get the repeat coordinates {repeat : [[chromo, start, end, orientation], [chromo, start, end, orientation]]} 
    repeat_coord = get_repeats_coord(repeatmasker_output, True)
    
    # create dict
    TE_targets = {}    
    
    
    # convert set to list to use indexing and slicing
    for repname in repeat_coord:
        repeat_coord[repname] = list(repeat_coord[repname])
    
    # search the antisense matches in TE
    # loop over pirnas
    for seqname in pirnas:
        print(seqname, len(TE_targets))
        # take the reverse complement of the prina sequence
        rev_seq = reverse_complement(pirnas[seqname]).upper()
        # loop over all TE, find a match
        for repname in repeat_coord:
            # get the repeat family
            fam = famrep[repname]
            # loop over all instances of TE
            for i in range(len(repeat_coord[repname])):
                # get chromo
                chromo = repeat_coord[repname][i][0]
                # get start and end position
                TE_start = repeat_coord[repname][i][1]
                TE_end = repeat_coord[repname][i][2]
                # get orientation
                orientation = repeat_coord[repname][i][3]
                # get the repeat sequence
                TE_seq = genome[chromo][TE_start:TE_end].upper()
                # check orientation
                if orientation == '-':
                    TE_seq = reverse_complement(TE_seq)                 
                # check how many mismatches are allowed
                if mismatch == 0:
                    # count the number of times pirna matches in TE seq
                    N_found = TE_seq.count(rev_seq)
                    # if N_found != 0: get all the indices with a sequence match
                    # set up index to look for a match
                    j = -1
                    while N_found != 0:
                        # get the index of pirna match
                        j = TE_seq.index(rev_seq, j+1)
                        site_start = j
                        site_end = j + len(rev_seq)
                        # get pirna target seq
                        pirna_target_seq = TE_seq[site_start:site_end]
                        # compute genomic coordinates 
                        genomic_start, genomic_end = convert_seq_coord_genome_coord(site_start, site_end, TE_start, TE_end, orientation, len(chromo))
                        # check that coordinates are coorrect
                        genome_target = genome[chromo][genomic_start:genomic_end]
                        # check orientation
                        if orientation == '-':
                            genome_target = reverse_complement(genome_target)
                        assert genome_target.upper() == pirna_target_seq.upper(), 'genomic coordinates seem incorrect'
                        # check if repeat in dict
                        if repname in TE_targets:
                            TE_targets[repname].append([fam, chromo, genomic_start, genomic_end, orientation, seqname, pirna_target_seq])
                        else:
                            TE_targets[repname] = [[fam, chromo, genomic_start, genomic_end, orientation, seqname, pirna_target_seq]]
                        N_found -= 1
                elif mismatch != 0:
                    # scan the entire sequence
                    for k in range(len(TE_seq) - len(rev_seq) +1):
                        # compare the reverse seq of the pirna and the potential target
                        pirna_target_seq = TE_seq[k:k+ len(rev_seq)]
                        diff = match_diff(rev_seq, pirna_target_seq)
                        # up to N mismatch included are tolerated 
                        if diff <= mismatch:
                            # recored target site
                            # compute genomic coordinates 
                            genomic_start, genomic_end = convert_seq_coord_genome_coord(k, k+len(rev_seq), TE_start, TE_end, orientation, len(chromo))
                            genome_target = genome[chromo][genomic_start:genomic_end]                            
                            # check orientation
                            if orientation == '-':
                                genome_target = reverse_complement(genome_target)
                            assert genome_target.upper() == pirna_target_seq.upper(), 'genomic coordinates seem incorrect'
                            # check if repeat in dict
                            if repname in TE_targets:
                                TE_targets[repname].append([fam, chromo, k, k+len(rev_seq), orientation, seqname, pirna_target_seq])
                            else:
                                TE_targets[repname] = [[fam, chromo, k, k+len(rev_seq), orientation, seqname, pirna_target_seq]]
                                
    return TE_targets                
    
    
# use this function to find potential piRNA targets
def find_pirna_CDS_targets(pirna_coord_file, CDS_fasta, genome_fasta, mismatch, unique_transcripts):
    '''
    (file, file, file, int, file, file) -> dict
    
    
    return a dictionary with target gene as key and a list with start, end 
    position in the CDS and the target sequence
    '''
    
    # make a set of valid transcripts
    valid_transcripts = get_valid_transcripts(unique_transcripts)    
        
    # get the pirna sequences
    pirnas = grab_piRNA_sequences(pirna_coord_file, genome_fasta)
    
    # get the CDS sequences
    CDS = convert_fasta(CDS_fasta)
    # remove non-valid transcripts
    to_remove = [i for i in CDS if i not in valid_transcripts]
    if len(to_remove) != 0:
        for i in to_remove:
            del CDS[i]
       
    # create dicts to store the targets {gene: [[start, end, pirna name, target sequence]]}
    # start and end are relative to the the CDS sequence, not the genomic coordinates
    CDS_targets = {}
                
    # search the antisense matches in CDS
    # loop over pirnas
    for seqname in pirnas:
        # take the reverse complement of the prina sequence
        rev_seq = reverse_complement(pirnas[seqname]).upper()
        # loop over all CDS, find a match
        for gene in CDS:
            # check how many mismatch are allowed
            if mismatch == 0:
                # count the number of times pirna matches in CDS seq
                N_found = CDS[gene].count(rev_seq)
                # if N_found != 0: get all the indices with a sequence match
                # set up index to look for a match
                j = -1
                while N_found != 0:
                    # get the index of pirna match
                    j = CDS[gene].index(rev_seq, j+1)
                    start = j
                    end = j + len(rev_seq)
                    # get pirna target seq
                    pirna_target_seq = CDS[gene][start:end]
                    # check if gene in dict
                    if gene in CDS_targets:
                        CDS_targets[gene].append([start, end, seqname, pirna_target_seq])
                    else:
                        CDS_targets[gene] = [[start, end, seqname, pirna_target_seq]]
                    N_found -= 1
            elif mismatch != 0:
                # scan the entire sequence
                for k in range(len(CDS[gene]) - len(rev_seq) +1):
                    # compare the reverse seq of the pirna and the potential target
                    pirna_target_seq = CDS[gene][k:k+ len(rev_seq)]
                    diff = match_diff(rev_seq, pirna_target_seq)
                    # up to N mismatch included are tolerated 
                    if diff <= mismatch:
                        # recored target site
                        # check if gene in dict
                        if gene in CDS_targets:
                            CDS_targets[gene].append([k, k+len(rev_seq), seqname, pirna_target_seq])
                        else:
                            CDS_targets[gene] = [[k, k+len(rev_seq), seqname, pirna_target_seq]]
                                
    return CDS_targets                
     


    

# use this function to count piRNAs in different locations
def find_pirna_locations(genome, pirna_coord, CDS_pos, UTR_pos, intron_pos, intergenic_pos):
    '''
    (dict, dict, dict, dict, dict, dict, dict) -> dict
    Take the dictionary of genome sequences, the  dictionary of piRNA coordinates,
    a dictionary with the CDS positions, a dictionary with the predicted UTR positions
    a dictionary with the intron positions, a dictionary with intergenic positions 
    and return a dict of counts of piRNAs located in CDS, predicted UTR, introns,
    intergenic and piRNAs that are not fully included in any of these categories
    (ie. overlapping different categories)
    Precondition: pirna_coord is in the form {chromo: [list of pirna coordinates]},
    the other dicts are in the form {chromo: set(positions on chromo)}
    '''
    
    # create a lit of counts
    locations = {}
    
    # create dict with sets of pirna already recorded
    already_recorded = {}
    # inititate sets
    for chromo in pirna_coord:
        already_recorded[chromo] = []
    
    # search pirnas in intergenic   
    # loop over chromo in pirna coord
    for chromo in pirna_coord:
        # check that chromo in intergenic_coord
        if chromo in intergenic_pos:
            # loop over pirna on that chromo
            for i in range(len(pirna_coord[chromo])):
                # check that pirna is not already recorded
                if pirna_coord[chromo][i] not in already_recorded[chromo]:
                    # search for pirna location
                    coordinate_pirna = set(range(pirna_coord[chromo][i][0], pirna_coord[chromo][i][1]))
                    # compare coordinates
                    if coordinate_pirna.issubset(intergenic_pos[chromo]):
                        # check that intergenic is key in dict
                        if 'intergenic' in locations:
                            # update count
                            locations['intergenic'] += 1
                        else:
                            # add key to dict
                            locations['intergenic'] = 1
                        # add pirna to doct already recorded
                        already_recorded[chromo].append(pirna_coord[chromo][i])

    # search pirnas in intron    
    # loop over chromo in pirna coord
    for chromo in pirna_coord:
        # check that chromo in intron_coord
        if chromo in intron_pos:
            # loop over pirna on that chromo
            for i in range(len(pirna_coord[chromo])):
                # check that pirna is not already recorded
                if pirna_coord[chromo][i] not in already_recorded[chromo]:
                    #search for pirna location
                    coordinate_pirna = set(range(pirna_coord[chromo][i][0], pirna_coord[chromo][i][1]))
                    # compare coordinates
                    if coordinate_pirna.issubset(intron_pos[chromo]):
                        # check that intron is key in dict
                        if 'intron' in locations:
                            # update count
                            locations['intron'] += 1
                        else:
                            # add key to dict
                            locations['intron'] = 1
                        # add pirna to doct already recorded
                        already_recorded[chromo].append(pirna_coord[chromo][i])

    # search pirnas in UTR    
    # loop over chromo in pirna coord
    for chromo in pirna_coord:
        # check that chromo in UTR_coord
        if chromo in UTR_pos:
            # loop over pirna on that chromo
            for i in range(len(pirna_coord[chromo])):
                # check that pirna is not already recorded
                if pirna_coord[chromo][i] not in already_recorded[chromo]:
                    # search for pirna location
                    coordinate_pirna = set(range(pirna_coord[chromo][i][0], pirna_coord[chromo][i][1]))
                    # compare coordinates
                    if coordinate_pirna.issubset(UTR_pos[chromo]):
                        # check that utr is key in dict
                        if 'UTR' in locations:
                            # update count
                            locations['UTR'] += 1
                        else:
                            # add key to dict
                            locations['UTR'] = 1
                        # add pirna to doct already recorded
                        already_recorded[chromo].append(pirna_coord[chromo][i])
    
    # search pirnas on CDS    
    # loop over chromo in pirna coord
    for chromo in pirna_coord:
        # check that chromo in CDS_coord
        if chromo in CDS_pos:
            # loop over pirna on that chromo
            for i in range(len(pirna_coord[chromo])):
                # check that pirna is not already recorded
                if pirna_coord[chromo][i] not in already_recorded[chromo]:
                    # search for pirna location
                    coordinate_pirna = set(range(pirna_coord[chromo][i][0], pirna_coord[chromo][i][1]))
                    # compare coordinates
                    if coordinate_pirna.issubset(CDS_pos[chromo]):
                        # check that cds is key in dict
                        if 'CDS' in locations:
                            # update count
                            locations['CDS'] += 1
                            # add pirna to already recorded set
                        else:
                            # add key to dict
                            locations['CDS'] = 1
                        # add pirna to doct already recorded
                        already_recorded[chromo].append(pirna_coord[chromo][i])
                                
    # loop over pirna
    for chromo in pirna_coord:
        # loop over pirna on that chromo
        for i in range(len(pirna_coord[chromo])):
            # check if pirna is recorded
            if pirna_coord[chromo][i] not in already_recorded[chromo]:
                # pirna is not fully contained in a site category, may be overlapping different sites
                # check if overlapping is key in dict
                if 'overlapping' in locations:
                    locations['overlapping'] += 1
                else:
                    locations['overlapping'] = 1
                    
    return locations
        