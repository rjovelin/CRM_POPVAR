from manipulate_sequences import *



    
# use this function to get the positions of the coding sequence for each transcript
def get_CDS_positions(caeno_gff):
    '''
    (file) -> dict
    Returns a dictionnary with transcript name as key and list of tuples containing the start and end positions
    of each coding regions of the given transcript, plus the chromosome, and the orientation of the gene on chromosome
    '''

    # create dictionnary to store positions
    # {TS1: [chromo, sense, [(s1, end1), (s2, end2)]]}
    CDS = {}

    # open gff file for reading
    gff = open(caeno_gff, 'r')
    for line in gff:
        line = line.rstrip()
        if line != '':
            line = line.split()
            if len(line) > 8:
                if line[2] == 'mRNA':
                    chromo = line[0]
                    orientation = line[6]
                    transcript = line[8][line[8].index('ID=')+3 : line[8].index(';')]
                    CDS[transcript] = [chromo, orientation, []]
                elif line[2] == 'CDS':
                    start = int(line[3])
                    end = int(line[4])
                    name = line[8][line[8].index('ID=') + 3: line[8].index(':')]
                    if name in CDS:
                        CDS[name][2].append((start, end))

    # sort the coordinates relative to their order on chromosome
    for transcript in CDS:
        CDS[transcript][2].sort()

    gff.close()
    return CDS

# use this function to generate a fasta file with the transcripts' coding sequences
def grab_CDS_sequences(caeno_gff, genome_file, output_file):
    '''
    (file, file) -> file
    Extracts the positions of the coding sequences for each transcript from the gff file,
    grab the coding sequences using the coordinates from the genome sequence, and
    save sequences excluding the terminal stop codon to output_file in fasta format
    '''

    # make a set with stop codons:
    stop_codons = {'TAA', 'TAG', 'TGA', 'taa', 'tag', 'tga'}

    # convert the genome sequence into a dictionnary
    genome = convert_fasta(genome_file)
    
    # get the coordinates of the coding sequences for each transcript
    CDS_coordinates = get_CDS_positions(caeno_gff)
    
    # create a dictionnary to store the CDS sequences
    CDS_sequences = {}
    
    # use CDS coordinates to fetch the coding sequences
    for transcript in CDS_coordinates:
        chromo = CDS_coordinates[transcript][0]
        chromo_seq = genome[chromo]
        transcript_seq = ''
        # for each CDS cordinate start, end pair, grab the sequence
        for pair in CDS_coordinates[transcript][2]:
            start = pair[0] - 1 # indice of seq starts at 0
            end = pair[1] # end is not included
            exon = chromo_seq[start:end]
            transcript_seq += exon
        # if the orientation is anti-sense, take the reverse complement
        if CDS_coordinates[transcript][1] == '-':
            rev_compl_seq = reverse_complement(transcript_seq)
            CDS_sequences[transcript] = rev_compl_seq
        elif CDS_coordinates[transcript][1] == '+':
            CDS_sequences[transcript] = transcript_seq

    # make the sequences in upper case
    for transcript in CDS_sequences:
        CDS_sequences[transcript] = CDS_sequences[transcript].upper()

    # exclude terminal codon
    for transcript in CDS_sequences:
        if CDS_sequences[transcript][-3:] in stop_codons:
            CDS_sequences[transcript] = CDS_sequences[transcript][:-3]

    # open file for writing
    newfile = open(output_file, 'w')
    for transcript in CDS_sequences:
        newfile.write('>' + transcript + '\n')
        newfile.write(CDS_sequences[transcript] + '\n')

    newfile.close()
        


    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    