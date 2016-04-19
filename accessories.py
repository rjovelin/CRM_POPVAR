#!/usr/bin/env python3




def cds_translate(cds):
    '''
    (str) -> str
    Translate a coding sequence into a protein sequence according to the standard genetic code

    >>> cds_translate('ATGGCCATGGCGCCCAGAACTGAGATCAATAGTACCCGTATTAACGGGTGA')
    MAMAPRTEINSTRING*
    >>> cds_translate('ATGTACTAA')
    MY*
    '''

    genetic_code = {'TTT': 'F', 'CTT': 'L', 'ATT': 'I', 'GTT': 'V',
                   'TTC': 'F', 'CTC': 'L', 'ATC': 'I', 'GTC': 'V',
                   'TTA': 'L', 'CTA': 'L', 'ATA': 'I', 'GTA': 'V',
                   'TTG': 'L', 'CTG': 'L', 'ATG': 'M', 'GTG': 'V',
                   'TCT': 'S', 'CCT': 'P', 'ACT': 'T', 'GCT': 'A',
                   'TCC': 'S', 'CCC': 'P', 'ACC': 'T', 'GCC': 'A',
                   'TCA': 'S', 'CCA': 'P', 'ACA': 'T', 'GCA': 'A',
                   'TCG': 'S', 'CCG': 'P', 'ACG': 'T', 'GCG': 'A',
                   'TAT': 'Y', 'CAT': 'H', 'AAT': 'N', 'GAT': 'D',
                   'TAC': 'Y', 'CAC': 'H', 'AAC': 'N', 'GAC': 'D',
                   'TAA': '*', 'CAA': 'Q', 'AAA': 'K', 'GAA': 'E',
                   'TAG': '*', 'CAG': 'Q', 'AAG': 'K', 'GAG': 'E',
                   'TGT': 'C', 'CGT': 'R', 'AGT': 'S', 'GGT': 'G',
                   'TGC': 'C', 'CGC': 'R', 'AGC': 'S', 'GGC': 'G',
                   'TGA': '*', 'CGA': 'R', 'AGA': 'R', 'GGA': 'G',
                   'TGG': 'W', 'CGG': 'R', 'AGG': 'R', 'GGG': 'G'}
   

    CDS = cds.upper()
    protein = ''

    for i in range(0, len(CDS), 3):
        codon = CDS[i:i+3]
        if codon not in genetic_code:
            protein += 'X'
        else:
            protein += genetic_code[codon]

    return protein


def convert_fasta(fasta):
    '''
    (file) -> dict
    Take a file with fasta sequences and return a dictionnary with
    sequence ID as key and single string sequence as value
    '''
    # convert nematode genomes into a dictionnary
    genome = {}
    infile = open(fasta, 'r')
    for line in infile:
        line = line.rstrip()
        if line != '':
            if line.startswith('>'):
                genome[line[1:]] = ""
                seq_name = line[1:]
            else:
                genome[seq_name] += line
    infile.close
    return genome


def reverse_complement(dna):
    '''
    (str) -> (str)
    Return the reverse complementary sequence of string dna

    >>> reverse_complement('atcg')
    'cgat'
    '''

    valid_bases = {'A', 'T', 'C', 'G'}

    dna2 = dna.upper()
    dna_comp = ''
    for i in dna2:
        if i == 'A':
            dna_comp += 'T'
        elif i == 'T':
            dna_comp += 'A'
        elif i == 'C':
            dna_comp += 'G'
        elif i == 'G':
            dna_comp += 'C'
        elif i not in valid_bases:
            dna_comp += 'N'

    reverse_comp_dna = ''
    for i in reversed(dna_comp):
        reverse_comp_dna += i

    if dna.islower():
        reverse_comp_dna = reverse_comp_dna.lower()
        
    return reverse_comp_dna


# use this function to get the complement of a DNA sequence
def seq_complement(dna):
    '''
    (str) -> str
    Take a DNA sequence and return its complement in upper case
    Precondition: accepts only N ambiguities, and treats any ambiguities as Ns
    
    >>> seq_complement('atcg')
    'TAGC'
    >>> seq_complement('CTGagTg')
    'GACTCAC'    
    >>> seq_complement('CTGARGTC')
    'GACTNCAG'
    '''

    valid_bases = {'A', 'T', 'C', 'G'}
        
    dna2 = dna.upper()
    dna_comp = ''
    for i in dna2:
        if i == 'A':
            dna_comp += 'T'
        elif i == 'T':
            dna_comp += 'A'
        elif i == 'C':
            dna_comp += 'G'
        elif i == 'G':
            dna_comp += 'C'
        elif i not in valid_bases:
            dna_comp += 'N'

    return dna_comp    


# use this function to reverse a sequence
def reverse_seq(seq):
    '''
    (str) -> str
    Take a sequence and return its reverse sequence in upper case
    
    >>> reverse_seq('actat')
    'TATCA'
    >>> reverse_seq('CGTCGTATCG')
    'GCTATGCTGC'
    '''
    
    seq_rv = ''
    
    for i in reversed(seq):
        seq_rv += i
    
    return seq_rv.upper()



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
    


