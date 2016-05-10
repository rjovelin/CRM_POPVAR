from Genomic_Coordinates import *

# use this function to translate a coding sequence
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

# use this function to convert a fasta file to a dictionary
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

# use this function to reverse complement a DNA sequence
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



# use this function to clean up the file of 1:1 orthologs from transcripts that belong to a same gene 
def clean_cel_crem_orthologs_file(best_blast_hits, crem_gff, cel_gff, outputfile):
    '''
    (file, file, file) -> file
    Clean up the file with best blat hits from transcripts that belong to a same
    gene. Write to the outputfile the orthologous pairs between remanei and elegans
    by keeping a single transcript per gene in each species
    '''
    
#    # create dicts with {gene: [transcript1, transcript2]} in remanei and elegans
#    remanei = parent_gene(crem_gff)
#    elegans = celegans_gene_to_transcripts(cel_gff)
#    
    # create sets of genes for which transcripts are already used in ortholog pairs
    rem_already_used = set()
    cel_already_used = set()
    
    # create dicts of transcript : gene pairs {transcript : gene}
    remanei_ts = transcript_to_gene(crem_gff)
    elegans_ts = celegans_transcript_to_gene(cel_gff)
    
    # create a dict to store the ortholog pairs {crem_TS : cel_TS}
    orthologs = {}    
    
    # open best blast hits file
    infile = open(best_blast_hits, 'r')
    # skip header
    infile.readline()
    
    # loop over file
    for line in infile:
        line = line.rstrip()
        line = line.split()
        cel_TS = line[0]
        crem_TS = line[1]
        # check that the remanei gene has not been used
        if remanei_ts[crem_TS] not in rem_already_used:
            # check that the latens gene has not been used
            if elegans_ts[cel_TS] not in cel_already_used:
                # populate ortholog dict
                orthologs[crem_TS] = cel_TS
                # add corresponding genes to sets
                rem_already_used.add(remanei_ts[crem_TS])
                cel_already_used.add(elegans_ts[cel_TS])
    
    
    # open file for writing
    newfile = open(outputfile, 'w')
    # write header
    newfile.write('# 1:1 orthologs between remanei and elegans obtained by reciprocal BLAST hits\n')
    newfile.write('# genes are represented by a single transcript\n')
    newfile.write('# ambiguous sites in the reference genome are fixed\n')
    newfile.write('Cremanei' + '\t' + 'Celegans' + '\n')
    
    # write content to file
    for gene in orthologs:
        newfile.write(gene + '\t' + orthologs[gene] + '\n')
    
    # close files
    infile.close()
    newfile.close()
    
    

    
    
    
    
    
    
    
    
    
    