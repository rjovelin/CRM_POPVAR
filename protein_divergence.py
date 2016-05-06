
import os
from accessories import *



# use this function to get the relationships between genes and transcripts
def parent_gene(gff):
    '''
    (file) -> dict
    Return a dictionnary with parent gene as key adn a list of transcripts
    as value from the genome annotation GFF file
    Precondition: parses the in-house remanei and latens GFF files
    '''
    # open file for reading
    infile = open(gff, 'r')
    
    # create dict to store the genes and transcripts
    # {gene: [transcript1, transcript2]}
    genes_ts = {}
    
    # loop over file
    for line in infile:
        line = line.rstrip()
        if line != '':
            line = line.split()
            if len(line) > 8:
                # get the gene and transcript names
                if line[2] == 'mRNA':
                    # get transcript name
                    transcript = line[8][line[8].index('ID=')+3 : line[8].index(';')]
                    # get gene name
                    parent = line[8][line[8].index('Parent=')+7: line[8].index(';', line[8].index('Parent'))]
                    # populate dict 
                    if parent in genes_ts:
                        genes_ts[parent].append(transcript)
                    else:
                        genes_ts[parent] = [transcript]
    infile.close()
    return genes_ts


# use this function to clean up the file of 1:1 orthologs from transcripts that belong to a same gene 
def clean_orthologs_file(best_blast_hits, crem_gff, clat_gff, outputfile):
    '''
    (file, file, file) -> file
    Clean up the file with best blat hits from transcripts that belong to a same
    gene. Write to the outputfile the orthologous pairs between remanei and latens
    by keeping a single transcript per gene in each species
    '''
    
    # create dicts with {gene: [transcript1, transcript2]} in remanei and latens
    remanei = parent_gene(crem_gff)
    latens = parent_gene(clat_gff)
    
    # create sets of genes for which transcripts are already used in ortholog pairs
    rem_already_used = set()
    lat_already_used = set()
    
    # reverse dicts with transcript : gene pairs {transcript : gene}
    remanei_ts = {}
    for gene in remanei:
        for ts in remanei[gene]:
            remanei_ts[ts] = gene
    latens_ts = {}
    for gene in latens:
        for ts in latens[gene]:
            latens_ts[ts] = gene
    
    # create a dict to store the ortholog pairs {crem_TS : clat_TS}
    orthologs = {}    
    
    # open best blast hits file
    infile = open(best_blast_hits, 'r')
    # skip header
    infile.readline()
    
    # loop over file
    for line in infile:
        line = line.rstrip()
        line = line.split()
        clat_TS = line[0]
        crem_TS = line[1]
        # check that the remanei gene has not been used
        if remanei_ts[crem_TS] not in rem_already_used:
            # check that the latens gene has not been used
            if latens_ts[clat_TS] not in lat_already_used:
                # populate ortholog dict
                orthologs[crem_TS] = clat_TS
                # add corresponding genes to sets
                rem_already_used.add(remanei_ts[crem_TS])
                lat_already_used.add(latens_ts[clat_TS])
    
    
    # open file for writing
    newfile = open(outputfile, 'w')
    # write header
    newfile.write('# 1:1 orthologs between remanei and latens obtained by reciprocal BLAST hits\n')
    newfile.write('# genes are represented by a single transcript\n')
    newfile.write('# ambiguous sites in the reference genome are fixed\n')
    newfile.write('Cremanei' + '\t' + 'Clatens' + '\n')
    
    # write content to file
    for gene in orthologs:
        newfile.write(gene + '\t' + orthologs[gene] + '\n')
    
    # close files
    infile.close()
    newfile.close()
    
# use this function to get the pairs of 1:1 orthologs    
def orthologous_pairs(ortholog_file):
    '''
    (file) -> dict
    Return a dictionnary with the pairs of orthologs from the ortholog_file
    '''
    
    # create a dict with gene names of orthologous pairs
    # { species1: species2}
    orthologs = {}
    
    # open file for reading
    infile = open(ortholog_file, 'r')
    # skip comment lines
    line = infile.readline()
    while line.startswith('#'):
        line = infile.readline()
    # skip header
    if not line.startswith('#'):
        header = line
    
    # loop over file, update dict
    for line in infile:
        line = line.rstrip()
        if line !='':
            line = line.split()
            crem = line[0]
            clat = line[1]
            orthologs[crem] = clat
       
    infile.close()
    return orthologs
    
# use this function to write separate fasta files contenaing each the sequences of the orthologs  
def save_orthologous_seq_pairs_to_file(ortholog_file, crem_CDS_file, cla_CDS_file, destination):
    '''
    (file, file, file) -> file
    Save the CDS sequences of the orthologous pairs in separate fasta files,
    by extracting the orthologous relationships from the ortholog file and the 
    CDS sequences from the CDS fasta files
    '''
    
    # convert CDS files into dictionnaries
    crem_CDS = convert_fasta(crem_CDS_file)
    cla_CDS = convert_fasta(cla_CDS_file)
    
    # create a dict with gene names of orthologous pairs
    # { species1: species2}
    orthologs = orthologous_pairs(ortholog_file)
    
    # loop over the pairs, grab the corresponding sequence, open file for writing and save sequences
    for gene in orthologs:
        # fetch the sequences of the cremanei and clatens transcripts
        crem_seq = crem_CDS[gene]
        cla_seq = cla_CDS[orthologs[gene]]
        # open file for writing
        outputfile = open(destination + gene + '_' + orthologs[gene] + '.tfa', 'w')
        outputfile.write('>' + gene + '\n')
        outputfile.write(crem_seq + '\n')
        outputfile.write('>' + orthologs[gene] + '\n')
        outputfile.write(cla_seq + '\n')
        outputfile.close()

# use this function to convert the fasta alignments to a format readbale by codeml
def alignment_file_to_PAML_format(filename, destination):
    '''
    Open alignment file generated by t-coffee, copy content into a new file
    with format compatible with PAML 
    '''
    # check that file is the t-coffee DNA alignment file
    if filename[-8:] == '_aln.tfa':
        # convert sequences to dictionnary
        alignment = convert_fasta(filename)
        # get the name of the file without '_aln.tfa'
        name = filename[:-8]
        # make a list with sequence names    
        seq_names = [gene for gene in alignment]
        # reverse order so that the remanei sequence is in the first
        seq_names.reverse()
        # open newfile for writing
        newfile = open(destination + name + '.txt', 'w')
        # add the number of sequences and the length of the alignment in the first line
        newfile.write(str(len(alignment)) + ' ' + str(len(alignment[seq_names[0]])) + '\n')
        # add orthologous sequences in fasta format
        for gene in seq_names[:-1]:
            newfile.write('>' + gene + '\n')
            newfile.write(alignment[gene] + '\n')
        newfile.write('>' + seq_names[-1] + '\n')
        newfile.write(alignment[seq_names[-1]] + '\n')
        newfile.close()
            
# use this function to create a codeml control file for a given alignment file
def generate_codeml_control_file(alignment_file, control_file, destination): 
    '''
    (file, file, str) -> file
    Take a PAML alignment file and template codeml control file and generate
    a new control file, with parameters seqfile and outfile modified to 
    correspond to the alignment file
    '''

    # check that alignment_file is a text file
    if alignment_file[-4:] == '.txt':
        # get the file name
        name = alignment_file[:-4]
        
    # open control file for reading
    ctlfile = open(control_file, 'r')
        
    # open outputfile file for writing
    newfile = open(destination + name + '_codeml.ctl', 'w')
    
    # go through the control file, copy to new file with new input file and output file
    for line in ctlfile:
        if 'seqfile' not in line and 'outfile' not in line:
            newfile.write(line)
        elif 'seqfile' in line:
            line = line.rstrip().split()
            # replace the file name with the alignment file name
            line[2] = alignment_file
            for item in line[:-1]:
                newfile.write(item + ' ')
            newfile.write(line[-1] + '\n')
        elif 'outfile' in line:
            line = line.rstrip().split()
            # replace the outputfile name with alignment_file.out
            line[2] = name + '.out.txt'
            for item in line[:-1]:
                newfile.write(item + ' ')
            newfile.write(line[-1] + '\n')
    ctlfile.close()
    newfile.close()


# use this function to parse the codeml output files and combine results into a single table
def save_divergence_results_to_file(crem_gff, cla_gff, outputfile, source):
    '''
    (file, file, file) -> file
    Extract the dN and dS values from the codml output files, 
    and save a table to outputfile with columns including the gene names
    and transcript names of each sequence using the GFF annotation files
    '''
    
    # make a list of files in current directory
    files = os.listdir(source)

    # open outfile for writing
    newfile = open(source + outputfile, 'w')
    newfile.write('Cremanei_transcript' + '\t'  + 'Cremanei_gene' + '\t' + 'Clatens_transcript' + '\t' + 'Clatens_gene' + '\t' + 'dN' + '\t' + 'dS' + '\t' + 'dN/dS' + '\n')

    # get latens transcript : gene pairs
    cla_transcripts = transcript_to_gene(cla_gff)

    # get remanei transcript : gene pairs
    crem_transcripts = transcript_to_gene(crem_gff)
    
    # loop over the files
    for filename in files:
        if '.out.txt' in filename:
            infile = open(filename, 'r')
            crem_TS = filename[:filename.index('_CLA')]
            cla_TS = filename[filename.index('CLA'): filename.index('.out')]
            for line in infile:
                if 'tree length for dN' in line:
                    line = line.split()
                    dN = float(line[-1])
                elif 'tree length for dS' in line:
                    line = line.split()
                    dS = float(line[-1])
            crem_gene = crem_transcripts[crem_TS]
            cla_gene = cla_transcripts[cla_TS]
            if dS != 0:
                omega = dN / dS
            elif dS == 0:
                omega = 'NA'
            newfile.write(crem_TS + '\t' + crem_gene + '\t' + cla_TS + '\t' + cla_gene + '\t' + str(dN) + '\t' + str(dS) + '\t' + str(omega) + '\n')
            infile.close()

    newfile.close()






