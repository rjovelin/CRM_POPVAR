
import os

# use this function to align sequence pairs located in separate files in directory
def run_tcoffee_noncoding(directory):
    '''
    (str) -> files
    Align sequences in each of the files contained in directory with t-coffee.
    Save alignments in fasta format
    '''
    
    files = os.listdir(directory)
    
    # fasta files have a .fas
    for filename in files:
        if filename[-4:] == '.fas':
            # align sequences (outpufile ends with .aln)
            alignCommand = "t_coffee " + directory + filename
            print(alignCommand)
            os.system(alignCommand)

# use this fucntion to convert the tcoffee format to fasta        
def convert_tcoffee_to_fasta(filename):
    '''
    (str) -> None
    Read the t-coffee alignment files and save the alignment in fasta format in a text file
    '''

    tcoffee = open(filename, 'r')
    tcoffee.readline()
    tcoffee.readline()

    ali = {}

    for line in tcoffee:
        line = line.rstrip()
        if line != '':
            line = line.split()
            if '*' in line[0]:
                continue
            elif line[0] in ali:
                ali[line[0]] += line[1]
            else:
                ali[line[0]] = line[1]

    alignment = open(filename[:-4] + '.txt', 'w')
    for gene in ali:
        alignment.write('>' + gene + '\n')
        alignment.write(ali[gene] + '\n')

    tcoffee.close()
    alignment.close()

# use this function to convert all tcoffee alignment files located in directory to fasta
def generate_fasta_from_tcoffee(directory):
    '''
    (str) -> None
    Convert each tcoffee files in the directory into a text file in a fasta format
    '''
    # create a list of files in directory
    files = os.listdir(directory)
    # loop over files and convert file to fasta
    for filename in files:
        if filename[-4:] == '.aln':
            convert_tcoffee_to_fasta(directory + filename)

        
# use this function to convert the output of protein alignments from t-coffee to fasta format
def convert_tcoffee_prot_to_fasta(filename):
    '''
    (file) -> file
    Take a t-coffee alignment file for protein sequences and save alignment
    to a new file in fasta format
    Precondition: the alignment contains only 2 sequences
    '''
    
    # open file for reading
    infile = open(filename, 'r')
    # skip 2 first lines
    infile.readline()
    infile.readline()
    
    # create a dict to store the 
    sequences = {}      
    
    # loop over file
    for line in infile:
        if '.' in line or ':' in line or '*' in line:
            continue
        else:
            if line.rstrip() != '':
                line = line.rstrip().split()
                if line[0] in sequences:
                    sequences[line[0]] += line[1]
                else:
                   sequences[line[0]] = line[1]
                   
    infile.close()
    
    alignment = open(filename[:-4] + '.txt', 'w')
    for gene in sequences:
        alignment.write('>' + gene + '\n')
        alignment.write(sequences[gene] + '\n')

    infile.close()
    alignment.close()
    
    
