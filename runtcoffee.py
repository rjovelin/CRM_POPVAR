
# run to-coffee command lines to align protein sequences and
# get the codon-based DNA alignment
# usage: have this script in the same directory as directory pairs
# which contains the input fasta files

import os

files = os.listdir('./pairs')

# fasta files have a .tfa extension
for filename in files:
    # get a outputfile name for the protein translation
    protein = filename[:-3] + '_p.tfa'
    # translate dna sequence into protein sequence
    dnaToPepCommand = "t_coffee -other_pg seq_reformat -in " + './pairs/' + filename + " -action +translate -output fasta_seq > " + protein
    print(dnaToPepcommand)
    os.system(dnaToPepcommand)
    # when sequences have different length, t-coffee padds the short sequence with 'o'
    # replace 'o' with '-'
    os.system("sed -i 's/o/-/g '" + protein)

    # align protein sequences
    alignPepCommand = "t_coffee " + protein
    print(alignPepCommand)
    os.system(alignPepCommand)

    # get name for aligned proteins
    protein_ali = protein[:-3] + '.aln'
    # get name for aligned DNA
    DNA_ali = filename[:-3] +'.aln'
    backtransCommand = 't_coffee -other_pg seq_reformat -in ' + './pairs/' + filename + ' -in2 ' + './pairs/' + protein_ali + ' -action +thread_dna_on_prot_aln' +\
                       ' -output clustalw > ' + DNA_ali
    print(backtransCommand)
    os.system(backtransCommand)



    # convert clustal format to fasta format
    DNA_ali_fasta = filename[:-3] + '_aln.tfa'
    convertDNAalignToFastaCommand = 't_coffee -other_pg seq_reformat -in' + './pairs'/ + DNA_ali + '-output fasta_aln > ' + DNA_ali_fasta
    print(convertDNAalignToFastaCommand)
    os.system(convertDNAalignToFastaCommand)
    
     
    
    
    
