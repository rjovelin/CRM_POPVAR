# -*- coding: utf-8 -*-
"""
Created on Mon Aug 10 09:59:13 2015

@author: RJovelin
"""


from piRNAs import *
from accessories import *
from get_coding_sequences import *
from Cel_UTR import *
from miRNA_target import *
from cel_UTR_length import *

# convert genome fasta to dict
genome = convert_fasta('../CREM_CLA_protein_divergence/noamb_356_v1_4.txt')

# get piRNA coordinates
pirna_coord = get_pirna_loci('PX356_piRNA_coord.txt') 


print('got piRNA coord')

# get CDS_coord {TS1: [chromo, sense, [(s1, end1), (s2, end2)]]}
CDS_coord = get_CDS_positions('../CREM_CLA_protein_divergence/356_10172014.gff3')
# create new dict in the form {chromo: {set of positions}}
CDS_pos = {}
for gene in CDS_coord:
    # get chromo
    chromo = CDS_coord[gene][0]
    # loop over cds coordinates 
    for i in range(len(CDS_coord[gene][2])):
        # get start and end positions in a list
        # convert to 0-based index
        start  = CDS_coord[gene][2][i][0] - 1
        end = CDS_coord[gene][2][i][1]
        # check if chromo in CDS_pos
        if chromo in CDS_pos:
            for j in range(start, end):
                CDS_pos[chromo].add(j)
        else:
            CDS_pos[chromo] = set()
            for j in range(start, end):
                CDS_pos[chromo].add(j)

print('got CDS coord')

# get intron coordinates
# create dict
intron_pos = {}
# loop over genes in CDS_pos
for gene in CDS_coord:
    # get chromo
    chromo = CDS_coord[gene][0]
    # loop over cds coordinates for that gene
    for i in range(0, len(CDS_coord[gene][2]) -1 ):
        # convert intron coordinates to 0-based index
        start = CDS_coord[gene][2][i][1]
        end = CDS_coord[gene][2][i+1][0]
        # check if chromo in intron pos
        if chromo in intron_pos:
            for j in range(start, end):
                intron_pos[chromo].add(j)
        else:
            intron_pos[chromo] = set()
            for j in range(start, end):
                intron_pos[chromo].add(j)
                
   
print('got intron coord')

# get UTR coordinates

# compute threshold based on the distribution of elegans UTR length
UTR_length = celegans_three_prime_UTR_length('../miRNA_Target_sites/c_elegans.PRJNA13758.WS248.annotations.gff3')
threshold = get_percentile(UTR_length, 99)

print('determined threshold based on Celegans UTRs: ', threshold)

# get UTR coord {TS1 : [chromo, start, end, orientation]}
three_prime = get_three_prime_UTR_positions('../CREM_CLA_protein_divergence/356_10172014.gff3', '../CREM_CLA_protein_divergence/noamb_356_v1_4.txt', threshold)
# create a dict UTR_pos
UTR_pos = {}
# loop over genes in three_prime
for gene in three_prime:
    # get chromo
    chromo = three_prime[gene][0]
    # convert to 0-based
    start = three_prime[gene][1] -1
    end = three_prime[gene][2]
    # check if chromo in UTR_pos
    if chromo in UTR_pos:
        for j in range(start, end):
            UTR_pos[chromo].add(j)
    else:
        UTR_pos[chromo] = set()
        for j in range(start, end):
            UTR_pos[chromo].add(j)
            
print('got UTR coord')

# find pirna locations
locations = find_pirna_locations(genome, pirna_coord, CDS_pos, UTR_pos, intron_pos)

print('got the piRNA locations')

# open file for writing
newfile = open('piRNA_site_locations.txt', 'w')

# get the number of pirna loci
N_pirnas = 0
for chromo in pirna_coord:
    for i in range(len(pirna_coord[chromo])):
        N_pirnas += 1
print('number of pirnas loci: {0}'.format(N_pirnas))

newfile.write('Number of piRNA loci in C. remanei:' + '\t' + str(N_pirnas) + '\n')
newfile.write('\n')
newfile.write('Number of piRNAs in different site categories:\n')
newfile.write('-' * 47 + '\n')
# loop over keys in locations dict
for key in locations:
    newfile.write(key+ ':' + '\t' + str(locations[key]) + '\n')
newfile.write('\n')

# verify the pirna counts
total = 0
for key in locations:
    total += locations[key]
    
print('number of pirnas in sites: {0}'.format(total))
    
newfile.write('total piRNAs in site categories:' + '\t' + str(total) + '\n')
missing = N_pirnas - total
newfile.write('number of pirnas not allocated to sites:' + '\t' + str(missing) + '\n')




