# -*- coding: utf-8 -*-
"""
Created on Fri Jun 12 22:50:32 2015

@author: Richard
"""

#!/usr/bin/env python3

from indel_coordinates import *

# use this small function to print the number of indel in the indel dict
def print_N_indels(indel_coord):
    total = 0
    for chromo in indel_coord:
        for i in indel_coord[chromo]:
            total += 1
    print(total)


# open file for writing
newfile = open('summary_indel_counts.txt', 'w')

newfile.write('Number of genes with indels in CDS\n')
newfile.write('-' * 35 + '\n')

# find transcrtips with deletions in KSR
KSR_del = identify_genes_with_indels('../CREM_CLA_protein_divergence/356_10172014.gff3', 'KSR_D.compiled_pindel.txt')
print('identified deletions in KSR')
newfile.write('genes with CDS deletions in KSR:' + '\t' + str(len(KSR_del)) + '\n')

# find transcripts with insertions in KSR
KSR_ins = identify_genes_with_indels('../CREM_CLA_protein_divergence/356_10172014.gff3', 'KSR_SI.compiled_pindel.txt')
print('identified insertions in KSR')
newfile.write('genes with CDS insertions in KSR:' + '\t' + str(len(KSR_ins)) + '\n')

# find transcripts with deletions in PX
PX_del = identify_genes_with_indels('../CREM_CLA_protein_divergence/356_10172014.gff3', 'PX_D.compiled_pindel.txt')
print('identified deletions in PX')
newfile.write('genes with CDS deletions in PX:' + '\t' + str(len(PX_del)) + '\n')

# find transcripts with insertions in PX
PX_ins = identify_genes_with_indels('../CREM_CLA_protein_divergence/356_10172014.gff3', 'PX_SI.compiled_pindel.txt')
print('identified insertions in PX')
newfile.write('genes with CDS insertions in PX:' + '\t' + str(len(PX_ins)) + '\n')

# find transcripts with deletions in PB
PB_del = identify_genes_with_indels('../CREM_CLA_protein_divergence/356_10172014.gff3', 'PB_D.compiled_pindel.txt')
print('identified deletions in PB')
newfile.write('genes with CDS deletions in PB:' + '\t' + str(len(PB_del)) + '\n')

# find transcripts with insertions in PB
PB_ins = identify_genes_with_indels('../CREM_CLA_protein_divergence/356_10172014.gff3', 'PB_SI.compiled_pindel.txt')
print('identified insertions in PB')
newfile.write('genes with CDS instertions in PB:' + '\t' + str(len(PB_ins)) + '\n')


# create a set of transcripts with indels in KSR
KSR_indels = KSR_del.union(KSR_ins)
# create a set of transcripts with indels in PX
PX_indels = PX_del.union(PX_ins)
# create a set of transcripts with indels in KSR and PX
KSR_PX_indels = KSR_indels.union(PX_indels)
# create a set of transcripts with indels in PB
PB_indels = PB_del.union(PB_ins)
# create a set of transcripts with indels in all strains
crm_indels = PB_indels.union(KSR_PX_indels)

# record number of genes with indels in KSR and PX 
newfile.write('number of genes with CDS indels in KSR and PX combined:' + '\t' + str(len(KSR_PX_indels)) + '\n')
# record the number of genes with indels in remanei
newfile.write('total number of genes with CDS indels in remanei:' + '\t' + str(len(crm_indels)) + '\n')


# create a dict of repeat name coordinate
repeat_coord = get_repeats_coord('../piRNAs/356_v1_4.fasta.out', False)
# create a dict of chromo : set of repeat indices
repeat_pos = get_repeat_positions(repeat_coord)

# get the coordinates of deletions in KSR
KSR_del_coord = get_short_indel_coordinates('KSR_D.compiled_pindel.txt')
# remove indels in repeat
KSR_del_coord = remove_indels_in_repeats(KSR_del_coord, repeat_pos)
print(len(KSR_del_coord))
print_N_indels(KSR_del_coord) 

# get the coordinates of instertions in KSR
KSR_ins_coord = get_short_indel_coordinates('KSR_SI.compiled_pindel.txt')
# remove indels in repeat
KSR_ins_coord = remove_indels_in_repeats(KSR_ins_coord, repeat_pos)
print(len(KSR_del_coord))
print_N_indels(KSR_ins_coord)

# get the coordinates of deletions in PX
PX_del_coord = get_short_indel_coordinates('PX_D.compiled_pindel.txt')
# remove indels in repeat
PX_del_coord = remove_indels_in_repeats(PX_del_coord, repeat_pos)
print(len(PX_del_coord))
print_N_indels(PX_del_coord)

# get the coordinates of instertions in PX
PX_ins_coord = get_short_indel_coordinates('PX_SI.compiled_pindel.txt')
# remove indels in repeat
PX_ins_coord = remove_indels_in_repeats(PX_ins_coord, repeat_pos)
print(len(PX_ins_coord))
print_N_indels(PX_ins_coord)

# get the coordinates of deletions in PB
PB_del_coord = get_short_indel_coordinates('PB_D.compiled_pindel.txt')
# remove indels in repeat
PB_del_coord = remove_indels_in_repeats(PB_del_coord, repeat_pos)
print(len(PB_del_coord))
print_N_indels(PB_del_coord)

# get the coordinates of instertions in PB
PB_ins_coord = get_short_indel_coordinates('PB_SI.compiled_pindel.txt')
# remove indels in repeat
PB_ins_coord = remove_indels_in_repeats(PB_ins_coord, repeat_pos)
print(len(PB_ins_coord))
print_N_indels(PB_ins_coord)

# record the number of insertions and deletions in file
newfile.write('\n')
newfile.write('total number of indels, including overlapping indels:\n')
newfile.write('-' * 54 + '\n')
newfile.write('number of deletions in KSR:' + '\t' + str(count_short_indels(KSR_del_coord)) + '\n')
newfile.write('number of insertions in KSR:' + '\t' + str(count_short_indels(KSR_ins_coord)) + '\n')
newfile.write('number of indels in KSR:' + '\t' + str(count_short_indels(KSR_del_coord) + count_short_indels(KSR_ins_coord)) + '\n')
newfile.write('number of deletions in PX:' + '\t' + str(count_short_indels(PX_del_coord)) + '\n')
newfile.write('number of insertions in PX:' + '\t' + str(count_short_indels(PX_ins_coord)) + '\n')
newfile.write('number of indels in PX:' + '\t' + str(count_short_indels(PX_del_coord) + count_short_indels(PX_ins_coord)) + '\n')
newfile.write('number of deletions in PB:' + '\t' + str(count_short_indels(PB_del_coord)) + '\n')
newfile.write('number of insertions in PB:' + '\t' + str(count_short_indels(PB_ins_coord)) + '\n')
newfile.write('number of indels in PB:' + '\t' + str(count_short_indels(PB_del_coord) + count_short_indels(PB_ins_coord)) + '\n')

# create a dict with chromo as key and combined unique deletions in KSR and PX
KSRPX_del_coord = combine_indels(KSR_del_coord, PX_del_coord)
# create a dict with chromo as key and combined unique insertions in KSR and PX
KSRPX_ins_coord = combine_indels(KSR_ins_coord, PX_ins_coord)
# create a dict with chromo as key and combined unique indels in KSR and PX
KSRPX_indel_coord = combine_indels(KSRPX_del_coord, KSRPX_ins_coord)
# create a dict by with chromo as key and combined unique deletions in KSR_PX and PB
rem_del_coord = combine_indels(KSRPX_del_coord, PB_del_coord)
# create a dict with chromo as key and combined unique insertions in KSR_PX and PB
rem_ins_coord = combine_indels(KSRPX_ins_coord, PX_ins_coord)
# create a dict with chromo as key and combined unique indels in remanei
rem_indel_coord = combine_indels(rem_del_coord, rem_ins_coord)

# get the number of deleletions in KSR and PX combined
newfile.write('number of unique deletions in KSR and PX:' + '\t' + str(count_short_indels(KSRPX_del_coord)) + '\n')

# get the number of insertions in KSR and PX combined
newfile.write('number of unique insertions in KSR and PX:' + '\t' + str(count_short_indels(KSRPX_ins_coord)) + '\n')

# get the number indels in KSR and PX combined
newfile.write('number of unique indels in KSR and PX:' + '\t' + str(count_short_indels(KSRPX_indel_coord)) + '\n')

# get the number of deletions in remanei
newfile.write('number of unique deletions in remanei:' + '\t' + str(count_short_indels(rem_del_coord)) + '\n')

# get the number of insertions in remanei
newfile.write('number of unique insertions in remanei:' + '\t' + str(count_short_indels(rem_ins_coord)) + '\n')

# get the number of indels in remanei
newfile.write('number of unique indels in remanei:' + '\t' + str(count_short_indels(rem_indel_coord)) + '\n')

newfile.write('\n')
newfile.write('indel length\n')
newfile.write('-' * 13 + '\n')

# get indel length for indels in KSR and PX
# first need to convert set to list
for chromo in KSRPX_indel_coord:
    KSRPX_indel_coord[chromo] = list(KSRPX_indel_coord[chromo])
KSRPX_indel_length = indel_length(KSRPX_indel_coord)
newfile.write('indel length range in KSR and PX:' + '\t' + str(min(KSRPX_indel_length)) + '\t' + str(max(KSRPX_indel_length)) + '\n')

# get indel length for indels in remanei
# first need to convert set to list
for chromo in rem_indel_coord:
    rem_indel_coord[chromo] = list(rem_indel_coord[chromo])
rem_indel_length = indel_length(rem_indel_coord)
newfile.write('indel length range in remanei:' + '\t' + str(min(rem_indel_length)) + '\t' + str(max(rem_indel_length)) + '\n')

newfile.write('\n')
newfile.write('indels in set of transcripts used in diversity analyses:\n')
newfile.write('-' * 57 + '\n')

# get the number of transcripts with indels, from the set of transcripts used in diversity analyses
# create a dict of transcripts: list of indel coordinates for the KSR+PX group
KSRPX_indel_CDS = indels_in_CDS('../CREM_CLA_protein_divergence/356_10172014.gff3', '../CREM_CLA_protein_divergence/unique_transcripts.txt', KSRPX_indel_coord)

newfile.write('number of genes with indels in CDS in KSR+PX combined:' + '\t' + str(len(KSRPX_indel_CDS)) + '\n')
newfile.write('number of indels in CDS in KSR+PX combined:' + '\t' + str(count_short_indels(KSRPX_indel_CDS)) + '\n')

# create a dict of transcripts : list of indel coordinates for remanei
rem_indel_CDS = indels_in_CDS('../CREM_CLA_protein_divergence/356_10172014.gff3', '../CREM_CLA_protein_divergence/unique_transcripts.txt', rem_indel_coord)

newfile.write('number of genes with indels in CDS in remanei:' + '\t' + str(len(rem_indel_CDS)) + '\n')
newfile.write('number of indels in CDS in remanei:' + '\t' + str(count_short_indels(rem_indel_CDS)) + '\n')

# close file after writing
newfile.close()

