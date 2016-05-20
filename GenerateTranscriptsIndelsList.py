# -*- coding: utf-8 -*-
"""
Created on Fri Jun 12 17:00:04 2015

@author: Richard
"""

# use this script to generate a set a transcripts with indels in coding sequences in KSR and PX strains
# save all transcript names to a text file


from indel_coordinates import *

# find transcrtips with deletions in KSR
KSR_del = identify_genes_with_indels('356_10172014.gff3', 'KSR_D.compiled_pindel.txt')
print('identified deletions in KSR')
print(len(KSR_del))

# find transcripts with insertions in KSR
KSR_ins = identify_genes_with_indels('356_10172014.gff3', 'KSR_SI.compiled_pindel.txt')
print('identified insertions in KSR')
print(len(KSR_ins))

# find transcripts with deletions in PX
PX_del = identify_genes_with_indels('356_10172014.gff3', 'PX_D.compiled_pindel.txt')
print('identified deletions in PX')
print(len(PX_del))

# find transcripts with insertions in PX
PX_ins = identify_genes_with_indels('356_10172014.gff3', 'PX_SI.compiled_pindel.txt')
print('identified insertions in PX')
print(len(PX_ins))

# create a set of transcripts with indels in KSR
KSR_indels = KSR_del.union(KSR_ins)
# create a set of transcripts with indels in PX
PX_indels = PX_del.union(PX_ins)
# create a set of transcripts with indels in KSR and PX
KSR_PX_indels = KSR_indels.union(PX_indels)

# open file for writing
newfile = open('transcripts_indels_CDS.txt', 'w')
for gene in KSR_PX_indels:
    newfile.write(gene + '\n')
    
# close file
newfile.close()