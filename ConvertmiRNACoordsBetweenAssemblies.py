# -*- coding: utf-8 -*-
"""
Created on Sat May  7 22:16:13 2016

@author: Richard
"""








infile = open('JW_miRNACoords.txt')


code 17
header = infile.readline().rstrip().split()


code 18
mirnas = {}


code 19
for line in infile:
    if line.rstrip() != '':
        line = line.rstrip().split()
        mirnas[line[0]] = [line[1].upper(), line[2], int(line[3])]
        


code 20
mirnas = {}


code 21
infile.close()


code 22
infile = open('JW_miRNACoords.txt')


code 23
header = infile.readline().rstrip().split()


code 24
for line in infile:
    if line.rstrip() != '':
        line = line.rstrip().split('\t')
        mirnas[line[0]] = [line[1].upper(), line[2], int(line[3])]
        


code 25
infile.close()


code 26
genome = convert_fasta('noamb_356_v1_4.txt')


code 27
sequences = {}


code 28
for mir in mirnas:
    LG = mirnas[mir][1]
    pos = mirnas[mir][3] -1
    seq = genome[LG][pos:pos + len(mirnas[mir][0])]
    revseq = reverse_complement(seq)
    sequences[mir] = [mir, seq, revseq]
    


code 29
sequences = {}


code 30
for mir in mirnas:
    LG = mirnas[mir][1]
    pos = mirnas[mir][2] -1
    seq = genome[LG][pos:pos + len(mirnas[mir][0])]
    revseq = reverse_complement(seq)
    sequences[mir] = [mir, seq, revseq]
    


code 31
len(sequences)


code 32
names = [i for i in sequences]


code 33
sequences[names[0]]


code 34
sequences = {}


code 35
for mir in mirnas:
    LG = mirnas[mir][1]
    pos = mirnas[mir][2] -1
    seq = genome[LG][pos:pos + len(mirnas[mir][0])]
    revseq = reverse_complement(seq)
    sequences[mir] = [mirnas[mir][0], seq, revseq]
    


code 36
names = [i for i in sequences]


code 37
sequences[names[0]]


code 38
differences = {}


code 39
from divergence import *


code 40
from divergence import *


code 41
from divergence import *


code 42
from divergence import *


code 43
for i in mirnas:
    diff1 = match_diff(sequences[i][0], sequences[i][1])
    diff2 = match_diff(sequences[i][0], sequences[i][2])
    differences[i] = [diff1, diff2]
    


code 44
for i in differences:
    print(i, differences[i][0], differences[i][1], sep = '\t')
    


code 45
for i in differences:
    if differences[i][0] > 10 and differences[i][1] > 10:
        print(i)
        print(sequences[i][0], sequences[i][1], sequences[i][2], end = '', sep = '\n')
        print('\n')
        


code 46
total = 0


code 47
for i in differences:
    if differences[i][0] > 10 and differences[i][1] > 10:
        total += 1
        print(i)
        print(sequences[i][0], sequences[i][1], sequences[i][2], end = '', sep = '\n')
        print('\n')
        


code 48
total


code 49
print(1)


code 50
print([1])




