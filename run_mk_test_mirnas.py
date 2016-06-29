# -*- coding: utf-8 -*-
"""
Created on Wed Jun 29 16:07:21 2016

@author: RJovelin
"""

# use this script to compute the MK test for miRNAs



# need a script to count P and D for non-coding sites



# extract from Lyu et al Plos genetics 2014
# New MicroRNAs in Drosophila—Birth, Death and Cycles of Adaptive Evolution

'''
We used the McDonald-Kreitman test (MK test) [32] framework
to detect positive selection in miRNAs from each age group
based on the polymorphisms within D. melanogaster and the
divergence between D. melanogaster and D. simulans. Precursor or
mature sequences of each miRNA group were combined and
treated as the functional category, while the 4-fold degenerate sites
in the whole genome were used as the neutral control. The
divergence is calculated by counting the number of changed
nucleotide sites between D. melanogaster (dm3) and D. simulans
(droSim1) based on the UCSC whole genome alignment.
Polymorphism data was retrieved from Drosophila Population
Genomics Project (DPGP, http://www.dpgp.org/, release 1.0).

SNPs that were detected on more than thirty individuals and
exhibited a derived allele frequency (DAF) > 5% were used for the
MK test.
The proportion of adaptively fixed mutations (a) was estimated
as previously described [75]. To estimate the evolutionary fate of
each miRNA, we first screened for adaptive miRNAs among the
238 candidates by using each miRNA’s precursor together with
the 50 bp flanking sequences on both sides as the functional sites.
The p-values of multiple MK tests were adjusted by the BenjaminiHochberg
method [76] and the adaptive significance of each
candidate is re-validated by using the precursor alone in the MK
test. We then identified the conservative miRNAs by comparing
the number of substitutions in the miRNA precursors (KmiR) with
the number of substitutions in the synonymous sites (KS) between
D.melanogaster and D.simulans. miRNAs with KmiR/KS < 0.5 were
considered to be conservatively evolving. Kimura’s 2-parameter
model [72] and the Nei-Gojobori model [77] were used to
calculate KmiR and KS, respectively. Finally, excluding the
adaptive and conservative miRNAs, the remaining were considered
to be in transition between adaptive to conservative/death
'''