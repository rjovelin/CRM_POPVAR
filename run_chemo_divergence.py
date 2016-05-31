# -*- coding: utf-8 -*-
"""
Created on Mon Jun 22 15:36:59 2015

@author: Richard
"""

# use this script to compute nucleotide divergence in membrane paritions of chemoreceptor genes


from protein_divergence import *
from chemoreceptors import *
import os


# get the alignment files in differnt partitions
TM_files = [filename for filename in os.listdir('./Partitions/Membrane/') if 'TM' in filename]
Out_files = [filename for filename in os.listdir('./Partitions/Outside/') if 'outside' in filename]
In_files = [filename for filename in os.listdir('./Partitions/Inside/') if 'inside' in filename]
Extra_files =[filename for filename in os.listdir('./Partitions/Extra_membrane/') if 'ExtraTM' in filename]

# generate codeml control files for each alignment files in each partition folder
for alignment_file in TM_files:
    generate_codeml_control_file(alignment_file, '../CREM_CLA_protein_divergence/codeml.ctl', './Partitions/Membrane/')
for alignment_file in Out_files:
    generate_codeml_control_file(alignment_file, '../CREM_CLA_protein_divergence/codeml.ctl', './Partitions/Outside/')
for alignment_file in In_files:
    generate_codeml_control_file(alignment_file, '../CREM_CLA_protein_divergence/codeml.ctl', './Partitions/Inside/')
for alignment_file in Extra_files:
    generate_codeml_control_file(alignment_file, '../CREM_CLA_protein_divergence/codeml.ctl', './Partitions/Extra_membrane/')

# run codeml on each partition
# get the the list of control files
TM_ctl = [filename for filename in os.listdir('./Partitions/Membrane/') if 'ctl' in filename]
Out_ctl = [filename for filename in os.listdir('./Partitions/Outside/') if 'ctl' in filename]
In_ctl = [filename for filename in os.listdir('./Partitions/Inside/') if 'ctl' in filename]
Extra_ctl = [filename for filename in os.listdir('./Partitions/Extra_membrane/') if 'ctl' in filename]

# copy tree files to directories
os.system('cp CREMCLATREE.tre.txt ./Partitions/Membrane/')
os.system('cp CREMCLATREE.tre.txt ./Partitions/Outside/')
os.system('cp CREMCLATREE.tre.txt ./Partitions/Inside/')
os.system('cp CREMCLATREE.tre.txt ./Partitions/Extra_membrane/')

# change directory to run codeml
os.chdir('./Partitions/Membrane/')
# loop over codml control files, run codeml
for filename in TM_ctl:
    os.system('codeml ' + filename)
# remove tree file
os.system('rm CREMCLATREE.tre.txt')

# change directory
os.chdir('../Outside')
# loop over control files, run codeml
for filename in Out_ctl:
    os.system('codeml ' + filename)
# remove tree file
os.system('rm CREMCLATREE.tre.txt')

# change directory
os.chdir('../Inside')
# loop over control files, run codeml
for filename in In_ctl:
    os.system('codeml ' + filename)
# remove tree file
os.system('rm CREMCLATREE.tre.txt')
    
# change directory
os.chdir('../Extra_membrane/')
# loop over control files, run codeml
for filename in Extra_ctl:
    os.system('codeml ' + filename)
# remove tree file
os.system('rm CREMCLATREE.tre.txt')


