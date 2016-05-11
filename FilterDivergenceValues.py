# -*- coding: utf-8 -*-
"""
Created on Wed May 11 11:29:46 2016

@author: RJovelin
"""

import os
import sys


divergFile = sys.argv[1]


# use this script to filter weird divergence values in the remanei-latens codeml divergence file

# set threshold for dN, dS and omega, remove genes instead of removing values
# genes are removed if any divergence estimate is greater than following thresholds
dN = 2
dS = 1.5
omega = 5

# open file for reading
infile = open(divergFile)
infile.readline()




  