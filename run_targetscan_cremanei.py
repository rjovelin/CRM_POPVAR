# -*- coding: utf-8 -*-
"""
Created on Fri May 29 12:45:56 2015

@author: Richard
"""

import os


# run TargetScan on Cremanei transcripts only
os.system('perl targetscan_60.pl Cremanei_miRfam_mirbase21_info.txt Crm_UTR_seq_targetscan.txt Cremanei_predicted_sites_targetscan.txt')

# run TargetScan on Cremanei and Clatens aligned orthologous UTRs
os.system('perl targetscan_60.pl Cremanei_miRfam_mirbase21_info.txt Crm_Cla_UTR_seq_targetscan.txt Crm_Cla_predicted_sites_targetscan.txt')

# run TargetScan on Cremanei and Celegans aligned orthologous files 
os.system('perl targetscan_60.pl Cremanei_miRfam_mirbase21_info.txt Crm_Cel_UTR_seq_targetscan.txt Crm_Cel_predicted_sites_targetscan.txt')





