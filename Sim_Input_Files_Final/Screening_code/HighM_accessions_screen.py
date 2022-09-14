import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.ticker as mtick
from scipy import stats
import scipy.stats
#import random
import pandas as pd
from timeit import default_timer as timer
import os
import csv
import bisect # to calculate local CG density
import dask.dataframe as dd
import argparse
import sys


#initial states folder:
# RL
initial_states_folder = "m1001_states_Relict"
file_name_output_Accessions = "Accessions_RelictsFilt_MinCov_0p7.tsv"
required_accessions_codes_file = os.path.join('..',"Relicts_filtered_IDs.txt" )
required_accessions_codes_df = pd.read_csv(required_accessions_codes_file,sep='\t')

# NS
# initial_states_folder = "m1001_states_NS"
# file_name_output_Accessions = "Accessions_NorthenSwedishFilt_MinCov_0p7.tsv"
# required_accessions_codes_file = os.path.join('..',"Northern_Swedish_filtered_IDs.csv" )
# required_accessions_codes_df = pd.read_csv(required_accessions_codes_file, usecols=['SRA_Accession','Genome_coverage'])


# create file
open(file_name_output_Accessions,'w').close


min_coverage = 0.7



# load in list of accessions codes: accessions_codes_file
required_accessions_codes_df = required_accessions_codes_df.set_index('SRA_Accession')
required_accessions_codes_list = required_accessions_codes_df.index.tolist()
#print(required_accessions_codes_df.head())


# Make list of codes for available initial state files
available_initial_state_codes_list = []
# synchronised list to store Decile number
available_initial_state_deciles_list = []

# sample filename in folder: m1001_states_all
# m1001_samples_meth_status_mCG_density_0_SRX1664742_3.tsv
filename_start = 'm1001_samples_meth_status_mCG_density_'
filename_end = '.tsv'

# # loop over all available initial state files
# for filename in os.listdir(os.path.join('..',"m1001_states_all")):
#     # trim down to accessions code only:
#     filename_CodeOnly = filename[40:]
#     filename_CodeOnly = filename_CodeOnly[:-6]
#     if '_' in filename_CodeOnly:
#         filename_CodeOnly = filename_CodeOnly[1:]
#         available_initial_state_codes_list.append(filename_CodeOnly)
#         available_initial_state_deciles_list.append(filename[38:40])
#         print(filename)
#         print(filename_CodeOnly)
#         print(filename[38:40])
#         print()
#     else:
#         available_initial_state_codes_list.append(filename_CodeOnly)
#         available_initial_state_deciles_list.append(filename[38])
#         print(filename)
#         print(filename_CodeOnly)
#         print(filename[38])
#         print()

# m1001_samples_meth_status_mCG_density_NS_SRX445925.tsv
# m1001_samples_meth_status_mCG_density_RL_SRX1665031.tsv
# loop over all available initial state files
for filename in os.listdir(os.path.join('..',initial_states_folder)):
    # trim down to accessions code only:
    filename_CodeOnly = filename[40:]
    filename_CodeOnly = filename_CodeOnly[:-4]
    if '_' in filename_CodeOnly:
        filename_CodeOnly = filename_CodeOnly[1:]
        available_initial_state_codes_list.append(filename_CodeOnly)
        available_initial_state_deciles_list.append(filename[38:40])
        print(filename)
        print(filename_CodeOnly)
        print(filename[38:40])
        print()
    else:
        available_initial_state_codes_list.append(filename_CodeOnly)
        available_initial_state_deciles_list.append(filename[38:40])
        print(filename)
        print(filename_CodeOnly)
        print(filename[38:40])
        print()

# make dataframe of available accessions
available_accesions_df = pd.DataFrame({'Code':available_initial_state_codes_list,'Decile':available_initial_state_deciles_list}, columns=['Code','Decile'])
available_accesions_df = available_accesions_df.set_index('Code')
print(available_accesions_df)


# make lists of:
required_codes_available = []
required_filenames_available = []
required_codes_missing = []

for i_ in range(len(required_accessions_codes_list)):
    test_code = required_accessions_codes_list[i_]
    #print(test_code)
    #print(required_accessions_codes_df.loc[test_code,'Genome_coverage'])
    if required_accessions_codes_df.loc[test_code,'Genome_coverage'] > min_coverage:
        if test_code in available_initial_state_codes_list:
            required_codes_available.append(test_code)
            required_filenames_available.append(filename_start+str( available_accesions_df.loc[test_code,'Decile'] )+'_'+test_code+filename_end)
            print(i_,'yes!!!!', available_initial_state_deciles_list[i_], str(available_accesions_df.loc[test_code,'Decile']))
        else:
            required_codes_missing.append(test_code)



print('already have')
print(available_initial_state_codes_list)
print()
print('will use')
print(required_filenames_available)
print()
print('missing')
print(required_codes_missing)
print()
print(available_initial_state_deciles_list)

# write out list of required filenames
for i_ in range(len(required_filenames_available)):
    with open(file_name_output_Accessions,'a',newline='') as output_file:
        line_writer = csv.writer(output_file,delimiter='\t')
        line_writer.writerow([required_filenames_available[i_]])
        print([required_filenames_available[i_]])


# read file back in to test

temp_input_df = pd.read_csv(file_name_output_Accessions,sep='\t',header=None)
print()
print()
print(temp_input_df.head())
temp_input_list = temp_input_df[0].tolist()
print(temp_input_list)
print(len(temp_input_list))

