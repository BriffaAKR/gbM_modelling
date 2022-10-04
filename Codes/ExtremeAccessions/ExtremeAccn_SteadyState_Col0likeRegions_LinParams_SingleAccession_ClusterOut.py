# Uses 'Col0' variables to analysise data from a single accession. 
# Also use this code to look at accessions that are not Col0 however. 

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

import input_params as params # import my params from separate python file

# set up argument read-ins
parser=argparse.ArgumentParser()
parser.add_argument('batch_no_in')
parser.add_argument('N_batch_in')
input_args = parser.parse_args()

i_batch = int(input_args.batch_no_in)
N_batch = int(input_args.N_batch_in) # plan to use N_batch = 250
N_total_IDs = params.N_total_IDs
batch_length = np.ceil(N_total_IDs/N_batch)
N_reps = params.N_reps
N_corrns_bins = params.N_corrns_bins

# define paramter coding for filenames
filename_params_code = params.filename_params_code

## define annoation files
InputFiles_folder = params.InputFiles_folder
# file_in_accessions_codes = params.file_in_accessions_codes

file_in_annotation = params.file_in_annotation
file_in_anno_filt_1 = params.file_in_anno_filt_1
file_in_anno_filt_2 = params.file_in_anno_filt_2
file_in_exclude_filt_1 = params.file_in_exclude_filt_1

chosen_accession_file = params.chosen_accession_file

# build up input path string
path_string = os.path.join( os.getcwd(), '..')
path_string = os.path.join( path_string, '..')
#path_string = os.path.join( path_string, '..')
#path_string = os.path.join( path_string, '..')
path_string = os.path.join( path_string, InputFiles_folder)

# accessions_codes_file = os.path.join( path_string, file_in_accessions_codes )

annotaion_file = os.path.join( path_string, file_in_annotation )
filt_1_file = os.path.join( path_string, file_in_anno_filt_1)
filt_2_file = os.path.join( path_string, file_in_anno_filt_2) 
exclude_filt_1_file = os.path.join( path_string, file_in_exclude_filt_1)
#file_in_Combined_Data_per_genes_meth_statuses = os.path.join( path_string, 'SeqRun9All_CG_2021-01-18_WTconsensus_meth_status_allgenes_withID.tsv' )
file_in_Combined_Data_per_genes_meth_statuses = chosen_accession_file

# make paths for output files
locus_type = params.locus_type
TE_filt = params.TE_filt
filename_start = params.filename_start

filename_ending = 'batch_'+str(i_batch)+'.tsv'

file_name_output_TEST = os.path.join('Output_batch_files/Output_TEST', filename_start + 'TEST_'+filename_ending)
# file_name_output_SimulatedStates_ = os.path.join('Output_batch_files/SimulatedStates', filename_start + 'SimState_')
# file_name_output_PairSeparationsMM_Sim = os.path.join('Output_batch_files/Correlations_Sim/PairSeparationsMM_Sim', filename_start + 'PairSeparationsMM_Sim_'+filename_ending)
# file_name_output_PairSeparationsMU_Sim = os.path.join('Output_batch_files/Correlations_Sim/PairSeparationsMU_Sim', filename_start + 'PairSeparationsMU_Sim_'+filename_ending)
# file_name_output_PairSeparationsUU_Sim = os.path.join('Output_batch_files/Correlations_Sim/PairSeparationsUU_Sim', filename_start + 'PairSeparationsUU_Sim_'+filename_ending)
# file_name_output_PairSeparationsXX_Sim = os.path.join('Output_batch_files/Correlations_Sim/PairSeparationsXX_Sim', filename_start + 'PairSeparationsXX_Sim_'+filename_ending)
# file_name_output_PairSeparationsMM_Col0 = os.path.join('Output_batch_files/Correlations_Col0/PairSeparationsMM_Col0', filename_start + 'PairSeparationsMM_Col0_'+filename_ending)
# file_name_output_PairSeparationsMU_Col0 = os.path.join('Output_batch_files/Correlations_Col0/PairSeparationsMU_Col0', filename_start + 'PairSeparationsMU_Col0_'+filename_ending)
# file_name_output_PairSeparationsUU_Col0 = os.path.join('Output_batch_files/Correlations_Col0/PairSeparationsUU_Col0', filename_start + 'PairSeparationsUU_Col0_'+filename_ending)
# file_name_output_PairSeparationsXX_Col0 = os.path.join('Output_batch_files/Correlations_Col0/PairSeparationsXX_Col0', filename_start + 'PairSeparationsXX_Col0_'+filename_ending)
# file_name_output_PairSeparationsMM_Data = os.path.join('Output_batch_files/Correlations_Data/PairSeparationsMM_Data', filename_start + 'PairSeparationsMM_Data_'+filename_ending)
# file_name_output_PairSeparationsMU_Data = os.path.join('Output_batch_files/Correlations_Data/PairSeparationsMU_Data', filename_start + 'PairSeparationsMU_Data_'+filename_ending)
# file_name_output_PairSeparationsUU_Data = os.path.join('Output_batch_files/Correlations_Data/PairSeparationsUU_Data', filename_start + 'PairSeparationsUU_Data_'+filename_ending)
# file_name_output_PairSeparationsXX_Data = os.path.join('Output_batch_files/Correlations_Data/PairSeparationsXX_Data', filename_start + 'PairSeparationsXX_Data_'+filename_ending)
# file_name_output_PairSeparationsMM_Dset2 = os.path.join('Output_batch_files/Correlations_Dset2/PairSeparationsMM_Dset2', filename_start + 'PairSeparationsMM_Dset2_'+filename_ending)
# file_name_output_PairSeparationsMU_Dset2 = os.path.join('Output_batch_files/Correlations_Dset2/PairSeparationsMU_Dset2', filename_start + 'PairSeparationsMU_Dset2_'+filename_ending)
# file_name_output_PairSeparationsUU_Dset2 = os.path.join('Output_batch_files/Correlations_Dset2/PairSeparationsUU_Dset2', filename_start + 'PairSeparationsUU_Dset2_'+filename_ending)
# file_name_output_PairSeparationsXX_Dset2 = os.path.join('Output_batch_files/Correlations_Dset2/PairSeparationsXX_Dset2', filename_start + 'PairSeparationsXX_Dset2_'+filename_ending)
# file_name_output_PairSeparationsMM_D5 = os.path.join('Output_batch_files/Correlations_D5/PairSeparationsMM_D5', filename_start + 'PairSeparationsMM_D5_'+filename_ending)
# file_name_output_PairSeparationsMU_D5 = os.path.join('Output_batch_files/Correlations_D5/PairSeparationsMU_D5', filename_start + 'PairSeparationsMU_D5_'+filename_ending)
# file_name_output_PairSeparationsUU_D5 = os.path.join('Output_batch_files/Correlations_D5/PairSeparationsUU_D5', filename_start + 'PairSeparationsUU_D5_'+filename_ending)
# file_name_output_PairSeparationsXX_D5 = os.path.join('Output_batch_files/Correlations_D5/PairSeparationsXX_D5', filename_start + 'PairSeparationsXX_D5_'+filename_ending)
# file_name_output_NearestPairSeparationsMM_Sim = os.path.join('Output_batch_files/Correlations_Sim/NearestPairSeparationsMM_Sim', filename_start + 'NearestPairSeparationsMM_Sim_'+filename_ending)
# file_name_output_NearestPairSeparationsMU_Sim = os.path.join('Output_batch_files/Correlations_Sim/NearestPairSeparationsMU_Sim', filename_start + 'NearestPairSeparationsMU_Sim_'+filename_ending)
# file_name_output_NearestPairSeparationsUU_Sim = os.path.join('Output_batch_files/Correlations_Sim/NearestPairSeparationsUU_Sim', filename_start + 'NearestPairSeparationsUU_Sim_'+filename_ending)
# file_name_output_NearestPairSeparationsXX_Sim = os.path.join('Output_batch_files/Correlations_Sim/NearestPairSeparationsXX_Sim', filename_start + 'NearestPairSeparationsXX_Sim_'+filename_ending)
# file_name_output_NearestPairSeparationsMM_Col0 = os.path.join('Output_batch_files/Correlations_Col0/NearestPairSeparationsMM_Col0', filename_start + 'NearestPairSeparationsMM_Col0_'+filename_ending)
# file_name_output_NearestPairSeparationsMU_Col0 = os.path.join('Output_batch_files/Correlations_Col0/NearestPairSeparationsMU_Col0', filename_start + 'NearestPairSeparationsMU_Col0_'+filename_ending)
# file_name_output_NearestPairSeparationsUU_Col0 = os.path.join('Output_batch_files/Correlations_Col0/NearestPairSeparationsUU_Col0', filename_start + 'NearestPairSeparationsUU_Col0_'+filename_ending)
# file_name_output_NearestPairSeparationsXX_Col0 = os.path.join('Output_batch_files/Correlations_Col0/NearestPairSeparationsXX_Col0', filename_start + 'NearestPairSeparationsXX_Col0_'+filename_ending)
# file_name_output_NearestPairSeparationsMM_Data = os.path.join('Output_batch_files/Correlations_Data/NearestPairSeparationsMM_Data', filename_start + 'NearestPairSeparationsMM_Data_'+filename_ending)
# file_name_output_NearestPairSeparationsMU_Data = os.path.join('Output_batch_files/Correlations_Data/NearestPairSeparationsMU_Data', filename_start + 'NearestPairSeparationsMU_Data_'+filename_ending)
# file_name_output_NearestPairSeparationsUU_Data = os.path.join('Output_batch_files/Correlations_Data/NearestPairSeparationsUU_Data', filename_start + 'NearestPairSeparationsUU_Data_'+filename_ending)
# file_name_output_NearestPairSeparationsXX_Data = os.path.join('Output_batch_files/Correlations_Data/NearestPairSeparationsXX_Data', filename_start + 'NearestPairSeparationsXX_Data_'+filename_ending)
# file_name_output_NearestPairSeparationsMM_Dset2 = os.path.join('Output_batch_files/Correlations_Dset2/NearestPairSeparationsMM_Dset2', filename_start + 'NearestPairSeparationsMM_Dset2_'+filename_ending)
# file_name_output_NearestPairSeparationsMU_Dset2 = os.path.join('Output_batch_files/Correlations_Dset2/NearestPairSeparationsMU_Dset2', filename_start + 'NearestPairSeparationsMU_Dset2_'+filename_ending)
# file_name_output_NearestPairSeparationsUU_Dset2 = os.path.join('Output_batch_files/Correlations_Dset2/NearestPairSeparationsUU_Dset2', filename_start + 'NearestPairSeparationsUU_Dset2_'+filename_ending)
# file_name_output_NearestPairSeparationsXX_Dset2 = os.path.join('Output_batch_files/Correlations_Dset2/NearestPairSeparationsXX_Dset2', filename_start + 'NearestPairSeparationsXX_Dset2_'+filename_ending)
# file_name_output_NearestPairSeparationsMM_D5 = os.path.join('Output_batch_files/Correlations_D5/NearestPairSeparationsMM_D5', filename_start + 'NearestPairSeparationsMM_D5_'+filename_ending)
# file_name_output_NearestPairSeparationsMU_D5 = os.path.join('Output_batch_files/Correlations_D5/NearestPairSeparationsMU_D5', filename_start + 'NearestPairSeparationsMU_D5_'+filename_ending)
# file_name_output_NearestPairSeparationsUU_D5 = os.path.join('Output_batch_files/Correlations_D5/NearestPairSeparationsUU_D5', filename_start + 'NearestPairSeparationsUU_D5_'+filename_ending)
# file_name_output_NearestPairSeparationsXX_D5 = os.path.join('Output_batch_files/Correlations_D5/NearestPairSeparationsXX_D5', filename_start + 'NearestPairSeparationsXX_D5_'+filename_ending)
file_name_output_LocusProperties = os.path.join('Output_batch_files/LocusProperties', filename_start + 'LocusProperties_'+filename_ending)
file_name_output_AllRepsMethFracs = os.path.join('Output_batch_files/AllRepsMethFracs', filename_start + 'AllRepsMethFracs_'+filename_ending)
file_name_output_Skipped_IDs = os.path.join('Output_batch_files/Skipped_IDs', filename_start + 'Skipped_IDs_'+filename_ending)
file_name_output_CodeProgress = os.path.join('Output_batch_files/CodeProgress', filename_start + 'CodeProgress_'+filename_ending)
#file_name_output_AllRepsXFracs = os.path.join('Output_batch_files/AllRepsXFracs', filename_start + 'AllRepsXFracs_'+filename_ending)

# define column headings for LocusProperties file
LocusProperties_headings = ['gene_ID', 'segment_number', 'N_CG', 'L_locus', 'CG_density',
                                        #    'Sim_mu_meth_frac', 'Sim_sigma_meth_frac', 
                                        #    'Sim10_mu_meth_frac', 'Sim10_sigma_meth_frac',
                                        #    'Dset1_mu_meth_frac', 'Dset1_sigma_meth_frac', 
                                        #   'Dset2_mu_meth_frac', 'Dset2_sigma_meth_frac', 
#                                            'D5_mu_meth_frac', 'D5_sigma_meth_frac',
#                                           'DataD4D5_mu_meth_frac', 'DataD4D5_sigma_meth_frac',
                                        #     'Sim_number_0M_frac', 'Sim10_number_0M_frac',
                                        #     'Dset1_number_0M_frac', 
                                        #    'Dset2_number_0M_frac', 
#                                            'D5_number_0M_frac',
#                                            'DataD4D5_number_0M_frac',
                                           'Col0_meth_frac', 
                                            'Col0_X_frac', 
                                        #     'Sim_meth_frac_0', 'Sim_meth_frac_1', 'Sim_meth_frac_2',
                                        #     'Dset1_meth_frac_0', 'Dset1_meth_frac_1', 'Dset1_meth_frac_2', 
                                        #    'Dset2_meth_frac_0', 'Dset2_meth_frac_1', 'Dset2_meth_frac_2', 
#                                            'D5_meth_frac_0', 'D5_meth_frac_1', 'D5_meth_frac_2',
                                        #     'Dset1_mu_meth_diff', 
                                        #    'Dset2_mu_meth_diff', 'D5_mu_meth_diff',
#                                            'DataD4D5_mu_meth_diff', 
                                            'CG_spacing', 'mean_h2az', 'max_h2az'
                                           ]
# define column headings for AllRepsMethFracs files
# Need to update manually to get correct number!! 
# AllRepsMethFracs_headings = ['gene_ID', 'segment_number',
#                             'Sim_meth_frac_0', 'Sim_meth_frac_1', 'Sim_meth_frac_2', 'Sim_meth_frac_3', 'Sim_meth_frac_4', 
#                             'Sim_meth_frac_5', 'Sim_meth_frac_6', 'Sim_meth_frac_7', 'Sim_meth_frac_8', 'Sim_meth_frac_9', 
#                             'Sim_meth_frac_10', 'Sim_meth_frac_11', 'Sim_meth_frac_12', 'Sim_meth_frac_13', 'Sim_meth_frac_14', 
#                             'Sim_meth_frac_15', 'Sim_meth_frac_16', 'Sim_meth_frac_17', 'Sim_meth_frac_18', 'Sim_meth_frac_19', 
#                             'Sim_meth_frac_20', 'Sim_meth_frac_21', 'Sim_meth_frac_22', 'Sim_meth_frac_23', 'Sim_meth_frac_24', 
#                             'Sim_meth_frac_25', 'Sim_meth_frac_26', 'Sim_meth_frac_27', 'Sim_meth_frac_28', 'Sim_meth_frac_29', 
#                             'Dset1_meth_frac_0','Dset1_meth_frac_1','Dset1_meth_frac_2','Dset1_meth_frac_3','Dset1_meth_frac_4',
#                             'Dset1_meth_frac_5','Dset1_meth_frac_6','Dset1_meth_frac_7','Dset1_meth_frac_8','Dset1_meth_frac_9',
#                             'Dset1_meth_frac_10','Dset1_meth_frac_11','Dset1_meth_frac_12','Dset1_meth_frac_13','Dset1_meth_frac_14',
#                             'Dset1_meth_frac_15','Dset1_meth_frac_16','Dset1_meth_frac_17','Dset1_meth_frac_18','Dset1_meth_frac_19',
#                             'Dset1_meth_frac_20','Dset1_meth_frac_21','Dset1_meth_frac_22','Dset1_meth_frac_23','Dset1_meth_frac_24',
#                             'Dset1_meth_frac_25','Dset1_meth_frac_26','Dset1_meth_frac_27','Dset1_meth_frac_28','Dset1_meth_frac_29']
# AllRepsXFracs_headings = ['gene_ID', 'segment_number',
#                             'Dset1_meth_frac_0','Dset1_meth_frac_1','Dset1_meth_frac_2','Dset1_meth_frac_3','Dset1_meth_frac_4',
#                             'Dset1_meth_frac_5','Dset1_meth_frac_6','Dset1_meth_frac_7','Dset1_meth_frac_8','Dset1_meth_frac_9',
#                             'Dset1_meth_frac_10','Dset1_meth_frac_11','Dset1_meth_frac_12','Dset1_meth_frac_13','Dset1_meth_frac_14',
#                             'Dset1_meth_frac_15','Dset1_meth_frac_16','Dset1_meth_frac_17','Dset1_meth_frac_18','Dset1_meth_frac_19',
#                             'Dset1_meth_frac_20','Dset1_meth_frac_21','Dset1_meth_frac_22','Dset1_meth_frac_23','Dset1_meth_frac_24',
#                             'Dset1_meth_frac_25','Dset1_meth_frac_26','Dset1_meth_frac_27','Dset1_meth_frac_28','Dset1_meth_frac_29']
# define column headings for Skipped_IDs files
Skipped_IDs_headings = ['gene_ID', 'N_CG', 'CG_desnity']
# define column headings for TEST files
TEST_headings = ['segment_no', 'gene_ID']
# # define column headings for SimStates
# SimState_headings = ['gene_ID','Chromosome','Start','End','Status']
# define column headings for CodeProgress
CodeProgress_headings = ['batch', 'N_gene', 'time_mins']

# if i_batch = 0 create output files and exit code using sys.exit()
if i_batch == 0:
    # write out column heading
    with open(file_name_output_TEST,'w',newline='') as output_file:
        line_writer = csv.writer(output_file,delimiter='\t')
        line_writer.writerow(TEST_headings)
    with open(file_name_output_LocusProperties,'w',newline='') as output_file:
        line_writer = csv.writer(output_file,delimiter='\t')
        line_writer.writerow(LocusProperties_headings)
    # with open(file_name_output_AllRepsMethFracs,'w',newline='') as output_file:
    #     line_writer = csv.writer(output_file,delimiter='\t')
    #     line_writer.writerow(AllRepsMethFracs_headings)
    # with open(file_name_output_AllRepsXFracs,'w',newline='') as output_file:
    #     line_writer = csv.writer(output_file,delimiter='\t')
    #     line_writer.writerow(AllRepsXFracs_headings)
    with open(file_name_output_Skipped_IDs,'w',newline='') as output_file:
        line_writer = csv.writer(output_file,delimiter='\t')
        line_writer.writerow(Skipped_IDs_headings)
    with open(file_name_output_CodeProgress,'w',newline='') as output_file:
        line_writer = csv.writer(output_file,delimiter='\t')
        line_writer.writerow(CodeProgress_headings)
    # # loop over output files for Simulated states
    # for i_ in range(N_reps):
    #    with open(file_name_output_SimulatedStates_+str(i_)+'_'+filename_ending,'w',newline='') as output_file:
    #        line_writer = csv.writer(output_file,delimiter='\t')
    #        line_writer.writerow(SimState_headings)
 #   create empty files for rest....
    # open(file_name_output_PairSeparationsMM_Sim,'w').close
    # open(file_name_output_PairSeparationsMU_Sim,'w').close
    # open(file_name_output_PairSeparationsUU_Sim,'w').close
    # open(file_name_output_PairSeparationsXX_Sim,'w').close
    # open(file_name_output_PairSeparationsMM_Col0,'w').close
    # open(file_name_output_PairSeparationsMU_Col0,'w').close
    # open(file_name_output_PairSeparationsUU_Col0,'w').close
    # open(file_name_output_PairSeparationsXX_Col0,'w').close
    # open(file_name_output_PairSeparationsMM_Data,'w').close
    # open(file_name_output_PairSeparationsMU_Data,'w').close
    # open(file_name_output_PairSeparationsUU_Data,'w').close
    # open(file_name_output_PairSeparationsXX_Data,'w').close
    # open(file_name_output_PairSeparationsMM_Dset2,'w').close
    # open(file_name_output_PairSeparationsMU_Dset2,'w').close
    # open(file_name_output_PairSeparationsUU_Dset2,'w').close
    # open(file_name_output_PairSeparationsXX_Dset2,'w').close
    # open(file_name_output_PairSeparationsMM_D5,'w').close
    # open(file_name_output_PairSeparationsMU_D5,'w').close
    # open(file_name_output_PairSeparationsUU_D5,'w').close
    # open(file_name_output_PairSeparationsXX_D5,'w').close
    # open(file_name_output_NearestPairSeparationsMM_Sim,'w').close
    # open(file_name_output_NearestPairSeparationsMU_Sim,'w').close
    # open(file_name_output_NearestPairSeparationsUU_Sim,'w').close
    # open(file_name_output_NearestPairSeparationsXX_Sim,'w').close
    # open(file_name_output_NearestPairSeparationsMM_Col0,'w').close
    # open(file_name_output_NearestPairSeparationsMU_Col0,'w').close
    # open(file_name_output_NearestPairSeparationsUU_Col0,'w').close
    # open(file_name_output_NearestPairSeparationsXX_Col0,'w').close
    # open(file_name_output_NearestPairSeparationsMM_Data,'w').close
    # open(file_name_output_NearestPairSeparationsMU_Data,'w').close
    # open(file_name_output_NearestPairSeparationsUU_Data,'w').close
    # open(file_name_output_NearestPairSeparationsXX_Data,'w').close
    # open(file_name_output_NearestPairSeparationsMM_Dset2,'w').close
    # open(file_name_output_NearestPairSeparationsMU_Dset2,'w').close
    # open(file_name_output_NearestPairSeparationsUU_Dset2,'w').close
    # open(file_name_output_NearestPairSeparationsXX_Dset2,'w').close
    # open(file_name_output_NearestPairSeparationsMM_D5,'w').close
    # open(file_name_output_NearestPairSeparationsMU_D5,'w').close
    # open(file_name_output_NearestPairSeparationsUU_D5,'w').close
    # open(file_name_output_NearestPairSeparationsXX_D5,'w').close
    sys.exit()
# otherwise create empty verions of ALL files
else: 
    open(file_name_output_TEST,'w').close
    open(file_name_output_LocusProperties,'w').close
    # open(file_name_output_AllRepsMethFracs,'w').close
    # open(file_name_output_AllRepsXFracs,'w').close
    open(file_name_output_Skipped_IDs,'w').close
    open(file_name_output_CodeProgress,'w').close
    # loop over output files for Simulated states
    # for i_ in range(N_reps):
    #    open(file_name_output_SimulatedStates_+str(i_)+'_'+filename_ending,'w').close
    # open(file_name_output_PairSeparationsMM_Sim,'w').close
    # open(file_name_output_PairSeparationsMU_Sim,'w').close
    # open(file_name_output_PairSeparationsUU_Sim,'w').close
    # open(file_name_output_PairSeparationsXX_Sim,'w').close
    # open(file_name_output_PairSeparationsMM_Col0,'w').close
    # open(file_name_output_PairSeparationsMU_Col0,'w').close
    # open(file_name_output_PairSeparationsUU_Col0,'w').close
    # open(file_name_output_PairSeparationsXX_Col0,'w').close
    # open(file_name_output_PairSeparationsMM_Data,'w').close
    # open(file_name_output_PairSeparationsMU_Data,'w').close
    # open(file_name_output_PairSeparationsUU_Data,'w').close
    # open(file_name_output_PairSeparationsXX_Data,'w').close
    # open(file_name_output_PairSeparationsMM_D4,'w').close
    # open(file_name_output_PairSeparationsMU_D4,'w').close
    # open(file_name_output_PairSeparationsUU_D4,'w').close
    # open(file_name_output_PairSeparationsXX_D4,'w').close
    # open(file_name_output_PairSeparationsMM_D5,'w').close
    # open(file_name_output_PairSeparationsMU_D5,'w').close
    # open(file_name_output_PairSeparationsUU_D5,'w').close
    # open(file_name_output_PairSeparationsXX_D5,'w').close
    # open(file_name_output_NearestPairSeparationsMM_Sim,'w').close
    # open(file_name_output_NearestPairSeparationsMU_Sim,'w').close
    # open(file_name_output_NearestPairSeparationsUU_Sim,'w').close
    # open(file_name_output_NearestPairSeparationsXX_Sim,'w').close
    # open(file_name_output_NearestPairSeparationsMM_Col0,'w').close
    # open(file_name_output_NearestPairSeparationsMU_Col0,'w').close
    # open(file_name_output_NearestPairSeparationsUU_Col0,'w').close
    # open(file_name_output_NearestPairSeparationsXX_Col0,'w').close
    # open(file_name_output_NearestPairSeparationsMM_Data,'w').close
    # open(file_name_output_NearestPairSeparationsMU_Data,'w').close
    # open(file_name_output_NearestPairSeparationsUU_Data,'w').close
    # open(file_name_output_NearestPairSeparationsXX_Data,'w').close
    # open(file_name_output_NearestPairSeparationsMM_D4,'w').close
    # open(file_name_output_NearestPairSeparationsMU_D4,'w').close
    # open(file_name_output_NearestPairSeparationsUU_D4,'w').close
    # open(file_name_output_NearestPairSeparationsXX_D4,'w').close
    # open(file_name_output_NearestPairSeparationsMM_D5,'w').close
    # open(file_name_output_NearestPairSeparationsMU_D5,'w').close
    # open(file_name_output_NearestPairSeparationsUU_D5,'w').close
    # open(file_name_output_NearestPairSeparationsXX_D5,'w').close



# work out which IDs to load in
segment_no_start = (i_batch-1)*batch_length
segment_no_end = i_batch*batch_length
segment_no_list = np.arange(segment_no_start+1,segment_no_end+1)

# # load in list of required data files
# required_accessions_filenames_df = pd.read_csv(accessions_codes_file,sep='\t',header=None)
# required_accessions_filenames_list = required_accessions_filenames_df[0].tolist()

# load in required segments to model from annotaion file: makes dataframe: AllLoci_df
# CG site at 'start' to be included: x >= start
# 'end' is past last CG site: x < end

# create lazy reader object
temp_input_df = dd.read_csv(annotaion_file,sep='\t',
                            usecols=['segment_no', 'gene_ID','start','end', 'mean_h2az', 'max_h2az'],
                            dtype={'segment_no':'int64', 'gene_ID':'str', 'start':'int64', 'end':'int64', 
                            'mean_h2az':'float', 'max_h2az':'float'})
# filter specific IDs for this core:
temp_input_df = temp_input_df[temp_input_df['segment_no'].isin(segment_no_list)  ]

# apply annofilt_1: list of non-TE contaminated genes UNLESS filename == 'NONE'
if file_in_anno_filt_1 != 'NONE': 
    with open(filt_1_file ) as file_in:
        anno_filt_1_IDs_array = np.array(file_in.read().splitlines())
    # do filtering
    temp_input_df = temp_input_df[temp_input_df['gene_ID'].isin(anno_filt_1_IDs_array)  ]

# apply annofilt_2: list of non-TE contaminated genes UNLESS filename == 'NONE'
if file_in_anno_filt_2 != 'NONE': 
    with open(filt_2_file ) as file_in:
        anno_filt_2_IDs_array = np.array(file_in.read().splitlines())
    # do filtering
    temp_input_df = temp_input_df[temp_input_df['gene_ID'].isin(anno_filt_2_IDs_array)  ]

# apply exclude_filt_1: Exclude filter! # need to read this one in as dataframe
if file_in_exclude_filt_1 != 'NONE': 
    exclude_filt_1_df = pd.read_csv(exclude_filt_1_file, sep='\t', usecols=['gene_ID'])
    exclude_filt_1_IDs_array = np.array( exclude_filt_1_df['gene_ID'].tolist() )
    # do filtering
    temp_input_df = temp_input_df[temp_input_df['gene_ID'].isin(exclude_filt_1_IDs_array) == False ]

# set index to gene_ID
temp_input_df = temp_input_df.set_index('gene_ID')
# apply filtering logic and convert to pandas dataframe
temp_input_df = temp_input_df.compute()
AllLoci_df = temp_input_df.copy()
# make list of IDs
genes_with_AllLoci_ID_list = AllLoci_df.index.drop_duplicates().tolist()

# write out segment_no and gene_ID to test file 
AllLoci_df.to_csv(file_name_output_TEST, sep='\t', columns=['segment_no'], header = False)

# load in reqired IDs from Lizzie's consensus state: makes dataframe: Col0_CG_site_data_df
# create lazy reader object
temp_input_df = dd.read_csv(file_in_Combined_Data_per_genes_meth_statuses,sep='\t',
                            usecols=['gene_ID','Start','Status'],
                            dtype={'gene_ID': 'str', 'Start': 'int64', 'Status': 'str'})
# define filtering logic
temp_input_df = temp_input_df[temp_input_df['gene_ID'].isin(genes_with_AllLoci_ID_list)  ] 
temp_input_df = temp_input_df.set_index('gene_ID')
# apply filtering logic and convert to pandas dataframe
temp_input_df = temp_input_df.compute()
Col0_CG_site_data_df = temp_input_df.copy()

# set input genes list to genes_with_AllLoci_ID_list
input_gene_ID_list = genes_with_AllLoci_ID_list

# start loading in 1001 genomes data files

# load in selected m1001 States
# # required_accessions_filenames_list
# ExptState_Dset1_df_list = []
# # build up path
# temp_path = os.path.join(path_string, 'm1001_states_NS')

# #print(temp_path)
# #print()
# for filename in os.listdir(temp_path):
#     if filename in required_accessions_filenames_list: 
#         #print(filename)
#         # create lazy reader object
#         temp_input_df = dd.read_csv(os.path.join(temp_path,filename),sep='\t',
#                                     usecols=['gene_ID','Start','Status'],
#                                 dtype={'gene_ID': 'str', 'Start': 'int64', 'Status': 'str'})
#         # define filtering logic
#         temp_input_df = temp_input_df[temp_input_df['gene_ID'].isin(input_gene_ID_list)  ] 
#         temp_input_df = temp_input_df.set_index('gene_ID')
#         # apply filtering logic and convert to pandas dataframe
#         temp_input_df = temp_input_df.compute()
#         ExptState_Dset1_df_list.append(temp_input_df.copy())
# N_ExptState_Dset1 = len(ExptState_Dset1_df_list)
# #print(N_ExptState_Dset1)
# #print(ExptState_Dset1_df_list[0].shape,ExptState_Dset1_df_list[1].shape)

# # load in Decile4 m1001 States
# Decile_number = 4
# ExptState_Decile4_df_list = []
# # build up path
# temp_path = os.path.join(path_string, 'm1001_states','Decile'+str(Decile_number)+'_3')

# #print(temp_path)
# #print()
# for filename in os.listdir(temp_path):
#     #print(filename)
#     # create lazy reader object
#     temp_input_df = dd.read_csv(os.path.join(temp_path,filename),sep='\t',
#                                 usecols=['gene_ID','Start','Status'],
#                                dtype={'gene_ID': 'str', 'Start': 'int64', 'Status': 'str'})
#     # define filtering logic
#     temp_input_df = temp_input_df[temp_input_df['gene_ID'].isin(input_gene_ID_list)  ] 
#     temp_input_df = temp_input_df.set_index('gene_ID')
#     # apply filtering logic and convert to pandas dataframe
#     temp_input_df = temp_input_df.compute()
#     ExptState_Decile4_df_list.append(temp_input_df.copy())
# N_ExptState_Decile4 = len(ExptState_Decile4_df_list)
# #print(N_ExptState_Decile4)
# #print(ExptState_Decile4_df_list[0].shape,ExptState_Decile4_df_list[1].shape)

# # load in Decile5 m1001 States
# Decile_number = 5
# ExptState_Decile5_df_list = []
# # build up path
# temp_path = os.path.join(path_string, 'm1001_states','Decile'+str(Decile_number)+'_3')

# #print(temp_path)
# #print()
# for filename in os.listdir(temp_path):
#     #print(filename)
#     # create lazy reader object
#     temp_input_df = dd.read_csv(os.path.join(temp_path,filename),sep='\t',
#                                 usecols=['gene_ID','Start','Status'],
#                                dtype={'gene_ID': 'str', 'Start': 'int64', 'Status': 'str'})
#     # define filtering logic
#     temp_input_df = temp_input_df[temp_input_df['gene_ID'].isin(input_gene_ID_list)  ] 
#     temp_input_df = temp_input_df.set_index('gene_ID')
#     # apply filtering logic and convert to pandas dataframe
#     temp_input_df = temp_input_df.compute()
#     ExptState_Decile5_df_list.append(temp_input_df.copy())
# N_ExptState_Decile5 = len(ExptState_Decile5_df_list)
# #print(N_ExptState_Decile5)
# #print(ExptState_Decile5_df_list[0].shape,ExptState_Decile5_df_list[1].shape)

# finish reading in 1001 methylome files


# Eqbrm initial state
def create_initial_state_eqbrm_Col0like_regions(ID_, initial_state_choice_, P_choice_,
                                           current_segment_start_, current_segment_end_):
    # possible values for initial_state_choice_
    # "100U"
    # "100M"
    # 'Expt'
    # '50M'

    # possilbe values for P_choice_
    # # 'UseU' 'Excl'

    # filter out CG sites for current ID_
    temp_df_ = current_CG_site_data_df[current_CG_site_data_df.index == ID_].copy()
    # filter out CG sites that lie within the trimmed locus
    temp_df_ = temp_df_.loc[ (temp_df_["Start"] >= current_segment_start_) & 
                            (temp_df_["Start"] < current_segment_end_) ]
    
    temp_df_.sort_values("Start", axis=0, ascending=True, inplace=True)
    
    CG_positions_ = temp_df_.Start.tolist()
    N_CG_ = len(CG_positions_)
    state_text_ = temp_df_.Status.tolist()
    state_expt_ = []

    for i_ in range(N_CG_):
        if state_text_[i_] == 'U':
            state_expt_.append(0)
        elif state_text_[i_] == 'M':
            state_expt_.append(2)
        elif ((state_text_[i_] == 'P') or (state_text_[i_] == 'I') ):
            if P_choice_ == "UseU":
                state_expt_.append(0)
            elif P_choice_ == "Excl":
                state_expt_.append(1)
        # elif (state_text_[i_] == 'I'):
        #     state_expt_.append(1)
        # equivalent to 'U'
        elif state_text_[i_] == '4':
            state_expt_.append(0)
        # equivalent to 'M'
        elif state_text_[i_] == '2':
            state_expt_.append(2)
        # eqivalent to 'P' and equivalent to 'I'
        elif ( (state_text_[i_] == '3') or (state_text_[i_] == '1') ):
            if P_choice_ == "UseU":
                state_expt_.append(0)
            elif P_choice_ == "Excl":
                state_expt_.append(1)
        # # equivalent to 'I'
        # elif (state_text_[i_] == '1'):
        #     state_expt_.append(1)
        elif (state_text_[i_] == 'X'):
            state_expt_.append(1)
        elif ( np.isnan(state_text_[i_]) ):
            state_expt_.append(1)
        else:
            state_expt_.append(1)
            print('fasd', state_text_[i_], type(state_text_[i_]), ID_, i_expts)
            print(N_CG_,state_text_)
            print(CG_positions_)
    
    if initial_state_choice_ == '100U':
        state_sim_input_ = np.full(N_CG_,0)
    elif initial_state_choice_ == '100M':
        state_sim_input_ = np.full(N_CG_,2)
    elif initial_state_choice_ == 'Expt':
        state_sim_input_ = state_expt_.copy()
        for i_ in range(N_CG_):
            # check for unknown sites an change to 'U' so they are also included in simulation
            if state_sim_input_[i_] == 1:
                state_sim_input_[i_] = 0

    return N_CG_, CG_positions_, state_expt_, state_sim_input_




def initial_propensities(N_CG_, state_in_, CG_positions_, spont_gain_params_, 
                         coop_gain_params_, coop_maint_params_):
    
    # unpack background gain parameters
    delta_ = spont_gain_params_[0]
    
    # unpack coop_gain parameters
    u_M_plus_PL_ = coop_gain_params_[0]
    lambda_coop_ = coop_gain_params_[1]
    r_div_ = coop_gain_params_[2]
    r_plat_ = coop_gain_params_[3]
    
    u_M_plus_PL_LR1_ = coop_gain_params_[4]
    lambda_coopIn_LR1_ = coop_gain_params_[5]
    lambda_coopOut_LR1_ = coop_gain_params_[6]
    r_div_LR1_ = coop_gain_params_[7]
    r_plat_LR1_ = coop_gain_params_[8]
    r_LR1_ = coop_gain_params_[9]
    
    # unpack coop_maint parameters
    epsilon_ = coop_maint_params_[0]
    gamma0_ = coop_maint_params_[1]
    lambda_gamma_ = coop_maint_params_[2]
    r_div_gamma_ = coop_maint_params_[3]
    r_gamma_ = coop_maint_params_[4]
    
    # create empty coop_gain propensity matrix
    coop_gains_propensities_array_ = np.full((N_CG_,N_CG_),0.)
    spont_gains_propensities_vector_ = np.full(N_CG_,0.)
    coop_maint_propensities_array_ = np.full((N_CG_,N_CG_),0.)
    # index array using (i_target_, i_passive_)
    # then sum over i_passive_ to produce vector of propensities
    
    for i_target_ in range(N_CG_):
    # loop over target sites. 
        
        # identify U sites for gains
        if state_in_[i_target_] == 0:
            
            # add in spontaneous gain propensity 
            # cal. background gain effective propensity
            # a_0_eff_ = (a_0_)(1/2) = 2delta_ /2      # note 2delta_ = U->H probability
            
            # add a_0_eff_ to existing propensity from coop. gains. 
            spont_gains_propensities_vector_[i_target_] = delta_
            
            #loop over i_passive_ (for coop. gains)
            for i_passive_ in range(i_target_):

            # 1. i_passive_ < i_target: r_LR_temp = abs(CG_positions_[i_passive_] - CG_positions_[i_target_] + r_LR1_)
                
                # identify it M site (to promote cooperative gain)
                if state_in_[i_passive_] == 2:
                    
                    r_temp = abs(CG_positions_[i_passive_] - CG_positions_[i_target_])
                    r_LR_temp = abs(CG_positions_[i_passive_] - CG_positions_[i_target_] + r_LR1_)
                    
                    # cal. short-range-coop. effective propensity
                    # a_coop_sr_eff_ = (a_coop_sr_)(1/2)
                    # a_coop_sr_eff_ = (2*u_M_plus_PL_)(1/2) = u_M_plus_PL_    for r_temp < r_plat_
                    # a_coop_sr_eff_ = (2.*u_M_plus_PL_*abs(r_plat_-r_div_)**lambda_coop_
                    #                   *abs(r_temp-r_div_)**(-lambda_coop_))*(1/2)    for r_temp >= r_plat_
                    #                = u_M_plus_PL_*abs(r_plat_-r_div_)**lambda_coop_
                    #                   *abs(r_temp-r_div_)**(-lambda_coop_)     for r_temp >= r_plat_
                    if r_temp < r_plat_:
                        a_coop_sr_eff_ = u_M_plus_PL_
                    elif r_temp >= r_plat_:
                        a_coop_sr_eff_ = u_M_plus_PL_*abs(r_plat_-r_div_)**lambda_coop_ \
                                            *abs(r_temp-r_div_)**(-lambda_coop_)
                            
                    # cal. long-range-coop. effective propensity
                    # a_coop_LR_eff_ = a_coop_LR_ / 2
                    # a_coop_LR_ = 2.*u_M_plus_PL_LR1_    for r_LR_temp >= r_plat_
                    # a_coop_LR_ = 2.*u_M_plus_PL_LR1_*abs(r_plat_LR1_-r_div_LR1_)**lambda_coopIn_LR1_     
                    #               *abs(r_LR_temp-r_div_LR1_)**(-lambda_coopIn_LR1_)*(1/2)     for r_LR_temp >= r_plat_LR1_ AND r_temp < r_LR1_
                    # a_coop_LR_ = 2.*u_M_plus_PL_LR1_*abs(r_plat_LR1_-r_div_LR1_)**lambda_coopOut_LR1_     
                    #               *abs(r_LR_temp-r_div_LR1_)**(-lambda_coopOut_LR1_)*(1/2)     for r_LR_temp >= r_plat_LR1_ AND r_temp > r_LR1_

                    if r_LR_temp < r_plat_LR1_:
                        a_coop_LR_eff_ = u_M_plus_PL_LR1_
                    elif r_LR_temp >= r_plat_LR1_:
                        # Use separate powerlaws depending if a) jumping gap 'IN towards nucleosome'
                        if r_temp < r_LR1_:
                            a_coop_LR_eff_ = u_M_plus_PL_LR1_*abs(r_plat_LR1_-r_div_LR1_)**lambda_coopIn_LR1_ \
                                               *abs(r_LR_temp-r_div_LR1_)**(-lambda_coopIn_LR1_)
                        # Use separate powerlaws depending if b) jumping gap 'OUT away from nucleosome'
                        elif r_temp > r_LR1_:
                            a_coop_LR_eff_ = u_M_plus_PL_LR1_*abs(r_plat_LR1_-r_div_LR1_)**lambda_coopOut_LR1_ \
                                               *abs(r_LR_temp-r_div_LR1_)**(-lambda_coopOut_LR1_)

                        
                    # update coop_gain propensity matrix-element 
                    coop_gains_propensities_array_[i_target_,i_passive_] = a_coop_sr_eff_ + a_coop_LR_eff_
                      
            for i_passive_ in range(i_target_+1, N_CG_):
            #loop over i_passive_
            # 1. i_passive_ > i_target: r_LR_temp = abs(CG_positions_[i_passive_] - CG_positions_[i_target_] - r_LR1_)
                
                # identify it M site (to promote cooperative gain)
                if state_in_[i_passive_] == 2:
                    
                    r_temp = abs(CG_positions_[i_passive_] - CG_positions_[i_target_])
                    r_LR_temp = abs(CG_positions_[i_passive_] - CG_positions_[i_target_] - r_LR1_)
                    
                    # cal. short-range-coop. effective propensity
                    # a_coop_sr_eff_ = (a_coop_sr_)(1/2)
                    # a_coop_sr_eff_ = (2*u_M_plus_PL_)(1/2) = u_M_plus_PL_    for r_temp < r_plat_
                    # a_coop_sr_eff_ = (2.*u_M_plus_PL_*abs(r_plat_-r_div_)**lambda_coop_
                    #                   *abs(r_temp-r_div_)**(-lambda_coop_))*(1/2)    for r_temp >= r_plat_
                    #                = u_M_plus_PL_*abs(r_plat_-r_div_)**lambda_coop_
                    #                   *abs(r_temp-r_div_)**(-lambda_coop_)     for r_temp >= r_plat_
                    if r_temp < r_plat_:
                        a_coop_sr_eff_ = u_M_plus_PL_
                    elif r_temp >= r_plat_:
                        a_coop_sr_eff_ = u_M_plus_PL_*abs(r_plat_-r_div_)**lambda_coop_ \
                                            *abs(r_temp-r_div_)**(-lambda_coop_)
                            
                    # cal. long-range-coop. effective propensity
                    # a_coop_LR_eff_ = a_coop_LR_ / 2
                    # a_coop_LR_ = 2.*u_M_plus_PL_LR1_    for r_LR_temp >= r_plat_
                    # a_coop_LR_ = 2.*u_M_plus_PL_LR1_*abs(r_plat_LR1_-r_div_LR1_)**lambda_coopIn_LR1_     
                    #               *abs(r_LR_temp-r_div_LR1_)**(-lambda_coopIn_LR1_)*(1/2)     for r_LR_temp >= r_plat_LR1_ AND r_temp < r_LR1_
                    # a_coop_LR_ = 2.*u_M_plus_PL_LR1_*abs(r_plat_LR1_-r_div_LR1_)**lambda_coopOut_LR1_     
                    #               *abs(r_LR_temp-r_div_LR1_)**(-lambda_coopOut_LR1_)*(1/2)     for r_LR_temp >= r_plat_LR1_ AND r_temp > r_LR1_

                    if r_LR_temp < r_plat_LR1_:
                        a_coop_LR_eff_ = u_M_plus_PL_LR1_
                    elif r_LR_temp >= r_plat_LR1_:
                        # Use separate powerlaws depending if a) jumping gap 'IN towards nucleosome'
                        if r_temp < r_LR1_:
                            a_coop_LR_eff_ = u_M_plus_PL_LR1_*abs(r_plat_LR1_-r_div_LR1_)**lambda_coopIn_LR1_ \
                                               *abs(r_LR_temp-r_div_LR1_)**(-lambda_coopIn_LR1_)
                        # Use separate powerlaws depending if b) jumping gap 'OUT away from nucleosome'
                        elif r_temp > r_LR1_:
                            a_coop_LR_eff_ = u_M_plus_PL_LR1_*abs(r_plat_LR1_-r_div_LR1_)**lambda_coopOut_LR1_ \
                                               *abs(r_LR_temp-r_div_LR1_)**(-lambda_coopOut_LR1_)
                            
                        
                    # update coop_gain propensity matrix-element 
                    coop_gains_propensities_array_[i_target_,i_passive_] = a_coop_sr_eff_ + a_coop_LR_eff_
                            
        # identify M sites (for coop. maint.)
        elif state_in_[i_target_] == 2:

            #loop over i_passive_
            for i_passive_ in range(i_target_):

                # identify it M site (to promote coop. maint.)
                if state_in_[i_passive_] == 2:
                    
                    r_temp = abs(CG_positions_[i_passive_] - CG_positions_[i_target_])

                    # cal. coop. maint. failure propensity
                    # a_MaintFail_eff_ = Prob(M->H)Prob(H->H)Prob(H->U)   (all probs at replication)
                    #                  = a_MaintFail/2    (if here a_MaintFail == Prob(H->H) at replication)
                    # a_MaintFail = epsilon_*PI_i{1. - gamma[r_i]}
                    # gamma[r_i] = gamma0_     for r_temp < r_gamma_
                    # gamma[r_i] = gamma0_*abs(r_gamma_-r_div_gamma_)**(lambda_gamma_)*
                    #              abs(r_temp - r_div_gamma_)**(-lambda_gamma_)     for r_temp >= r_gamma_
                    
                    if r_temp < r_gamma_:
                        a_MaintFail_eff_ = 1.0 - gamma0_
                    elif r_temp >= r_gamma_:
                        a_MaintFail_eff_ = 1.0 - gamma0_*abs(r_gamma_-r_div_gamma_)**lambda_gamma_ \
                                            *abs(r_temp-r_div_gamma_)**(-lambda_gamma_)
                        
                    # update coop_maint propensity matrix-element 
                    coop_maint_propensities_array_[i_target_,i_passive_] = a_MaintFail_eff_
                    
                # otherwise, no coop. maint. so set to '1'
                else:
                    coop_maint_propensities_array_[i_target_,i_passive_] = 1.

            # set self-interaction element to 1.
            coop_maint_propensities_array_[i_target_,i_target_] = 1.
                    
            for i_passive_ in range(i_target_+1, N_CG_):
            #loop over i_passive_
                
                # identify if M site (to promote coop. maint. )
                if state_in_[i_passive_] == 2:
                    
                    r_temp = abs(CG_positions_[i_passive_] - CG_positions_[i_target_])
                    
                    # cal. coop. maint. failure propensity
                    # a_MaintFail_eff_ = Prob(M->H)Prob(H->H)Prob(H->U)   (all probs at replication)
                    #                  = a_MaintFail/2    (if here a_MaintFail == Prob(H->H) at replication)
                    # a_MaintFail = epsilon_*PI_i{1. - gamma[r_i]}
                    # gamma[r_i] = gamma0_     for r_temp < r_gamma_
                    # gamma[r_i] = gamma0_*abs(r_gamma_-r_div_gamma_)**(lambda_gamma_)*
                    #              abs(r_temp - r_div_gamma_)**(-lambda_gamma_)     for r_temp >= r_gamma_
                    
                    if r_temp < r_gamma_:
                        a_MaintFail_eff_ = 1.0 - gamma0_
                    elif r_temp >= r_gamma_:
                        a_MaintFail_eff_ = 1.0 - gamma0_*abs(r_gamma_-r_div_gamma_)**lambda_gamma_ \
                                            *abs(r_temp-r_div_gamma_)**(-lambda_gamma_)
                            
                    # update coop_maint propensity matrix-element 
                    coop_maint_propensities_array_[i_target_,i_passive_] = a_MaintFail_eff_
                    
                # otherwise, no coop. maint. so set to '1'
                else:
                    coop_maint_propensities_array_[i_target_,i_passive_] = 1.
                    
                    # end of all loops
                    
    # calculate composite propensities vector (one propensity for each site in gene)
    # sum over i_passive_ (ie. along rows) to create coop_gains_propensities_vector_
    coop_gains_propensities_vector_ = np.sum(coop_gains_propensities_array_, axis=1)
    # find product over i_passive_ (i.e. along rows) to create coop_maint_propensities_vector
    coop_maint_propensities_vector_ = (epsilon_/2.)*np.prod( coop_maint_propensities_array_,axis=1)

    total_propensities_vector_ = coop_gains_propensities_vector_ + \
                                    spont_gains_propensities_vector_ + coop_maint_propensities_vector_

    
    return coop_gains_propensities_array_, spont_gains_propensities_vector_, \
            coop_maint_propensities_array_, total_propensities_vector_




def update_propensities(N_CG_, state_in_, CG_positions_, spont_gain_params_, coop_gain_params_, 
                        coop_maint_params_, coop_gains_propensities_array_, 
                        spont_gains_propensities_vector_, coop_maint_propensities_array_, change_index_):
    
    # unpack background gain parameters
    delta_ = spont_gain_params_[0]
    
    # unpack coop_gain parameters
    u_M_plus_PL_ = coop_gain_params_[0]
    lambda_coop_ = coop_gain_params_[1]
    r_div_ = coop_gain_params_[2]
    r_plat_ = coop_gain_params_[3]
    
    u_M_plus_PL_LR1_ = coop_gain_params_[4]
    lambda_coopIn_LR1_ = coop_gain_params_[5]
    lambda_coopOut_LR1_ = coop_gain_params_[6]
    r_div_LR1_ = coop_gain_params_[7]
    r_plat_LR1_ = coop_gain_params_[8]
    r_LR1_ = coop_gain_params_[9]
    
    # unpack coop_maint parameters
    epsilon_ = coop_maint_params_[0]
    gamma0_ = coop_maint_params_[1]
    lambda_gamma_ = coop_maint_params_[2]
    r_div_gamma_ = coop_maint_params_[3]
    r_gamma_ = coop_maint_params_[4]
    
    
    
    # update propensities
    # need to update one rown and one column
    
    # a) use i_target_ = change_index_
    i_target_ = change_index_
    
    
    # if i_target is M site:
    if state_in_[i_target_] == 2:
        #set row of coop. gain. propensities to zero:
        coop_gains_propensities_array_[i_target_,:] = 0
        
        # calc. coop. maint. propensities for row:
        #loop over i_passive_
        for i_passive_ in range(i_target_):
            # identify if M site (to promote coop. maint.)
            if state_in_[i_passive_] == 2:

                r_temp = abs(CG_positions_[i_passive_] - CG_positions_[i_target_])

                # cal. coop. maint. failure propensity
                # a_MaintFail_eff_ = Prob(M->H)Prob(H->H)Prob(H->U)   (all probs at replication)
                #                  = a_MaintFail/2    (if here a_MaintFail == Prob(H->H) at replication)
                # a_MaintFail = epsilon_*PI_i{1. - gamma[r_i]}
                # gamma[r_i] = gamma0_     for r_temp < r_gamma_
                # gamma[r_i] = gamma0_*abs(r_gamma_-r_div_gamma_)**(lambda_gamma_)*
                #              abs(r_temp - r_div_gamma_)**(-lambda_gamma_)     for r_temp >= r_gamma_

                if r_temp < r_gamma_:
                    a_MaintFail_eff_ = 1.0 - gamma0_
                elif r_temp >= r_gamma_:
                    a_MaintFail_eff_ = 1.0 - gamma0_*abs(r_gamma_-r_div_gamma_)**lambda_gamma_ \
                                        *abs(r_temp-r_div_gamma_)**(-lambda_gamma_)

                # update coop_maint propensity matrix-element 
                coop_maint_propensities_array_[i_target_,i_passive_] = a_MaintFail_eff_
                
            # otherwise, no coop. maint. so set to '1'
            else:
                coop_maint_propensities_array_[i_target_,i_passive_] = 1.
                
        # set self-interaction element to 1.
        coop_maint_propensities_array_[i_target_,i_target_] = 1.

        for i_passive_ in range(i_target_+1, N_CG_):
        #loop over i_passive_

            # identify if M site (to promote coop. maint. )
            if state_in_[i_passive_] == 2:

                r_temp = abs(CG_positions_[i_passive_] - CG_positions_[i_target_])

                # cal. coop. maint. failure propensity
                # a_MaintFail_eff_ = Prob(M->H)Prob(H->H)Prob(H->U)   (all probs at replication)
                #                  = a_MaintFail/2    (if here a_MaintFail == Prob(H->H) at replication)
                # a_MaintFail = epsilon_*PI_i{1. - gamma[r_i]}
                # gamma[r_i] = gamma0_     for r_temp < r_gamma_
                # gamma[r_i] = gamma0_*abs(r_gamma_-r_div_gamma_)**(lambda_gamma_)*
                #              abs(r_temp - r_div_gamma_)**(-lambda_gamma_)     for r_temp >= r_gamma_

                if r_temp < r_gamma_:
                    a_MaintFail_eff_ = 1.0 - gamma0_
                elif r_temp >= r_gamma_:
                    a_MaintFail_eff_ = 1.0 - gamma0_*abs(r_gamma_-r_div_gamma_)**lambda_gamma_ \
                                        *abs(r_temp-r_div_gamma_)**(-lambda_gamma_)

                # update coop_maint propensity matrix-element 
                coop_maint_propensities_array_[i_target_,i_passive_] = a_MaintFail_eff_
        
            # otherwise, no coop. maint. so set to '1'
            else:
                coop_maint_propensities_array_[i_target_,i_passive_] = 1.

    
    # if i_target is U site:
    elif state_in_[i_target_] == 0:
        
        #set row of coop. maint. propensities to zero:
        coop_maint_propensities_array_[i_target_,:] = 0
        
        # calc. coop. gain propensities for row:
        #loop over i_passive_
        for i_passive_ in range(i_target_):

        # 1. i_passive_ < i_target: r_LR_temp = abs(CG_positions_[i_passive_] - CG_positions_[i_target_] + r_LR1_)

            # identify it M site (to promote cooperative gain)
            if state_in_[i_passive_] == 2:

                r_temp = abs(CG_positions_[i_passive_] - CG_positions_[i_target_])
                r_LR_temp = abs(CG_positions_[i_passive_] - CG_positions_[i_target_] + r_LR1_)

                # cal. short-range-coop. effective propensity
                # a_coop_sr_eff_ = (a_coop_sr_)(1/2)
                # a_coop_sr_eff_ = (2*u_M_plus_PL_)(1/2) = u_M_plus_PL_    for r_temp < r_plat_
                # a_coop_sr_eff_ = (2.*u_M_plus_PL_*abs(r_plat_-r_div_)**lambda_coop_
                #                   *abs(r_temp-r_div_)**(-lambda_coop_))*(1/2)    for r_temp >= r_plat_
                #                = u_M_plus_PL_*abs(r_plat_-r_div_)**lambda_coop_
                #                   *abs(r_temp-r_div_)**(-lambda_coop_)     for r_temp >= r_plat_
                if r_temp < r_plat_:
                    a_coop_sr_eff_ = u_M_plus_PL_
                elif r_temp >= r_plat_:
                    a_coop_sr_eff_ = u_M_plus_PL_*abs(r_plat_-r_div_)**lambda_coop_ \
                                        *abs(r_temp-r_div_)**(-lambda_coop_)

                # cal. long-range-coop. effective propensity
                # a_coop_LR_eff_ = a_coop_LR_ / 2
                # a_coop_LR_ = 2.*u_M_plus_PL_LR1_    for r_LR_temp >= r_plat_
                # a_coop_LR_ = 2.*u_M_plus_PL_LR1_*abs(r_plat_LR1_-r_div_LR1_)**lambda_coopIn_LR1_     
                #               *abs(r_LR_temp-r_div_LR1_)**(-lambda_coopIn_LR1_)*(1/2)     for r_LR_temp >= r_plat_LR1_ AND r_temp < r_LR1_
                # a_coop_LR_ = 2.*u_M_plus_PL_LR1_*abs(r_plat_LR1_-r_div_LR1_)**lambda_coopOut_LR1_     
                #               *abs(r_LR_temp-r_div_LR1_)**(-lambda_coopOut_LR1_)*(1/2)     for r_LR_temp >= r_plat_LR1_ AND r_temp > r_LR1_

                if r_LR_temp < r_plat_LR1_:
                    a_coop_LR_eff_ = u_M_plus_PL_LR1_
                elif r_LR_temp >= r_plat_LR1_:
                    # Use separate powerlaws depending if a) jumping gap 'IN towards nucleosome'
                    if r_temp < r_LR1_:
                        a_coop_LR_eff_ = u_M_plus_PL_LR1_*abs(r_plat_LR1_-r_div_LR1_)**lambda_coopIn_LR1_ \
                                           *abs(r_LR_temp-r_div_LR1_)**(-lambda_coopIn_LR1_)
                    # Use separate powerlaws depending if b) jumping gap 'OUT away from nucleosome'
                    elif r_temp > r_LR1_:
                        a_coop_LR_eff_ = u_M_plus_PL_LR1_*abs(r_plat_LR1_-r_div_LR1_)**lambda_coopOut_LR1_ \
                                           *abs(r_LR_temp-r_div_LR1_)**(-lambda_coopOut_LR1_)


                # update coop_gain propensity matrix-element 
                coop_gains_propensities_array_[i_target_,i_passive_] = a_coop_sr_eff_ + a_coop_LR_eff_

        # previously i_target_ = M so coop. gain was imposible. Diagonal matrix element already 0.
                
        for i_passive_ in range(i_target_+1, N_CG_):
        #loop over i_passive_
        # 1. i_passive_ > i_target: r_LR_temp = abs(CG_positions_[i_passive_] - CG_positions_[i_target_] - r_LR1_)

            # identify it M site (to promote cooperative gain)
            if state_in_[i_passive_] == 2:

                r_temp = abs(CG_positions_[i_passive_] - CG_positions_[i_target_])
                r_LR_temp = abs(CG_positions_[i_passive_] - CG_positions_[i_target_] - r_LR1_)

                # cal. short-range-coop. effective propensity
                # a_coop_sr_eff_ = (a_coop_sr_)(1/2)
                # a_coop_sr_eff_ = (2*u_M_plus_PL_)(1/2) = u_M_plus_PL_    for r_temp < r_plat_
                # a_coop_sr_eff_ = (2.*u_M_plus_PL_*abs(r_plat_-r_div_)**lambda_coop_
                #                   *abs(r_temp-r_div_)**(-lambda_coop_))*(1/2)    for r_temp >= r_plat_
                #                = u_M_plus_PL_*abs(r_plat_-r_div_)**lambda_coop_
                #                   *abs(r_temp-r_div_)**(-lambda_coop_)     for r_temp >= r_plat_
                if r_temp < r_plat_:
                    a_coop_sr_eff_ = u_M_plus_PL_
                elif r_temp >= r_plat_:
                    a_coop_sr_eff_ = u_M_plus_PL_*abs(r_plat_-r_div_)**lambda_coop_ \
                                        *abs(r_temp-r_div_)**(-lambda_coop_)

                # cal. long-range-coop. effective propensity
                # a_coop_LR_eff_ = a_coop_LR_ / 2
                # a_coop_LR_ = 2.*u_M_plus_PL_LR1_    for r_LR_temp >= r_plat_
                # a_coop_LR_ = 2.*u_M_plus_PL_LR1_*abs(r_plat_LR1_-r_div_LR1_)**lambda_coopIn_LR1_     
                #               *abs(r_LR_temp-r_div_LR1_)**(-lambda_coopIn_LR1_)*(1/2)     for r_LR_temp >= r_plat_LR1_ AND r_temp < r_LR1_
                # a_coop_LR_ = 2.*u_M_plus_PL_LR1_*abs(r_plat_LR1_-r_div_LR1_)**lambda_coopOut_LR1_     
                #               *abs(r_LR_temp-r_div_LR1_)**(-lambda_coopOut_LR1_)*(1/2)     for r_LR_temp >= r_plat_LR1_ AND r_temp > r_LR1_

                if r_LR_temp < r_plat_LR1_:
                    a_coop_LR_eff_ = u_M_plus_PL_LR1_
                elif r_LR_temp >= r_plat_LR1_:
                    # Use separate powerlaws depending if a) jumping gap 'IN towards nucleosome'
                    if r_temp < r_LR1_:
                        a_coop_LR_eff_ = u_M_plus_PL_LR1_*abs(r_plat_LR1_-r_div_LR1_)**lambda_coopIn_LR1_ \
                                           *abs(r_LR_temp-r_div_LR1_)**(-lambda_coopIn_LR1_)
                    # Use separate powerlaws depending if b) jumping gap 'OUT away from nucleosome'
                    elif r_temp > r_LR1_:
                        a_coop_LR_eff_ = u_M_plus_PL_LR1_*abs(r_plat_LR1_-r_div_LR1_)**lambda_coopOut_LR1_ \
                                           *abs(r_LR_temp-r_div_LR1_)**(-lambda_coopOut_LR1_)


                # update coop_gain propensity matrix-element 
                coop_gains_propensities_array_[i_target_,i_passive_] = a_coop_sr_eff_ + a_coop_LR_eff_

                    
                    
                # end of loops
                
                
    # b) use i_passive_ = change_index_
    i_passive_ = change_index_
    
    # if i_passive is U site:
    if state_in_[i_passive_] == 0:
        # for coop. gain set whole column to zero: 
        coop_gains_propensities_array_[:,i_passive_] = 0
        
        # calc. coop. maint. propensities for column:
        # no coop. maint. so propensity = 0 if i_target = 0 and propensity = 1 if i_targent = 2
        for i_target_ in range(N_CG_):
            if state_in_[i_target_] == 0:
                # U sites: loss not possilbe
                coop_maint_propensities_array_[i_target_,i_passive_] = 0
            elif state_in_[i_target_] == 2:
                # M site: so could loose spontaneously
                coop_maint_propensities_array_[i_target_,i_passive_] = 1.
     
        
    # if i_passive is M:
    elif state_in_[i_passive_] == 2:
        
        # calc. coop. maint. propensities for column:
        # here i_target_ state unchanged, so only need to update for i_target_ = M. 
        # When i_target_ = U matrix element for coop. maint. already 0
        for i_target_ in range(i_passive_):
            
            # identify if M site (to promote coop. maint.)
            if state_in_[i_target_] == 2:

                r_temp = abs(CG_positions_[i_passive_] - CG_positions_[i_target_])

                # cal. coop. maint. failure propensity
                # a_MaintFail_eff_ = Prob(M->H)Prob(H->H)Prob(H->U)   (all probs at replication)
                #                  = a_MaintFail/2    (if here a_MaintFail == Prob(H->H) at replication)
                # a_MaintFail = epsilon_*PI_i{1. - gamma[r_i]}
                # gamma[r_i] = gamma0_     for r_temp < r_gamma_
                # gamma[r_i] = gamma0_*abs(r_gamma_-r_div_gamma_)**(lambda_gamma_)*
                #              abs(r_temp - r_div_gamma_)**(-lambda_gamma_)     for r_temp >= r_gamma_

                if r_temp < r_gamma_:
                    a_MaintFail_eff_ = 1.0 - gamma0_
                elif r_temp >= r_gamma_:
                    a_MaintFail_eff_ = 1.0 - gamma0_*abs(r_gamma_-r_div_gamma_)**lambda_gamma_ \
                                        *abs(r_temp-r_div_gamma_)**(-lambda_gamma_)

                # update coop_maint propensity matrix-element 
                coop_maint_propensities_array_[i_target_,i_passive_] = a_MaintFail_eff_
                
        # site at i_target_ = i_passive_ was set in section (a)
        
        for i_target_ in range(i_passive_+1, N_CG_):
            
            # identify if M site (to promote coop. maint.)
            if state_in_[i_target_] == 2:

                r_temp = abs(CG_positions_[i_passive_] - CG_positions_[i_target_])

                # cal. coop. maint. failure propensity
                # a_MaintFail_eff_ = Prob(M->H)Prob(H->H)Prob(H->U)   (all probs at replication)
                #                  = a_MaintFail/2    (if here a_MaintFail == Prob(H->H) at replication)
                # a_MaintFail = epsilon_*PI_i{1. - gamma[r_i]}
                # gamma[r_i] = gamma0_     for r_temp < r_gamma_
                # gamma[r_i] = gamma0_*abs(r_gamma_-r_div_gamma_)**(lambda_gamma_)*
                #              abs(r_temp - r_div_gamma_)**(-lambda_gamma_)     for r_temp >= r_gamma_

                if r_temp < r_gamma_:
                    a_MaintFail_eff_ = 1.0 - gamma0_
                elif r_temp >= r_gamma_:
                    a_MaintFail_eff_ = 1.0 - gamma0_*abs(r_gamma_-r_div_gamma_)**lambda_gamma_ \
                                        *abs(r_temp-r_div_gamma_)**(-lambda_gamma_)

                # update coop_maint propensity matrix-element 
                coop_maint_propensities_array_[i_target_,i_passive_] = a_MaintFail_eff_
        
        

        # calc. coop. gain propensities for column:
        # here i_target_ state unchanged, so only need to update for i_target_ = U. 
        # When i_target_ = M matrix element for coop. maint. already 0
        for i_target_ in range(i_passive_):
        #loop over i_target_
        # 1. i_passive_ > i_target: r_LR_temp = abs(CG_positions_[i_passive_] - CG_positions_[i_target_] - r_LR1_)

            # identify it M site (to promote cooperative gain)
            if state_in_[i_target_] == 0:

                r_temp = abs(CG_positions_[i_passive_] - CG_positions_[i_target_])
                r_LR_temp = abs(CG_positions_[i_passive_] - CG_positions_[i_target_] - r_LR1_)

                # cal. short-range-coop. effective propensity
                # a_coop_sr_eff_ = (a_coop_sr_)(1/2)
                # a_coop_sr_eff_ = (2*u_M_plus_PL_)(1/2) = u_M_plus_PL_    for r_temp < r_plat_
                # a_coop_sr_eff_ = (2.*u_M_plus_PL_*abs(r_plat_-r_div_)**lambda_coop_
                #                   *abs(r_temp-r_div_)**(-lambda_coop_))*(1/2)    for r_temp >= r_plat_
                #                = u_M_plus_PL_*abs(r_plat_-r_div_)**lambda_coop_
                #                   *abs(r_temp-r_div_)**(-lambda_coop_)     for r_temp >= r_plat_
                if r_temp < r_plat_:
                    a_coop_sr_eff_ = u_M_plus_PL_
                elif r_temp >= r_plat_:
                    a_coop_sr_eff_ = u_M_plus_PL_*abs(r_plat_-r_div_)**lambda_coop_ \
                                        *abs(r_temp-r_div_)**(-lambda_coop_)

                # cal. long-range-coop. effective propensity
                # a_coop_LR_eff_ = a_coop_LR_ / 2
                # a_coop_LR_ = 2.*u_M_plus_PL_LR1_    for r_LR_temp >= r_plat_
                # a_coop_LR_ = 2.*u_M_plus_PL_LR1_*abs(r_plat_LR1_-r_div_LR1_)**lambda_coopIn_LR1_     
                #               *abs(r_LR_temp-r_div_LR1_)**(-lambda_coopIn_LR1_)*(1/2)     for r_LR_temp >= r_plat_LR1_ AND r_temp < r_LR1_
                # a_coop_LR_ = 2.*u_M_plus_PL_LR1_*abs(r_plat_LR1_-r_div_LR1_)**lambda_coopOut_LR1_     
                #               *abs(r_LR_temp-r_div_LR1_)**(-lambda_coopOut_LR1_)*(1/2)     for r_LR_temp >= r_plat_LR1_ AND r_temp > r_LR1_

                if r_LR_temp < r_plat_LR1_:
                    a_coop_LR_eff_ = u_M_plus_PL_LR1_
                elif r_LR_temp >= r_plat_LR1_:
                    # Use separate powerlaws depending if a) jumping gap 'IN towards nucleosome'
                    if r_temp < r_LR1_:
                        a_coop_LR_eff_ = u_M_plus_PL_LR1_*abs(r_plat_LR1_-r_div_LR1_)**lambda_coopIn_LR1_ \
                                           *abs(r_LR_temp-r_div_LR1_)**(-lambda_coopIn_LR1_)
                    # Use separate powerlaws depending if b) jumping gap 'OUT away from nucleosome'
                    elif r_temp > r_LR1_:
                        a_coop_LR_eff_ = u_M_plus_PL_LR1_*abs(r_plat_LR1_-r_div_LR1_)**lambda_coopOut_LR1_ \
                                           *abs(r_LR_temp-r_div_LR1_)**(-lambda_coopOut_LR1_)


                # update coop_gain propensity matrix-element 
                coop_gains_propensities_array_[i_target_,i_passive_] = a_coop_sr_eff_ + a_coop_LR_eff_



        for i_target_ in range(i_passive_+1, N_CG_):
        #loop over i_target_
        # 1. i_passive_ < i_target: r_LR_temp = abs(CG_positions_[i_passive_] - CG_positions_[i_target_] + r_LR1_)

            # identify it M site (to promote cooperative gain)
            if state_in_[i_target_] == 0:

                r_temp = abs(CG_positions_[i_passive_] - CG_positions_[i_target_])
                r_LR_temp = abs(CG_positions_[i_passive_] - CG_positions_[i_target_] + r_LR1_)

                # cal. short-range-coop. effective propensity
                # a_coop_sr_eff_ = (a_coop_sr_)(1/2)
                # a_coop_sr_eff_ = (2*u_M_plus_PL_)(1/2) = u_M_plus_PL_    for r_temp < r_plat_
                # a_coop_sr_eff_ = (2.*u_M_plus_PL_*abs(r_plat_-r_div_)**lambda_coop_
                #                   *abs(r_temp-r_div_)**(-lambda_coop_))*(1/2)    for r_temp >= r_plat_
                #                = u_M_plus_PL_*abs(r_plat_-r_div_)**lambda_coop_
                #                   *abs(r_temp-r_div_)**(-lambda_coop_)     for r_temp >= r_plat_
                if r_temp < r_plat_:
                    a_coop_sr_eff_ = u_M_plus_PL_
                elif r_temp >= r_plat_:
                    a_coop_sr_eff_ = u_M_plus_PL_*abs(r_plat_-r_div_)**lambda_coop_ \
                                        *abs(r_temp-r_div_)**(-lambda_coop_)

                # cal. long-range-coop. effective propensity
                # a_coop_LR_eff_ = a_coop_LR_ / 2
                # a_coop_LR_ = 2.*u_M_plus_PL_LR1_    for r_LR_temp >= r_plat_
                # a_coop_LR_ = 2.*u_M_plus_PL_LR1_*abs(r_plat_LR1_-r_div_LR1_)**lambda_coopIn_LR1_     
                #               *abs(r_LR_temp-r_div_LR1_)**(-lambda_coopIn_LR1_)*(1/2)     for r_LR_temp >= r_plat_LR1_ AND r_temp < r_LR1_
                # a_coop_LR_ = 2.*u_M_plus_PL_LR1_*abs(r_plat_LR1_-r_div_LR1_)**lambda_coopOut_LR1_     
                #               *abs(r_LR_temp-r_div_LR1_)**(-lambda_coopOut_LR1_)*(1/2)     for r_LR_temp >= r_plat_LR1_ AND r_temp > r_LR1_

                if r_LR_temp < r_plat_LR1_:
                    a_coop_LR_eff_ = u_M_plus_PL_LR1_
                elif r_LR_temp >= r_plat_LR1_:
                    # Use separate powerlaws depending if a) jumping gap 'IN towards nucleosome'
                    if r_temp < r_LR1_:
                        a_coop_LR_eff_ = u_M_plus_PL_LR1_*abs(r_plat_LR1_-r_div_LR1_)**lambda_coopIn_LR1_ \
                                           *abs(r_LR_temp-r_div_LR1_)**(-lambda_coopIn_LR1_)
                    # Use separate powerlaws depending if b) jumping gap 'OUT away from nucleosome'
                    elif r_temp > r_LR1_:
                        a_coop_LR_eff_ = u_M_plus_PL_LR1_*abs(r_plat_LR1_-r_div_LR1_)**lambda_coopOut_LR1_ \
                                           *abs(r_LR_temp-r_div_LR1_)**(-lambda_coopOut_LR1_)


                # update coop_gain propensity matrix-element 
                coop_gains_propensities_array_[i_target_,i_passive_] = a_coop_sr_eff_ + a_coop_LR_eff_
                    
                    
                # end of loops
                

    # update spontaneous gains
    i_target_ = change_index_

    # if M site set this element to zero
    if state_in_[i_target_] == 2:
        spont_gains_propensities_vector_[i_target_] = 0
    
    # if U site add in spontaneous gain propensity 
    elif state_in_[i_target_] == 0:

        # cal. background gain effective propensity
        # a_0_eff_ = (a_0_)(1/2) = 2delta_ /2      # note 2delta_ = U->H probability

        # add a_0_eff_ to existing propensity from coop. gains. 
        spont_gains_propensities_vector_[i_target_] = delta_
                    
                  
                
    # calculate composite propensities vector (one propensity for each site in gene)
    # sum over i_passive_ (ie. along rows) to create coop_gains_propensities_vector_
    coop_gains_propensities_vector_ = np.sum(coop_gains_propensities_array_, axis=1)
    # find product over i_passive_ (i.e. along rows) to create coop_maint_propensities_vector
    coop_maint_propensities_vector_ = (epsilon_/2.)*np.prod( coop_maint_propensities_array_,axis=1)

    total_propensities_vector_ = coop_gains_propensities_vector_ + \
                                    spont_gains_propensities_vector_ + coop_maint_propensities_vector_

    
    return coop_gains_propensities_array_, spont_gains_propensities_vector_, \
            coop_maint_propensities_array_, total_propensities_vector_




# function to count methylation levels in state
def count_state(state_,N_CG_):
# update so that n_m = N_m/(N_m+N_u)
# keep n_u and n_h normalised to N_CG,
# so could still reconstruct N_M (if knew n_u though)

    # count number of sites in each state
    N_u_ = 0.
    N_h_ = 0.
    N_m_ = 0.
    for i_ in range(N_CG_):
        if state_[i_] == 0:
            N_u_ += 1.
        elif state_[i_] == 2:
            N_m_ += 1.
        else:
            N_h_ += 1.

    n_h_ = N_h_/N_CG_
    if (N_m_+N_u_) > 0:
        n_m_ = N_m_/(N_m_ + N_u_)
        n_u_ = N_u_/(N_m_ + N_u_)
    else:
        n_m_ = 0
        n_u_ = 0

    return (n_u_,n_h_,n_m_)



def whole_locus_gaps_analysis(N_, state_in_, CG_positions_):
    UU_gaps_ = []
    UM_gaps_ = []
    MM_gaps_ = []
    XX_gaps_ = []
    
    for i_ in range(0,N_-1):
        
        if (state_in_[i_] == 0 and state_in_[i_+1] == 0):
            UU_gaps_.append( CG_positions_[i_+1] - CG_positions_[i_] )
        elif (state_in_[i_] == 2 and state_in_[i_+1] == 0):
            UM_gaps_.append( CG_positions_[i_+1] - CG_positions_[i_] )
        elif (state_in_[i_] == 0 and state_in_[i_+1] == 2):
            UM_gaps_.append( CG_positions_[i_+1] - CG_positions_[i_] )
        elif (state_in_[i_] == 2 and state_in_[i_+1] ==2):
            MM_gaps_.append( CG_positions_[i_+1] - CG_positions_[i_] )
        else:
            XX_gaps_.append( CG_positions_[i_+1] - CG_positions_[i_] )
            #print(i_, state_in_[i_], i_+1, state_in_[i_+1])
    
    return UU_gaps_, UM_gaps_, MM_gaps_, XX_gaps_    




def whole_locus_nearest_gaps_analysis(N_, state_in_, CG_positions_):
    UU_gaps_ = []
    UM_gaps_ = []
    MM_gaps_ = []
    XX_gaps_ = []
    
    for i_ in range(1,N_-1):
        
        gap_left_temp = CG_positions_[i_] - CG_positions_[i_-1]
        gap_right_temp = CG_positions_[i_+1] - CG_positions_[i_]
        left_temp = 'XX'
        right_temp = 'XX'
        if state_in_[i_] == 0:
            if state_in_[i_-1] == 0:
                left_temp = 'UU'
            elif state_in_[i_-1] == 2:
                left_temp = 'UM'
            if state_in_[i_+1] == 0:
                right_temp = 'UU'
            elif state_in_[i_+1] == 2:
                right_temp = 'UM'
        if state_in_[i_] == 2:
            if state_in_[i_-1] == 0:
                left_temp = 'UM'
            elif state_in_[i_-1] == 2:
                left_temp = 'MM'
            if state_in_[i_+1] == 0:
                right_temp = 'UM'
            elif state_in_[i_+1] == 2:
                right_temp = 'MM'
        
        # if both gaps equal.. defult to the left gap
        if gap_left_temp > gap_right_temp:
            gap_nearest_temp = gap_right_temp
            nearest_temp = right_temp
        else:
            gap_nearest_temp = gap_left_temp
            nearest_temp = left_temp
            
            
        if nearest_temp == 'UU':
            UU_gaps_.append( gap_nearest_temp )
        elif nearest_temp == 'UM':
            UM_gaps_.append( gap_nearest_temp)
        elif nearest_temp == 'MM':
            MM_gaps_.append( gap_nearest_temp)
        elif gap_nearest_temp:
            XX_gaps_.append( gap_nearest_temp)
            #print(i_, state_in_[i_], i_+1, state_in_[i_+1])
    
    return UU_gaps_, UM_gaps_, MM_gaps_, XX_gaps_




def single_gene_sim_eqbrm_single_time_stats(N_CG_, state_in_, CG_positions_,
                       spont_gain_params_, coop_gain_params_, coop_maint_params_):
    # unnecessary inputs? : data_number_state_, data_number_mean_properties_, N_CG_, N_gen_
    

    #####    # N_gen_cc_max = maximum number of cell cycles to simulate for
    # burn_in_time = when to start recording statistics
    #####    # N_output = number of times to write out data (unless N_gen_cc_max exceeded)
    #####    # min_output_time_interval = space data write outs by at least this timegap
    
    #####    # next_output_time_ = variable to decide when to trigger next data write out
    #####    # i_output_ = varialbe to count how many outputs have been made so far

    state_ = state_in_.copy()
    
    # variables to store statistics
    # whole gene quantaties
    UU_gaps_ = []
    UM_gaps_ = []
    MM_gaps_ = []
    XX_gaps_ = []

    UU_nearest_gaps_ = []
    UM_nearest_gaps_ = []
    MM_nearest_gaps_ = []
    XX_nearest_gaps_ = []
    
    single_whole_gene_mean_M_level_ = []
    single_whole_gene_mean_M_time_ = []
    

    # initilise system
    time_ = 0.0
    
    
    # calculate initial propensities
    coop_gains_propensities_array_, spont_gains_propensities_vector_, \
            coop_maint_propensities_array_, total_propensities_vector_ = \
                initial_propensities(N_CG_, state_, CG_positions_, spont_gain_params_, 
                                         coop_gain_params_, coop_maint_params_)
    
    
    # churn through the simulation
    
    while time_ < burn_in_time:
        
        # implement timestep using Gillespies direct method

        # find indices that sort propensities into ascending order
        propensities_sorted_indicies_ = np.argsort(total_propensities_vector_)

        cumulative_propensities_ = np.cumsum( np.array(total_propensities_vector_)[propensities_sorted_indicies_] )

        # make array listing site indicies to apply propensities-ordering to so that the site 
        # which changes can be identified
        ordered_site_list_ = np.array( np.arange(N_CG_)[propensities_sorted_indicies_] )

        # get a pair of random numbers
        random_pair_ = np.random.rand(2)

        # calc. time incrememnt delta_t_
        delta_t_ = -np.log(random_pair_[0])*(1/cumulative_propensities_[-1])

        
        # check whether have overshot simulation end and if so break
        if (time_ + delta_t_) > burn_in_time:

            # set time to burn-in time
            time_ = burn_in_time
                
            # output stats. 
            # calc. whole gene states
            UU_gaps_temp, UM_gaps_temp, MM_gaps_temp, XX_gaps_temp = whole_locus_gaps_analysis(N_CG_, state_, 
                                                                                                 CG_positions_)
            UU_gaps_nearest_temp, UM_gaps_nearest_temp, MM_gaps_nearest_temp, XX_gaps_nearest_temp = \
            whole_locus_nearest_gaps_analysis(N_CG_, state_, CG_positions_)

            UU_gaps_.extend(UU_gaps_temp)
            UM_gaps_.extend(UM_gaps_temp)
            MM_gaps_.extend(MM_gaps_temp)
            XX_gaps_.extend(XX_gaps_temp)

            UU_nearest_gaps_.extend(UU_gaps_nearest_temp)
            UM_nearest_gaps_.extend(UM_gaps_nearest_temp)
            MM_nearest_gaps_.extend(MM_gaps_nearest_temp)
            XX_nearest_gaps_.extend(XX_gaps_nearest_temp)

            # analyse state
            n_u_, n_h_, n_m_ = count_state(state_,N_CG_)

            single_whole_gene_mean_M_level_.append(n_m_)
            single_whole_gene_mean_M_time_.append(time_)

            # end of stats output
            # exit loop without updating state
            break
            
        # if haven't overshot burn_in_time, will not have broken loop, so, complete the timestep updates

        # cal. which site will change (specified by change_index_)
        propensities_cut_off_ = random_pair_[1]*cumulative_propensities_[-1]

        change_index_ = ordered_site_list_[ np.searchsorted(cumulative_propensities_, propensities_cut_off_) ]


        # update the time
        time_ += delta_t_

        # update state
        if state_[change_index_] == 0:
            state_[change_index_] += 2
        elif state_[change_index_] == 2:
            state_[change_index_] -= 2
        else:
            print('ahdla')

        # update the propensiteis
        coop_gains_propensities_array_, spont_gains_propensities_vector_, \
            coop_maint_propensities_array_, total_propensities_vector_ = \
            update_propensities(N_CG_, state_, CG_positions_, spont_gain_params_, 
                                         coop_gain_params_, coop_maint_params_, coop_gains_propensities_array_, 
                            spont_gains_propensities_vector_, coop_maint_propensities_array_, change_index_)


        # end of timestep
        
    
            
        
    # end of: while time_ < burn_in_time: loop
    
    
    #    # calc. summary statistics
    #    mu_whole_gene_mean_M_level_ = np.mean(single_whole_gene_mean_M_level_)
    #    sigma_whole_gene_mean_M_level_ = np.std(single_whole_gene_mean_M_level_)
    
    # package sim_outputs
    # 1. 'eqbrm_stats_whole_gene_gaps'
    # 2. 'eqbrm_stats_whole_gene_meth_level'
    # 3. 'eqbrm_stats_summary_whole_gene_'
    
    eqbrm_stats_whole_gene_gaps_ = [UU_gaps_, UM_gaps_, MM_gaps_, XX_gaps_, 
                                  UU_nearest_gaps_, UM_nearest_gaps_, MM_nearest_gaps_, XX_nearest_gaps_]
    eqbrm_stats_whole_gene_meth_level_ = [single_whole_gene_mean_M_time_, single_whole_gene_mean_M_level_]
    #    eqbrm_stats_summary_whole_gene_ = [mu_whole_gene_mean_M_level_, sigma_whole_gene_mean_M_level_]
    
    #    sim_outputs_ = [eqbrm_stats_whole_gene_gaps_, eqbrm_stats_whole_gene_meth_level_, eqbrm_stats_summary_whole_gene_]
    sim_outputs_ = [eqbrm_stats_whole_gene_gaps_, eqbrm_stats_whole_gene_meth_level_, state_]
            
    return sim_outputs_



# function to update seed deterministically for each simulation replicate
def calc_seed(N_reps_, initial_seed_, current_segment_no_, i_reps_):
    current_seed_ = initial_seed_ + (current_segment_no_-1)*N_reps_ + i_reps_
    return current_seed_



# main code body!!
# loop to simulate using: single_gene_sim_eqbrm_stats

# define parameters

initial_seed = params.initial_seed

# possible values for P_choice
# "100% U"
# "50% U & 50% M"
# "100% M"

P_choice = params.P_choice

rep_time = params.rep_time
n_cc = params.n_cc
off_rate = params.off_rate

initial_state_choice = params.initial_state_choice

N_gen_burn_in = params.N_gen_burn_in

burn_in_time = N_gen_burn_in*n_cc

# define linear relationships:
u_scale_val = params.u_scale_val

u_m_val = params.u_m_val
u_c_val = params.u_c_val
u_cap_val = params.u_cap_val

e_m_val = params.e_m_val
e_c_val = params.e_c_val
e_cap_val = params.e_cap_val

g_m_val = params.g_m_val
g_c_val = params.g_c_val
g_cap_val = params.g_cap_val

# 0 < epsilon < 1
# 0 < gamma0 < 1
# Jay's parameters
#####epsilon = 3.2e-4
delta =  params.delta

#####gamma0 = 0.43
lambda_gamma = params.lambda_gamma
r_div_gamma = params.r_div_gamma
r_gamma = params.r_gamma

# Jay's params
#####u_M_plus_PL = 0.1E-04
lambda_coop = params.lambda_coop
r_div = params.r_div
r_plat = params.r_plat

# LR (long-range) parameters
#####u_M_plus_PL_LR1 = u_M_plus_PL*0.08
lambda_coopIn_LR1 = params.lambda_coopIn_LR1
lambda_coopOut_LR1 = params.lambda_coopOut_LR1
r_div_LR1 = params.r_div_LR1
r_plat_LR1 = params.r_plat_LR1

r_LR1 = params.r_LR1


# # Define arrays to store histogram count totals
# Col0_UU_gaps_counts = np.zeros(N_corrns_bins)
# Col0_MU_gaps_counts = np.zeros(N_corrns_bins)
# Col0_MM_gaps_counts = np.zeros(N_corrns_bins)
# Col0_XX_gaps_counts = np.zeros(N_corrns_bins)

# Col0_UU_nearest_gaps_counts = np.zeros(N_corrns_bins)
# Col0_MU_nearest_gaps_counts = np.zeros(N_corrns_bins)
# Col0_MM_nearest_gaps_counts = np.zeros(N_corrns_bins)
# Col0_XX_nearest_gaps_counts = np.zeros(N_corrns_bins)

# Data_UU_gaps_counts = np.zeros(N_corrns_bins)
# Data_MU_gaps_counts = np.zeros(N_corrns_bins)
# Data_MM_gaps_counts = np.zeros(N_corrns_bins)
# Data_XX_gaps_counts = np.zeros(N_corrns_bins)

# Data_UU_nearest_gaps_counts = np.zeros(N_corrns_bins)
# Data_MU_nearest_gaps_counts = np.zeros(N_corrns_bins)
# Data_MM_nearest_gaps_counts = np.zeros(N_corrns_bins)
# Data_XX_nearest_gaps_counts = np.zeros(N_corrns_bins)

# D4_UU_gaps_counts = np.zeros(N_corrns_bins)
# D4_MU_gaps_counts = np.zeros(N_corrns_bins)
# D4_MM_gaps_counts = np.zeros(N_corrns_bins)
# D4_XX_gaps_counts = np.zeros(N_corrns_bins)

# D4_UU_nearest_gaps_counts = np.zeros(N_corrns_bins)
# D4_MU_nearest_gaps_counts = np.zeros(N_corrns_bins)
# D4_MM_nearest_gaps_counts = np.zeros(N_corrns_bins)
# D4_XX_nearest_gaps_counts = np.zeros(N_corrns_bins)

# D5_UU_gaps_counts = np.zeros(N_corrns_bins)
# D5_MU_gaps_counts = np.zeros(N_corrns_bins)
# D5_MM_gaps_counts = np.zeros(N_corrns_bins)
# D5_XX_gaps_counts = np.zeros(N_corrns_bins)

# D5_UU_nearest_gaps_counts = np.zeros(N_corrns_bins)
# D5_MU_nearest_gaps_counts = np.zeros(N_corrns_bins)
# D5_MM_nearest_gaps_counts = np.zeros(N_corrns_bins)
# D5_XX_nearest_gaps_counts = np.zeros(N_corrns_bins)

# Sim_UU_gaps_counts = np.zeros(N_corrns_bins)
# Sim_MU_gaps_counts = np.zeros(N_corrns_bins)
# Sim_MM_gaps_counts = np.zeros(N_corrns_bins)
# Sim_XX_gaps_counts = np.zeros(N_corrns_bins)

# Sim_UU_nearest_gaps_counts = np.zeros(N_corrns_bins)
# Sim_MU_nearest_gaps_counts = np.zeros(N_corrns_bins)
# Sim_MM_nearest_gaps_counts = np.zeros(N_corrns_bins)
# Sim_XX_nearest_gaps_counts = np.zeros(N_corrns_bins)



timer_start_total = timer()
timer_start = timer()


i_gene_start = 0
i_gene_end = len(input_gene_ID_list)
i_gene = 0 # set this before loop incase gene list is actually empty on this core, in which case loop below won't activate


# counter to record how many segments modeled in total
i_seg_total = 0
i_reps = 0

N_CG_min = params.N_CG_min
N_CG_max = params.N_CG_max
N_CG_density_min = params.N_CG_density_min
N_CG_density_max = params.N_CG_density_max

# invalidate_current_ID
# if = 0 proceed
# if = 1 skip 


for i_gene in range(i_gene_start,i_gene_end):
#for i_gene in range(2000):

    # initialise validatation variable
    invalidate_current_ID = 0
    
    current_gene_ID = input_gene_ID_list[i_gene]
    # find how many segments exist for this gene
    current_segments_df = AllLoci_df.loc[[current_gene_ID]]
    current_number_of_segments = len(current_segments_df) # assume that this will be 1
    current_segment_no = current_segments_df.loc[current_gene_ID,'segment_no']
    mean_h2az = current_segments_df.loc[current_gene_ID,'mean_h2az']
    max_h2az = current_segments_df.loc[current_gene_ID,'max_h2az']

    # loop over number of segments in this gene
    for i_seg_current_gene in range(0,current_number_of_segments):
        
        # find start and end of current segment
        region_start_index = current_segments_df.columns.get_loc('start')
        region_end_index = current_segments_df.columns.get_loc('end')
        current_segment_start = current_segments_df.iloc[i_seg_current_gene,region_start_index]
        current_segment_end = current_segments_df.iloc[i_seg_current_gene,region_end_index]
        
        # calculate Col0 expt. state        
        # set current dataframe to first Col0 state
        current_CG_site_data_df = Col0_CG_site_data_df

        N_CG, CG_positions_gene, state_expt, state_gene = create_initial_state_eqbrm_Col0like_regions(current_gene_ID, 
                                    initial_state_choice, P_choice, current_segment_start, current_segment_end)
        
        
        if N_CG < N_CG_min:
            invalidate_current_ID = 1
            temp_output_list = [current_gene_ID, N_CG, np.nan]
            temp_output_file = file_name_output_Skipped_IDs
            with open(temp_output_file,'a',newline='') as output_file:
                line_writer = csv.writer(output_file,delimiter='\t')
                line_writer.writerow(temp_output_list)
            continue
        if N_CG >= N_CG_max:
            invalidate_current_ID = 1
            temp_output_list = [current_gene_ID, N_CG, np.nan]
            temp_output_file = file_name_output_Skipped_IDs
            with open(temp_output_file,'a',newline='') as output_file:
                line_writer = csv.writer(output_file,delimiter='\t')
                line_writer.writerow(temp_output_list)
            continue
            
        N_CG_density = (N_CG-1)/(CG_positions_gene[-1]-CG_positions_gene[0])
        L_locus = CG_positions_gene[-1]-CG_positions_gene[0]
            
        if N_CG_density < N_CG_density_min:
            invalidate_current_ID = 1
            temp_output_list = [current_gene_ID, N_CG, N_CG_density]
            temp_output_file = file_name_output_Skipped_IDs
            with open(temp_output_file,'a',newline='') as output_file:
                line_writer = csv.writer(output_file,delimiter='\t')
                line_writer.writerow(temp_output_list)
            continue
        if N_CG_density >= N_CG_density_max:
            invalidate_current_ID = 1
            temp_output_list = [current_gene_ID, N_CG, N_CG_density]
            temp_output_file = file_name_output_Skipped_IDs
            with open(temp_output_file,'a',newline='') as output_file:
                line_writer = csv.writer(output_file,delimiter='\t')
                line_writer.writerow(temp_output_list)
            continue
            
            
        # calculate variable paramters, setting to zero of linear function goes negative. 
        u_M_plus_PL = max(1.E-4*(u_m_val*N_CG_density  + u_c_val),1.E-4*u_cap_val)
        u_M_plus_PL_LR1 = u_scale_val*u_M_plus_PL
        epsilon = min(1.E-4*e_cap_val,1.E-4*(e_m_val*N_CG_density + e_c_val))
        gamma0 = max(1.E-4*g_cap_val,g_m_val*N_CG_density + g_c_val)
        ##print(u_M_plus_PL, epsilon, gamma0)

        # package spontaneous background gain parameters into list
        spont_gain_params = [delta]

        # package replication parameters into list
        coop_maint_params = [epsilon, gamma0, lambda_gamma, r_div_gamma, r_gamma]

        # package coop_gain_parameters into list:
        coop_gain_params = [u_M_plus_PL, lambda_coop, r_div, r_plat, 
                             u_M_plus_PL_LR1, lambda_coopIn_LR1, lambda_coopOut_LR1,r_div_LR1, r_plat_LR1, r_LR1]
            
            
        # # gaps analysis on expt. state
        # Col0_UU_gaps_temp, Col0_MU_gaps_temp, Col0_MM_gaps_temp, Col0_XX_gaps_temp = whole_locus_gaps_analysis(N_CG, 
        #                                                                                state_expt, CG_positions_gene)
        # Col0_UU_gaps_nearest_temp, Col0_MU_gaps_nearest_temp, Col0_MM_gaps_nearest_temp, Col0_XX_gaps_nearest_temp = \
        # whole_locus_nearest_gaps_analysis(N_CG, state_expt, CG_positions_gene)

        #find mean meth. level for experimental state
        n_u_expt_temp, n_h_expt_temp, n_m_expt_temp = count_state(state_expt,N_CG)
        Col0_locus_mean_M_level_temp = n_m_expt_temp
        Col0_locus_mean_X_level_temp = n_h_expt_temp

        # # no loop here so set pending values too... 
        # Col0_UU_gaps_pending = Col0_UU_gaps_temp
        # Col0_MU_gaps_pending = Col0_MU_gaps_temp
        # Col0_MM_gaps_pending = Col0_MM_gaps_temp
        # Col0_XX_gaps_pending = Col0_XX_gaps_temp
        
        # Col0_UU_nearest_gaps_pending = Col0_UU_gaps_nearest_temp
        # Col0_MU_nearest_gaps_pending = Col0_MU_gaps_nearest_temp
        # Col0_MM_nearest_gaps_pending = Col0_MM_gaps_nearest_temp
        # Col0_XX_nearest_gaps_pending = Col0_XX_gaps_nearest_temp

        Col0_locus_mean_M_level_pending = Col0_locus_mean_M_level_temp
        Col0_locus_mean_X_level_pending = Col0_locus_mean_X_level_temp

        
        # #loop over analysis of selected experimental states
        # # define pending oupt-put variables before start replicates
        # # Data_UU_gaps_pending = []
        # # Data_MU_gaps_pending = []
        # # Data_MM_gaps_pending = []
        # # Data_XX_gaps_pending = []
        
        # # Data_UU_nearest_gaps_pending = []
        # # Data_MU_nearest_gaps_pending = []
        # # Data_MM_nearest_gaps_pending = []
        # # Data_XX_nearest_gaps_pending = []
        
        # Dset1_locus_mean_M_level_pending = []
        # Dset1_locus_mean_X_level_pending = [] # For sites not U or M... will show up as H in count sites
        
        # for i_expts in range(0,N_ExptState_Dset1):
        # # set current dataframe to first expt. state
        #     current_CG_site_data_df = ExptState_Dset1_df_list[i_expts]

        #     N_CG_temp, CG_positions_gene_temp, state_expt_temp, state_gene_temp = \
        #                 create_initial_state_eqbrm_Col0like_regions(current_gene_ID, 
        #                         initial_state_choice, P_choice, current_segment_start, current_segment_end)
                
        #     # if N_CG different to Col0 values skip segment
        #     if N_CG_temp != N_CG:
        #         invalidate_current_ID = 1
        #         temp_output_list = [current_gene_ID, N_CG_temp, np.nan]
        #         temp_output_file = file_name_output_Skipped_IDs
        #         with open(temp_output_file,'a',newline='') as output_file:
        #             line_writer = csv.writer(output_file,delimiter='\t')
        #             line_writer.writerow(temp_output_list)
        #         break
            
        #     N_CG_density_temp = (N_CG_temp-1)/(CG_positions_gene_temp[-1]-CG_positions_gene_temp[0])

        #     # if N_CG_density different to Col0 values skip segment
        #     if N_CG_density_temp != N_CG_density:
        #         invalidate_current_ID = 1
        #         temp_output_list = [current_gene_ID, N_CG_temp, N_CG_density_temp]
        #         temp_output_file = file_name_output_Skipped_IDs
        #         with open(temp_output_file,'a',newline='') as output_file:
        #             line_writer = csv.writer(output_file,delimiter='\t')
        #             line_writer.writerow(temp_output_list)
        #         break

        #     # # gaps analysis on first expt. state
        #     # expt_UU_gaps_temp, expt_MU_gaps_temp, expt_MM_gaps_temp, expt_XX_gaps_temp = whole_locus_gaps_analysis(N_CG_temp, 
        #     #                                                                            state_expt_temp, CG_positions_gene_temp)
        #     # expt_UU_gaps_nearest_temp, expt_MU_gaps_nearest_temp, expt_MM_gaps_nearest_temp, expt_XX_gaps_nearest_temp = \
        #     # whole_locus_nearest_gaps_analysis(N_CG_temp, state_expt_temp, CG_positions_gene_temp)

        #     # Data_UU_gaps_pending.extend(expt_UU_gaps_temp)
        #     # Data_MU_gaps_pending.extend(expt_MU_gaps_temp)
        #     # Data_MM_gaps_pending.extend(expt_MM_gaps_temp)
        #     # Data_XX_gaps_pending.extend(expt_XX_gaps_temp)

        #     # Data_UU_nearest_gaps_pending.extend(expt_UU_gaps_nearest_temp)
        #     # Data_MU_nearest_gaps_pending.extend(expt_MU_gaps_nearest_temp)
        #     # Data_MM_nearest_gaps_pending.extend(expt_MM_gaps_nearest_temp)
        #     # Data_XX_nearest_gaps_pending.extend(expt_XX_gaps_nearest_temp)

        #     #find mean meth. level for experimental state
        #     n_u_expt_temp, n_h_expt_temp, n_m_expt_temp = count_state(state_expt_temp,N_CG_temp)
        #     Dset1_locus_mean_M_level_pending.append(n_m_expt_temp)
        #     Dset1_locus_mean_X_level_pending.append(n_h_expt_temp)
            
        #     # end of loop over expt. selected states (i_expts)

  

            
        # # check whether ID has been invalidated and if so move to next segment
        # if invalidate_current_ID == 1:
        #     continue
        
        # # loop over analysis of Decile4 experimental states
        # # define pending oupt-put variables before start replicates
        # # D4_UU_gaps_pending = []
        # # D4_MU_gaps_pending = []
        # # D4_MM_gaps_pending = []
        # # D4_XX_gaps_pending = []

        # # D4_UU_nearest_gaps_pending = []
        # # D4_MU_nearest_gaps_pending = []
        # # D4_MM_nearest_gaps_pending = []
        # # D4_XX_nearest_gaps_pending = []
        
        # Dset2_locus_mean_M_level_pending = []
        
        # for i_expts in range(0,N_ExptState_Dset2):
        # # set current dataframe to first expt. state
        #     current_CG_site_data_df = ExptState_Dset2_df_list[i_expts]

        #     N_CG_temp, CG_positions_gene_temp, state_expt_temp, state_gene_temp = \
        #                 create_initial_state_eqbrm_Col0like_regions(current_gene_ID, 
        #                         initial_state_choice, P_choice, current_segment_start, current_segment_end)
                
        #     # if N_CG different to Col0 values skip segment
        #     if N_CG_temp != N_CG:
        #         invalidate_current_ID = 1
        #         temp_output_list = [current_gene_ID, N_CG_temp, np.nan]
        #         temp_output_file = file_name_output_Skipped_IDs
        #         with open(temp_output_file,'a',newline='') as output_file:
        #             line_writer = csv.writer(output_file,delimiter='\t')
        #             line_writer.writerow(temp_output_list)
        #         break
            
        #     N_CG_density_temp = (N_CG_temp-1)/(CG_positions_gene_temp[-1]-CG_positions_gene_temp[0])

        #     # if N_CG_density different to Col0 values skip segment
        #     if N_CG_density_temp != N_CG_density:
        #         invalidate_current_ID = 1
        #         temp_output_list = [current_gene_ID, N_CG_temp, N_CG_densitytemp]
        #         temp_output_file = file_name_output_Skipped_IDs
        #         with open(temp_output_file,'a',newline='') as output_file:
        #             line_writer = csv.writer(output_file,delimiter='\t')
        #             line_writer.writerow(temp_output_list)
        #         break

        #     # # gaps analysis on first expt. state
        #     # expt_UU_gaps_temp, expt_MU_gaps_temp, expt_MM_gaps_temp, expt_XX_gaps_temp = whole_locus_gaps_analysis(N_CG_temp, 
        #     #                                                                            state_expt_temp, CG_positions_gene_temp)
        #     # expt_UU_gaps_nearest_temp, expt_MU_gaps_nearest_temp, expt_MM_gaps_nearest_temp, expt_XX_gaps_nearest_temp = \
        #     # whole_locus_nearest_gaps_analysis(N_CG_temp, state_expt_temp, CG_positions_gene_temp)

        #     # D4_UU_gaps_pending.extend(expt_UU_gaps_temp)
        #     # D4_MU_gaps_pending.extend(expt_MU_gaps_temp)
        #     # D4_MM_gaps_pending.extend(expt_MM_gaps_temp)
        #     # D4_XX_gaps_pending.extend(expt_XX_gaps_temp)

        #     # D4_UU_nearest_gaps_pending.extend(expt_UU_gaps_nearest_temp)
        #     # D4_MU_nearest_gaps_pending.extend(expt_MU_gaps_nearest_temp)
        #     # D4_MM_nearest_gaps_pending.extend(expt_MM_gaps_nearest_temp)
        #     # D4_XX_nearest_gaps_pending.extend(expt_XX_gaps_nearest_temp)

        #     #find mean meth. level for experimental state
        #     n_u_expt_temp, n_h_expt_temp, n_m_expt_temp = count_state(state_expt_temp,N_CG_temp)
        #     Dset2_locus_mean_M_level_pending.append(n_m_expt_temp)
            
        #     # end of loop over expt. Decile4 states (i_expts)
            
        # # check whether ID has been invalidated and if so move to next segment
        # if invalidate_current_ID == 1:
        #     continue
        
        # # # loop over analysis of Decile5 experimental states
        # # # define pending oupt-put variables before start replicates
        # # D5_UU_gaps_pending = []
        # # D5_MU_gaps_pending = []
        # # D5_MM_gaps_pending = []
        # # D5_XX_gaps_pending = []
        
        # # D5_UU_nearest_gaps_pending = []
        # # D5_MU_nearest_gaps_pending = []
        # # D5_MM_nearest_gaps_pending = []
        # # D5_XX_nearest_gaps_pending = []
        
        # # D5_locus_mean_M_level_pending = []
        
        # # for i_expts in range(0,N_ExptState_Decile5):
        # # # set current dataframe to first expt. state
        # #     current_CG_site_data_df = ExptState_Decile5_df_list[i_expts]

        # #     N_CG_temp, CG_positions_gene_temp, state_expt_temp, state_gene_temp = \
        # #                 create_initial_state_eqbrm_Col0like_regions(current_gene_ID, 
        # #                         initial_state_choice, P_choice, current_segment_start, current_segment_end)
                
        # #     # if N_CG different to Col0 values skip segment
        # #     if N_CG_temp != N_CG:
        # #         invalidate_current_ID = 1
        # #         temp_output_list = [current_gene_ID, N_CG_temp, np.nan]
        # #         temp_output_file = file_name_output_Skipped_IDs
        # #         with open(temp_output_file,'a',newline='') as output_file:
        # #             line_writer = csv.writer(output_file,delimiter='\t')
        # #             line_writer.writerow(temp_output_list)
        # #         break
            
        # #     N_CG_density_temp = (N_CG_temp-1)/(CG_positions_gene_temp[-1]-CG_positions_gene_temp[0])

        # #     # if N_CG_density different to Col0 values skip segment
        # #     if N_CG_density_temp != N_CG_density:
        # #         invalidate_current_ID = 1
        # #         temp_output_list = [current_gene_ID, N_CG_temp, N_CG_density_temp]
        # #         temp_output_file = file_name_output_Skipped_IDs
        # #         with open(temp_output_file,'a',newline='') as output_file:
        # #             line_writer = csv.writer(output_file,delimiter='\t')
        # #             line_writer.writerow(temp_output_list)
        # #         break

        # #     # gaps analysis on first expt. state
        # #     expt_UU_gaps_temp, expt_MU_gaps_temp, expt_MM_gaps_temp, expt_XX_gaps_temp = whole_locus_gaps_analysis(N_CG_temp, 
        # #                                                                                state_expt_temp, CG_positions_gene_temp)
        # #     expt_UU_gaps_nearest_temp, expt_MU_gaps_nearest_temp, expt_MM_gaps_nearest_temp, expt_XX_gaps_nearest_temp = \
        # #     whole_locus_nearest_gaps_analysis(N_CG_temp, state_expt_temp, CG_positions_gene_temp)

        # #     D5_UU_gaps_pending.extend(expt_UU_gaps_temp)
        # #     D5_MU_gaps_pending.extend(expt_MU_gaps_temp)
        # #     D5_MM_gaps_pending.extend(expt_MM_gaps_temp)
        # #     D5_XX_gaps_pending.extend(expt_XX_gaps_temp)

        # #     D5_UU_nearest_gaps_pending.extend(expt_UU_gaps_nearest_temp)
        # #     D5_MU_nearest_gaps_pending.extend(expt_MU_gaps_nearest_temp)
        # #     D5_MM_nearest_gaps_pending.extend(expt_MM_gaps_nearest_temp)
        # #     D5_XX_nearest_gaps_pending.extend(expt_XX_gaps_nearest_temp)

        # #     #find mean meth. level for experimental state
        # #     n_u_expt_temp, n_h_expt_temp, n_m_expt_temp = count_state(state_expt_temp,N_CG_temp)
        # #     D5_locus_mean_M_level_pending.append(n_m_expt_temp)
            
        # #     # end of loop over expt. Decile5 states (i_expts)
            
        # # # check whether ID has been invalidated and if so move to next segment
        # # if invalidate_current_ID == 1:
        # #     continue
        
        
        # # simulation
        # # set sim. oupt-put variables to zero before start replicates
        # Sim_locus_mean_M_time_pending = []
        # Sim_locus_mean_M_level_pending = []

        # # Sim_UU_gaps_pending = []
        # # Sim_MU_gaps_pending = []
        # # Sim_MM_gaps_pending = []
        # # Sim_XX_gaps_pending = []
        
        # # Sim_UU_nearest_gaps_pending = []
        # # Sim_MU_nearest_gaps_pending = []
        # # Sim_MM_nearest_gaps_pending = []
        # # Sim_XX_nearest_gaps_pending = []

        # for i_reps in range(0,N_reps):
        #     # calculate the seed
        #     current_seed = calc_seed(N_reps, initial_seed, current_segment_no, i_reps)
        #     # seed random number generator
        #     np.random.seed(current_seed)

        #     #print()
        #     #print(current_gene_ID,N_CG,'1/',1./N_CG_density)
        #     #print(state_gene)
        #     #print(CG_positions_gene)
        #     #print(spont_gain_params)
        #     #print(coop_gain_params)
        #     #print(coop_maint_params)
            
        #     # run simulation          
        #     sim_outputs_temp = single_gene_sim_eqbrm_single_time_stats(N_CG, 
        #                             state_gene, CG_positions_gene, spont_gain_params, coop_gain_params, coop_maint_params)



        #     # unpack simulation_outputs
        #     # eqbrm_stats_locus_gaps_temp = sim_outputs_temp[0]
        #     eqbrm_stats_locus_meth_level_temp = sim_outputs_temp[1]
        #     eqbrm_state_output_temp = sim_outputs_temp[2]

        #     # # eqbrm_stats_locus_gaps_temp
        #     # Sim_UU_gaps_temp = eqbrm_stats_locus_gaps_temp[0]
        #     # Sim_MU_gaps_temp = eqbrm_stats_locus_gaps_temp[1]
        #     # Sim_MM_gaps_temp = eqbrm_stats_locus_gaps_temp[2]
        #     # Sim_XX_gaps_temp = eqbrm_stats_locus_gaps_temp[3]
        #     # Sim_UU_nearest_gaps_temp = eqbrm_stats_locus_gaps_temp[4]
        #     # Sim_MU_nearest_gaps_temp = eqbrm_stats_locus_gaps_temp[5]
        #     # Sim_MM_nearest_gaps_temp = eqbrm_stats_locus_gaps_temp[6]
        #     # Sim_XX_nearest_gaps_temp = eqbrm_stats_locus_gaps_temp[7]

        #     # Sim_UU_gaps_pending.extend(Sim_UU_gaps_temp)
        #     # Sim_MU_gaps_pending.extend(Sim_MU_gaps_temp)
        #     # Sim_MM_gaps_pending.extend(Sim_MM_gaps_temp)
        #     # Sim_XX_gaps_pending.extend(Sim_XX_gaps_temp)
        #     # Sim_UU_nearest_gaps_pending.extend(Sim_UU_nearest_gaps_temp)
        #     # Sim_MU_nearest_gaps_pending.extend(Sim_MU_nearest_gaps_temp)
        #     # Sim_MM_nearest_gaps_pending.extend(Sim_MM_nearest_gaps_temp)
        #     # Sim_XX_nearest_gaps_pending.extend(Sim_XX_nearest_gaps_temp)

        #     # eqbrm_stats_locus_meth_level_temp
        #     Sim_locus_mean_M_time_pending.extend(eqbrm_stats_locus_meth_level_temp[0])
        #     Sim_locus_mean_M_level_pending.extend(eqbrm_stats_locus_meth_level_temp[1])
            

        #     # # write out simulated state to file: 
        #     # with open(file_name_output_SimulatedStates_+str(i_reps)+'_'+filename_ending,'a',newline='') as output_file:
        #     #    for i_ in range(len(CG_positions_gene)):
        #     #        # transform number-code back to U or M
        #     #        if eqbrm_state_output_temp[i_] == 0:
        #     #            site_state_temp = 'U'
        #     #        elif eqbrm_state_output_temp[i_] == 2:
        #     #            site_state_temp = 'M'
        #     #        else:
        #     #            site_state_temp = 'oh_no!!!!_what_is_this_site_if_not_U_or_M???'
        #     #        line_writer = csv.writer(output_file,delimiter='\t')
        #     #        line_writer.writerow([current_gene_ID,current_gene_ID[2],CG_positions_gene[i_],CG_positions_gene[i_]+1,site_state_temp])
        #     # # finished write out simulated state to file

            
        #     # end of replicates

        # # calculate histogram counts for gaps-analysis using base pair bins
        # Col0_UU_gaps_counts_pending, bins_dummy = np.histogram(Col0_UU_gaps_pending, bins = N_corrns_bins,range=(-0.5,N_corrns_bins-0.5))
        # Col0_MU_gaps_counts_pending, bins_dummy = np.histogram(Col0_MU_gaps_pending, bins = N_corrns_bins,range=(-0.5,N_corrns_bins-0.5))
        # Col0_MM_gaps_counts_pending, bins_dummy = np.histogram(Col0_MM_gaps_pending, bins = N_corrns_bins,range=(-0.5,N_corrns_bins-0.5))
        # Col0_XX_gaps_counts_pending, bins_dummy = np.histogram(Col0_XX_gaps_pending, bins = N_corrns_bins,range=(-0.5,N_corrns_bins-0.5))

        # Col0_UU_nearest_gaps_counts_pending, bins_dummy = np.histogram(Col0_UU_nearest_gaps_pending, bins = N_corrns_bins,range=(-0.5,N_corrns_bins-0.5))
        # Col0_MU_nearest_gaps_counts_pending, bins_dummy = np.histogram(Col0_MU_nearest_gaps_pending, bins = N_corrns_bins,range=(-0.5,N_corrns_bins-0.5))
        # Col0_MM_nearest_gaps_counts_pending, bins_dummy = np.histogram(Col0_MM_nearest_gaps_pending, bins = N_corrns_bins,range=(-0.5,N_corrns_bins-0.5))
        # Col0_XX_nearest_gaps_counts_pending, bins_dummy = np.histogram(Col0_XX_nearest_gaps_pending, bins = N_corrns_bins,range=(-0.5,N_corrns_bins-0.5))

        # Data_UU_gaps_counts_pending, bins_dummy = np.histogram(Data_UU_gaps_pending, bins = N_corrns_bins,range=(-0.5,N_corrns_bins-0.5))
        # Data_MU_gaps_counts_pending, bins_dummy = np.histogram(Data_MU_gaps_pending, bins = N_corrns_bins,range=(-0.5,N_corrns_bins-0.5))
        # Data_MM_gaps_counts_pending, bins_dummy = np.histogram(Data_MM_gaps_pending, bins = N_corrns_bins,range=(-0.5,N_corrns_bins-0.5))
        # Data_XX_gaps_counts_pending, bins_dummy = np.histogram(Data_XX_gaps_pending, bins = N_corrns_bins,range=(-0.5,N_corrns_bins-0.5))

        # Data_UU_nearest_gaps_counts_pending, bins_dummy = np.histogram(Data_UU_nearest_gaps_pending, bins = N_corrns_bins,range=(-0.5,N_corrns_bins-0.5))
        # Data_MU_nearest_gaps_counts_pending, bins_dummy = np.histogram(Data_MU_nearest_gaps_pending, bins = N_corrns_bins,range=(-0.5,N_corrns_bins-0.5))
        # Data_MM_nearest_gaps_counts_pending, bins_dummy = np.histogram(Data_MM_nearest_gaps_pending, bins = N_corrns_bins,range=(-0.5,N_corrns_bins-0.5))
        # Data_XX_nearest_gaps_counts_pending, bins_dummy = np.histogram(Data_XX_nearest_gaps_pending, bins = N_corrns_bins,range=(-0.5,N_corrns_bins-0.5))

        # D4_UU_gaps_counts_pending, bins_dummy = np.histogram(D4_UU_gaps_pending, bins = N_corrns_bins,range=(-0.5,N_corrns_bins-0.5))
        # D4_MU_gaps_counts_pending, bins_dummy = np.histogram(D4_MU_gaps_pending, bins = N_corrns_bins,range=(-0.5,N_corrns_bins-0.5))
        # D4_MM_gaps_counts_pending, bins_dummy = np.histogram(D4_MM_gaps_pending, bins = N_corrns_bins,range=(-0.5,N_corrns_bins-0.5))
        # D4_XX_gaps_counts_pending, bins_dummy = np.histogram(D4_XX_gaps_pending, bins = N_corrns_bins,range=(-0.5,N_corrns_bins-0.5))

        # D4_UU_nearest_gaps_counts_pending, bins_dummy = np.histogram(D4_UU_nearest_gaps_pending, bins = N_corrns_bins,range=(-0.5,N_corrns_bins-0.5))
        # D4_MU_nearest_gaps_counts_pending, bins_dummy = np.histogram(D4_MU_nearest_gaps_pending, bins = N_corrns_bins,range=(-0.5,N_corrns_bins-0.5))
        # D4_MM_nearest_gaps_counts_pending, bins_dummy = np.histogram(D4_MM_nearest_gaps_pending, bins = N_corrns_bins,range=(-0.5,N_corrns_bins-0.5))
        # D4_XX_nearest_gaps_counts_pending, bins_dummy = np.histogram(D4_XX_nearest_gaps_pending, bins = N_corrns_bins,range=(-0.5,N_corrns_bins-0.5))

        # D5_UU_gaps_counts_pending, bins_dummy = np.histogram(D5_UU_gaps_pending, bins = N_corrns_bins,range=(-0.5,N_corrns_bins-0.5))
        # D5_MU_gaps_counts_pending, bins_dummy = np.histogram(D5_MU_gaps_pending, bins = N_corrns_bins,range=(-0.5,N_corrns_bins-0.5))
        # D5_MM_gaps_counts_pending, bins_dummy = np.histogram(D5_MM_gaps_pending, bins = N_corrns_bins,range=(-0.5,N_corrns_bins-0.5))
        # D5_XX_gaps_counts_pending, bins_dummy = np.histogram(D5_XX_gaps_pending, bins = N_corrns_bins,range=(-0.5,N_corrns_bins-0.5))

        # D5_UU_nearest_gaps_counts_pending, bins_dummy = np.histogram(D5_UU_nearest_gaps_pending, bins = N_corrns_bins,range=(-0.5,N_corrns_bins-0.5))
        # D5_MU_nearest_gaps_counts_pending, bins_dummy = np.histogram(D5_MU_nearest_gaps_pending, bins = N_corrns_bins,range=(-0.5,N_corrns_bins-0.5))
        # D5_MM_nearest_gaps_counts_pending, bins_dummy = np.histogram(D5_MM_nearest_gaps_pending, bins = N_corrns_bins,range=(-0.5,N_corrns_bins-0.5))
        # D5_XX_nearest_gaps_counts_pending, bins_dummy = np.histogram(D5_XX_nearest_gaps_pending, bins = N_corrns_bins,range=(-0.5,N_corrns_bins-0.5))

        # Sim_UU_gaps_counts_pending, bins_dummy = np.histogram(Sim_UU_gaps_pending, bins = N_corrns_bins,range=(-0.5,N_corrns_bins-0.5))
        # Sim_MU_gaps_counts_pending, bins_dummy = np.histogram(Sim_MU_gaps_pending, bins = N_corrns_bins,range=(-0.5,N_corrns_bins-0.5))
        # Sim_MM_gaps_counts_pending, bins_dummy = np.histogram(Sim_MM_gaps_pending, bins = N_corrns_bins,range=(-0.5,N_corrns_bins-0.5))
        # Sim_XX_gaps_counts_pending, bins_dummy = np.histogram(Sim_XX_gaps_pending, bins = N_corrns_bins,range=(-0.5,N_corrns_bins-0.5))

        # Sim_UU_nearest_gaps_counts_pending, bins_dummy = np.histogram(Sim_UU_nearest_gaps_pending, bins = N_corrns_bins,range=(-0.5,N_corrns_bins-0.5))
        # Sim_MU_nearest_gaps_counts_pending, bins_dummy = np.histogram(Sim_MU_nearest_gaps_pending, bins = N_corrns_bins,range=(-0.5,N_corrns_bins-0.5))
        # Sim_MM_nearest_gaps_counts_pending, bins_dummy = np.histogram(Sim_MM_nearest_gaps_pending, bins = N_corrns_bins,range=(-0.5,N_corrns_bins-0.5))
        # Sim_XX_nearest_gaps_counts_pending, bins_dummy = np.histogram(Sim_XX_nearest_gaps_pending, bins = N_corrns_bins,range=(-0.5,N_corrns_bins-0.5))

        # # And pending count values to totals
        # Col0_UU_gaps_counts = np.add(Col0_UU_gaps_counts,Col0_UU_gaps_counts_pending)
        # Col0_MU_gaps_counts = np.add(Col0_MU_gaps_counts,Col0_MU_gaps_counts_pending)
        # Col0_MM_gaps_counts = np.add(Col0_MM_gaps_counts,Col0_MM_gaps_counts_pending)
        # Col0_XX_gaps_counts = np.add(Col0_XX_gaps_counts,Col0_XX_gaps_counts_pending)

        # Col0_UU_nearest_gaps_counts = np.add(Col0_UU_nearest_gaps_counts,Col0_UU_nearest_gaps_counts_pending)
        # Col0_MU_nearest_gaps_counts = np.add(Col0_MU_nearest_gaps_counts,Col0_MU_nearest_gaps_counts_pending)
        # Col0_MM_nearest_gaps_counts = np.add(Col0_MM_nearest_gaps_counts,Col0_MM_nearest_gaps_counts_pending)
        # Col0_XX_nearest_gaps_counts = np.add(Col0_XX_nearest_gaps_counts,Col0_XX_nearest_gaps_counts_pending)

        # Data_UU_gaps_counts = np.add(Data_UU_gaps_counts,Data_UU_gaps_counts_pending)
        # Data_MU_gaps_counts = np.add(Data_MU_gaps_counts,Data_MU_gaps_counts_pending)
        # Data_MM_gaps_counts = np.add(Data_MM_gaps_counts,Data_MM_gaps_counts_pending)
        # Data_XX_gaps_counts = np.add(Data_XX_gaps_counts,Data_XX_gaps_counts_pending)

        # Data_UU_nearest_gaps_counts = np.add(Data_UU_nearest_gaps_counts,Data_UU_nearest_gaps_counts_pending)
        # Data_MU_nearest_gaps_counts = np.add(Data_MU_nearest_gaps_counts,Data_MU_nearest_gaps_counts_pending)
        # Data_MM_nearest_gaps_counts = np.add(Data_MM_nearest_gaps_counts,Data_MM_nearest_gaps_counts_pending)
        # Data_XX_nearest_gaps_counts = np.add(Data_XX_nearest_gaps_counts,Data_XX_nearest_gaps_counts_pending)

        # D4_UU_gaps_counts = np.add(D4_UU_gaps_counts,D4_UU_gaps_counts_pending)
        # D4_MU_gaps_counts = np.add(D4_MU_gaps_counts,D4_MU_gaps_counts_pending)
        # D4_MM_gaps_counts = np.add(D4_MM_gaps_counts,D4_MM_gaps_counts_pending)
        # D4_XX_gaps_counts = np.add(D4_XX_gaps_counts,D4_XX_gaps_counts_pending)

        # D4_UU_nearest_gaps_counts = np.add(D4_UU_nearest_gaps_counts,D4_UU_nearest_gaps_counts_pending)
        # D4_MU_nearest_gaps_counts = np.add(D4_MU_nearest_gaps_counts,D4_MU_nearest_gaps_counts_pending)
        # D4_MM_nearest_gaps_counts = np.add(D4_MM_nearest_gaps_counts,D4_MM_nearest_gaps_counts_pending)
        # D4_XX_nearest_gaps_counts = np.add(D4_XX_nearest_gaps_counts,D4_XX_nearest_gaps_counts_pending)

        # D5_UU_gaps_counts = np.add(D5_UU_gaps_counts,D5_UU_gaps_counts_pending)
        # D5_MU_gaps_counts = np.add(D5_MU_gaps_counts,D5_MU_gaps_counts_pending)
        # D5_MM_gaps_counts = np.add(D5_MM_gaps_counts,D5_MM_gaps_counts_pending)
        # D5_XX_gaps_counts = np.add(D5_XX_gaps_counts,D5_XX_gaps_counts_pending)

        # D5_UU_nearest_gaps_counts = np.add(D5_UU_nearest_gaps_counts,D5_UU_nearest_gaps_counts_pending)
        # D5_MU_nearest_gaps_counts = np.add(D5_MU_nearest_gaps_counts,D5_MU_nearest_gaps_counts_pending)
        # D5_MM_nearest_gaps_counts = np.add(D5_MM_nearest_gaps_counts,D5_MM_nearest_gaps_counts_pending)
        # D5_XX_nearest_gaps_counts = np.add(D5_XX_nearest_gaps_counts,D5_XX_nearest_gaps_counts_pending)

        # Sim_UU_gaps_counts = np.add(Sim_UU_gaps_counts,Sim_UU_gaps_counts_pending)
        # Sim_MU_gaps_counts = np.add(Sim_MU_gaps_counts,Sim_MU_gaps_counts_pending)
        # Sim_MM_gaps_counts = np.add(Sim_MM_gaps_counts,Sim_MM_gaps_counts_pending)
        # Sim_XX_gaps_counts = np.add(Sim_XX_gaps_counts,Sim_XX_gaps_counts_pending)

        # Sim_UU_nearest_gaps_counts = np.add(Sim_UU_nearest_gaps_counts,Sim_UU_nearest_gaps_counts_pending)
        # Sim_MU_nearest_gaps_counts = np.add(Sim_MU_nearest_gaps_counts,Sim_MU_nearest_gaps_counts_pending)
        # Sim_MM_nearest_gaps_counts = np.add(Sim_MM_nearest_gaps_counts,Sim_MM_nearest_gaps_counts_pending)
        # Sim_XX_nearest_gaps_counts = np.add(Sim_XX_nearest_gaps_counts,Sim_XX_nearest_gaps_counts_pending)

        # calculate additional individual locus properties 

        # Sim_mu_methylation_fraction = np.mean(Sim_locus_mean_M_level_pending)
        # Sim_sigma_methylation_fraction = np.std(Sim_locus_mean_M_level_pending)
        # Sim10_mu_methylation_fraction = np.mean(Sim_locus_mean_M_level_pending[:10])
        # Sim10_sigma_methylation_fraction = np.std(Sim_locus_mean_M_level_pending[:10])
        # Dset1_mu_methylation_fraction = np.mean(Dset1_locus_mean_M_level_pending)
        # Dset1_sigma_methylation_fraction = np.std(Dset1_locus_mean_M_level_pending)
        # Dset2_mu_methylation_fraction = np.mean(Dset2_locus_mean_M_level_pending)
        # Dset2_sigma_methylation_fraction = np.std(Dset2_locus_mean_M_level_pending)
        # D5_mu_methylation_fraction = np.mean(D5_locus_mean_M_level_pending)
        # D5_sigma_methylation_fraction = np.std(D5_locus_mean_M_level_pending)
        # DataD4D5_locus_mean_M_level_pending = np.hstack((Data_locus_mean_M_level_pending,
        #                                        Dset2_locus_mean_M_level_pending,
        #                                         D5_locus_mean_M_level_pending))
        # DataD4D5_mu_methylation_fraction = np.mean(DataD4D5_locus_mean_M_level_pending)
        # DataD4D5_sigma_methylation_fraction = np.std(DataD4D5_locus_mean_M_level_pending)

        # Sim_number_0M_frac = np.count_nonzero(np.array(Sim_locus_mean_M_level_pending) == 0.0)
        # Sim10_number_0M_frac = np.count_nonzero(np.array(Sim_locus_mean_M_level_pending[:10]) == 0.0)
        # Dset1_number_0M_frac = np.count_nonzero(np.array(Dset1_locus_mean_M_level_pending) == 0.0)
        #  Dset2_number_0M_frac = np.count_nonzero(np.array(Dset2_locus_mean_M_level_pending) == 0.0)
        # D5_number_0M_frac = np.count_nonzero(np.array(D5_locus_mean_M_level_pending) == 0.0)
        # DataD4D5_number_0M_frac = np.count_nonzero(np.array(DataD4D5_locus_mean_M_level_pending) == 0.0)

        # Dset1_mu_meth_diff = Sim10_mu_methylation_fraction - Dset1_mu_methylation_fraction
        # Dset2_mu_meth_diff = Sim10_mu_methylation_fraction - Dset2_mu_methylation_fraction
        # D5_mu_meth_diff = Sim10_mu_methylation_fraction - D5_mu_methylation_fraction
        # DataD4D5_mu_meth_diff = Sim_mu_methylation_fraction -  DataD4D5_mu_methylation_fraction

        LocusProperties_temp = [current_gene_ID, current_segment_no,N_CG, L_locus, N_CG_density,
                                        #    Sim_mu_methylation_fraction, Sim_sigma_methylation_fraction, 
                                        #    Sim10_mu_methylation_fraction, Sim10_sigma_methylation_fraction,
                                        #    Dset1_mu_methylation_fraction, Dset1_sigma_methylation_fraction, 
                                        #     Dset2_mu_methylation_fraction, Dset2_sigma_methylation_fraction, 
                                        #     D5_mu_methylation_fraction, D5_sigma_methylation_fraction,
                                        #    DataD4D5_mu_methylation_fraction, DataD4D5_sigma_methylation_fraction,
                                            # Sim_number_0M_frac, Sim10_number_0M_frac,
                                            # Dset1_number_0M_frac, 
                                            # Dset2_number_0M_frac, 
                                            # D5_number_0M_frac,
                                            # DataD4D5_number_0M_frac,
                                           Col0_locus_mean_M_level_pending, 
                                           Col0_locus_mean_X_level_pending, 
                                            # Sim_locus_mean_M_level_pending[0], Sim_locus_mean_M_level_pending[1], Sim_locus_mean_M_level_pending[2],
                                            # Dset1_locus_mean_M_level_pending[0], Dset1_locus_mean_M_level_pending[1], Dset1_locus_mean_M_level_pending[2], 
                                            #  Dset2_locus_mean_M_level_pending[0], Dset2_locus_mean_M_level_pending[1], Dset2_locus_mean_M_level_pending[2], 
                                            # D5_locus_mean_M_level_pending[0], D5_locus_mean_M_level_pending[1], D5_locus_mean_M_level_pending[2],
                                            # Dset1_mu_meth_diff, 
                                            # Dset2_mu_meth_diff, D5_mu_meth_diff,
                                            # DataD4D5_mu_meth_diff, 
                                            1./N_CG_density, mean_h2az, max_h2az
                                           ]

        temp_output_list = LocusProperties_temp
        temp_output_file = file_name_output_LocusProperties
        with open(temp_output_file,'a',newline='') as output_file:
            line_writer = csv.writer(output_file,delimiter='\t')
            line_writer.writerow(temp_output_list)

        # AllRepsMethFracs_temp = [current_gene_ID, current_segment_no,
        #                         Sim_locus_mean_M_level_pending[0], Sim_locus_mean_M_level_pending[1], Sim_locus_mean_M_level_pending[2],
        #                          Sim_locus_mean_M_level_pending[3], Sim_locus_mean_M_level_pending[4], Sim_locus_mean_M_level_pending[5],
        #                          Sim_locus_mean_M_level_pending[6], Sim_locus_mean_M_level_pending[7], Sim_locus_mean_M_level_pending[8],
        #                          Sim_locus_mean_M_level_pending[9], 
        #                          Sim_locus_mean_M_level_pending[10], Sim_locus_mean_M_level_pending[11], Sim_locus_mean_M_level_pending[12],
        #                          Sim_locus_mean_M_level_pending[13], Sim_locus_mean_M_level_pending[14], Sim_locus_mean_M_level_pending[15],
        #                          Sim_locus_mean_M_level_pending[16], Sim_locus_mean_M_level_pending[17], Sim_locus_mean_M_level_pending[18],
        #                          Sim_locus_mean_M_level_pending[19], 
        #                          Sim_locus_mean_M_level_pending[20], Sim_locus_mean_M_level_pending[21], Sim_locus_mean_M_level_pending[22],
        #                          Sim_locus_mean_M_level_pending[23], Sim_locus_mean_M_level_pending[24], Sim_locus_mean_M_level_pending[25],
        #                          Sim_locus_mean_M_level_pending[26], Sim_locus_mean_M_level_pending[27], Sim_locus_mean_M_level_pending[28],
        #                          Sim_locus_mean_M_level_pending[29], 
        #                          Dset1_locus_mean_M_level_pending[0], Dset1_locus_mean_M_level_pending[1], Dset1_locus_mean_M_level_pending[2],
        #                          Dset1_locus_mean_M_level_pending[3], Dset1_locus_mean_M_level_pending[4], Dset1_locus_mean_M_level_pending[5],
        #                          Dset1_locus_mean_M_level_pending[6], Dset1_locus_mean_M_level_pending[7], Dset1_locus_mean_M_level_pending[8],
        #                          Dset1_locus_mean_M_level_pending[9],
        #                          Dset1_locus_mean_M_level_pending[10], Dset1_locus_mean_M_level_pending[11], Dset1_locus_mean_M_level_pending[12],
        #                          Dset1_locus_mean_M_level_pending[13], Dset1_locus_mean_M_level_pending[14], Dset1_locus_mean_M_level_pending[15],
        #                          Dset1_locus_mean_M_level_pending[16], Dset1_locus_mean_M_level_pending[17], Dset1_locus_mean_M_level_pending[18],
        #                          Dset1_locus_mean_M_level_pending[19],
        #                          Dset1_locus_mean_M_level_pending[20], Dset1_locus_mean_M_level_pending[21], Dset1_locus_mean_M_level_pending[22],
        #                          Dset1_locus_mean_M_level_pending[23], Dset1_locus_mean_M_level_pending[24], Dset1_locus_mean_M_level_pending[25],
        #                          Dset1_locus_mean_M_level_pending[26], Dset1_locus_mean_M_level_pending[27], Dset1_locus_mean_M_level_pending[28],
        #                          Dset1_locus_mean_M_level_pending[29] 
        #                          ]

        # AllRepsXFracs_temp = [current_gene_ID, current_segment_no,
        #                          Dset1_locus_mean_X_level_pending[0], Dset1_locus_mean_X_level_pending[1], Dset1_locus_mean_X_level_pending[2],
        #                          Dset1_locus_mean_X_level_pending[3], Dset1_locus_mean_X_level_pending[4], Dset1_locus_mean_X_level_pending[5],
        #                          Dset1_locus_mean_X_level_pending[6], Dset1_locus_mean_X_level_pending[7], Dset1_locus_mean_X_level_pending[8],
        #                          Dset1_locus_mean_X_level_pending[9],
        #                          Dset1_locus_mean_X_level_pending[10], Dset1_locus_mean_X_level_pending[11], Dset1_locus_mean_X_level_pending[12],
        #                          Dset1_locus_mean_X_level_pending[13], Dset1_locus_mean_X_level_pending[14], Dset1_locus_mean_X_level_pending[15],
        #                          Dset1_locus_mean_X_level_pending[16], Dset1_locus_mean_X_level_pending[17], Dset1_locus_mean_X_level_pending[18],
        #                          Dset1_locus_mean_X_level_pending[19],
        #                          Dset1_locus_mean_X_level_pending[20], Dset1_locus_mean_X_level_pending[21], Dset1_locus_mean_X_level_pending[22],
        #                          Dset1_locus_mean_X_level_pending[23], Dset1_locus_mean_X_level_pending[24], Dset1_locus_mean_X_level_pending[25],
        #                          Dset1_locus_mean_X_level_pending[26], Dset1_locus_mean_X_level_pending[27], Dset1_locus_mean_X_level_pending[28],
        #                          Dset1_locus_mean_X_level_pending[29] 
        #                          ]

        # temp_output_list = AllRepsMethFracs_temp
        # temp_output_file = file_name_output_AllRepsMethFracs
        # with open(temp_output_file,'a',newline='') as output_file:
        #     line_writer = csv.writer(output_file,delimiter='\t')
        #     line_writer.writerow(temp_output_list)

        # temp_output_list = AllRepsXFracs_temp
        # temp_output_file = file_name_output_AllRepsXFracs
        # with open(temp_output_file,'a',newline='') as output_file:
        #     line_writer = csv.writer(output_file,delimiter='\t')
        #     line_writer.writerow(temp_output_list)
        
        i_seg_total +=1
        
        # end of loop over segments
    
    #if (i_gene%10 == 0):
    #    timer_end = timer()
    #    #print('finished', i_gene, i_reps,'time', timer_end - timer_start, flush=True)
    #    
    #    with open(file_name_output+'CodeProgress'+'.tsv','a',newline='') as output_file:##~~++~~##
    #        line_writer = csv.writer(output_file,delimiter='\t')
    #        line_writer.writerow(['finished', i_gene, i_reps,'time', timer_end - timer_start])
    #    
    #    timer_start = timer()
    
    # end of current gene loop

# # write out Pair correlations Histograms to file
# temp_output_list = Sim_MM_gaps_counts
# temp_output_file = file_name_output_PairSeparationsMM_Sim
# with open(temp_output_file,'a',newline='') as output_file:
#    line_writer = csv.writer(output_file,delimiter='\t')
#    line_writer.writerow(temp_output_list)

# temp_output_list = Sim_MU_gaps_counts
# temp_output_file = file_name_output_PairSeparationsMU_Sim
# with open(temp_output_file,'a',newline='') as output_file:
#    line_writer = csv.writer(output_file,delimiter='\t')
#    line_writer.writerow(temp_output_list)

# temp_output_list = Sim_UU_gaps_counts
# temp_output_file = file_name_output_PairSeparationsUU_Sim
# with open(temp_output_file,'a',newline='') as output_file:
#    line_writer = csv.writer(output_file,delimiter='\t')
#    line_writer.writerow(temp_output_list)

# temp_output_list = Sim_XX_gaps_counts
# temp_output_file = file_name_output_PairSeparationsXX_Sim
# with open(temp_output_file,'a',newline='') as output_file:
#    line_writer = csv.writer(output_file,delimiter='\t')
#    line_writer.writerow(temp_output_list)

# temp_output_list = Sim_MM_nearest_gaps_counts
# temp_output_file = file_name_output_NearestPairSeparationsMM_Sim
# with open(temp_output_file,'a',newline='') as output_file:
#    line_writer = csv.writer(output_file,delimiter='\t')
#    line_writer.writerow(temp_output_list)

# temp_output_list = Sim_MU_nearest_gaps_counts
# temp_output_file = file_name_output_NearestPairSeparationsMU_Sim
# with open(temp_output_file,'a',newline='') as output_file:
#    line_writer = csv.writer(output_file,delimiter='\t')
#    line_writer.writerow(temp_output_list)

# temp_output_list = Sim_UU_nearest_gaps_counts
# temp_output_file = file_name_output_NearestPairSeparationsUU_Sim
# with open(temp_output_file,'a',newline='') as output_file:
#    line_writer = csv.writer(output_file,delimiter='\t')
#    line_writer.writerow(temp_output_list)

# temp_output_list = Sim_XX_nearest_gaps_counts
# temp_output_file = file_name_output_NearestPairSeparationsXX_Sim
# with open(temp_output_file,'a',newline='') as output_file:
#    line_writer = csv.writer(output_file,delimiter='\t')
#    line_writer.writerow(temp_output_list)


# temp_output_list = Col0_MM_gaps_counts
# temp_output_file = file_name_output_PairSeparationsMM_Col0
# with open(temp_output_file,'a',newline='') as output_file:
#    line_writer = csv.writer(output_file,delimiter='\t')
#    line_writer.writerow(temp_output_list)

# temp_output_list = Col0_MU_gaps_counts
# temp_output_file = file_name_output_PairSeparationsMU_Col0
# with open(temp_output_file,'a',newline='') as output_file:
#    line_writer = csv.writer(output_file,delimiter='\t')
#    line_writer.writerow(temp_output_list)

# temp_output_list = Col0_UU_gaps_counts
# temp_output_file = file_name_output_PairSeparationsUU_Col0
# with open(temp_output_file,'a',newline='') as output_file:
#    line_writer = csv.writer(output_file,delimiter='\t')
#    line_writer.writerow(temp_output_list)

# temp_output_list = Col0_XX_gaps_counts
# temp_output_file = file_name_output_PairSeparationsXX_Col0
# with open(temp_output_file,'a',newline='') as output_file:
#    line_writer = csv.writer(output_file,delimiter='\t')
#    line_writer.writerow(temp_output_list)

# temp_output_list = Col0_MM_nearest_gaps_counts
# temp_output_file = file_name_output_NearestPairSeparationsMM_Col0
# with open(temp_output_file,'a',newline='') as output_file:
#    line_writer = csv.writer(output_file,delimiter='\t')
#    line_writer.writerow(temp_output_list)

# temp_output_list = Col0_MU_nearest_gaps_counts
# temp_output_file = file_name_output_NearestPairSeparationsMU_Col0
# with open(temp_output_file,'a',newline='') as output_file:
#    line_writer = csv.writer(output_file,delimiter='\t')
#    line_writer.writerow(temp_output_list)

# temp_output_list = Col0_UU_nearest_gaps_counts
# temp_output_file = file_name_output_NearestPairSeparationsUU_Col0
# with open(temp_output_file,'a',newline='') as output_file:
#    line_writer = csv.writer(output_file,delimiter='\t')
#    line_writer.writerow(temp_output_list)

# temp_output_list = Col0_XX_nearest_gaps_counts
# temp_output_file = file_name_output_NearestPairSeparationsXX_Col0
# with open(temp_output_file,'a',newline='') as output_file:
#    line_writer = csv.writer(output_file,delimiter='\t')
#    line_writer.writerow(temp_output_list)


# temp_output_list = Data_MM_gaps_counts
# temp_output_file = file_name_output_PairSeparationsMM_Data
# with open(temp_output_file,'a',newline='') as output_file:
#    line_writer = csv.writer(output_file,delimiter='\t')
#    line_writer.writerow(temp_output_list)

# temp_output_list = Data_MU_gaps_counts
# temp_output_file = file_name_output_PairSeparationsMU_Data
# with open(temp_output_file,'a',newline='') as output_file:
#    line_writer = csv.writer(output_file,delimiter='\t')
#    line_writer.writerow(temp_output_list)

# temp_output_list = Data_UU_gaps_counts
# temp_output_file = file_name_output_PairSeparationsUU_Data
# with open(temp_output_file,'a',newline='') as output_file:
#    line_writer = csv.writer(output_file,delimiter='\t')
#    line_writer.writerow(temp_output_list)

# temp_output_list = Data_XX_gaps_counts
# temp_output_file = file_name_output_PairSeparationsXX_Data
# with open(temp_output_file,'a',newline='') as output_file:
#    line_writer = csv.writer(output_file,delimiter='\t')
#    line_writer.writerow(temp_output_list)

# temp_output_list = Data_MM_nearest_gaps_counts
# temp_output_file = file_name_output_NearestPairSeparationsMM_Data
# with open(temp_output_file,'a',newline='') as output_file:
#    line_writer = csv.writer(output_file,delimiter='\t')
#    line_writer.writerow(temp_output_list)

# temp_output_list = Data_MU_nearest_gaps_counts
# temp_output_file = file_name_output_NearestPairSeparationsMU_Data
# with open(temp_output_file,'a',newline='') as output_file:
#    line_writer = csv.writer(output_file,delimiter='\t')
#    line_writer.writerow(temp_output_list)

# temp_output_list = Data_UU_nearest_gaps_counts
# temp_output_file = file_name_output_NearestPairSeparationsUU_Data
# with open(temp_output_file,'a',newline='') as output_file:
#    line_writer = csv.writer(output_file,delimiter='\t')
#    line_writer.writerow(temp_output_list)

# temp_output_list = Data_XX_nearest_gaps_counts
# temp_output_file = file_name_output_NearestPairSeparationsXX_Data
# with open(temp_output_file,'a',newline='') as output_file:
#    line_writer = csv.writer(output_file,delimiter='\t')
#    line_writer.writerow(temp_output_list)


# temp_output_list = D4_MM_gaps_counts
# temp_output_file = file_name_output_PairSeparationsMM_D4
# with open(temp_output_file,'a',newline='') as output_file:
#    line_writer = csv.writer(output_file,delimiter='\t')
#    line_writer.writerow(temp_output_list)

# temp_output_list = D4_MU_gaps_counts
# temp_output_file = file_name_output_PairSeparationsMU_D4
# with open(temp_output_file,'a',newline='') as output_file:
#    line_writer = csv.writer(output_file,delimiter='\t')
#    line_writer.writerow(temp_output_list)

# temp_output_list = D4_UU_gaps_counts
# temp_output_file = file_name_output_PairSeparationsUU_D4
# with open(temp_output_file,'a',newline='') as output_file:
#    line_writer = csv.writer(output_file,delimiter='\t')
#    line_writer.writerow(temp_output_list)

# temp_output_list = D4_XX_gaps_counts
# temp_output_file = file_name_output_PairSeparationsXX_D4
# with open(temp_output_file,'a',newline='') as output_file:
#    line_writer = csv.writer(output_file,delimiter='\t')
#    line_writer.writerow(temp_output_list)

# temp_output_list = D4_MM_nearest_gaps_counts
# temp_output_file = file_name_output_NearestPairSeparationsMM_D4
# with open(temp_output_file,'a',newline='') as output_file:
#    line_writer = csv.writer(output_file,delimiter='\t')
#    line_writer.writerow(temp_output_list)

# temp_output_list = D4_MU_nearest_gaps_counts
# temp_output_file = file_name_output_NearestPairSeparationsMU_D4
# with open(temp_output_file,'a',newline='') as output_file:
#    line_writer = csv.writer(output_file,delimiter='\t')
#    line_writer.writerow(temp_output_list)

# temp_output_list = D4_UU_nearest_gaps_counts
# temp_output_file = file_name_output_NearestPairSeparationsUU_D4
# with open(temp_output_file,'a',newline='') as output_file:
#    line_writer = csv.writer(output_file,delimiter='\t')
#    line_writer.writerow(temp_output_list)

# temp_output_list = D4_XX_nearest_gaps_counts
# temp_output_file = file_name_output_NearestPairSeparationsXX_D4
# with open(temp_output_file,'a',newline='') as output_file:
#    line_writer = csv.writer(output_file,delimiter='\t')
#    line_writer.writerow(temp_output_list)


# temp_output_list = D5_MM_gaps_counts
# temp_output_file = file_name_output_PairSeparationsMM_D5
# with open(temp_output_file,'a',newline='') as output_file:
#    line_writer = csv.writer(output_file,delimiter='\t')
#    line_writer.writerow(temp_output_list)

# temp_output_list = D5_MU_gaps_counts
# temp_output_file = file_name_output_PairSeparationsMU_D5
# with open(temp_output_file,'a',newline='') as output_file:
#    line_writer = csv.writer(output_file,delimiter='\t')
#    line_writer.writerow(temp_output_list)

# temp_output_list = D5_UU_gaps_counts
# temp_output_file = file_name_output_PairSeparationsUU_D5
# with open(temp_output_file,'a',newline='') as output_file:
#    line_writer = csv.writer(output_file,delimiter='\t')
#    line_writer.writerow(temp_output_list)

# temp_output_list = D5_XX_gaps_counts
# temp_output_file = file_name_output_PairSeparationsXX_D5
# with open(temp_output_file,'a',newline='') as output_file:
#    line_writer = csv.writer(output_file,delimiter='\t')
#    line_writer.writerow(temp_output_list)

# temp_output_list = D5_MM_nearest_gaps_counts
# temp_output_file = file_name_output_NearestPairSeparationsMM_D5
# with open(temp_output_file,'a',newline='') as output_file:
#    line_writer = csv.writer(output_file,delimiter='\t')
#    line_writer.writerow(temp_output_list)

# temp_output_list = D5_MU_nearest_gaps_counts
# temp_output_file = file_name_output_NearestPairSeparationsMU_D5
# with open(temp_output_file,'a',newline='') as output_file:
#    line_writer = csv.writer(output_file,delimiter='\t')
#    line_writer.writerow(temp_output_list)

# temp_output_list = D5_UU_nearest_gaps_counts
# temp_output_file = file_name_output_NearestPairSeparationsUU_D5
# with open(temp_output_file,'a',newline='') as output_file:
#    line_writer = csv.writer(output_file,delimiter='\t')
#    line_writer.writerow(temp_output_list)

# temp_output_list = D5_XX_nearest_gaps_counts
# temp_output_file = file_name_output_NearestPairSeparationsXX_D5
# with open(temp_output_file,'a',newline='') as output_file:
#    line_writer = csv.writer(output_file,delimiter='\t')
#    line_writer.writerow(temp_output_list)



timer_end_total = timer()

with open(file_name_output_CodeProgress,'a',newline='') as output_file:
    line_writer = csv.writer(output_file,delimiter='\t')
    line_writer.writerow([i_batch, i_gene+1, (timer_end_total - timer_start_total)/60.])


#print('end')
#print('time', (timer_end_total - timer_start_total)/60., 'mins')
