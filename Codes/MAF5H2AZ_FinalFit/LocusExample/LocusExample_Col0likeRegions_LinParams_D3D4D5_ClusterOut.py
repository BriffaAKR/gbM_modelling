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
#import argparse
import sys

import input_params as params # import my params from separate python file

current_gene_ID = params.current_gene_ID

# define paramter coding for filenames
filename_params_code = params.filename_params_code

## define annoation files
InputFiles_folder = params.InputFiles_folder
file_in_annotation = params.file_in_annotation
file_in_anno_filt_1 = params.file_in_anno_filt_1
file_in_anno_filt_2 = params.file_in_anno_filt_2
file_in_exclude_filt_1 = params.file_in_exclude_filt_1

# build up input path string
path_string = os.path.join( os.getcwd(), '..')
path_string = os.path.join( path_string, '..')
path_string = os.path.join( path_string, '..')
#path_string = os.path.join( path_string, '..')
path_string = os.path.join( path_string, InputFiles_folder)

annotaion_file = os.path.join( path_string, file_in_annotation )
filt_1_file = os.path.join( path_string, file_in_anno_filt_1)
filt_2_file = os.path.join( path_string, file_in_anno_filt_2) 
exclude_filt_1_file = os.path.join( path_string, file_in_exclude_filt_1)
file_in_Combined_Data_per_genes_meth_statuses = os.path.join( path_string, 'SeqRun9All_CG_2021-01-18_WTconsensus_meth_status_allgenes_withID.tsv' )

# make paths for output files
locus_type = params.locus_type
TE_filt = params.TE_filt
filename_start = params.filename_start
filename_ending = '.tsv'

file_name_output_TEST = os.path.join('Output_files', filename_start + 'TEST'+filename_ending)
file_name_output_MethLevel = os.path.join('Output_files', filename_start + 'MethLevel'+filename_ending)
file_name_output_MethStateTime = os.path.join('Output_files', filename_start + 'MethStateTime'+filename_ending)
file_name_output_MethStateMatrix = os.path.join('Output_files', filename_start + 'MethStateMatrix'+filename_ending)
file_name_output_ExptStateCol0 = os.path.join('Output_files', filename_start + 'ExptStateCol0'+filename_ending)
file_name_output_CGpositions = os.path.join('Output_files', filename_start + 'CGpositions'+filename_ending)
file_name_output_ExptStatesD3 = os.path.join('Output_files', filename_start + 'MethStateExptStatesD3'+filename_ending)
file_name_output_ExptStatesD4 = os.path.join('Output_files', filename_start + 'MethStateExptStatesD4'+filename_ending)
file_name_output_ExptStatesD5 = os.path.join('Output_files', filename_start + 'MethStateExptStatesD5'+filename_ending)
file_name_output_ExptStatesD3D4D5 = os.path.join('Output_files', filename_start + 'MethStateExptStatesD3D4D5'+filename_ending)
file_name_output_SRXcodesD3 = os.path.join('Output_files', filename_start + 'SRXcodesD3'+filename_ending)
file_name_output_SRXcodesD4 = os.path.join('Output_files', filename_start + 'SRXcodesD4'+filename_ending)
file_name_output_SRXcodesD5 = os.path.join('Output_files', filename_start + 'SRXcodesD5'+filename_ending)
file_name_output_Skipped_IDs = os.path.join('Output_files', filename_start + 'Skipped_IDs'+filename_ending)
file_name_output_CodeProgress = os.path.join('Output_files', filename_start + 'CodeProgress'+filename_ending)

# define column headings for LocusProperties file
MethLevel_headings = ['Time', 'M_Frac','U_Frac']
# define column headings for Skipped_IDs files
Skipped_IDs_headings = ['gene_ID', 'N_CG', 'CG_desnity']
# define column headings for TEST files
TEST_headings = ['segment_no', 'gene_ID']
# define column headings for CodeProgress
CodeProgress_headings = ['batch', 'N_gene', 'time_mins']

# write out column heading
with open(file_name_output_TEST,'w',newline='') as output_file:
    line_writer = csv.writer(output_file,delimiter='\t')
    line_writer.writerow(TEST_headings)
with open(file_name_output_MethLevel,'w',newline='') as output_file:
    line_writer = csv.writer(output_file,delimiter='\t')
    line_writer.writerow(MethLevel_headings)
with open(file_name_output_Skipped_IDs,'w',newline='') as output_file:
    line_writer = csv.writer(output_file,delimiter='\t')
    line_writer.writerow(Skipped_IDs_headings)
with open(file_name_output_CodeProgress,'w',newline='') as output_file:
    line_writer = csv.writer(output_file,delimiter='\t')
    line_writer.writerow(CodeProgress_headings)
#create empty files for rest....
open(file_name_output_MethStateTime,'w').close
open(file_name_output_MethStateMatrix,'w').close
open(file_name_output_ExptStateCol0,'w').close
open(file_name_output_CGpositions,'w').close
open(file_name_output_ExptStatesD3,'w').close
open(file_name_output_ExptStatesD4,'w').close
open(file_name_output_ExptStatesD5,'w').close
open(file_name_output_ExptStatesD3D4D5,'w').close
open(file_name_output_SRXcodesD3,'w').close
open(file_name_output_SRXcodesD4,'w').close
open(file_name_output_SRXcodesD5,'w').close


# load in required segments to model from annotaion file: makes dataframe: AllLoci_df
# CG site at 'start' to be included: x >= start
# 'end' is past last CG site: x < end

# create lazy reader object
temp_input_df = dd.read_csv(annotaion_file,sep='\t',
                            usecols=['segment_no', 'gene_ID','start','end', 'mean_h2az', 'max_h2az'],
                            dtype={'segment_no':'int64', 'gene_ID':'str', 'start':'int64', 'end':'int64', 
                            'mean_h2az':'float', 'max_h2az':'float'})
# filter specific ID for this example:
# temp_input_df = temp_input_df[temp_input_df['segment_no'].isin(segment_no_list)  ]
temp_input_df = temp_input_df[temp_input_df['gene_ID'].isin([current_gene_ID])  ]


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

# load in Decile3 m1001 States
Decile_number = 3
ExptState_Decile3_df_list = []
ExptState_Decile3_SRXcodes_list = []
# build up path
temp_path = os.path.join(path_string, 'm1001_states','Decile'+str(Decile_number)+'_3')

#print(temp_path)
#print()
for filename in os.listdir(temp_path):
    #print(filename)
    # isolate SRX code from filename
    temp_SRX_code = filename[40:]
    temp_SRX_code = temp_SRX_code[:-6]
    ExptState_Decile3_SRXcodes_list.append(temp_SRX_code)
    # create lazy reader object
    temp_input_df = dd.read_csv(os.path.join(temp_path,filename),sep='\t',
                                usecols=['gene_ID','Start','Status'],
                               dtype={'gene_ID': 'str', 'Start': 'int64', 'Status': 'str'})
    # define filtering logic
    temp_input_df = temp_input_df[temp_input_df['gene_ID'].isin(input_gene_ID_list)  ] 
    temp_input_df = temp_input_df.set_index('gene_ID')
    # apply filtering logic and convert to pandas dataframe
    temp_input_df = temp_input_df.compute()
    ExptState_Decile3_df_list.append(temp_input_df.copy())
N_ExptState_Decile3 = len(ExptState_Decile3_df_list)
#print(N_ExptState_Decile3)
#print(ExptState_Decile3_df_list[0].shape,ExptState_Decile3_df_list[1].shape)

# load in Decile4 m1001 States
Decile_number = 4
ExptState_Decile4_df_list = []
ExptState_Decile4_SRXcodes_list = []
# build up path
temp_path = os.path.join(path_string, 'm1001_states','Decile'+str(Decile_number)+'_3')

#print(temp_path)
#print()
for filename in os.listdir(temp_path):
    # isolate SRX code from filename
    temp_SRX_code = filename[40:]
    temp_SRX_code = temp_SRX_code[:-6]
    ExptState_Decile4_SRXcodes_list.append(temp_SRX_code)
    #print(filename)
    # create lazy reader object
    temp_input_df = dd.read_csv(os.path.join(temp_path,filename),sep='\t',
                                usecols=['gene_ID','Start','Status'],
                               dtype={'gene_ID': 'str', 'Start': 'int64', 'Status': 'str'})
    # define filtering logic
    temp_input_df = temp_input_df[temp_input_df['gene_ID'].isin(input_gene_ID_list)  ] 
    temp_input_df = temp_input_df.set_index('gene_ID')
    # apply filtering logic and convert to pandas dataframe
    temp_input_df = temp_input_df.compute()
    ExptState_Decile4_df_list.append(temp_input_df.copy())
N_ExptState_Decile4 = len(ExptState_Decile4_df_list)
#print(N_ExptState_Decile4)
#print(ExptState_Decile4_df_list[0].shape,ExptState_Decile4_df_list[1].shape)

# load in Decile5 m1001 States
Decile_number = 5
ExptState_Decile5_df_list = []
ExptState_Decile5_SRXcodes_list = []
# build up path
temp_path = os.path.join(path_string, 'm1001_states','Decile'+str(Decile_number)+'_3')

#print(temp_path)
#print()
for filename in os.listdir(temp_path):
    # isolate SRX code from filename
    temp_SRX_code = filename[40:]
    temp_SRX_code = temp_SRX_code[:-6]
    ExptState_Decile5_SRXcodes_list.append(temp_SRX_code)
    #print(filename)
    # create lazy reader object
    temp_input_df = dd.read_csv(os.path.join(temp_path,filename),sep='\t',
                                usecols=['gene_ID','Start','Status'],
                               dtype={'gene_ID': 'str', 'Start': 'int64', 'Status': 'str'})
    # define filtering logic
    temp_input_df = temp_input_df[temp_input_df['gene_ID'].isin(input_gene_ID_list)  ] 
    temp_input_df = temp_input_df.set_index('gene_ID')
    # apply filtering logic and convert to pandas dataframe
    temp_input_df = temp_input_df.compute()
    ExptState_Decile5_df_list.append(temp_input_df.copy())
N_ExptState_Decile5 = len(ExptState_Decile5_df_list)
#print(N_ExptState_Decile5)
#print(ExptState_Decile5_df_list[0].shape,ExptState_Decile5_df_list[1].shape)

print('D3')
print(ExptState_Decile3_SRXcodes_list)
print()
print('D4')
print(ExptState_Decile4_SRXcodes_list)
print()
print('D5')
print(ExptState_Decile5_SRXcodes_list)
print()

# finish reading in 1001 methylome files


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

    elif initial_state_choice_ == '50M':
        state_sim_input_ = []
        for i_ in range(N_CG_):
            if np.random.rand(1) <= 0.5:
                state_sim_input_.append(2)
            else:
                state_sim_input_.append(0)
        print('initial value', np.sum(state_sim_input_)/(2.*N_CG_))
    else:
        print('alsaasa')

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
#        print(N_,len(state_in_), len(CG_positions_))
#        print(i_,i_+1)
#        print(state_in_[i_],state_in_[i_+1])
#          print()
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



def single_gene_sim_transients(N_gen_cc_, N_CG_, state_, CG_positions_,
                       spont_gain_params_, coop_gain_params_, coop_maint_params_):

    
    # N_gen_cc_ = number of cell cycles to simulate for

    # write out state data at every update

    
    
    
    state_time_ = [] # record states at specified timepoints
    state_matrix_ = []
    
    
    dat_time_ = [] # for simulated methylation levels
    dat_n_u_ = []
    dat_n_m_ = []

    # initilise system
    time_ = 0.0
    
    # analyse state
    n_u_, n_h_, n_m_ = count_state(state_,N_CG_)

    # update data lists
    dat_time_.append(time_)
    dat_n_u_.append(n_u_)
    dat_n_m_.append(n_m_)
    
    state_time_.append(0)
    state_matrix_ = np.array(state_)

    
    # calculate initial propensities
    coop_gains_propensities_array_, spont_gains_propensities_vector_, \
            coop_maint_propensities_array_, total_propensities_vector_ = \
                initial_propensities(N_CG_, state_, CG_positions_, spont_gain_params_, 
                                         coop_gain_params_, coop_maint_params_)
    
    
    # churn through the simulation
    # total simulation time = no. of gens * no. of cell cycles per gen.
    while True:
        
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
        if (time_ + delta_t_) > N_gen_cc_:
            # record final state at N_gen_cc_ by duplicating last state before updated
            state_time_.append(N_gen_cc_)
            state_matrix_ = np.vstack((state_matrix_,state_matrix_[-1,:])) 
            dat_time_.append(N_gen_cc_)
            dat_n_u_.append(dat_n_u_[-1])
            dat_n_m_.append(dat_n_m_[-1])
            break
            
        # if haven't overshot simulation end will not have broken loop, so, complete the timestep updates

        # cal. which site will change (specified by change_index_)
        propensities_cut_off_ = random_pair_[1]*cumulative_propensities_[-1]

        change_index_ = ordered_site_list_[ np.searchsorted(cumulative_propensities_, propensities_cut_off_) ]


        # update the time
        time_ += delta_t_

        # first record existing state at this timepoint
        # write out state data        
        state_time_.append(time_)
        state_matrix_ = np.vstack((state_matrix_,state_))            
        # end of 'write out state data'
            
        # write out mean properties of state
        # whole gene stats.
        # analyse state
        n_u_, n_h_, n_m_ = count_state(state_,N_CG_)
        # update data lists
        dat_time_.append(time_)
        dat_n_u_.append(n_u_)
        dat_n_m_.append(n_m_)

        
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
        
        
        
        # write out state data        
        state_time_.append(time_)
        state_matrix_ = np.vstack((state_matrix_,state_))            
        # end of 'write out state data'
            
        # write out mean properties of state
        # whole gene stats.
        # analyse state
        n_u_, n_h_, n_m_ = count_state(state_,N_CG_)
        # update data lists
        dat_time_.append(time_)
        dat_n_u_.append(n_u_)
        dat_n_m_.append(n_m_)


    
    # package sim_outputs
    # 1. 'transients_state'
    # 2. 'transients_meth_level'

    
    transients_state_ = [dat_time_, dat_n_u_, dat_n_m_]
    transients_meth_level_ = [state_ ,  state_time_, state_matrix_]

    
    sim_outputs_ = [transients_state_, transients_meth_level_]
            
    return sim_outputs_




# function to update seed deterministically for each simulation replicate
def calc_seed(N_reps_, initial_seed_, current_segment_no_, i_reps_):
    current_seed_ = initial_seed_ + (current_segment_no_-1)*N_reps_ + i_reps_
    return current_seed_



# main code body!!
# loop to simulate using: single_gene_sim_eqbrm_stats

# define parameters

initial_seed = params.initial_seed

# possilbe values for P_choice_
# # 'UseU' 'Excl'

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
# initialise validatation variable
invalidate_current_ID = 0

   
# find how many segments exist for this gene
current_segments_df = AllLoci_df.loc[[current_gene_ID]]
current_number_of_segments = len(current_segments_df) # assume that this will be 1
# check only a single locus!
if current_number_of_segments > 1:
    print('error multpiple segments',current_number_of_segments)

print(current_segments_df)
current_segment_no = current_segments_df.loc[current_gene_ID,'segment_no']
mean_h2az = current_segments_df.loc[current_gene_ID,'mean_h2az']
max_h2az = current_segments_df.loc[current_gene_ID,'max_h2az']
# instead of looping over segments:
i_seg_current_gene = 0

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
    sys.exit()
if N_CG >= N_CG_max:
    invalidate_current_ID = 1
    temp_output_list = [current_gene_ID, N_CG, np.nan]
    temp_output_file = file_name_output_Skipped_IDs
    with open(temp_output_file,'a',newline='') as output_file:
        line_writer = csv.writer(output_file,delimiter='\t')
        line_writer.writerow(temp_output_list)
    sys.exit()
    
N_CG_density = (N_CG-1)/(CG_positions_gene[-1]-CG_positions_gene[0])
L_locus = CG_positions_gene[-1]-CG_positions_gene[0]
    
if N_CG_density < N_CG_density_min:
    invalidate_current_ID = 1
    temp_output_list = [current_gene_ID, N_CG, N_CG_density]
    temp_output_file = file_name_output_Skipped_IDs
    with open(temp_output_file,'a',newline='') as output_file:
        line_writer = csv.writer(output_file,delimiter='\t')
        line_writer.writerow(temp_output_list)
    sys.exit()
if N_CG_density >= N_CG_density_max:
    invalidate_current_ID = 1
    temp_output_list = [current_gene_ID, N_CG, N_CG_density]
    temp_output_file = file_name_output_Skipped_IDs
    with open(temp_output_file,'a',newline='') as output_file:
        line_writer = csv.writer(output_file,delimiter='\t')
        line_writer.writerow(temp_output_list)
    sys.exit()
    
    
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
    
    
#find mean meth. level for experimental state
n_u_expt_temp, n_h_expt_temp, n_m_expt_temp = count_state(state_expt,N_CG)
Col0_locus_mean_M_level_temp = n_m_expt_temp

# no loop here so set pending values too... 
Col0_locus_mean_M_level_pending = Col0_locus_mean_M_level_temp

# write out Col0 state
temp_output_list = state_expt
temp_output_file = file_name_output_ExptStateCol0
with open(temp_output_file,'a',newline='') as output_file:
    line_writer = csv.writer(output_file,delimiter='\t')
    line_writer.writerow(temp_output_list)
# write out CGpositions state
temp_output_list = CG_positions_gene
temp_output_file = file_name_output_CGpositions
with open(temp_output_file,'a',newline='') as output_file:
    line_writer = csv.writer(output_file,delimiter='\t')
    line_writer.writerow(temp_output_list)


# loop over analysis of Decile3 experimental states
# define pending oupt-put variables before start replicates
D3_locus_mean_M_level_pending = []

for i_expts in range(0,N_ExptState_Decile3):
# set current dataframe to first expt. state
    current_CG_site_data_df = ExptState_Decile3_df_list[i_expts]

    N_CG_temp, CG_positions_gene_temp, state_expt_temp, state_gene_temp = \
                create_initial_state_eqbrm_Col0like_regions(current_gene_ID, 
                        initial_state_choice, P_choice, current_segment_start, current_segment_end)
        
    # if N_CG different to Col0 values skip segment
    if N_CG_temp != N_CG:
        invalidate_current_ID = 1
        temp_output_list = [current_gene_ID, N_CG_temp, np.nan]
        temp_output_file = file_name_output_Skipped_IDs
        with open(temp_output_file,'a',newline='') as output_file:
            line_writer = csv.writer(output_file,delimiter='\t')
            line_writer.writerow(temp_output_list)
        break
    
    N_CG_density_temp = (N_CG_temp-1)/(CG_positions_gene_temp[-1]-CG_positions_gene_temp[0])

    # if N_CG_density different to Col0 values skip segment
    if N_CG_density_temp != N_CG_density:
        invalidate_current_ID = 1
        temp_output_list = [current_gene_ID, N_CG_temp, N_CG_density_temp]
        temp_output_file = file_name_output_Skipped_IDs
        with open(temp_output_file,'a',newline='') as output_file:
            line_writer = csv.writer(output_file,delimiter='\t')
            line_writer.writerow(temp_output_list)
        break

    #find mean meth. level for experimental state
    n_u_expt_temp, n_h_expt_temp, n_m_expt_temp = count_state(state_expt_temp,N_CG_temp)
    D3_locus_mean_M_level_pending.append(n_m_expt_temp)

    # write out state
    temp_output_list = state_expt_temp
    temp_output_file = file_name_output_ExptStatesD3
    with open(temp_output_file,'a',newline='') as output_file:
        line_writer = csv.writer(output_file,delimiter='\t')
        line_writer.writerow(temp_output_list)
    # write out SRX
    temp_output_list = [ExptState_Decile3_SRXcodes_list[i_expts]]
    temp_output_file = file_name_output_SRXcodesD3
    with open(temp_output_file,'a',newline='') as output_file:
        line_writer = csv.writer(output_file,delimiter='\t')
        line_writer.writerow(temp_output_list)
    # write out (all D3D4D5) state + SRX code
    temp_output_list = [ExptState_Decile3_SRXcodes_list[i_expts]]+state_expt_temp
    temp_output_file = file_name_output_ExptStatesD3D4D5
    with open(temp_output_file,'a',newline='') as output_file:
        line_writer = csv.writer(output_file,delimiter='\t')
        line_writer.writerow(temp_output_list)




    
    # end of loop over expt. Decile3 states (i_expts)



    
# check whether ID has been invalidated and if so move to next segment
if invalidate_current_ID == 1:
    sys.exit()

# loop over analysis of Decile4 experimental states
# define pending oupt-put variables before start replicates
D4_locus_mean_M_level_pending = []

for i_expts in range(0,N_ExptState_Decile4):
# set current dataframe to first expt. state
    current_CG_site_data_df = ExptState_Decile4_df_list[i_expts]

    N_CG_temp, CG_positions_gene_temp, state_expt_temp, state_gene_temp = \
                create_initial_state_eqbrm_Col0like_regions(current_gene_ID, 
                        initial_state_choice, P_choice, current_segment_start, current_segment_end)
        
    # if N_CG different to Col0 values skip segment
    if N_CG_temp != N_CG:
        invalidate_current_ID = 1
        temp_output_list = [current_gene_ID, N_CG_temp, np.nan]
        temp_output_file = file_name_output_Skipped_IDs
        with open(temp_output_file,'a',newline='') as output_file:
            line_writer = csv.writer(output_file,delimiter='\t')
            line_writer.writerow(temp_output_list)
        break
    
    N_CG_density_temp = (N_CG_temp-1)/(CG_positions_gene_temp[-1]-CG_positions_gene_temp[0])

    # if N_CG_density different to Col0 values skip segment
    if N_CG_density_temp != N_CG_density:
        invalidate_current_ID = 1
        temp_output_list = [current_gene_ID, N_CG_temp, N_CG_densitytemp]
        temp_output_file = file_name_output_Skipped_IDs
        with open(temp_output_file,'a',newline='') as output_file:
            line_writer = csv.writer(output_file,delimiter='\t')
            line_writer.writerow(temp_output_list)
        break

    #find mean meth. level for experimental state
    n_u_expt_temp, n_h_expt_temp, n_m_expt_temp = count_state(state_expt_temp,N_CG_temp)
    D4_locus_mean_M_level_pending.append(n_m_expt_temp)

    # write out state
    temp_output_list = state_expt_temp
    temp_output_file = file_name_output_ExptStatesD4
    with open(temp_output_file,'a',newline='') as output_file:
        line_writer = csv.writer(output_file,delimiter='\t')
        line_writer.writerow(temp_output_list)
    # write out SRX
    temp_output_list = [ExptState_Decile4_SRXcodes_list[i_expts]]
    temp_output_file = file_name_output_SRXcodesD4
    with open(temp_output_file,'a',newline='') as output_file:
        line_writer = csv.writer(output_file,delimiter='\t')
        line_writer.writerow(temp_output_list)
    # write out (all D3D4D5) state + SRX code
    temp_output_list = [ExptState_Decile4_SRXcodes_list[i_expts]]+state_expt_temp
    temp_output_file = file_name_output_ExptStatesD3D4D5
    with open(temp_output_file,'a',newline='') as output_file:
        line_writer = csv.writer(output_file,delimiter='\t')
        line_writer.writerow(temp_output_list)


    
    # end of loop over expt. Decile4 states (i_expts)
    
# check whether ID has been invalidated and if so move to next segment
if invalidate_current_ID == 1:
    sys.exit()

# loop over analysis of Decile5 experimental states
# define pending oupt-put variables before start replicates
D5_locus_mean_M_level_pending = []

for i_expts in range(0,N_ExptState_Decile5):
# set current dataframe to first expt. state
    current_CG_site_data_df = ExptState_Decile5_df_list[i_expts]

    N_CG_temp, CG_positions_gene_temp, state_expt_temp, state_gene_temp = \
                create_initial_state_eqbrm_Col0like_regions(current_gene_ID, 
                        initial_state_choice, P_choice, current_segment_start, current_segment_end)
        
    # if N_CG different to Col0 values skip segment
    if N_CG_temp != N_CG:
        invalidate_current_ID = 1
        temp_output_list = [current_gene_ID, N_CG_temp, np.nan]
        temp_output_file = file_name_output_Skipped_IDs
        with open(temp_output_file,'a',newline='') as output_file:
            line_writer = csv.writer(output_file,delimiter='\t')
            line_writer.writerow(temp_output_list)
        break
    
    N_CG_density_temp = (N_CG_temp-1)/(CG_positions_gene_temp[-1]-CG_positions_gene_temp[0])

    # if N_CG_density different to Col0 values skip segment
    if N_CG_density_temp != N_CG_density:
        invalidate_current_ID = 1
        temp_output_list = [current_gene_ID, N_CG_temp, N_CG_density_temp]
        temp_output_file = file_name_output_Skipped_IDs
        with open(temp_output_file,'a',newline='') as output_file:
            line_writer = csv.writer(output_file,delimiter='\t')
            line_writer.writerow(temp_output_list)
        break

    #find mean meth. level for experimental state
    n_u_expt_temp, n_h_expt_temp, n_m_expt_temp = count_state(state_expt_temp,N_CG_temp)
    D5_locus_mean_M_level_pending.append(n_m_expt_temp)
    
    # write out state
    temp_output_list = state_expt_temp
    temp_output_file = file_name_output_ExptStatesD5
    with open(temp_output_file,'a',newline='') as output_file:
        line_writer = csv.writer(output_file,delimiter='\t')
        line_writer.writerow(temp_output_list)
    # write out SRX
    temp_output_list = [ExptState_Decile5_SRXcodes_list[i_expts]]
    temp_output_file = file_name_output_SRXcodesD5
    with open(temp_output_file,'a',newline='') as output_file:
        line_writer = csv.writer(output_file,delimiter='\t')
        line_writer.writerow(temp_output_list)
    # write out (all D3D4D5) state + SRX code
    temp_output_list = [ExptState_Decile5_SRXcodes_list[i_expts]]+state_expt_temp
    temp_output_file = file_name_output_ExptStatesD3D4D5
    with open(temp_output_file,'a',newline='') as output_file:
        line_writer = csv.writer(output_file,delimiter='\t')
        line_writer.writerow(temp_output_list)


    # end of loop over expt. Decile5 states (i_expts)
    
# check whether ID has been invalidated and if so move to next segment
if invalidate_current_ID == 1:
    sys.exit()


# simulation

# seed random number generator
np.random.seed(initial_seed)

# set sim. oupt-put variables to zero before start replicates
Sim_locus_mean_M_time_pending = []
Sim_locus_mean_M_level_pending = []

#print()
#print(current_gene_ID,N_CG,'1/',1./N_CG_density)
#print(state_gene)
#print(CG_positions_gene)
#print(spont_gain_params)
#print(coop_gain_params)
#print(coop_maint_params)

# run simulation          
sim_outputs_temp = single_gene_sim_transients(burn_in_time, N_CG, 
                state_gene, CG_positions_gene, spont_gain_params, coop_gain_params, coop_maint_params)

# unpack sim_outputs_temp
transients_state_temp = sim_outputs_temp[0]
transients_meth_level_temp = sim_outputs_temp[1]

dat_time_out = transients_state_temp[0]
dat_n_u_out = transients_state_temp[1]
dat_n_m_out = transients_state_temp[2]


state_out = transients_meth_level_temp[0]
state_time_out = transients_meth_level_temp[1]
state_matrix_out = transients_meth_level_temp[2]            


# write out MethLevel
#MethLevel_headings = ['Time', 'M_Frac','U_Frac']
temp_output_file = file_name_output_MethLevel
with open(temp_output_file,'a',newline='') as output_file:
    for i_ in range(len(dat_time_out)):
        temp_output_list = [dat_time_out[i_],dat_n_m_out[i_],dat_n_u_out[i_]]
        line_writer = csv.writer(output_file,delimiter='\t')
        line_writer.writerow(temp_output_list)

# write out simulated state to file: 
temp_output_file = file_name_output_MethStateTime
with open(temp_output_file,'a',newline='') as output_file:
    for i_ in range(len(state_time_out)):
        line_writer = csv.writer(output_file,delimiter='\t')
        line_writer.writerow([state_time_out[i_]])
temp_output_file = file_name_output_MethStateMatrix
with open(temp_output_file,'a',newline='') as output_file:
    for i_ in range(len(state_time_out)):
        line_writer = csv.writer(output_file,delimiter='\t')
        line_writer.writerow(state_matrix_out[i_,:].tolist())
# finished write out simulated state to file

    
# finished simulation


#if (i_gene%10 == 0):
#    timer_end = timer()
#    #print('finished', i_gene, i_reps,'time', timer_end - timer_start, flush=True)
#    
#    with open(file_name_output+'CodeProgress'+'.tsv','a',newline='') as output_file:##~~++~~##
#        line_writer = csv.writer(output_file,delimiter='\t')
#        line_writer.writerow(['finished', i_gene, i_reps,'time', timer_end - timer_start])
#    
#    timer_start = timer()
    


timer_end_total = timer()

with open(file_name_output_CodeProgress,'a',newline='') as output_file:
    line_writer = csv.writer(output_file,delimiter='\t')
    line_writer.writerow([(timer_end_total - timer_start_total)/60.])


#print('end')
#print('time', (timer_end_total - timer_start_total)/60., 'mins')

print('Completed')
