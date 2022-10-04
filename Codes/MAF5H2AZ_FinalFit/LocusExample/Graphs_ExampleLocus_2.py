import numpy as np
import matplotlib as mpl
mpl.use('Agg') # need this to run on cluster?
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
from scipy import stats
#import math
#import sys
import scipy.stats
#import random
import pandas as pd
#from timeit import default_timer as timer
import os
import csv
#import bisect # to calculate local CG density
#import dask.dataframe as dd

import input_params as params

# fetch parameters

# define annoation files
file_in_annotation = params.file_in_annotation
file_in_anno_filt_1 = params.file_in_anno_filt_1
file_in_anno_filt_2 = params.file_in_anno_filt_2
file_in_exclude_filt_1 = params.file_in_exclude_filt_1

filename_params_code = params.filename_params_code

locus_type = params.locus_type
TE_filt = params.TE_filt


P_choice = params.P_choice
rep_time = params.rep_time
initial_state_choice = params.initial_state_choice
N_gen_burn_in = params.N_gen_burn_in

off_rate = params.off_rate

N_CG_min = params.N_CG_min
N_CG_max = params.N_CG_max
N_CG_density_min = params.N_CG_density_min
N_CG_density_max = params.N_CG_density_max

# define linear relationships:
u_scale_val = params.u_scale_val
spacing_cap = params.spacing_cap

u_m_val = params.u_m_val
u_c_val = params.u_c_val
u_cap_val = params.u_cap_val

e_m_val = params.e_m_val
e_c_val = params.e_c_val
e_cap_val = params.e_cap_val

g_m_val = params.g_m_val
g_c_val = params.g_c_val
g_cap_val = params.g_cap_val

delta = params.delta

lambda_gamma = params.lambda_gamma
r_div_gamma = params.r_div_gamma
r_gamma = params.r_gamma

# SR (short-range) parameters
lambda_coop = params.lambda_coop
r_div = params.r_div
r_plat = params.r_plat

# LR (long-range) parameters
lambda_coopIn_LR1 = params.lambda_coopIn_LR1
lambda_coopOut_LR1 = params.lambda_coopOut_LR1
r_div_LR1 = params.r_div_LR1
r_plat_LR1 = params.r_plat_LR1

r_LR1 = params.r_LR1

# unlikely to edit
N_total_IDs = params.N_total_IDs
n_cc = params.n_cc

N_ExptState_Decile3 = 10
N_ExptState_Decile4 = 10
N_ExptState_Decile5 = 10

TimeAverage_StartGen = 50*1000

GraphFolder = "Graphs"
filename_ending = '.tsv'


# Graph inputs #
fig_label = 'F5D'

current_gene_ID_gene1 = 'AT3G22425'
current_gene_ID_gene2 = 'AT1G01230'
current_gene_ID_gene3 = 'AT2G20540'


initial_seed_gene1 = 0
initial_seed_gene2 = 0
initial_seed_gene3 = 0

initial_state_choice_gene1 = '100U'
initial_state_choice_gene2 = '100U'
initial_state_choice_gene3 = '100U'

filename_start_gene1 = 'MAF5_H2AZWT_EqbrmOutput_FinalFit_Pchoice_Excl_InitialState_100U_AT3G22425_Seed_0_'
filename_start_gene2 = 'MAF5_H2AZWT_EqbrmOutput_FinalFit_Pchoice_Excl_InitialState_100U_AT1G01230_Seed_0_'
filename_start_gene3 = 'MAF5_H2AZWT_EqbrmOutput_FinalFit_Pchoice_Excl_InitialState_100U_AT2G20540_Seed_0_'

# windows to average over local CG_density
base_pair_window1 = 61 # 41.


# Load in data files gene1
file_in_MethLevel_gene1 = os.path.join('Output_files', filename_start_gene1 + 'MethLevel'+filename_ending)
file_in_ExptStateCol0_gene1 = os.path.join('Output_files', filename_start_gene1 + 'ExptStateCol0'+filename_ending)
file_in_CGpositions_gene1 = os.path.join('Output_files', filename_start_gene1 + 'CGpositions'+filename_ending)
file_in_ExptStatesD3D4D5_gene1 = os.path.join('Output_files', filename_start_gene1 + 'MethStateExptStatesD3D4D5'+filename_ending)
file_in_MethStateMatrix_gene1 = os.path.join('Output_files', filename_start_gene1 + 'MethStateMatrix'+filename_ending)
file_in_MethStateTime_gene1 = os.path.join('Output_files', filename_start_gene1 + 'MethStateTime'+filename_ending)

MethLevel_gene1_df = pd.read_csv(file_in_MethLevel_gene1, sep="\t{1}",engine='python')
ExptStateCol0_gene1_df = pd.read_csv(file_in_ExptStateCol0_gene1, sep="\t{1}",engine='python',header=None)
CGpositions_gene1_df = pd.read_csv(file_in_CGpositions_gene1, sep="\t{1}",engine='python',header=None)
ExptStatesD3D4D5_gene1_df = pd.read_csv(file_in_ExptStatesD3D4D5_gene1, sep="\t{1}",engine='python',header=None,index_col=0)
MethStateMatrix_gene1_df = pd.read_csv(file_in_MethStateMatrix_gene1, sep="\t{1}",engine='python',header=None)
MethStateTime_gene1_df = pd.read_csv(file_in_MethStateTime_gene1, sep="\t{1}",engine='python',header=None)

StateExptCol0_gene1 = ExptStateCol0_gene1_df.values.flatten()
CG_positions_gene1 = CGpositions_gene1_df.values.flatten()
state_matrix_gene1 = MethStateMatrix_gene1_df.values
state_time_gene1 = MethStateTime_gene1_df.values.flatten()
N_CG_gene1 = len(CG_positions_gene1)
L_CG_gene1 = CG_positions_gene1[-1]-CG_positions_gene1[0]
CG_density_gene1 = (N_CG_gene1-1)/L_CG_gene1
StateExptD3D4D5_gene1 = ExptStatesD3D4D5_gene1_df.values
print(current_gene_ID_gene1)
print('N_CG', N_CG_gene1, 'L_locus', L_CG_gene1, 'CG_density', 1./CG_density_gene1)

# Find mean meth. level of D3D4D5
import collections
MethCountsD3D4D5_gene1 = collections.Counter( StateExptD3D4D5_gene1.flatten() )
MeanMethFracD3D4D5_gene1 = MethCountsD3D4D5_gene1[2]/(MethCountsD3D4D5_gene1[0]+MethCountsD3D4D5_gene1[2])
#print(MethCountsD3D4D5_gene1)
#print(MeanMethFracD3D4D5_gene1)
#print(TimeAverage_StartGen*n_cc)

# Load in data files gene2
file_in_MethLevel_gene2 = os.path.join('Output_files', filename_start_gene2 + 'MethLevel'+filename_ending)
file_in_ExptStateCol0_gene2 = os.path.join('Output_files', filename_start_gene2 + 'ExptStateCol0'+filename_ending)
file_in_CGpositions_gene2 = os.path.join('Output_files', filename_start_gene2 + 'CGpositions'+filename_ending)
file_in_ExptStatesD3D4D5_gene2 = os.path.join('Output_files', filename_start_gene2 + 'MethStateExptStatesD3D4D5'+filename_ending)
file_in_MethStateMatrix_gene2 = os.path.join('Output_files', filename_start_gene2 + 'MethStateMatrix'+filename_ending)
file_in_MethStateTime_gene2 = os.path.join('Output_files', filename_start_gene2 + 'MethStateTime'+filename_ending)

MethLevel_gene2_df = pd.read_csv(file_in_MethLevel_gene2, sep="\t{1}",engine='python')
ExptStateCol0_gene2_df = pd.read_csv(file_in_ExptStateCol0_gene2, sep="\t{1}",engine='python',header=None)
CGpositions_gene2_df = pd.read_csv(file_in_CGpositions_gene2, sep="\t{1}",engine='python',header=None)
ExptStatesD3D4D5_gene2_df = pd.read_csv(file_in_ExptStatesD3D4D5_gene2, sep="\t{1}",engine='python',header=None,index_col=0)
MethStateMatrix_gene2_df = pd.read_csv(file_in_MethStateMatrix_gene2, sep="\t{1}",engine='python',header=None)
MethStateTime_gene2_df = pd.read_csv(file_in_MethStateTime_gene2, sep="\t{1}",engine='python',header=None)

StateExptCol0_gene2 = ExptStateCol0_gene2_df.values.flatten()
CG_positions_gene2 = CGpositions_gene2_df.values.flatten()
state_matrix_gene2 = MethStateMatrix_gene2_df.values
state_time_gene2 = MethStateTime_gene2_df.values.flatten()
N_CG_gene2 = len(CG_positions_gene2)
L_CG_gene2 = CG_positions_gene2[-1]-CG_positions_gene2[0]
CG_density_gene2 = (N_CG_gene2-1)/L_CG_gene2
StateExptD3D4D5_gene2 = ExptStatesD3D4D5_gene2_df.values
print(current_gene_ID_gene2)
print('N_CG', N_CG_gene2, 'L_locus', L_CG_gene2, 'CG_density', 1./CG_density_gene2)

MethCountsD3D4D5_gene2 = collections.Counter( StateExptD3D4D5_gene2.flatten() )
MeanMethFracD3D4D5_gene2 = MethCountsD3D4D5_gene2[2]/(MethCountsD3D4D5_gene2[0]+MethCountsD3D4D5_gene2[2])

# Load in data files gene3
file_in_MethLevel_gene3 = os.path.join('Output_files', filename_start_gene3 + 'MethLevel'+filename_ending)
file_in_ExptStateCol0_gene3 = os.path.join('Output_files', filename_start_gene3 + 'ExptStateCol0'+filename_ending)
file_in_CGpositions_gene3 = os.path.join('Output_files', filename_start_gene3 + 'CGpositions'+filename_ending)
file_in_ExptStatesD3D4D5_gene3 = os.path.join('Output_files', filename_start_gene3 + 'MethStateExptStatesD3D4D5'+filename_ending)
file_in_MethStateMatrix_gene3 = os.path.join('Output_files', filename_start_gene3 + 'MethStateMatrix'+filename_ending)
file_in_MethStateTime_gene3 = os.path.join('Output_files', filename_start_gene3 + 'MethStateTime'+filename_ending)

MethLevel_gene3_df = pd.read_csv(file_in_MethLevel_gene3, sep="\t{1}",engine='python')
ExptStateCol0_gene3_df = pd.read_csv(file_in_ExptStateCol0_gene3, sep="\t{1}",engine='python',header=None)
CGpositions_gene3_df = pd.read_csv(file_in_CGpositions_gene3, sep="\t{1}",engine='python',header=None)
ExptStatesD3D4D5_gene3_df = pd.read_csv(file_in_ExptStatesD3D4D5_gene3, sep="\t{1}",engine='python',header=None,index_col=0)
MethStateMatrix_gene3_df = pd.read_csv(file_in_MethStateMatrix_gene3, sep="\t{1}",engine='python',header=None)
MethStateTime_gene3_df = pd.read_csv(file_in_MethStateTime_gene3, sep="\t{1}",engine='python',header=None)

StateExptCol0_gene3 = ExptStateCol0_gene3_df.values.flatten()
CG_positions_gene3 = CGpositions_gene3_df.values.flatten()
state_matrix_gene3 = MethStateMatrix_gene3_df.values
state_time_gene3 = MethStateTime_gene3_df.values.flatten()
N_CG_gene3 = len(CG_positions_gene3)
L_CG_gene3 = CG_positions_gene3[-1]-CG_positions_gene3[0]
CG_density_gene3 = (N_CG_gene3-1)/L_CG_gene3
StateExptD3D4D5_gene3 = ExptStatesD3D4D5_gene3_df.values
print(current_gene_ID_gene3)
print('N_CG', N_CG_gene3, 'L_locus', L_CG_gene3, 'CG_density', 1./CG_density_gene3)
print()

MethCountsD3D4D5_gene3 = collections.Counter( StateExptD3D4D5_gene3.flatten() )
MeanMethFracD3D4D5_gene3 = MethCountsD3D4D5_gene3[2]/(MethCountsD3D4D5_gene3[0]+MethCountsD3D4D5_gene3[2])
SRX_labels_list = ExptStatesD3D4D5_gene3_df.index.tolist()


#print()
#print(ExptStatesD3D4D5_gene3_df)
#print(StateExptD3D4D5_gene3)
#print(SRX_labels_list)

def TimeSeries_stats(time_vals_, M_frac_vals_, TimeAverage_StartGen_, N_gen_burn_in_, initial_state_choice_):
    # assume each timepoint sampled twice: before and after state updated!

    TimeInterval_Start_cc_ = TimeAverage_StartGen_*n_cc
    TimeInterval_End_cc_ = N_gen_burn_in_*n_cc
    TimeInterval_total_ = (TimeInterval_End_cc_-TimeInterval_Start_cc_)
    
    if time_vals_[-1] != TimeInterval_End_cc_:
        print('TimeAveraging Error1', time_vals_[-1], TimeInterval_End_cc_)

    # record time of first change
    FirstTimeChangedSite_ = time_vals_[1]

    # select for values above threshold. 
    start_index_ = np.argmax(time_vals_ >= TimeInterval_Start_cc_)
    # check first timepoint is correct
    if time_vals_[start_index_] < TimeInterval_Start_cc_:
        print('TimeAveraging Error2', time_vals_[0], TimeInterval_Start_cc_)
    time_vals_SteadyState_ = time_vals_[start_index_-1:].copy()
    M_frac_vals_SteadyState_ = M_frac_vals_[start_index_-1:].copy()
    # set first element of time_vals to TimeInterval_Start_cc_
    time_vals_SteadyState_[0] = TimeInterval_Start_cc_

    #print('Time','M_Frac')
    #for i_ in range(len(time_vals_)):
    #    print(time_vals_[i_],M_frac_vals_[i_])

    mu_ = 0
    sigma_sq_ = 0

    # loop over every second timepoint. 
    for i_ in range(0,len(time_vals_SteadyState_),2):
        time_interval_ = time_vals_SteadyState_[i_+1]-time_vals_SteadyState_[i_]
        mu_ += M_frac_vals_SteadyState_[i_]*time_interval_
        sigma_sq_ += M_frac_vals_SteadyState_[i_]**2*time_interval_

    mu_ = mu_/TimeInterval_total_
    sigma_sq_ = sigma_sq_/TimeInterval_total_
    sigma_sq_ = sigma_sq_ - mu_**2

    # find extreme values of M_frac within averaging windodw
    # and indicies where it occurs
    M_frac_max_ = np.max(M_frac_vals_SteadyState_)
    M_frac_max_indicies_ = np.argwhere(M_frac_vals_SteadyState_ >= M_frac_max_)
    # print('M_frac_Max')
    # print(M_frac_max_,M_frac_max_indicies_)
    M_frac_min_ = np.min(M_frac_vals_SteadyState_)
    M_frac_min_indicies_ = np.argwhere(M_frac_vals_SteadyState_ <= M_frac_min_)
    # print('M_frac_Min')
    # print(M_frac_min_,M_frac_min_indicies_)
    # print()

    # find time M_Frac first crosses steady state value
    if initial_state_choice_ == '100U':
        FirstTimeSteadyState_index_= np.argmax(M_frac_vals_ >= mu_)
        FirstTimeSteadyState_ = time_vals_[FirstTimeSteadyState_index_]
    elif initial_state_choice_ == '100M':
        FirstTimeSteadyState_index_= np.argmax(M_frac_vals_ <= mu_)
        FirstTimeSteadyState_ = time_vals_[FirstTimeSteadyState_index_]
    else:
        FirstTimeSteadyState_ = np.nan

    # find magnitude and duration of greatest fluctuation
    # bounded by lengh of averaging window
    MaxFluct_Mag_ = np.nan
    MaxFluct_Time_ = np.nan
    
    # print('abs(M_frac_max_ - mu_) > abs(M_frac_min_ - mu_)')
    # print((M_frac_max_ - mu_), (M_frac_min_ - mu_))
    if abs(M_frac_max_ - mu_) > abs(M_frac_min_ - mu_):
        # fluctuation will have positive amplitude! 
        # print('positive amplitude')
        MaxFluct_Mag_ = abs(M_frac_max_ - mu_)
        # check fluctuation has correct sign
        if M_frac_max_ - mu_ <= 0:
            print('MaxFluct_Mag error', M_frac_max_ - mu_)

        # loop over all occurances of extreme fluctuation
        MaxFluct_Time_ = 0.0
        # i_1_ loops over values of M_frac_max_indicies_
        # i_2_ loops over selected elements of M_frac_vals_SteadyState_
        for i_1_ in M_frac_max_indicies_:
            i_1_ = np.int(i_1_)
            #print(type(i_1_),i_1_)
            Fluct_start_time_ = TimeInterval_Start_cc_
            Fluct_end_time_ = TimeInterval_End_cc_
            # search back in time for first occurance of negative amplitude
            for i_2_ in range(i_1_-1,-1,-1):
                if (M_frac_vals_SteadyState_[i_2_] - mu_) <=0:
                    # update Fluct_start_time_ and break loop
                    Fluct_start_time_ = time_vals_SteadyState_[i_2_]
                    break
            # search forward in time for first occurance of negative amplitude
            for i_2_ in range(i_1_+1,len(M_frac_vals_SteadyState_)):
                if (M_frac_vals_SteadyState_[i_2_] - mu_) <=0:
                    # update Fluct_end_time_ and break loop
                    Fluct_end_time_ = time_vals_SteadyState_[i_2_]
                    break
            # test whether exceeds current MaxFluct_Time_ and if so update
            if (Fluct_end_time_ - Fluct_start_time_) > MaxFluct_Time_:
                MaxFluct_Time_ = Fluct_end_time_ - Fluct_start_time_
    else:
        # print('negative amplitude')
        # neglect unlikely posibility +ve and -ve extreme fluctuations have identical magnitudes. 
        # fluctuation will have negative amplitude! 
        MaxFluct_Mag_ = abs(M_frac_min_ - mu_)
        # check fluctuation has correct sign
        if M_frac_min_ - mu_ >= 0:
            print('MaxFluct_Mag error2', M_frac_min_ - mu_)

        # loop over all occurances of extreme fluctuation
        MaxFluct_Time_ = 0.0
        # i_1_ loops over values of M_frac_min_indicies_
        # i_2_ loops over selected elements of M_frac_vals_SteadyState_
        for i_1_ in M_frac_min_indicies_:
            i_1_ = np.int(i_1_)
            Fluct_start_time_ = TimeInterval_Start_cc_
            Fluct_end_time_ = TimeInterval_End_cc_
            # search back in time for first occurance of positive amplitude
            for i_2_ in range(i_1_-1,-1,-1):
                if (M_frac_vals_SteadyState_[i_2_] - mu_) >=0:
                    # update Fluct_start_time_ and break loop
                    Fluct_start_time_ = time_vals_SteadyState_[i_2_]
                    break
            # search forward in time for first occurance of positive amplitude
            for i_2_ in range(i_1_+1,len(M_frac_vals_SteadyState_)):
                if (M_frac_vals_SteadyState_[i_2_] - mu_) >=0:
                    # update Fluct_end_time_ and break loop
                    Fluct_end_time_ = time_vals_SteadyState_[i_2_]
                    break
            # test whether exceeds current MaxFluct_Time_ and if so update
            if (Fluct_end_time_ - Fluct_start_time_) > MaxFluct_Time_:
                MaxFluct_Time_ = Fluct_end_time_ - Fluct_start_time_

    return mu_, np.sqrt(sigma_sq_), M_frac_max_, M_frac_min_, FirstTimeSteadyState_, MaxFluct_Mag_, MaxFluct_Time_, FirstTimeChangedSite_

###

# calculate time averaged statistics

M_Frac_star_mu_gene1, M_Frac_star_simga_gene1, M_Frac_star_max_gene1, M_Frac_star_min_gene1, FirstTimeSteadyState_gene1, MaxFluct_Mag_gene1, MaxFluct_Time_gene1, FirstTimeChangedSite_gene1 = \
        TimeSeries_stats(MethLevel_gene1_df['Time'].values, MethLevel_gene1_df['M_Frac'].values, TimeAverage_StartGen, N_gen_burn_in, initial_state_choice_gene1)
M_Frac_star_mu_gene2, M_Frac_star_simga_gene2, M_Frac_star_max_gene2, M_Frac_star_min_gene2, FirstTimeSteadyState_gene2, MaxFluct_Mag_gene2, MaxFluct_Time_gene2, FirstTimeChangedSite_gene2 = \
        TimeSeries_stats(MethLevel_gene2_df['Time'].values, MethLevel_gene2_df['M_Frac'].values, TimeAverage_StartGen, N_gen_burn_in, initial_state_choice_gene2)
M_Frac_star_mu_gene3, M_Frac_star_simga_gene3, M_Frac_star_max_gene3, M_Frac_star_min_gene3, FirstTimeSteadyState_gene3, MaxFluct_Mag_gene3, MaxFluct_Time_gene3, FirstTimeChangedSite_gene3 = \
        TimeSeries_stats(MethLevel_gene3_df['Time'].values, MethLevel_gene3_df['M_Frac'].values, TimeAverage_StartGen, N_gen_burn_in, initial_state_choice_gene3)
print('M_Frac_star_mu_gene1', 'M_Frac_star_simga_gene1', 'M_Frac_star_max_gene1', 'M_Frac_star_min_gene1', 'FirstTimeSteadyState_gene1', 'MaxFluct_Mag_gene1', 'MaxFluct_Time_gene1', 'FirstTimeChangedSite_gene1')
print(M_Frac_star_mu_gene1, M_Frac_star_simga_gene1, M_Frac_star_max_gene1, M_Frac_star_min_gene1, FirstTimeSteadyState_gene1/n_cc, MaxFluct_Mag_gene1, MaxFluct_Time_gene1/n_cc, FirstTimeChangedSite_gene1/n_cc)
print('M_Frac_star_mu_gene2', 'M_Frac_star_simga_gene2', 'M_Frac_star_max_gene2', 'M_Frac_star_min_gene2', 'FirstTimeSteadyState_gene2', 'MaxFluct_Mag_gene2', 'MaxFluct_Time_gene2', 'FirstTimeChangedSite_gene2')
print(M_Frac_star_mu_gene2, M_Frac_star_simga_gene2, M_Frac_star_max_gene2, M_Frac_star_min_gene2, FirstTimeSteadyState_gene2/n_cc, MaxFluct_Mag_gene2, MaxFluct_Time_gene2/n_cc, FirstTimeChangedSite_gene2/n_cc)
print('M_Frac_star_mu_gene3', 'M_Frac_star_simga_gene3', 'M_Frac_star_max_gene3', 'M_Frac_star_min_gene3', 'FirstTimeSteadyState_gene3', 'MaxFluct_Mag_gene3', 'MaxFluct_Time_gene3', 'FirstTimeChangedSite_gene3')
print(M_Frac_star_mu_gene3, M_Frac_star_simga_gene3, M_Frac_star_max_gene3, M_Frac_star_min_gene3, FirstTimeSteadyState_gene3/n_cc, MaxFluct_Mag_gene3, MaxFluct_Time_gene3/n_cc, FirstTimeChangedSite_gene3/n_cc)
print()

###

# function to caclulate local CG_density at a given point: x_
# base_pair_window assumed to be an odd integer
import bisect

def local_CG_density_function(L_, CG_positions_,base_pair_window_,x_):
    x_start_ = max(0, x_ - (base_pair_window_-1)/2)
    x_end_ = min(L_, x_ + (base_pair_window_-1)/2)
    
    i_start_ = bisect.bisect_left(CG_positions_,x_start_)
    i_end_ = bisect.bisect_right(CG_positions_,x_end_)
    
    CG_count_ = len(CG_positions_[i_start_:i_end_])
    CG_local_density_ = CG_count_/base_pair_window_

    return CG_local_density_

###

# Calculate local CG density 
# rescale CG_positions to start at 0
CG_positions_origin_temp = CG_positions_gene1[0]
CG_positions_gene1 = [CG_positions_gene1[i_]-CG_positions_origin_temp for i_ in range(len(CG_positions_gene1))]
CG_positions_origin_temp = CG_positions_gene2[0]
CG_positions_gene2 = [CG_positions_gene2[i_]-CG_positions_origin_temp for i_ in range(len(CG_positions_gene2))]
CG_positions_origin_temp = CG_positions_gene3[0]
CG_positions_gene3 = [CG_positions_gene3[i_]-CG_positions_origin_temp for i_ in range(len(CG_positions_gene3))]

CG_site_values_1 = np.arange(1,len(CG_positions_gene1)+1)
CG_site_values_2 = np.arange(1,len(CG_positions_gene2)+1)
CG_site_values_3 = np.arange(1,len(CG_positions_gene3)+1)

CG_local_density_1 = [local_CG_density_function(CG_positions_gene1[-1], CG_positions_gene1,base_pair_window1,x_) 
                      for x_ in CG_positions_gene1]
CG_local_density_2 = [local_CG_density_function(CG_positions_gene2[-1], CG_positions_gene2,base_pair_window1,x_) 
                      for x_ in CG_positions_gene2]
CG_local_density_3 = [local_CG_density_function(CG_positions_gene3[-1], CG_positions_gene3,base_pair_window1,x_) 
                      for x_ in CG_positions_gene3]


local_CG_number_1 = [x_*base_pair_window1 for x_ in CG_local_density_1]
local_CG_number_2 = [x_*base_pair_window1 for x_ in CG_local_density_2]
local_CG_number_3 = [x_*base_pair_window1 for x_ in CG_local_density_3]

print( 'local_CG_number Min', 'local_CG_number Max')
print(min(local_CG_number_1),max(local_CG_number_1))
print(min(local_CG_number_2),max(local_CG_number_2))
print(min(local_CG_number_3),max(local_CG_number_3))
print()



# # graph format #
# for spine in ['left','right','top','bottom']:
#     ax.spines[spine].set_color('k')
#     ax.spines[spine].set_linewidth(0.8)
# ax.set_facecolor('white')

# #ax.grid(False)
# ax.grid(b=True, which='major', color='lightgrey', linestyle=':',linewidth=1)
# # end graph format #




# Start Graph
time_scaler = 1.e+3
time_scaler_label = ' ($\\times 10^{3}$)'
time_start = 0
time_end = 100000/time_scaler

from matplotlib import colors

# rescale axes from cell cycles to years
# assume n_cc cell cycles per plant generation in germline
state_data_in1 = state_matrix_gene1
time_gen_in1 = [x_ / (time_scaler*n_cc) for x_ in state_time_gene1]

state_data_in2 = state_matrix_gene2
time_gen_in2 = [x_ / (time_scaler*n_cc) for x_ in state_time_gene2]

state_data_in3 = state_matrix_gene3
time_gen_in3 = [x_ / (time_scaler*n_cc) for x_ in state_time_gene3]


# make mesh-grid
x_grid_1 = np.tile(np.arange(N_CG_gene1+1),(len(time_gen_in1),1))
y_grid_1 = np.tile(time_gen_in1,(N_CG_gene1+1,1)).T
x_grid_2 = np.tile(np.arange(N_CG_gene2+1),(len(time_gen_in2),1))
y_grid_2 = np.tile(time_gen_in2,(N_CG_gene2+1,1)).T
x_grid_3 = np.tile(np.arange(N_CG_gene3+1),(len(time_gen_in3),1))
y_grid_3 = np.tile(time_gen_in3,(N_CG_gene3+1,1)).T

# x_ExpGrid_1 = np.tile(np.arange(N_CG_gene1+1),(len(np.arange(0,11)),1))
# y_ExpGrid_1 = np.tile(np.arange(0,11),(N_CG_gene1+1,1)).T
# x_ExpGrid_2 = np.tile(np.arange(N_CG_gene2+1),(len(np.arange(0,11)),1))
# y_ExpGrid_2 = np.tile(np.arange(0,11),(N_CG_gene2+1,1)).T
# x_ExpGrid_3 = np.tile(np.arange(N_CG_gene3+1),(len(np.arange(0,11)),1))
# y_ExpGrid_3 = np.tile(np.arange(0,11),(N_CG_gene3+1,1)).T
x_ExpGrid_1 = np.tile(np.arange(N_CG_gene1+1),(len(np.arange(0,31)),1))
y_ExpGrid_1 = np.tile(np.arange(0,31),(N_CG_gene1+1,1)).T
x_ExpGrid_2 = np.tile(np.arange(N_CG_gene2+1),(len(np.arange(0,31)),1))
y_ExpGrid_2 = np.tile(np.arange(0,31),(N_CG_gene2+1,1)).T
x_ExpGrid_3 = np.tile(np.arange(N_CG_gene3+1),(len(np.arange(0,31)),1))
y_ExpGrid_3 = np.tile(np.arange(0,31),(N_CG_gene3+1,1)).T


# make a color map of fixed colors
my_cmap_sim = colors.ListedColormap(['b', 'r'])
cmap_bounds_sim=[-0.5,1e-100,2.5]
cmap_norm_sim = colors.BoundaryNorm(cmap_bounds_sim, my_cmap_sim.N)

my_cmap_data = colors.ListedColormap(['b', 'grey', 'r'])
cmap_bounds_data=[-0.5,0.5,1.5,2.5]
cmap_norm_data = colors.BoundaryNorm(cmap_bounds_data, my_cmap_data.N)


fig1, ax1 = plt.subplots(6,8,figsize=(8,11), gridspec_kw={"height_ratios":[1, 0.01, 0.04, 0.01, 0.04, 0.3], 
                                                           "width_ratios":[2*N_CG_gene1,30,4,2*N_CG_gene2,30,4,2*N_CG_gene3,30]})
# figsize is in inches!

ax1[0,2].set_visible(False)
ax1[0,5].set_visible(False)
ax1[1,0].set_visible(False)
ax1[1,2].set_visible(False)
ax1[1,5].set_visible(False)
ax1[2,2].set_visible(False)
ax1[2,5].set_visible(False)
ax1[3,0].set_visible(False)
ax1[3,2].set_visible(False)
ax1[3,5].set_visible(False)
ax1[4,2].set_visible(False)
ax1[4,5].set_visible(False)
ax1[5,2].set_visible(False)
ax1[5,5].set_visible(False)

ax1[1,3].set_visible(False)
ax1[1,6].set_visible(False)
ax1[3,3].set_visible(False)
ax1[3,6].set_visible(False)

ax1[1,1].set_visible(False)
ax1[1,4].set_visible(False)
ax1[1,7].set_visible(False)
ax1[2,1].set_visible(False)
ax1[2,4].set_visible(False)
ax1[2,7].set_visible(False)
ax1[3,1].set_visible(False)
ax1[3,4].set_visible(False)
ax1[3,7].set_visible(False)
ax1[4,1].set_visible(False)
ax1[4,4].set_visible(False)
ax1[4,7].set_visible(False)
ax1[5,1].set_visible(False)
ax1[5,4].set_visible(False)
ax1[5,7].set_visible(False)

# # merge three of the axis: # not working! 
# import matplotlib.gridspec as gridspec
# # find the gridspec of the 'first' axis
# gs_vlaues = ax1[1,7].get_gridspec()
# # remove the three axis to be combined
# ax1[1,7].remove()
# ax1[2,7].remove()
# ax1[3,7].remove()
# ax1_merged = fig1.add_subplot(gs[1:4,7])

ax1[0,0].set_title("                    %s\n\n" 
                    % (current_gene_ID_gene1,), fontsize=10)
ax1[0,3].set_title("                    %s\n\n"
                    % (current_gene_ID_gene2), fontsize=10)
ax1[0,6].set_title("                    %s\n\n"
                    % (current_gene_ID_gene3), fontsize=10)


im00_1 = ax1[0,0].pcolormesh(x_grid_1,y_grid_1,state_data_in1, cmap=my_cmap_sim, norm=cmap_norm_sim , vmin=0.,vmax=2.,shading='flat',snap=True)
im01_1 = ax1[0,3].pcolormesh(x_grid_2,y_grid_2,state_data_in2, cmap=my_cmap_sim, norm=cmap_norm_sim , vmin=0.,vmax=2.,shading='flat',snap=True)
im02_1 = ax1[0,6].pcolormesh(x_grid_3,y_grid_3,state_data_in3, cmap=my_cmap_sim, norm=cmap_norm_sim , vmin=0.,vmax=2.,shading='flat',snap=True)

ax1[0,0].set_facecolor('whitesmoke')
ax1[0,3].set_facecolor('whitesmoke')
ax1[0,6].set_facecolor('whitesmoke')

ax1[0,0].set_ylabel("Generation"+time_scaler_label, fontsize=10)
#ax1[0,3].set_xlabel("Site number", fontsize=10)

ax1[0,0].set_ylim(time_end,time_start)
ax1[0,3].set_ylim(time_end,time_start)
ax1[0,6].set_ylim(time_end,time_start)

ax1[0,3].set_yticks([])
ax1[0,6].set_yticks([])

ax1[0,0].xaxis.set_tick_params(labelsize=12)
ax1[0,0].yaxis.set_tick_params(labelsize=12)
ax1[0,3].xaxis.set_tick_params(labelsize=12)
ax1[0,3].yaxis.set_tick_params(labelsize=12)
ax1[0,6].xaxis.set_tick_params(labelsize=12)
ax1[0,6].yaxis.set_tick_params(labelsize=12)

ax1[0,0].yaxis.set_major_locator(plt.MaxNLocator(11))


linewidth_choice = 0.5

ax1[0,1].fill_betweenx(MethLevel_gene1_df['Time'].values/(time_scaler*n_cc), 
                        M_Frac_star_mu_gene1-M_Frac_star_simga_gene1,
                        M_Frac_star_mu_gene1+M_Frac_star_simga_gene1,
                        alpha=0.15,color='r')
ax1[0,1].axvline(M_Frac_star_mu_gene1, linestyle='--', color='r',linewidth=1)
ax1[0,1].axvline(MeanMethFracD3D4D5_gene1, linestyle=':', color='k',linewidth=2)
# ax1[0,1].plot(MethLevel_gene1_df['U_Frac'].values, MethLevel_gene1_df['Time'].values/(time_scaler*n_cc), 
#                 color='b',linestyle='-',linewidth=1)
ax1[0,1].plot(MethLevel_gene1_df['M_Frac'].values, MethLevel_gene1_df['Time'].values/(time_scaler*n_cc), 
                color='r',linestyle='-',linewidth=linewidth_choice)
ax1[0,1].set_ylim(time_end,time_start)
ax1[0,1].set_xlim(0,1)
ax1[0,1].xaxis.set_tick_params(labelsize=12)
ax1[0,1].set_yticks([])
ax1[0,1].set_xlabel("mCG", fontsize=10)
ax1[0,1].xaxis.tick_top()
ax1[0,1].xaxis.set_label_position('top') 
ax1[0,1].xaxis.set_major_locator(plt.MaxNLocator(1))

ax1[0,4].fill_betweenx(MethLevel_gene2_df['Time'].values/(time_scaler*n_cc), 
                        M_Frac_star_mu_gene2-M_Frac_star_simga_gene2,
                        M_Frac_star_mu_gene2+M_Frac_star_simga_gene2,
                        alpha=0.15,color='r')
ax1[0,4].axvline(M_Frac_star_mu_gene2, linestyle='--', color='r',linewidth=1)
ax1[0,4].axvline(MeanMethFracD3D4D5_gene2, linestyle=':', color='k',linewidth=2)
# ax1[0,4].plot(MethLevel_gene2_df['U_Frac'].values, MethLevel_gene2_df['Time'].values/(time_scaler*n_cc), 
#                 color='b',linestyle='-',linewidth=1)
ax1[0,4].plot(MethLevel_gene2_df['M_Frac'].values, MethLevel_gene2_df['Time'].values/(time_scaler*n_cc), 
                color='r',linestyle='-',linewidth=linewidth_choice)
ax1[0,4].set_ylim(time_end,time_start)
ax1[0,4].set_xlim(0,1)
ax1[0,4].xaxis.set_tick_params(labelsize=12)
ax1[0,4].set_yticks([])
ax1[0,4].set_xlabel("mCG", fontsize=10)
ax1[0,4].xaxis.tick_top()
ax1[0,4].xaxis.set_label_position('top') 
ax1[0,4].xaxis.set_major_locator(plt.MaxNLocator(1))

ax1[0,7].fill_betweenx(MethLevel_gene3_df['Time'].values/(time_scaler*n_cc), 
                        M_Frac_star_mu_gene3-M_Frac_star_simga_gene3,
                        M_Frac_star_mu_gene3+M_Frac_star_simga_gene3,
                        alpha=0.15,color='r')
ax1[0,7].axvline(M_Frac_star_mu_gene3, linestyle='--', color='r',linewidth=1)
ax1[0,7].axvline(MeanMethFracD3D4D5_gene3, linestyle=':', color='k',linewidth=2)
# ax1[0,7].plot(MethLevel_gene3_df['U_Frac'].values, MethLevel_gene3_df['Time'].values/(time_scaler*n_cc), 
#                 color='b',linestyle='-',linewidth=1)
ax1[0,7].plot(MethLevel_gene3_df['M_Frac'].values, MethLevel_gene3_df['Time'].values/(time_scaler*n_cc), 
                color='r',linestyle='-',linewidth=linewidth_choice)
ax1[0,7].set_ylim(time_end,time_start)
ax1[0,7].set_xlim(0,1)
ax1[0,7].xaxis.set_tick_params(labelsize=12)
ax1[0,7].set_yticks([])
ax1[0,7].set_xlabel("mCG", fontsize=10)
ax1[0,7].xaxis.tick_top()
ax1[0,7].xaxis.set_label_position('top') 
ax1[0,7].xaxis.set_major_locator(plt.MaxNLocator(1))

ax1[0,1].xaxis.labelpad = -8
ax1[0,4].xaxis.labelpad = -8
ax1[0,7].xaxis.labelpad = -8

####

print('max density',np.max(local_CG_number_1),np.max(local_CG_number_2),np.max(local_CG_number_3))
Max_CG_denisty1 = np.max(local_CG_number_1)
Max_CG_denisty2 = np.max(local_CG_number_2)
Max_CG_denisty3 = np.max(local_CG_number_3)
Max_CG_denstiy = np.max([Max_CG_denisty1,Max_CG_denisty2,Max_CG_denisty3])
print(Max_CG_denstiy)

from matplotlib import cm
my_cmap2 = cm.get_cmap('Oranges', Max_CG_denstiy)    # N discrete colors


im10 = ax1[2,0].imshow(np.array(local_CG_number_1).reshape(1,-1), interpolation='none', cmap=my_cmap2,
                      extent=[0,N_CG_gene1,0,1],vmin=1,vmax=Max_CG_denstiy)
ax1[2,0].set_aspect('auto')
ax1[2,0].set_yticks([])
ax1[2,0].grid(False)
#ax1[2,0].set_xlabel("Site number", fontsize=10)
#ax1[2,0].set_title("CG density", fontsize=10)

im13 = ax1[2,3].imshow(np.array(local_CG_number_2).reshape(1,-1), interpolation='none', cmap=my_cmap2,
                      extent=[0,N_CG_gene2,0,1],vmin=1,vmax=Max_CG_denstiy)
ax1[2,3].set_aspect('auto')
ax1[2,3].set_yticks([])
ax1[2,3].grid(False)
#ax1[2,3].set_xlabel("Site number", fontsize=10)
ax1[2,3].set_title("CG density", fontsize=10)

im16 = ax1[2,6].imshow(np.array(local_CG_number_3).reshape(1,-1), interpolation='none', cmap=my_cmap2,
                      extent=[0,N_CG_gene3,0,1],vmin=1,vmax=Max_CG_denstiy)
ax1[2,6].set_aspect('auto')
ax1[2,6].set_yticks([])
ax1[2,6].grid(False)
#ax1[2,6].set_xlabel("Site number", fontsize=10)
#ax1[2,6].set_title("CG density", fontsize=10)

# make colour-bar
cbar_0 = fig1.colorbar(im10,ax=[ax1[1,7],ax1[2,7],ax1[3,7]], orientation='vertical',extend='both',aspect=6,fraction=0.4,shrink=1.2,anchor=(-1.0, -0.0)) 
#cbar_0.set_label("\nNo. CG within\n$\pm$ %.0f bp"% (1.*(base_pair_window1-1.)/2.), fontsize=10)
cbar_0.ax.set_title("\nNo. CG within\n$\pm$ %.0f bp"% (1.*(base_pair_window1-1.)/2.), fontsize=10)


# define colorbar ticks
cbar_tick_locs = []
cbar_tick_label = []
for i_ in range(1,np.int(Max_CG_denstiy)+1):
    tick_loc_temp = 1. + (i_-0.5)*(Max_CG_denstiy-1.)/(Max_CG_denstiy)
    cbar_tick_locs.append(tick_loc_temp)
    cbar_tick_label.append(str(i_))
cbar_0.set_ticks(cbar_tick_locs)
cbar_0.set_ticklabels(cbar_tick_label)
cbar_0.ax.tick_params(labelsize=10)



####
# fig1.colorbar(im10,cax=ax1[3,0], orientation='horizontal',extend='both')
# #ax1[3,0].set_xlabel("No. CG sites\n within $\pm$ %.0f bp"% (1.*(base_pair_window1-1.)/2.), fontsize=10)
# #ax1[3,0].set_title("Color-bar" , fontsize=10)
# ax1[3,0].xaxis.set_major_locator(plt.MaxNLocator(3))
# ax1[3,0].set_yticks([])

# fig1.colorbar(im13,cax=ax1[3,3], orientation='horizontal',extend='both')
# #ax1[3,3].set_xlabel("Color-bar: No. CG sites within $\pm$ %.0f bp"% (1.*(base_pair_window1-1.)/2.), fontsize=10)
# ax1[3,3].set_title("Color-bar: No. CG sites within $\pm$ %.0f bp"% (1.*(base_pair_window1-1.)/2.) , fontsize=10)
# ax1[3,3].xaxis.set_major_locator(plt.MaxNLocator(4))
# ax1[3,3].set_yticks([])

# fig1.colorbar(im16,cax=ax1[3,6], orientation='horizontal',extend='both')
# #ax1[3,6].set_xlabel("No. CG sites\n within $\pm$ %.0f bp"% (1.*(base_pair_window1-1.)/2.), fontsize=10)
# #ax1[3,6].set_title("Color-bar" , fontsize=10)
# ax1[3,6].xaxis.set_major_locator(plt.MaxNLocator(4))
# ax1[3,6].set_yticks([])

####
im11 = ax1[4,0].imshow(np.expand_dims(StateExptCol0_gene1, axis=0), interpolation='none', cmap=my_cmap_data,
                      extent=[0,N_CG_gene1,0,1],norm=cmap_norm_data, aspect='auto')
#ax1[4,0].set_xlabel("Site number", fontsize=10)
#ax1[4,0].set_title("Col-0 state", fontsize=10)
ax1[4,0].set_yticks([])
ax1[4,0].grid(False)

im11 = ax1[4,3].imshow(np.expand_dims(StateExptCol0_gene2, axis=0), interpolation='none', cmap=my_cmap_data,
                      extent=[0,N_CG_gene2,0,1],norm=cmap_norm_data, aspect='auto')
#ax1[4,3].set_xlabel("Site number", fontsize=10)
ax1[4,3].set_title("Col-0 state", fontsize=10)
ax1[4,3].set_yticks([])
ax1[4,3].grid(False)

im11 = ax1[4,6].imshow(np.expand_dims(StateExptCol0_gene3, axis=0), interpolation='none', cmap=my_cmap_data,
                      extent=[0,N_CG_gene3,0,1],norm=cmap_norm_data, aspect='auto')
#ax1[4,6].set_xlabel("Site number", fontsize=10)
#ax1[4,6].set_title("Col-0 state", fontsize=10)
ax1[4,6].set_yticks([])
ax1[4,6].grid(False)






# ax1[5,0].pcolormesh(x_ExpGrid_1,y_ExpGrid_1,StateExptD4_gene1, cmap=my_cmap_data, norm=cmap_norm_data , vmin=0.,vmax=2.,shading='flat')
# ax1[5,3].pcolormesh(x_ExpGrid_2,y_ExpGrid_2,StateExptD4_gene2, cmap=my_cmap_data, norm=cmap_norm_data , vmin=0.,vmax=2.,shading='flat')
# ax1[5,6].pcolormesh(x_ExpGrid_3,y_ExpGrid_3,StateExptD4_gene3, cmap=my_cmap_data, norm=cmap_norm_data , vmin=0.,vmax=2.,shading='flat')
ax1[5,0].pcolormesh(x_ExpGrid_1,y_ExpGrid_1,StateExptD3D4D5_gene1
, cmap=my_cmap_data, norm=cmap_norm_data , vmin=0.,vmax=2.,shading='flat')
ax1[5,3].pcolormesh(x_ExpGrid_2,y_ExpGrid_2,StateExptD3D4D5_gene2
, cmap=my_cmap_data, norm=cmap_norm_data , vmin=0.,vmax=2.,shading='flat')
ax1[5,6].pcolormesh(x_ExpGrid_3,y_ExpGrid_3,StateExptD3D4D5_gene3
, cmap=my_cmap_data, norm=cmap_norm_data , vmin=0.,vmax=2.,shading='flat')

ax1[0,3].set_xlabel("Site number", fontsize=10)
ax1[2,3].set_xlabel("Site number", fontsize=10)
ax1[5,3].set_xlabel("Site number", fontsize=10)

ax1[5,0].set_ylabel("Accession", fontsize=10)
#ax1[5,0].set_title("Col-like\n accessions", fontsize=10)
ax1[5,3].set_title("\nCol-like accessions", fontsize=10)
#ax142,6].set_title("Col-like\n accessions", fontsize=10)

ax1[5,0].yaxis.tick_right()
ax1[5,3].yaxis.tick_right()
ax1[5,6].yaxis.tick_right()

ExptReps_tick_pos1 = [0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5]
ExptReps_tick_pos2 = [10.5,11.5,12.5,13.5,14.5,15.5,16.5,17.5,18.5,19.5]
ExptReps_tick_pos3 = [20.5,21.5,22.5,23.5,24.5,25.5,26.5,27.5,28.5,29.5]
ExptReps_tick_pos = ExptReps_tick_pos1+ExptReps_tick_pos2+ExptReps_tick_pos3

# ExptRep_tick_lables_D3 = ['Accn_1','Accn_2','Accn_3','Accn_4','Accn_5','Accn_6','Accn_7','Accn_8','Accn_9','Accn_10']
# ExptRep_tick_lables_D4 = ['Accn_11','Accn_12','Accn_13','Accn_14','Accn_15','Accn_16','Accn_17','Accn_18','Accn_19','Accn_20']
# ExptRep_tick_lables_D5 = ['Accn_21','Accn_22','Accn_23','Accn_24','Accn_25','Accn_26','Accn_27','Accn_28','Accn_29','Accn_30']
ExptRep_tick_lables_D3 = ['','','','','5','','','','','10']
ExptRep_tick_lables_D4 = ['','','','','15','','','','','20']
ExptRep_tick_lables_D5 = ['','','','','25','','','','','30']

ExptRep_tick_lables = ExptRep_tick_lables_D3 + ExptRep_tick_lables_D4 + ExptRep_tick_lables_D5
ax1[5,0].set_yticks(ExptReps_tick_pos)
ax1[5,0].set_yticklabels(ExptRep_tick_lables, minor=False, rotation=0,fontsize=10)
ax1[5,3].set_yticks(ExptReps_tick_pos)
ax1[5,3].set_yticklabels(ExptRep_tick_lables, minor=False, rotation=0,fontsize=10)
ax1[5,6].set_yticks(ExptReps_tick_pos)
ax1[5,6].set_yticklabels(ExptRep_tick_lables, minor=False, rotation=0,fontsize=10)

for i_ in range(1,30):
    ax1[5,0].axhline(y=i_,color='k',linewidth=0.4)
    ax1[5,3].axhline(y=i_,color='k',linewidth=0.4)
    ax1[5,6].axhline(y=i_,color='k',linewidth=0.4)




fig1.subplots_adjust(wspace=0.0)
#fig1.subplots_adjust(hspace=0.45)
fig1.subplots_adjust(hspace=0.4)


#fig1.subplots_adjust(top=graph_scale2)


fig1.savefig(os.path.join(GraphFolder,"ExampleLocus_2_"+fig_label+".png"),bbox_inches='tight')
# End Graph


