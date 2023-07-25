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

### fetch parameters  ###

# # define annoation files
file_in_annotation = params.file_in_annotation
# file_in_anno_filt_1 = params.file_in_anno_filt_1
# file_in_anno_filt_2 = params.file_in_anno_filt_2
# file_in_exclude_filt_1 = params.file_in_exclude_filt_1

# ##filename_params_code = params.filename_params_code

# locus_type = params.locus_type
# TE_filt = params.TE_filt

# N_corrns_bins = params.N_corrns_bins

N_reps = params.N_reps
# initial_seed = params.initial_seed
# P_choice = params.P_choice
# rep_time = params.rep_time
# initial_state_choice = params.initial_state_choice
# N_gen_burn_in = params.N_gen_burn_in

# off_rate = params.off_rate

# N_CG_min = params.N_CG_min
# N_CG_max = params.N_CG_max
# N_CG_density_min = params.N_CG_density_min
# N_CG_density_max = params.N_CG_density_max

# # define linear relationships:
# u_scale_val = params.u_scale_val
# spacing_cap = params.spacing_cap

# u_m_val = params.u_m_val
# u_c_val = params.u_c_val
# u_cap_val = params.u_cap_val

# e_m_val = params.e_m_val
# e_c_val = params.e_c_val
# e_cap_val = params.e_cap_val

# g_m_val = params.g_m_val
# g_c_val = params.g_c_val
# g_cap_val = params.g_cap_val

# delta = params.delta

# lambda_gamma = params.lambda_gamma
# r_div_gamma = params.r_div_gamma
# r_gamma = params.r_gamma

# # SR (short-range) parameters
# lambda_coop = params.lambda_coop
# r_div = params.r_div
# r_plat = params.r_plat

# # LR (long-range) parameters
# lambda_coopIn_LR1 = params.lambda_coopIn_LR1
# lambda_coopOut_LR1 = params.lambda_coopOut_LR1
# r_div_LR1 = params.r_div_LR1
# r_plat_LR1 = params.r_plat_LR1

# r_LR1 = params.r_LR1

# # unlikely to edit
# N_total_IDs = params.N_total_IDs
# n_cc = params.n_cc

# N_ExptState_Data = 30


GraphFolder = "Graphs"
filename_params_code = 'ExtremeAccns_output_'
filename_start = file_in_annotation[:-4] + '_' + filename_params_code

filename_ending = '.tsv'


# Drop this from simulatation dataframes
# And accessiosn: 
# Col0
# Can0
# Cvi0
# Dor10
# UKID116
missing_gene_in_data_ID = 'AT4G23000'


###  load in DatOnly files  ###

# RLFilt
N_ExptState_RLFilt = 21

filename_start_RLFilt = filename_start+'RLFilt_'
file_in_LocusProperties_RLFilt = os.path.join('Output_files', filename_start_RLFilt + 'LocusProperties'+filename_ending)
file_in_AllRepsMethFracs_RLFilt = os.path.join('Output_files', filename_start_RLFilt + 'AllRepsMethFracs'+filename_ending)
file_in_AllRepsXFracs_RLFilt = os.path.join('Output_files', filename_start_RLFilt + 'AllRepsXFracs'+filename_ending)

LocusProperties_RLFilt_df = pd.read_csv(file_in_LocusProperties_RLFilt, sep="\t{1}",engine='python')
LocusProperties_RLFilt_df.set_index("gene_ID", inplace = True)
N_gene = len(LocusProperties_RLFilt_df)
IDs_RLFilt_list = LocusProperties_RLFilt_df.index.values.tolist()

AllRepsMethFracs_RLFilt_df = pd.read_csv(file_in_AllRepsMethFracs_RLFilt, sep="\t{1}",engine='python')
AllRepsMethFracs_RLFilt_df.set_index("gene_ID", inplace = True)

AllRepsXFracs_RLFilt_df = pd.read_csv(file_in_AllRepsXFracs_RLFilt, sep="\t{1}",engine='python')
AllRepsXFracs_RLFilt_df.set_index("gene_ID", inplace = True)

print('RL', len(LocusProperties_RLFilt_df))


# NSFilt
N_ExptState_NSFilt = 34

filename_start_NSFilt = filename_start+'NSFilt_'
file_in_LocusProperties_NSFilt = os.path.join('Output_files', filename_start_NSFilt + 'LocusProperties'+filename_ending)
file_in_AllRepsMethFracs_NSFilt = os.path.join('Output_files', filename_start_NSFilt + 'AllRepsMethFracs'+filename_ending)
file_in_AllRepsXFracs_NSFilt = os.path.join('Output_files', filename_start_NSFilt + 'AllRepsXFracs'+filename_ending)

LocusProperties_NSFilt_df = pd.read_csv(file_in_LocusProperties_NSFilt, sep="\t{1}",engine='python')
LocusProperties_NSFilt_df.set_index("gene_ID", inplace = True)
N_gene = len(LocusProperties_NSFilt_df)
IDs_NSFilt_list = LocusProperties_NSFilt_df.index.values.tolist()

AllRepsMethFracs_NSFilt_df = pd.read_csv(file_in_AllRepsMethFracs_NSFilt, sep="\t{1}",engine='python')
AllRepsMethFracs_NSFilt_df.set_index("gene_ID", inplace = True)

AllRepsXFracs_NSFilt_df = pd.read_csv(file_in_AllRepsXFracs_NSFilt, sep="\t{1}",engine='python')
AllRepsXFracs_NSFilt_df.set_index("gene_ID", inplace = True)

print('NS', len(LocusProperties_NSFilt_df))


# Col0
N_ExptState_Col0 = 1

filename_start_Col0 = filename_start+'Col0_'
file_in_LocusProperties_Col0 = os.path.join('Output_files', filename_start_Col0 + 'LocusProperties'+filename_ending)

LocusProperties_Col0_df = pd.read_csv(file_in_LocusProperties_Col0, sep="\t{1}",engine='python')
LocusProperties_Col0_df.set_index("gene_ID", inplace = True)
N_gene = len(LocusProperties_Col0_df)
IDs_Col0_list = LocusProperties_Col0_df.index.values.tolist()

LocusProperties_Col0_df = LocusProperties_Col0_df.drop([missing_gene_in_data_ID])
print('Col0', len(LocusProperties_Col0_df))


# Can0
N_ExptState_Can0 = 1

filename_start_Can0 = filename_start+'Can0_'
file_in_LocusProperties_Can0 = os.path.join('Output_files', filename_start_Can0 + 'LocusProperties'+filename_ending)

LocusProperties_Can0_df = pd.read_csv(file_in_LocusProperties_Can0, sep="\t{1}",engine='python')
LocusProperties_Can0_df.set_index("gene_ID", inplace = True)
N_gene = len(LocusProperties_Can0_df)
IDs_Can0_list = LocusProperties_Can0_df.index.values.tolist()

LocusProperties_Can0_df = LocusProperties_Can0_df.drop([missing_gene_in_data_ID])
print('Can0', len(LocusProperties_Can0_df))


# Cvi0
N_ExptState_Cvi0 = 1

filename_start_Cvi0 = filename_start+'Cvi0_'
file_in_LocusProperties_Cvi0 = os.path.join('Output_files', filename_start_Cvi0 + 'LocusProperties'+filename_ending)

LocusProperties_Cvi0_df = pd.read_csv(file_in_LocusProperties_Cvi0, sep="\t{1}",engine='python')
LocusProperties_Cvi0_df.set_index("gene_ID", inplace = True)
N_gene = len(LocusProperties_Cvi0_df)
IDs_Cvi0_list = LocusProperties_Cvi0_df.index.values.tolist()

LocusProperties_Cvi0_df = LocusProperties_Cvi0_df.drop([missing_gene_in_data_ID])
print('Cvi0', len(LocusProperties_Cvi0_df))


# Dor10
N_ExptState_Dor10 = 1

filename_start_Dor10 = filename_start+'Dor10_'
file_in_LocusProperties_Dor10 = os.path.join('Output_files', filename_start_Dor10 + 'LocusProperties'+filename_ending)

LocusProperties_Dor10_df = pd.read_csv(file_in_LocusProperties_Dor10, sep="\t{1}",engine='python')
LocusProperties_Dor10_df.set_index("gene_ID", inplace = True)
N_gene = len(LocusProperties_Dor10_df)
IDs_Dor10_list = LocusProperties_Dor10_df.index.values.tolist()

LocusProperties_Dor10_df = LocusProperties_Dor10_df.drop([missing_gene_in_data_ID])
print('Dor10', len(LocusProperties_Dor10_df))


# UKID116
N_ExptState_UKID116 = 1

filename_start_UKID116 = filename_start+'UKID116_'
file_in_LocusProperties_UKID116 = os.path.join('Output_files', filename_start_UKID116 + 'LocusProperties'+filename_ending)

LocusProperties_UKID116_df = pd.read_csv(file_in_LocusProperties_UKID116, sep="\t{1}",engine='python')
LocusProperties_UKID116_df.set_index("gene_ID", inplace = True)
N_gene = len(LocusProperties_UKID116_df)
IDs_UKID116_list = LocusProperties_UKID116_df.index.values.tolist()

LocusProperties_UKID116_df = LocusProperties_UKID116_df.drop([missing_gene_in_data_ID])
print('UKID116', len(LocusProperties_UKID116_df))
print()

# ###  Load in simulations  ###

# CoopStrengthSweep
ParamCode_list_CoopStrengthSweep = ['Sim_Coop0p65_delta4p0_','Sim_Coop0p7_delta4p0_','Sim_Coop0p8_delta4p0_',
                                    'Sim_Coop0p9_delta4p0_','Sim_Coop1p0_delta4p0_','Sim_Coop1p1_delta4p0_',
                                    'Sim_Coop1p19_delta4p0_']
Labels_CoopStrengthSweep = ['$r^*$ = 0.65', '$r^*$ = 0.70', '$r^*$ = 0.80', 
                            '$r^*$ = 0.90', '$r^*$ = 1.00', '$r^*$ = 1.10', 
                            '$r^*$ = 1.19']
LocusProperties_CoopStrengthSweep_df_list = []
AllRepsMethFracs_CoopStrengthSweep_df_list = []

for i_ in range(len(ParamCode_list_CoopStrengthSweep)):
    filename_start_SimTemp = filename_start + ParamCode_list_CoopStrengthSweep[i_]
    file_in_LocusProperties_SimTemp = os.path.join('Output_files', filename_start_SimTemp + 'LocusProperties'+filename_ending)
    file_in_AllRepsMethFracs_SimTemp = os.path.join('Output_files', filename_start_SimTemp + 'AllRepsMethFracs'+filename_ending)

    LocusProperties_SimTemp_df = pd.read_csv(file_in_LocusProperties_SimTemp, sep="\t{1}",engine='python')
    LocusProperties_SimTemp_df.set_index("gene_ID", inplace = True)
    IDs_SimTemp_list = LocusProperties_SimTemp_df.index.values.tolist()

    AllRepsMethFracs_SimTemp_df = pd.read_csv(file_in_AllRepsMethFracs_SimTemp, sep="\t{1}",engine='python')
    AllRepsMethFracs_SimTemp_df.set_index("gene_ID", inplace = True)

    LocusProperties_SimTemp_df = LocusProperties_SimTemp_df.drop([missing_gene_in_data_ID])
    AllRepsMethFracs_SimTemp_df = AllRepsMethFracs_SimTemp_df.drop([missing_gene_in_data_ID])
    N_gene = len(LocusProperties_SimTemp_df)
    print(ParamCode_list_CoopStrengthSweep[i_], N_gene)

    LocusProperties_CoopStrengthSweep_df_list.append(LocusProperties_SimTemp_df)
    AllRepsMethFracs_CoopStrengthSweep_df_list.append(AllRepsMethFracs_SimTemp_df)
print()



# DeltaSweep
ParamCode_list_DeltaSweep = ['Sim_Coop1p0_delta0p01_','Sim_Coop1p0_delta0p1_','Sim_Coop1p0_delta1p0_',
                                    'Sim_Coop1p0_delta4p0_','Sim_Coop1p0_delta10p0_','Sim_Coop1p0_delta50p0_',
                                    'Sim_Coop1p0_delta100p0_']
# Labels_DeltaSweep = ['$r_0^+$ = 0.01', '$r_0^+$ = 0.1', '$r_0^+$ = 1.0', 
#                             '$r_0^+$ = 4.0', '$r_0^+$ = 10.0', '$r_0^+$ = 50.0', 
#                             '$r_0^+$ = 100']
Labels_DeltaSweep = ['$r_0^+$ = $1.0\\times10^{-8}$', '$r_0^+$ = $1.0\\times10^{-7}$', '$r_0^+$ = $1.0\\times10^{-6}$', 
                            '$r_0^+$ = $4.0\\times10^{-6}$', '$r_0^+$ = $1.0\\times10^{-5}$', '$r_0^+$ = $5.0\\times10^{-5}$', 
                            '$r_0^+$ = $1.0\\times10^{-4}$']


LocusProperties_DeltaSweep_df_list = []
AllRepsMethFracs_DeltaSweep_df_list = []

for i_ in range(len(ParamCode_list_DeltaSweep)):
    filename_start_SimTemp = filename_start + ParamCode_list_DeltaSweep[i_]
    file_in_LocusProperties_SimTemp = os.path.join('Output_files', filename_start_SimTemp + 'LocusProperties'+filename_ending)
    file_in_AllRepsMethFracs_SimTemp = os.path.join('Output_files', filename_start_SimTemp + 'AllRepsMethFracs'+filename_ending)

    LocusProperties_SimTemp_df = pd.read_csv(file_in_LocusProperties_SimTemp, sep="\t{1}",engine='python')
    LocusProperties_SimTemp_df.set_index("gene_ID", inplace = True)
    IDs_SimTemp_list = LocusProperties_SimTemp_df.index.values.tolist()

    AllRepsMethFracs_SimTemp_df = pd.read_csv(file_in_AllRepsMethFracs_SimTemp, sep="\t{1}",engine='python')
    AllRepsMethFracs_SimTemp_df.set_index("gene_ID", inplace = True)

    LocusProperties_SimTemp_df = LocusProperties_SimTemp_df.drop([missing_gene_in_data_ID])
    AllRepsMethFracs_SimTemp_df = AllRepsMethFracs_SimTemp_df.drop([missing_gene_in_data_ID])
    N_gene = len(LocusProperties_SimTemp_df)
    print(ParamCode_list_DeltaSweep[i_], N_gene)

    LocusProperties_DeltaSweep_df_list.append(LocusProperties_SimTemp_df)
    AllRepsMethFracs_DeltaSweep_df_list.append(AllRepsMethFracs_SimTemp_df)
print()



# Col0Fit
ParamCode_list_Col0Fit = ['Sim_Coop1p0_delta4p0_']
Labels_Col0Fit = ['Simulation']

LocusProperties_Col0Fit_df_list = []
AllRepsMethFracs_Col0Fit_df_list = []

for i_ in range(len(ParamCode_list_Col0Fit)):
    filename_start_SimTemp = filename_start + ParamCode_list_Col0Fit[i_]
    file_in_LocusProperties_SimTemp = os.path.join('Output_files', filename_start_SimTemp + 'LocusProperties'+filename_ending)
    file_in_AllRepsMethFracs_SimTemp = os.path.join('Output_files', filename_start_SimTemp + 'AllRepsMethFracs'+filename_ending)

    LocusProperties_SimTemp_df = pd.read_csv(file_in_LocusProperties_SimTemp, sep="\t{1}",engine='python')
    LocusProperties_SimTemp_df.set_index("gene_ID", inplace = True)
    IDs_SimTemp_list = LocusProperties_SimTemp_df.index.values.tolist()

    AllRepsMethFracs_SimTemp_df = pd.read_csv(file_in_AllRepsMethFracs_SimTemp, sep="\t{1}",engine='python')
    AllRepsMethFracs_SimTemp_df.set_index("gene_ID", inplace = True)

    LocusProperties_SimTemp_df = LocusProperties_SimTemp_df.drop([missing_gene_in_data_ID])
    AllRepsMethFracs_SimTemp_df = AllRepsMethFracs_SimTemp_df.drop([missing_gene_in_data_ID])
    N_gene = len(LocusProperties_SimTemp_df)
    print(ParamCode_list_Col0Fit[i_], N_gene)

    LocusProperties_Col0Fit_df_list.append(LocusProperties_SimTemp_df)
    AllRepsMethFracs_Col0Fit_df_list.append(AllRepsMethFracs_SimTemp_df)
print()


# Can0UKID116Cvi0_CoopStrength
ParamCode_list_Can0UKID116Cvi0_CoopStrength = ['Sim_Coop0p8_delta4p0_','Sim_Coop0p85_delta4p0_','Sim_Coop0p9_delta4p0_']
Labels_Can0UKID116Cvi0_CoopStrength = ['$r^*$ = 0.80', '$r^*$ = 0.85','$r^*$ = 0.90']
colors_list_Can0UKID116Cvi0_CoopStrength = ['m', 'b','lightseagreen']


LocusProperties_Can0UKID116Cvi0_CoopStrength_df_list = []
AllRepsMethFracs_Can0UKID116Cvi0_CoopStrength_df_list = []

for i_ in range(len(ParamCode_list_Can0UKID116Cvi0_CoopStrength)):
    filename_start_SimTemp = filename_start + ParamCode_list_Can0UKID116Cvi0_CoopStrength[i_]
    file_in_LocusProperties_SimTemp = os.path.join('Output_files', filename_start_SimTemp + 'LocusProperties'+filename_ending)
    file_in_AllRepsMethFracs_SimTemp = os.path.join('Output_files', filename_start_SimTemp + 'AllRepsMethFracs'+filename_ending)

    LocusProperties_SimTemp_df = pd.read_csv(file_in_LocusProperties_SimTemp, sep="\t{1}",engine='python')
    LocusProperties_SimTemp_df.set_index("gene_ID", inplace = True)
    IDs_SimTemp_list = LocusProperties_SimTemp_df.index.values.tolist()

    AllRepsMethFracs_SimTemp_df = pd.read_csv(file_in_AllRepsMethFracs_SimTemp, sep="\t{1}",engine='python')
    AllRepsMethFracs_SimTemp_df.set_index("gene_ID", inplace = True)

    LocusProperties_SimTemp_df = LocusProperties_SimTemp_df.drop([missing_gene_in_data_ID])
    AllRepsMethFracs_SimTemp_df = AllRepsMethFracs_SimTemp_df.drop([missing_gene_in_data_ID])
    N_gene = len(LocusProperties_SimTemp_df)
    print(ParamCode_list_Can0UKID116Cvi0_CoopStrength[i_], N_gene)

    LocusProperties_Can0UKID116Cvi0_CoopStrength_df_list.append(LocusProperties_SimTemp_df)
    AllRepsMethFracs_Can0UKID116Cvi0_CoopStrength_df_list.append(AllRepsMethFracs_SimTemp_df)
print()


# Can0UKID116Cvi0_Delta
ParamCode_list_Can0UKID116Cvi0_Delta = ['Sim_Coop1p0_delta0p1_','Sim_Coop1p0_delta1p0_']
Labels_Can0UKID116Cvi0_Delta = ['$r_0^+$ = 0.1','$r_0^+$ = 1.0']
colors_list_Can0UKID116Cvi0_Delta = ['m', 'b']

LocusProperties_Can0UKID116Cvi0_Delta_df_list = []
AllRepsMethFracs_Can0UKID116Cvi0_Delta_df_list = []

for i_ in range(len(ParamCode_list_Can0UKID116Cvi0_Delta)):
    filename_start_SimTemp = filename_start + ParamCode_list_Can0UKID116Cvi0_Delta[i_]
    file_in_LocusProperties_SimTemp = os.path.join('Output_files', filename_start_SimTemp + 'LocusProperties'+filename_ending)
    file_in_AllRepsMethFracs_SimTemp = os.path.join('Output_files', filename_start_SimTemp + 'AllRepsMethFracs'+filename_ending)

    LocusProperties_SimTemp_df = pd.read_csv(file_in_LocusProperties_SimTemp, sep="\t{1}",engine='python')
    LocusProperties_SimTemp_df.set_index("gene_ID", inplace = True)
    IDs_SimTemp_list = LocusProperties_SimTemp_df.index.values.tolist()

    AllRepsMethFracs_SimTemp_df = pd.read_csv(file_in_AllRepsMethFracs_SimTemp, sep="\t{1}",engine='python')
    AllRepsMethFracs_SimTemp_df.set_index("gene_ID", inplace = True)

    LocusProperties_SimTemp_df = LocusProperties_SimTemp_df.drop([missing_gene_in_data_ID])
    AllRepsMethFracs_SimTemp_df = AllRepsMethFracs_SimTemp_df.drop([missing_gene_in_data_ID])
    N_gene = len(LocusProperties_SimTemp_df)
    print(ParamCode_list_Can0UKID116Cvi0_Delta[i_], N_gene)

    LocusProperties_Can0UKID116Cvi0_Delta_df_list.append(LocusProperties_SimTemp_df)
    AllRepsMethFracs_Can0UKID116Cvi0_Delta_df_list.append(AllRepsMethFracs_SimTemp_df)
print()


# Can0UKID116Cvi0_Both
ParamCode_list_Can0UKID116Cvi0_Both = ['Sim_Coop0p85_delta4p0_','Sim_Coop0p92_delta1p0_','Sim_Coop1p0_delta0p1_']
# Labels_Can0UKID116Cvi0_Both = ['$r^*$ = 0.85   $r_0^+$ = 4.0','$r^*$ = 0.92   $r_0^+$ = 1.0','$r^*$ = 1.0   $r_0^+$ = 0.1']
Labels_Can0UKID116Cvi0_Both = ['$r^*$ = 0.85,   $r_0^+$ = $4.0\\times10^{-6}$','$r^*$ = 0.92,   $r_0^+$ = $1.0\\times10^{-6}$','$r^*$ = 1.00,   $r_0^+$ = $1.0\\times10^{-5}$']

colors_list_Can0UKID116Cvi0_Both = ['m', 'b','lightseagreen']

LocusProperties_Can0UKID116Cvi0_Both_df_list = []
AllRepsMethFracs_Can0UKID116Cvi0_Both_df_list = []

for i_ in range(len(ParamCode_list_Can0UKID116Cvi0_Both)):
    filename_start_SimTemp = filename_start + ParamCode_list_Can0UKID116Cvi0_Both[i_]
    file_in_LocusProperties_SimTemp = os.path.join('Output_files', filename_start_SimTemp + 'LocusProperties'+filename_ending)
    file_in_AllRepsMethFracs_SimTemp = os.path.join('Output_files', filename_start_SimTemp + 'AllRepsMethFracs'+filename_ending)

    LocusProperties_SimTemp_df = pd.read_csv(file_in_LocusProperties_SimTemp, sep="\t{1}",engine='python')
    LocusProperties_SimTemp_df.set_index("gene_ID", inplace = True)
    IDs_SimTemp_list = LocusProperties_SimTemp_df.index.values.tolist()

    AllRepsMethFracs_SimTemp_df = pd.read_csv(file_in_AllRepsMethFracs_SimTemp, sep="\t{1}",engine='python')
    AllRepsMethFracs_SimTemp_df.set_index("gene_ID", inplace = True)

    LocusProperties_SimTemp_df = LocusProperties_SimTemp_df.drop([missing_gene_in_data_ID])
    AllRepsMethFracs_SimTemp_df = AllRepsMethFracs_SimTemp_df.drop([missing_gene_in_data_ID])
    N_gene = len(LocusProperties_SimTemp_df)
    print(ParamCode_list_Can0UKID116Cvi0_Both[i_], N_gene)

    LocusProperties_Can0UKID116Cvi0_Both_df_list.append(LocusProperties_SimTemp_df)
    AllRepsMethFracs_Can0UKID116Cvi0_Both_df_list.append(AllRepsMethFracs_SimTemp_df)
print()


# Dor10_CoopStrength
ParamCode_list_Dor10_CoopStrength = ['Sim_Coop1p1_delta4p0_','Sim_Coop1p19_delta4p0_','Sim_Coop1p25_delta4p0_']
Labels_Dor10_CoopStrength = ['$r^*$ = 1.10', '$r^*$ = 1.19', '$r^{\dagger}$ = 1.25']
colors_list_Dor10_CoopStrength = ['limegreen', 'darkorange', 'r']

LocusProperties_Dor10_CoopStrength_df_list = []
AllRepsMethFracs_Dor10_CoopStrength_df_list = []

for i_ in range(len(ParamCode_list_Dor10_CoopStrength)):
    filename_start_SimTemp = filename_start + ParamCode_list_Dor10_CoopStrength[i_]
    file_in_LocusProperties_SimTemp = os.path.join('Output_files', filename_start_SimTemp + 'LocusProperties'+filename_ending)
    file_in_AllRepsMethFracs_SimTemp = os.path.join('Output_files', filename_start_SimTemp + 'AllRepsMethFracs'+filename_ending)

    LocusProperties_SimTemp_df = pd.read_csv(file_in_LocusProperties_SimTemp, sep="\t{1}",engine='python')
    LocusProperties_SimTemp_df.set_index("gene_ID", inplace = True)
    IDs_SimTemp_list = LocusProperties_SimTemp_df.index.values.tolist()

    AllRepsMethFracs_SimTemp_df = pd.read_csv(file_in_AllRepsMethFracs_SimTemp, sep="\t{1}",engine='python')
    AllRepsMethFracs_SimTemp_df.set_index("gene_ID", inplace = True)

    LocusProperties_SimTemp_df = LocusProperties_SimTemp_df.drop([missing_gene_in_data_ID])
    AllRepsMethFracs_SimTemp_df = AllRepsMethFracs_SimTemp_df.drop([missing_gene_in_data_ID])
    N_gene = len(LocusProperties_SimTemp_df)
    print(ParamCode_list_Dor10_CoopStrength[i_], N_gene)

    LocusProperties_Dor10_CoopStrength_df_list.append(LocusProperties_SimTemp_df)
    AllRepsMethFracs_Dor10_CoopStrength_df_list.append(AllRepsMethFracs_SimTemp_df)
print()


# Dor10_Delta
ParamCode_list_Dor10_Delta = ['Sim_Coop1p0_delta25p0_','Sim_Coop1p0_delta30p0_','Sim_Coop1p0_delta40p0_']
Labels_Dor10_Delta = ['$r_0^+$ = 25.0','$r_0^+$ = 30.0','$r_0^+$ = 40.0']
colors_list_Dor10_Delta = ['limegreen', 'darkorange', 'r']

LocusProperties_Dor10_Delta_df_list = []
AllRepsMethFracs_Dor10_Delta_df_list = []

for i_ in range(len(ParamCode_list_Dor10_Delta)):
    filename_start_SimTemp = filename_start + ParamCode_list_Dor10_Delta[i_]
    file_in_LocusProperties_SimTemp = os.path.join('Output_files', filename_start_SimTemp + 'LocusProperties'+filename_ending)
    file_in_AllRepsMethFracs_SimTemp = os.path.join('Output_files', filename_start_SimTemp + 'AllRepsMethFracs'+filename_ending)

    LocusProperties_SimTemp_df = pd.read_csv(file_in_LocusProperties_SimTemp, sep="\t{1}",engine='python')
    LocusProperties_SimTemp_df.set_index("gene_ID", inplace = True)
    IDs_SimTemp_list = LocusProperties_SimTemp_df.index.values.tolist()

    AllRepsMethFracs_SimTemp_df = pd.read_csv(file_in_AllRepsMethFracs_SimTemp, sep="\t{1}",engine='python')
    AllRepsMethFracs_SimTemp_df.set_index("gene_ID", inplace = True)

    LocusProperties_SimTemp_df = LocusProperties_SimTemp_df.drop([missing_gene_in_data_ID])
    AllRepsMethFracs_SimTemp_df = AllRepsMethFracs_SimTemp_df.drop([missing_gene_in_data_ID])
    N_gene = len(LocusProperties_SimTemp_df)
    print(ParamCode_list_Dor10_Delta[i_], N_gene)

    LocusProperties_Dor10_Delta_df_list.append(LocusProperties_SimTemp_df)
    AllRepsMethFracs_Dor10_Delta_df_list.append(AllRepsMethFracs_SimTemp_df)
print()

# Dor10_Both
ParamCode_list_Dor10_Both = ['Sim_Coop1p19_delta4p0_','Sim_Coop1p0_delta30p0_']
# Labels_Dor10_Both = ['$r^*$ = 1.19   $r_0^+$ = 4.0','$r^*$ = 1.00   $r_0^+$ = 30.0']
Labels_Dor10_Both = ['$r^*$ = 1.19,   $r_0^+$ = $4.0\\times10^{-6}$','$r^*$ = 1.00,   $r_0^+$ = $3.0\\times10^{-5}$']

colors_list_Dor10_Both = ['darkorange', 'r']

LocusProperties_Dor10_Both_df_list = []
AllRepsMethFracs_Dor10_Both_df_list = []

for i_ in range(len(ParamCode_list_Dor10_Both)):
    filename_start_SimTemp = filename_start + ParamCode_list_Dor10_Both[i_]
    file_in_LocusProperties_SimTemp = os.path.join('Output_files', filename_start_SimTemp + 'LocusProperties'+filename_ending)
    file_in_AllRepsMethFracs_SimTemp = os.path.join('Output_files', filename_start_SimTemp + 'AllRepsMethFracs'+filename_ending)

    LocusProperties_SimTemp_df = pd.read_csv(file_in_LocusProperties_SimTemp, sep="\t{1}",engine='python')
    LocusProperties_SimTemp_df.set_index("gene_ID", inplace = True)
    IDs_SimTemp_list = LocusProperties_SimTemp_df.index.values.tolist()

    AllRepsMethFracs_SimTemp_df = pd.read_csv(file_in_AllRepsMethFracs_SimTemp, sep="\t{1}",engine='python')
    AllRepsMethFracs_SimTemp_df.set_index("gene_ID", inplace = True)

    LocusProperties_SimTemp_df = LocusProperties_SimTemp_df.drop([missing_gene_in_data_ID])
    AllRepsMethFracs_SimTemp_df = AllRepsMethFracs_SimTemp_df.drop([missing_gene_in_data_ID])
    N_gene = len(LocusProperties_SimTemp_df)
    print(ParamCode_list_Dor10_Both[i_], N_gene)

    LocusProperties_Dor10_Both_df_list.append(LocusProperties_SimTemp_df)
    AllRepsMethFracs_Dor10_Both_df_list.append(AllRepsMethFracs_SimTemp_df)
print()



# NS_CoopStrength
ParamCode_list_NS_CoopStrength = ['Sim_Coop1p1_delta4p0_']
Labels_NS_CoopStrength = ['$r^*$ = 1.1']
colors_list_NS_CoopStrength = ['darkorange']

LocusProperties_NS_CoopStrength_df_list = []
AllRepsMethFracs_NS_CoopStrength_df_list = []

for i_ in range(len(ParamCode_list_NS_CoopStrength)):
    filename_start_SimTemp = filename_start + ParamCode_list_NS_CoopStrength[i_]
    file_in_LocusProperties_SimTemp = os.path.join('Output_files', filename_start_SimTemp + 'LocusProperties'+filename_ending)
    file_in_AllRepsMethFracs_SimTemp = os.path.join('Output_files', filename_start_SimTemp + 'AllRepsMethFracs'+filename_ending)

    LocusProperties_SimTemp_df = pd.read_csv(file_in_LocusProperties_SimTemp, sep="\t{1}",engine='python')
    LocusProperties_SimTemp_df.set_index("gene_ID", inplace = True)
    IDs_SimTemp_list = LocusProperties_SimTemp_df.index.values.tolist()

    AllRepsMethFracs_SimTemp_df = pd.read_csv(file_in_AllRepsMethFracs_SimTemp, sep="\t{1}",engine='python')
    AllRepsMethFracs_SimTemp_df.set_index("gene_ID", inplace = True)

    LocusProperties_SimTemp_df = LocusProperties_SimTemp_df.drop([missing_gene_in_data_ID])
    AllRepsMethFracs_SimTemp_df = AllRepsMethFracs_SimTemp_df.drop([missing_gene_in_data_ID])
    N_gene = len(LocusProperties_SimTemp_df)
    print(ParamCode_list_NS_CoopStrength[i_], N_gene)

    LocusProperties_NS_CoopStrength_df_list.append(LocusProperties_SimTemp_df)
    AllRepsMethFracs_NS_CoopStrength_df_list.append(AllRepsMethFracs_SimTemp_df)
print()

# NS_Delta
ParamCode_list_NS_Delta = ['Sim_Coop1p0_delta10p0_','Sim_Coop1p0_delta12p5_','Sim_Coop1p0_delta15p0_']
Labels_NS_Delta = ['$r_0^+$ = 10.0','$r_0^+$ = 12.5','$r_0^+$ = 15.0']
colors_list_NS_Delta = ['limegreen', 'darkorange', 'r']

LocusProperties_NS_Delta_df_list = []
AllRepsMethFracs_NS_Delta_df_list = []

for i_ in range(len(ParamCode_list_NS_Delta)):
    filename_start_SimTemp = filename_start + ParamCode_list_NS_Delta[i_]
    file_in_LocusProperties_SimTemp = os.path.join('Output_files', filename_start_SimTemp + 'LocusProperties'+filename_ending)
    file_in_AllRepsMethFracs_SimTemp = os.path.join('Output_files', filename_start_SimTemp + 'AllRepsMethFracs'+filename_ending)

    LocusProperties_SimTemp_df = pd.read_csv(file_in_LocusProperties_SimTemp, sep="\t{1}",engine='python')
    LocusProperties_SimTemp_df.set_index("gene_ID", inplace = True)
    IDs_SimTemp_list = LocusProperties_SimTemp_df.index.values.tolist()

    AllRepsMethFracs_SimTemp_df = pd.read_csv(file_in_AllRepsMethFracs_SimTemp, sep="\t{1}",engine='python')
    AllRepsMethFracs_SimTemp_df.set_index("gene_ID", inplace = True)

    LocusProperties_SimTemp_df = LocusProperties_SimTemp_df.drop([missing_gene_in_data_ID])
    AllRepsMethFracs_SimTemp_df = AllRepsMethFracs_SimTemp_df.drop([missing_gene_in_data_ID])
    N_gene = len(LocusProperties_SimTemp_df)
    print(ParamCode_list_NS_Delta[i_], N_gene)

    LocusProperties_NS_Delta_df_list.append(LocusProperties_SimTemp_df)
    AllRepsMethFracs_NS_Delta_df_list.append(AllRepsMethFracs_SimTemp_df)
print()

# NS_Both
ParamCode_list_NS_Both = ['Sim_Coop1p1_delta4p0_','Sim_Coop1p0_delta12p5_']
# Labels_NS_Both = ['$r^*$ = 1.1   $r_0^+$ = 4.0','$r^*$ = 1.0   $r_0^+$ = 12.5']
Labels_NS_Both = ['$r^*$ = 1.1,   $r_0^+$ = $4.00\\times10^{-6}$','$r^*$ = 1.0,   $r_0^+$ = $1.25\\times10^{-5}$']

colors_list_NS_Both = ['darkorange', 'r' ]

LocusProperties_NS_Both_df_list = []
AllRepsMethFracs_NS_Both_df_list = []

for i_ in range(len(ParamCode_list_NS_Both)):
    filename_start_SimTemp = filename_start + ParamCode_list_NS_Both[i_]
    file_in_LocusProperties_SimTemp = os.path.join('Output_files', filename_start_SimTemp + 'LocusProperties'+filename_ending)
    file_in_AllRepsMethFracs_SimTemp = os.path.join('Output_files', filename_start_SimTemp + 'AllRepsMethFracs'+filename_ending)

    LocusProperties_SimTemp_df = pd.read_csv(file_in_LocusProperties_SimTemp, sep="\t{1}",engine='python')
    LocusProperties_SimTemp_df.set_index("gene_ID", inplace = True)
    IDs_SimTemp_list = LocusProperties_SimTemp_df.index.values.tolist()

    AllRepsMethFracs_SimTemp_df = pd.read_csv(file_in_AllRepsMethFracs_SimTemp, sep="\t{1}",engine='python')
    AllRepsMethFracs_SimTemp_df.set_index("gene_ID", inplace = True)

    LocusProperties_SimTemp_df = LocusProperties_SimTemp_df.drop([missing_gene_in_data_ID])
    AllRepsMethFracs_SimTemp_df = AllRepsMethFracs_SimTemp_df.drop([missing_gene_in_data_ID])
    N_gene = len(LocusProperties_SimTemp_df)
    print(ParamCode_list_NS_Both[i_], N_gene)

    LocusProperties_NS_Both_df_list.append(LocusProperties_SimTemp_df)
    AllRepsMethFracs_NS_Both_df_list.append(AllRepsMethFracs_SimTemp_df)
print()


# RL_CoopStrength
ParamCode_list_RL_CoopStrength = ['Sim_Coop0p96_delta4p0_','Sim_Coop0p98_delta4p0_']
Labels_RL_CoopStrength = ['$r^*$ = 0.96','$r^*$ = 0.98']
colors_list_RL_CoopStrength = ['green', 'limegreen']

LocusProperties_RL_CoopStrength_df_list = []
AllRepsMethFracs_RL_CoopStrength_df_list = []

for i_ in range(len(ParamCode_list_RL_CoopStrength)):
    filename_start_SimTemp = filename_start + ParamCode_list_RL_CoopStrength[i_]
    file_in_LocusProperties_SimTemp = os.path.join('Output_files', filename_start_SimTemp + 'LocusProperties'+filename_ending)
    file_in_AllRepsMethFracs_SimTemp = os.path.join('Output_files', filename_start_SimTemp + 'AllRepsMethFracs'+filename_ending)

    LocusProperties_SimTemp_df = pd.read_csv(file_in_LocusProperties_SimTemp, sep="\t{1}",engine='python')
    LocusProperties_SimTemp_df.set_index("gene_ID", inplace = True)
    IDs_SimTemp_list = LocusProperties_SimTemp_df.index.values.tolist()

    AllRepsMethFracs_SimTemp_df = pd.read_csv(file_in_AllRepsMethFracs_SimTemp, sep="\t{1}",engine='python')
    AllRepsMethFracs_SimTemp_df.set_index("gene_ID", inplace = True)

    LocusProperties_SimTemp_df = LocusProperties_SimTemp_df.drop([missing_gene_in_data_ID])
    AllRepsMethFracs_SimTemp_df = AllRepsMethFracs_SimTemp_df.drop([missing_gene_in_data_ID])
    N_gene = len(LocusProperties_SimTemp_df)
    print(ParamCode_list_RL_CoopStrength[i_], N_gene)

    LocusProperties_RL_CoopStrength_df_list.append(LocusProperties_SimTemp_df)
    AllRepsMethFracs_RL_CoopStrength_df_list.append(AllRepsMethFracs_SimTemp_df)
print()

# RL_Delta
ParamCode_list_RL_Delta = ['Sim_Coop1p0_delta1p0_','Sim_Coop1p0_delta2p5_','Sim_Coop1p0_delta4p0_']
Labels_RL_Delta = ['$r_0^+$ = 1.0','$r_0^+$ = 2.5','$r_0^+$ = 4.0']
colors_list_RL_Delta = ['lightseagreen','green', 'limegreen']

LocusProperties_RL_Delta_df_list = []
AllRepsMethFracs_RL_Delta_df_list = []

for i_ in range(len(ParamCode_list_RL_Delta)):
    filename_start_SimTemp = filename_start + ParamCode_list_RL_Delta[i_]
    file_in_LocusProperties_SimTemp = os.path.join('Output_files', filename_start_SimTemp + 'LocusProperties'+filename_ending)
    file_in_AllRepsMethFracs_SimTemp = os.path.join('Output_files', filename_start_SimTemp + 'AllRepsMethFracs'+filename_ending)

    LocusProperties_SimTemp_df = pd.read_csv(file_in_LocusProperties_SimTemp, sep="\t{1}",engine='python')
    LocusProperties_SimTemp_df.set_index("gene_ID", inplace = True)
    IDs_SimTemp_list = LocusProperties_SimTemp_df.index.values.tolist()

    AllRepsMethFracs_SimTemp_df = pd.read_csv(file_in_AllRepsMethFracs_SimTemp, sep="\t{1}",engine='python')
    AllRepsMethFracs_SimTemp_df.set_index("gene_ID", inplace = True)

    LocusProperties_SimTemp_df = LocusProperties_SimTemp_df.drop([missing_gene_in_data_ID])
    AllRepsMethFracs_SimTemp_df = AllRepsMethFracs_SimTemp_df.drop([missing_gene_in_data_ID])
    N_gene = len(LocusProperties_SimTemp_df)
    print(ParamCode_list_RL_Delta[i_], N_gene)

    LocusProperties_RL_Delta_df_list.append(LocusProperties_SimTemp_df)
    AllRepsMethFracs_RL_Delta_df_list.append(AllRepsMethFracs_SimTemp_df)
print()

# RL_Both
ParamCode_list_RL_Both = ['Sim_Coop0p98_delta4p0_','Sim_Coop1p0_delta2p5_']
# Labels_RL_Both = ['$r^*$ = 0.98   $r_0^+$ = 4.0','$r^*$ = 1.00   $r_0^+$ = 2.5']
Labels_RL_Both = ['$r^*$ = 0.98,   $r_0^+$ = $4.0\\times10^{-6}$','$r^*$ = 1.00,   $r_0^+$ = $2.5\\times10^{-6}$']

colors_list_RL_Both = ['lightseagreen', 'm']

LocusProperties_RL_Both_df_list = []
AllRepsMethFracs_RL_Both_df_list = []

for i_ in range(len(ParamCode_list_RL_Both)):
    filename_start_SimTemp = filename_start + ParamCode_list_RL_Both[i_]
    file_in_LocusProperties_SimTemp = os.path.join('Output_files', filename_start_SimTemp + 'LocusProperties'+filename_ending)
    file_in_AllRepsMethFracs_SimTemp = os.path.join('Output_files', filename_start_SimTemp + 'AllRepsMethFracs'+filename_ending)

    LocusProperties_SimTemp_df = pd.read_csv(file_in_LocusProperties_SimTemp, sep="\t{1}",engine='python')
    LocusProperties_SimTemp_df.set_index("gene_ID", inplace = True)
    IDs_SimTemp_list = LocusProperties_SimTemp_df.index.values.tolist()

    AllRepsMethFracs_SimTemp_df = pd.read_csv(file_in_AllRepsMethFracs_SimTemp, sep="\t{1}",engine='python')
    AllRepsMethFracs_SimTemp_df.set_index("gene_ID", inplace = True)

    LocusProperties_SimTemp_df = LocusProperties_SimTemp_df.drop([missing_gene_in_data_ID])
    AllRepsMethFracs_SimTemp_df = AllRepsMethFracs_SimTemp_df.drop([missing_gene_in_data_ID])
    N_gene = len(LocusProperties_SimTemp_df)
    print(ParamCode_list_RL_Both[i_], N_gene)

    LocusProperties_RL_Both_df_list.append(LocusProperties_SimTemp_df)
    AllRepsMethFracs_RL_Both_df_list.append(AllRepsMethFracs_SimTemp_df)
print()


linewidth_1 = 1
linewidth_2_Dat = 2
linewidth_2_Sim = 2

colors_list_sweep = ['m', 'b', 'lightseagreen', 'green', 'limegreen', 'darkorange', 'r']



n_bin = 50
# convert to percentage
values_scaler = 100


### Summary Histograms ###


# Start Graph
fig, ax = plt.subplots(1,1,figsize=(4,4))

for spine in ['left','right','top','bottom']:
    ax.spines[spine].set_color('k')
    ax.spines[spine].set_linewidth(0.8)
ax.set_facecolor('white')

#ax.grid(False)
ax.grid(b=True, which='major', color='lightgrey', linestyle=':',linewidth=1)

x_start = 0.
x_end = 100. 
x_hist = np.arange(x_start+x_end/(2.*n_bin),x_end,x_end/n_bin)

# temp_cols_list = []
# for i_temp in range(N_reps):
#     temp_cols_list.append('Sim_meth_frac_'+str(i_temp))
# hist_Sim, bins_Sim = np.histogram(AllRepsMethFracs_df[temp_cols_list].values, bins=n_bin, range=(0,1))
# hist_Sim = hist_Sim/N_reps # normalise
# ax[0].axvline(np.nanmean(AllRepsMethFracs_df[temp_cols_list].values), linestyle='--', color='r')
# ax[1].axvline(np.nanmean(AllRepsMethFracs_df[temp_cols_list].values), linestyle='--', color='r')
# ax[0].plot(x_hist, hist_Sim,  linewidth = 2, color='r',
#                   label='N_bins = %d\nSim. $\mu$ = %.0f' % (n_bin,np.nanmean(AllRepsMethFracs_df[temp_cols_list].values)))
# ax[1].plot(x_hist, hist_Sim,  linewidth = 2, color='r',
#                   label='Normalised\nSim. $\mu$ = %.0f' % np.nanmean(AllRepsMethFracs_df[temp_cols_list].values))




# Cvi0
hist_Cvi0, bins_Cvi0 = np.histogram(LocusProperties_Cvi0_df[['Col0_meth_frac']].values*values_scaler, bins=n_bin, range=(x_start,x_end))
# ax.axvline(np.nanmean(LocusProperties_Cvi0_df[['Col0_meth_frac']].values*values_scaler), linestyle='--', color=colors_list_sweep[0],  linewidth = linewidth_1)
ax.plot(x_hist, hist_Cvi0,  linewidth = linewidth_1, color=colors_list_sweep[0],
                  label='Cvi0: Mean = %.0f %%' % np.nanmean(LocusProperties_Cvi0_df[['Col0_meth_frac']].values*values_scaler))


# UKID116
hist_UKID116, bins_UKID116 = np.histogram(LocusProperties_UKID116_df[['Col0_meth_frac']].values*values_scaler, bins=n_bin, range=(x_start,x_end))
# ax.axvline(np.nanmean(LocusProperties_UKID116_df[['Col0_meth_frac']].values*values_scaler), linestyle='--', color=colors_list_sweep[1],  linewidth = linewidth_1)
ax.plot(x_hist, hist_UKID116,  linewidth = linewidth_1, color=colors_list_sweep[1],
                  label='UKID116: Mean = %.0f %%' % np.nanmean(LocusProperties_UKID116_df[['Col0_meth_frac']].values*values_scaler))


# Can0
hist_Can0, bins_Can0 = np.histogram(LocusProperties_Can0_df[['Col0_meth_frac']].values*values_scaler, bins=n_bin, range=(x_start,x_end))
# ax.axvline(np.nanmean(LocusProperties_Can0_df[['Col0_meth_frac']].values*values_scaler), linestyle='--', color=colors_list_sweep[2],  linewidth = linewidth_1)
ax.plot(x_hist, hist_Can0,  linewidth = linewidth_1, color=colors_list_sweep[2],
                  label='Can0: Mean = %.0f %%' % np.nanmean(LocusProperties_Can0_df[['Col0_meth_frac']].values*values_scaler))


# RLFilt
temp_cols_list = []
for i_temp in range(N_ExptState_RLFilt):
    temp_cols_list.append('Dset1_meth_frac_'+str(i_temp))
hist_RLFilt, bins_RLFilt = np.histogram(AllRepsMethFracs_RLFilt_df[temp_cols_list].values*values_scaler, bins=n_bin, range=(x_start,x_end))
hist_RLFilt = hist_RLFilt/N_ExptState_RLFilt # normalise
# ax.axvline(np.nanmean(AllRepsMethFracs_RLFilt_df[temp_cols_list].values*values_scaler), linestyle='--', color=colors_list_sweep[3],  linewidth = linewidth_1)
ax.plot(x_hist, hist_RLFilt,  linewidth = linewidth_1, color=colors_list_sweep[3],
                  label='RL: Mean = %.0f %%' % np.nanmean(AllRepsMethFracs_RLFilt_df[temp_cols_list].values*values_scaler))


# Col0
hist_Col0, bins_Col0 = np.histogram(LocusProperties_Col0_df[['Col0_meth_frac']].values*values_scaler, bins=n_bin, range=(x_start,x_end))
# ax.axvline(np.nanmean(LocusProperties_Col0_df[['Col0_meth_frac']].values*values_scaler), linestyle='--', color=colors_list_sweep[4],  linewidth = linewidth_1)
ax.plot(x_hist, hist_Col0,  linewidth = linewidth_1, color=colors_list_sweep[4],
                  label='Col0: Mean = %.0f %%' % np.nanmean(LocusProperties_Col0_df[['Col0_meth_frac']].values*values_scaler))


# NSFilt
temp_cols_list = []
for i_temp in range(N_ExptState_NSFilt):
    temp_cols_list.append('Dset1_meth_frac_'+str(i_temp))
hist_NSFilt, bins_NSFilt = np.histogram(AllRepsMethFracs_NSFilt_df[temp_cols_list].values*values_scaler, bins=n_bin, range=(x_start,x_end))
hist_NSFilt = hist_NSFilt/N_ExptState_NSFilt # normalise
# ax.axvline(np.nanmean(AllRepsMethFracs_NSFilt_df[temp_cols_list].values*values_scaler), linestyle='--', color=colors_list_sweep[5],  linewidth = linewidth_1)
ax.plot(x_hist, hist_NSFilt,  linewidth = linewidth_1, color=colors_list_sweep[5],
                  label='NS: Mean = %.0f %%' % np.nanmean(AllRepsMethFracs_NSFilt_df[temp_cols_list].values*values_scaler))


# Dor10
hist_Dor10, bins_Dor10 = np.histogram(LocusProperties_Dor10_df[['Col0_meth_frac']].values*values_scaler, bins=n_bin, range=(x_start,x_end))
# ax.axvline(np.nanmean(LocusProperties_Dor10_df[['Col0_meth_frac']].values*values_scaler), linestyle='--', color=colors_list_sweep[6],  linewidth = linewidth_1)
ax.plot(x_hist, hist_Dor10,  linewidth = linewidth_1, color=colors_list_sweep[6],
                  label='Dor10: Mean = %.0f %%' % np.nanmean(LocusProperties_Dor10_df[['Col0_meth_frac']].values*values_scaler))



ax.legend(framealpha=1)

ax.set_xlim(-10,110)

ax.set_ylim(0,2500)

ax.xaxis.set_tick_params(labelsize=12)
ax.yaxis.set_tick_params(labelsize=12)


ax.set_xlabel("Methylation percentage", fontsize=12)
ax.set_ylabel("Number of loci", fontsize=12)

#fig.tight_layout()
# fig.subplots_adjust(top=0.72)

fig.savefig( os.path.join(GraphFolder,filename_start +"ExtemeAccns_AllData.png"), bbox_inches = 'tight')
fig.savefig( os.path.join(GraphFolder,filename_start +"ExtemeAccns_AllData.pdf"), bbox_inches = 'tight')

plt.show()
# End Graph

############################

# Start Graph
fig, ax = plt.subplots(1,1,figsize=(4,4))

for spine in ['left','right','top','bottom']:
    ax.spines[spine].set_color('k')
    ax.spines[spine].set_linewidth(0.8)
ax.set_facecolor('white')

#ax.grid(False)
ax.grid(b=True, which='major', color='lightgrey', linestyle=':',linewidth=1)

x_start = 0.
x_end = 100. 
x_hist = np.arange(x_start+x_end/(2.*n_bin),x_end,x_end/n_bin)

# temp_cols_list = []
# for i_temp in range(N_reps):
#     temp_cols_list.append('Sim_meth_frac_'+str(i_temp))
# hist_Sim, bins_Sim = np.histogram(AllRepsMethFracs_df[temp_cols_list].values, bins=n_bin, range=(0,1))
# hist_Sim = hist_Sim/N_reps # normalise
# ax[0].axvline(np.nanmean(AllRepsMethFracs_df[temp_cols_list].values), linestyle='--', color='r')
# ax[1].axvline(np.nanmean(AllRepsMethFracs_df[temp_cols_list].values), linestyle='--', color='r')
# ax[0].plot(x_hist, hist_Sim,  linewidth = 2, color='r',
#                   label='N_bins = %d\nSim. $\mu$ = %.0f' % (n_bin,np.nanmean(AllRepsMethFracs_df[temp_cols_list].values)))
# ax[1].plot(x_hist, hist_Sim,  linewidth = 2, color='r',
#                   label='Normalised\nSim. $\mu$ = %.0f' % np.nanmean(AllRepsMethFracs_df[temp_cols_list].values))


# Cvi0
hist_Cvi0, bins_Cvi0 = np.histogram(LocusProperties_Cvi0_df[['Col0_X_frac']].values*values_scaler, bins=n_bin, range=(x_start,x_end))
# ax.axvline(np.nanmean(LocusProperties_Cvi0_df[['Col0_X_frac']].values*values_scaler), linestyle='--', color=colors_list_sweep[0],  linewidth = linewidth_1)
ax.plot(x_hist, hist_Cvi0,  linewidth = linewidth_1, color=colors_list_sweep[0],
                  label='Cvi0: Mean = %.0f %%' % np.nanmean(LocusProperties_Cvi0_df[['Col0_X_frac']].values*values_scaler))


# UKID116
hist_UKID116, bins_UKID116 = np.histogram(LocusProperties_UKID116_df[['Col0_X_frac']].values*values_scaler, bins=n_bin, range=(x_start,x_end))
# ax.axvline(np.nanmean(LocusProperties_UKID116_df[['Col0_X_frac']].values*values_scaler), linestyle='--', color=colors_list_sweep[1],  linewidth = linewidth_1)
ax.plot(x_hist, hist_UKID116,  linewidth = linewidth_1, color=colors_list_sweep[1],
                  label='UKID116: Mean = %.0f %%' % np.nanmean(LocusProperties_UKID116_df[['Col0_X_frac']].values*values_scaler))


# Can0
hist_Can0, bins_Can0 = np.histogram(LocusProperties_Can0_df[['Col0_X_frac']].values*values_scaler, bins=n_bin, range=(x_start,x_end))
# ax.axvline(np.nanmean(LocusProperties_Can0_df[['Col0_X_frac']].values*values_scaler), linestyle='--', color=colors_list_sweep[2],  linewidth = linewidth_1)
ax.plot(x_hist, hist_Can0,  linewidth = linewidth_1, color=colors_list_sweep[2],
                  label='Can0: Mean = %.0f %%' % np.nanmean(LocusProperties_Can0_df[['Col0_X_frac']].values*values_scaler))


# RLFilt
temp_cols_list = []
for i_temp in range(N_ExptState_RLFilt):
    temp_cols_list.append('Dset1_meth_frac_'+str(i_temp))
hist_RLFilt, bins_RLFilt = np.histogram(AllRepsXFracs_RLFilt_df[temp_cols_list].values*values_scaler, bins=n_bin, range=(x_start,x_end))
hist_RLFilt = hist_RLFilt/N_ExptState_RLFilt # normalise
# ax.axvline(np.nanmean(AllRepsXFracs_RLFilt_df[temp_cols_list].values*values_scaler), linestyle='--', color=colors_list_sweep[3],  linewidth = linewidth_1)
ax.plot(x_hist, hist_RLFilt,  linewidth = linewidth_1, color=colors_list_sweep[3],
                  label='RL: Mean = %.0f %%' % np.nanmean(AllRepsXFracs_RLFilt_df[temp_cols_list].values*values_scaler))


# Col0
hist_Col0, bins_Col0 = np.histogram(LocusProperties_Col0_df[['Col0_X_frac']].values*values_scaler, bins=n_bin, range=(x_start,x_end))
# ax.axvline(np.nanmean(LocusProperties_Col0_df[['Col0_X_frac']].values*values_scaler), linestyle='--', color=colors_list_sweep[4],  linewidth = linewidth_1)
ax.plot(x_hist, hist_Col0,  linewidth = linewidth_1, color=colors_list_sweep[4],
                  label='Col0: Mean = %.0f %%' % np.nanmean(LocusProperties_Col0_df[['Col0_X_frac']].values*values_scaler))


# NSFilt
temp_cols_list = []
for i_temp in range(N_ExptState_NSFilt):
    temp_cols_list.append('Dset1_meth_frac_'+str(i_temp))
hist_NSFilt, bins_NSFilt = np.histogram(AllRepsXFracs_NSFilt_df[temp_cols_list].values*values_scaler, bins=n_bin, range=(x_start,x_end))
hist_NSFilt = hist_NSFilt/N_ExptState_NSFilt # normalise
# ax.axvline(np.nanmean(AllRepsXFracs_NSFilt_df[temp_cols_list].values*values_scaler), linestyle='--', color=colors_list_sweep[5],  linewidth = linewidth_1)
ax.plot(x_hist, hist_NSFilt,  linewidth = linewidth_1, color=colors_list_sweep[5],
                  label='NS: Mean = %.0f %%' % np.nanmean(AllRepsXFracs_NSFilt_df[temp_cols_list].values*values_scaler))


# Dor10
hist_Dor10, bins_Dor10 = np.histogram(LocusProperties_Dor10_df[['Col0_X_frac']].values*values_scaler, bins=n_bin, range=(x_start,x_end))
# ax.axvline(np.nanmean(LocusProperties_Dor10_df[['Col0_X_frac']].values*values_scaler), linestyle='--', color=colors_list_sweep[6],  linewidth = linewidth_1)
ax.plot(x_hist, hist_Dor10,  linewidth = linewidth_1, color=colors_list_sweep[6],
                  label='Dor10: Mean = %.0f %%' % np.nanmean(LocusProperties_Dor10_df[['Col0_X_frac']].values*values_scaler))



ax.legend(framealpha=1)

ax.set_xlim(-10,110)

ax.set_ylim(0,2500)

ax.xaxis.set_tick_params(labelsize=12)
ax.yaxis.set_tick_params(labelsize=12)


ax.set_xlabel("$X$-site percentage", fontsize=12)
ax.set_ylabel("Number of loci", fontsize=12)

#fig.tight_layout()
# fig.subplots_adjust(top=0.72)

fig.savefig( os.path.join(GraphFolder,filename_start +"ExtemeAccns_AllData_Xfrac.png"), bbox_inches = 'tight')

plt.show()
# End Graph



###################


# Start Graph
fig, ax = plt.subplots(1,1,figsize=(4,4))

for spine in ['left','right','top','bottom']:
    ax.spines[spine].set_color('k')
    ax.spines[spine].set_linewidth(0.8)
ax.set_facecolor('white')

#ax.grid(False)
ax.grid(b=True, which='major', color='lightgrey', linestyle=':',linewidth=1)

x_start = 0.
x_end = 100. 
x_hist = np.arange(x_start+x_end/(2.*n_bin),x_end,x_end/n_bin)


# loop over params. 
temp_cols_list = []
for i_temp in range(N_reps):
    temp_cols_list.append('Sim_meth_frac_'+str(i_temp))

for i_ in range(len(AllRepsMethFracs_CoopStrengthSweep_df_list)):
    temp_df = AllRepsMethFracs_CoopStrengthSweep_df_list[i_]
    hist_temp, bins_temp = np.histogram(temp_df[temp_cols_list].values*values_scaler, bins=n_bin, range=(x_start,x_end))
    hist_temp = hist_temp/N_reps # normalise
    # ax.axvline(np.nanmean(temp_df[temp_cols_list].values*values_scaler), linestyle='--', color=colors_list_sweep[i_],linewidth = linewidth_1)
    ax.plot(x_hist, hist_temp,  linewidth = linewidth_1, color=colors_list_sweep[i_],
                    label='%s: Mean = %.0f %%' % (Labels_CoopStrengthSweep[i_],np.nanmean(temp_df[temp_cols_list].values*values_scaler)) )

ax.legend(framealpha=1)

ax.set_xlim(-10,110)
ax.set_ylim(0,2500)

ax.xaxis.set_tick_params(labelsize=12)
ax.yaxis.set_tick_params(labelsize=12)

ax.set_xlabel("Methylation percentage", fontsize=12)
ax.set_ylabel("Number of loci", fontsize=12)

#fig.tight_layout()
# fig.subplots_adjust(top=0.72)

fig.savefig( os.path.join(GraphFolder,filename_start +"ExtemeAccns_CoopStrengthSweep.png"), bbox_inches = 'tight')
fig.savefig( os.path.join(GraphFolder,filename_start +"ExtemeAccns_CoopStrengthSweep.pdf"), bbox_inches = 'tight')

plt.show()
# End Graph

#########################

# Start Graph
fig, ax = plt.subplots(1,1,figsize=(4,4))

for spine in ['left','right','top','bottom']:
    ax.spines[spine].set_color('k')
    ax.spines[spine].set_linewidth(0.8)
ax.set_facecolor('white')

#ax.grid(False)
ax.grid(b=True, which='major', color='lightgrey', linestyle=':',linewidth=1)

x_start = 0.
x_end = 100. 
x_hist = np.arange(x_start+x_end/(2.*n_bin),x_end,x_end/n_bin)


# loop over params. 
temp_cols_list = []
for i_temp in range(N_reps):
    temp_cols_list.append('Sim_meth_frac_'+str(i_temp))

for i_ in range(len(AllRepsMethFracs_DeltaSweep_df_list)):
    temp_df = AllRepsMethFracs_DeltaSweep_df_list[i_]
    hist_temp, bins_temp = np.histogram(temp_df[temp_cols_list].values*values_scaler, bins=n_bin, range=(x_start,x_end))
    hist_temp = hist_temp/N_reps # normalise
    # ax.axvline(np.nanmean(temp_df[temp_cols_list].values*values_scaler), linestyle='--', color=colors_list_sweep[i_],linewidth = linewidth_1)
    ax.plot(x_hist, hist_temp,  linewidth = linewidth_1, color=colors_list_sweep[i_],
                    label='%s: Mean = %.0f %%' % (Labels_DeltaSweep[i_],np.nanmean(temp_df[temp_cols_list].values*values_scaler)) )

ax.legend(framealpha=1, fontsize=9)

ax.set_xlim(-10,110)
ax.set_ylim(0,2500)

ax.xaxis.set_tick_params(labelsize=12)
ax.yaxis.set_tick_params(labelsize=12)

ax.set_xlabel("Methylation percentage", fontsize=12)
ax.set_ylabel("Number of loci", fontsize=12)

#fig.tight_layout()
# fig.subplots_adjust(top=0.72)

fig.savefig( os.path.join(GraphFolder,filename_start +"ExtemeAccns_DeltaSweep.png"), bbox_inches = 'tight')
fig.savefig( os.path.join(GraphFolder,filename_start +"ExtemeAccns_DeltaSweep.pdf"), bbox_inches = 'tight')

plt.show()
# End Graph



##############################

# Start Graph
fig, ax = plt.subplots(1,1,figsize=(4,4))

for spine in ['left','right','top','bottom']:
    ax.spines[spine].set_color('k')
    ax.spines[spine].set_linewidth(0.8)
ax.set_facecolor('white')

#ax.grid(False)
ax.grid(b=True, which='major', color='lightgrey', linestyle=':',linewidth=1)

x_start = 0.
x_end = 100. 
x_hist = np.arange(x_start+x_end/(2.*n_bin),x_end,x_end/n_bin)

# Col0
hist_Col0, bins_Col0 = np.histogram(LocusProperties_Col0_df[['Col0_meth_frac']].values*values_scaler, bins=n_bin, range=(x_start,x_end))
# ax.axvline(np.nanmean(LocusProperties_Col0_df[['Col0_X_frac']].values*values_scaler), linestyle='--', color='k',  linewidth = linewidth_2_Dat)
ax.plot(x_hist, hist_Col0,  linewidth = linewidth_2_Dat, color='k',
                  label='Col0' )

# Single Sim
temp_cols_list = []
for i_temp in range(N_reps):
    temp_cols_list.append('Sim_meth_frac_'+str(i_temp))

for i_ in range(len(AllRepsMethFracs_Col0Fit_df_list)):
    temp_df = AllRepsMethFracs_Col0Fit_df_list[i_]
    hist_temp, bins_temp = np.histogram(temp_df[temp_cols_list].values*values_scaler, bins=n_bin, range=(x_start,x_end))
    hist_temp = hist_temp/N_reps # normalise
    # ax.axvline(np.nanmean(temp_df[temp_cols_list].values*values_scaler), linestyle='--', color='green',linewidth = linewidth_2_Sim)
    ax.plot(x_hist, hist_temp,  linewidth = linewidth_2_Sim, color='green',
                    label='%s' % (Labels_Col0Fit[i_]) )

ax.legend(framealpha=1)

ax.set_xlim(-10,110)
ax.set_ylim(0,1000)

ax.xaxis.set_tick_params(labelsize=12)
ax.yaxis.set_tick_params(labelsize=12)

ax.set_xlabel("Methylation percentage", fontsize=12)
ax.set_ylabel("Number of loci", fontsize=12)

#fig.tight_layout()
# fig.subplots_adjust(top=0.72)

fig.savefig( os.path.join(GraphFolder,filename_start +"ExtemeAccns_Col0Fit.png"), bbox_inches = 'tight')

plt.show()
# End Graph



##############################

# Can0UKID116Cvi0_CoopStrength
# Start Graph
fig, ax = plt.subplots(1,1,figsize=(4,4))

for spine in ['left','right','top','bottom']:
    ax.spines[spine].set_color('k')
    ax.spines[spine].set_linewidth(0.8)
ax.set_facecolor('white')

#ax.grid(False)
ax.grid(b=True, which='major', color='lightgrey', linestyle=':',linewidth=1)

x_start = 0.
x_end = 100. 
x_hist = np.arange(x_start+x_end/(2.*n_bin),x_end,x_end/n_bin)

# Cvi0
hist_Cvi0, bins_Cvi0 = np.histogram(LocusProperties_Cvi0_df[['Col0_meth_frac']].values*values_scaler, bins=n_bin, range=(x_start,x_end))
# ax.axvline(np.nanmean(LocusProperties_Cvi0_df[['Col0_meth_frac']].values*values_scaler), linestyle='--', color='k',  linewidth = linewidth_2_Dat)
ax.plot(x_hist, hist_Cvi0,  linewidth = linewidth_2_Dat, color='k', linestyle='-',
                  label='Cvi0: Mean = %.1f %%' % np.nanmean(LocusProperties_Cvi0_df[['Col0_meth_frac']].values*values_scaler))


# UKID116
hist_UKID116, bins_UKID116 = np.histogram(LocusProperties_UKID116_df[['Col0_meth_frac']].values*values_scaler, bins=n_bin, range=(x_start,x_end))
# ax.axvline(np.nanmean(LocusProperties_UKID116_df[['Col0_meth_frac']].values*values_scaler), linestyle='--', color='k',  linewidth = linewidth_2_Dat)
ax.plot(x_hist, hist_UKID116,  linewidth = linewidth_2_Dat, color='k', linestyle='--',
                  label='UKID116: Mean = %.1f %%' % np.nanmean(LocusProperties_UKID116_df[['Col0_meth_frac']].values*values_scaler))


# Can0
hist_Can0, bins_Can0 = np.histogram(LocusProperties_Can0_df[['Col0_meth_frac']].values*values_scaler, bins=n_bin, range=(x_start,x_end))
# ax.axvline(np.nanmean(LocusProperties_Can0_df[['Col0_meth_frac']].values*values_scaler), linestyle='--', color='k',  linewidth = linewidth_2_Dat)
ax.plot(x_hist, hist_Can0,  linewidth = linewidth_2_Dat, color='k', linestyle=':',
                  label='Can0: Mean = %.1f %%' % np.nanmean(LocusProperties_Can0_df[['Col0_meth_frac']].values*values_scaler))

# Single Sim
temp_cols_list = []
for i_temp in range(N_reps):
    temp_cols_list.append('Sim_meth_frac_'+str(i_temp))

for i_ in range(len(AllRepsMethFracs_Can0UKID116Cvi0_CoopStrength_df_list)):
    temp_df = AllRepsMethFracs_Can0UKID116Cvi0_CoopStrength_df_list[i_]
    hist_temp, bins_temp = np.histogram(temp_df[temp_cols_list].values*values_scaler, bins=n_bin, range=(x_start,x_end))
    hist_temp = hist_temp/N_reps # normalise
    # ax.axvline(np.nanmean(temp_df[temp_cols_list].values*values_scaler), linestyle='--', color=colors_list_Can0UKID116Cvi0_CoopStrength[i_],linewidth = linewidth_2_Sim)
    ax.plot(x_hist, hist_temp,  linewidth = linewidth_2_Sim, color=colors_list_Can0UKID116Cvi0_CoopStrength[i_],
                    label='%s: Mean = %.1f %%' % (Labels_Can0UKID116Cvi0_CoopStrength[i_],np.nanmean(temp_df[temp_cols_list].values*values_scaler)) )

ax.legend(framealpha=1)

ax.set_xlim(-10,110)
ax.set_ylim(0,2500)

ax.xaxis.set_tick_params(labelsize=12)
ax.yaxis.set_tick_params(labelsize=12)

ax.set_xlabel("Methylation percentage", fontsize=12)
ax.set_ylabel("Number of loci", fontsize=12)

#fig.tight_layout()
# fig.subplots_adjust(top=0.72)

fig.savefig( os.path.join(GraphFolder,filename_start +"ExtemeAccns_Can0UKID116Cvi0_CoopStrength.png"), bbox_inches = 'tight')

plt.show()
# End Graph



##############################

# Can0UKID116Cvi0_Delta
# Start Graph
fig, ax = plt.subplots(1,1,figsize=(4,4))

for spine in ['left','right','top','bottom']:
    ax.spines[spine].set_color('k')
    ax.spines[spine].set_linewidth(0.8)
ax.set_facecolor('white')

#ax.grid(False)
ax.grid(b=True, which='major', color='lightgrey', linestyle=':',linewidth=1)

x_start = 0.
x_end = 100. 
x_hist = np.arange(x_start+x_end/(2.*n_bin),x_end,x_end/n_bin)

# Cvi0
hist_Cvi0, bins_Cvi0 = np.histogram(LocusProperties_Cvi0_df[['Col0_meth_frac']].values*values_scaler, bins=n_bin, range=(x_start,x_end))
# ax.axvline(np.nanmean(LocusProperties_Cvi0_df[['Col0_meth_frac']].values*values_scaler), linestyle='--', color='k',  linewidth = linewidth_2_Dat)
ax.plot(x_hist, hist_Cvi0,  linewidth = linewidth_2_Dat, color='k', linestyle='-',
                  label='Cvi0: Mean = %.1f %%' % np.nanmean(LocusProperties_Cvi0_df[['Col0_meth_frac']].values*values_scaler))


# UKID116
hist_UKID116, bins_UKID116 = np.histogram(LocusProperties_UKID116_df[['Col0_meth_frac']].values*values_scaler, bins=n_bin, range=(x_start,x_end))
# ax.axvline(np.nanmean(LocusProperties_UKID116_df[['Col0_meth_frac']].values*values_scaler), linestyle='--', color='k',  linewidth = linewidth_2_Dat)
ax.plot(x_hist, hist_UKID116,  linewidth = linewidth_2_Dat, color='k', linestyle='--',
                  label='UKID116: Mean = %.1f %%' % np.nanmean(LocusProperties_UKID116_df[['Col0_meth_frac']].values*values_scaler))


# Can0
hist_Can0, bins_Can0 = np.histogram(LocusProperties_Can0_df[['Col0_meth_frac']].values*values_scaler, bins=n_bin, range=(x_start,x_end))
# ax.axvline(np.nanmean(LocusProperties_Can0_df[['Col0_meth_frac']].values*values_scaler), linestyle='--', color='k',  linewidth = linewidth_2_Dat)
ax.plot(x_hist, hist_Can0,  linewidth = linewidth_2_Dat, color='k', linestyle=':',
                  label='Can0: Mean = %.1f %%' % np.nanmean(LocusProperties_Can0_df[['Col0_meth_frac']].values*values_scaler))

# Single Sim
temp_cols_list = []
for i_temp in range(N_reps):
    temp_cols_list.append('Sim_meth_frac_'+str(i_temp))

for i_ in range(len(AllRepsMethFracs_Can0UKID116Cvi0_Delta_df_list)):
    temp_df = AllRepsMethFracs_Can0UKID116Cvi0_Delta_df_list[i_]
    hist_temp, bins_temp = np.histogram(temp_df[temp_cols_list].values*values_scaler, bins=n_bin, range=(x_start,x_end))
    hist_temp = hist_temp/N_reps # normalise
    # ax.axvline(np.nanmean(temp_df[temp_cols_list].values*values_scaler), linestyle='--', color=colors_list_Can0UKID116Cvi0_Delta[i_],linewidth = linewidth_2_Sim)
    ax.plot(x_hist, hist_temp,  linewidth = linewidth_2_Sim, color=colors_list_Can0UKID116Cvi0_Delta[i_],
                    label='%s: Mean = %.1f %%' % (Labels_Can0UKID116Cvi0_Delta[i_],np.nanmean(temp_df[temp_cols_list].values*values_scaler)) )

ax.legend(framealpha=1)

ax.set_xlim(-10,110)
ax.set_ylim(0,2500)

ax.xaxis.set_tick_params(labelsize=12)
ax.yaxis.set_tick_params(labelsize=12)

ax.set_xlabel("Methylation percentage", fontsize=12)
ax.set_ylabel("Number of loci", fontsize=12)

#fig.tight_layout()
# fig.subplots_adjust(top=0.72)

fig.savefig( os.path.join(GraphFolder,filename_start +"ExtemeAccns_Can0UKID116Cvi0_Delta.png"), bbox_inches = 'tight')

plt.show()
# End Graph



##############################

# Can0UKID116Cvi0_Both
# Start Graph
fig, ax = plt.subplots(1,1,figsize=(4,4))

for spine in ['left','right','top','bottom']:
    ax.spines[spine].set_color('k')
    ax.spines[spine].set_linewidth(0.8)
ax.set_facecolor('white')

#ax.grid(False)
ax.grid(b=True, which='major', color='lightgrey', linestyle=':',linewidth=1)

x_start = 0.
x_end = 100. 
x_hist = np.arange(x_start+x_end/(2.*n_bin),x_end,x_end/n_bin)

# Col0
hist_Col0, bins_Col0 = np.histogram(LocusProperties_Col0_df[['Col0_meth_frac']].values*values_scaler, bins=n_bin, range=(x_start,x_end))
# ax.axvline(np.nanmean(LocusProperties_Col0_df[['Col0_X_frac']].values*values_scaler), linestyle='--', color='k',  linewidth = linewidth_2_Dat)
ax.plot(x_hist, hist_Col0,  linewidth = linewidth_2_Dat, color='green',
                  label='Col0' )

# Cvi0
hist_Cvi0, bins_Cvi0 = np.histogram(LocusProperties_Cvi0_df[['Col0_meth_frac']].values*values_scaler, bins=n_bin, range=(x_start,x_end))
# ax.axvline(np.nanmean(LocusProperties_Cvi0_df[['Col0_meth_frac']].values*values_scaler), linestyle='--', color='k',  linewidth = linewidth_2_Dat)
ax.plot(x_hist, hist_Cvi0,  linewidth = linewidth_2_Dat, color='k', linestyle='-',
                  label='Cvi0')


# UKID116
hist_UKID116, bins_UKID116 = np.histogram(LocusProperties_UKID116_df[['Col0_meth_frac']].values*values_scaler, bins=n_bin, range=(x_start,x_end))
# ax.axvline(np.nanmean(LocusProperties_UKID116_df[['Col0_meth_frac']].values*values_scaler), linestyle='--', color='k',  linewidth = linewidth_2_Dat)
ax.plot(x_hist, hist_UKID116,  linewidth = linewidth_2_Dat, color='k', linestyle='--',
                  label='UKID116' )


# Can0
hist_Can0, bins_Can0 = np.histogram(LocusProperties_Can0_df[['Col0_meth_frac']].values*values_scaler, bins=n_bin, range=(x_start,x_end))
# ax.axvline(np.nanmean(LocusProperties_Can0_df[['Col0_meth_frac']].values*values_scaler), linestyle='--', color='k',  linewidth = linewidth_2_Dat)
ax.plot(x_hist, hist_Can0,  linewidth = linewidth_2_Dat, color='k', linestyle=':',
                  label='Can0' )

# Single Sim
temp_cols_list = []
for i_temp in range(N_reps):
    temp_cols_list.append('Sim_meth_frac_'+str(i_temp))

for i_ in range(len(AllRepsMethFracs_Can0UKID116Cvi0_Both_df_list)):
    temp_df = AllRepsMethFracs_Can0UKID116Cvi0_Both_df_list[i_]
    hist_temp, bins_temp = np.histogram(temp_df[temp_cols_list].values*values_scaler, bins=n_bin, range=(x_start,x_end))
    hist_temp = hist_temp/N_reps # normalise
    # ax.axvline(np.nanmean(temp_df[temp_cols_list].values*values_scaler), linestyle='--', color=colors_list_Can0UKID116Cvi0_Both[i_],linewidth = linewidth_2_Sim)
    ax.plot(x_hist, hist_temp,  linewidth = linewidth_2_Sim, color=colors_list_Can0UKID116Cvi0_Both[i_],
                    label='%s' % (Labels_Can0UKID116Cvi0_Both[i_]) )

ax.legend(framealpha=1)

ax.set_xlim(-10,110)
ax.set_ylim(0,2500)

ax.xaxis.set_tick_params(labelsize=12)
ax.yaxis.set_tick_params(labelsize=12)

ax.set_xlabel("Methylation percentage", fontsize=12)
ax.set_ylabel("Number of loci", fontsize=12)

#fig.tight_layout()
# fig.subplots_adjust(top=0.72)

fig.savefig( os.path.join(GraphFolder,filename_start +"ExtemeAccns_Can0UKID116Cvi0_Both_v1.png"), bbox_inches = 'tight')
fig.savefig( os.path.join(GraphFolder,filename_start +"ExtemeAccns_Can0UKID116Cvi0_Both_v1.pdf"), bbox_inches = 'tight')

plt.show()
# End Graph



##############################

# Dor10_CoopStrength
# Start Graph
fig, ax = plt.subplots(1,1,figsize=(4,4))

for spine in ['left','right','top','bottom']:
    ax.spines[spine].set_color('k')
    ax.spines[spine].set_linewidth(0.8)
ax.set_facecolor('white')

#ax.grid(False)
ax.grid(b=True, which='major', color='lightgrey', linestyle=':',linewidth=1)

x_start = 0.
x_end = 100. 
x_hist = np.arange(x_start+x_end/(2.*n_bin),x_end,x_end/n_bin)

# Dor10
hist_Dor10, bins_Dor10 = np.histogram(LocusProperties_Dor10_df[['Col0_meth_frac']].values*values_scaler, bins=n_bin, range=(x_start,x_end))
# ax.axvline(np.nanmean(LocusProperties_Dor10_df[['Col0_meth_frac']].values*values_scaler), linestyle='--', color='k',  linewidth = linewidth_2_Dat)
ax.plot(x_hist, hist_Dor10,  linewidth = linewidth_2_Dat, color='k',
                  label='Dor10: Mean = %.1f %%' % np.nanmean(LocusProperties_Dor10_df[['Col0_meth_frac']].values*values_scaler))

# Single Sim
temp_cols_list = []
for i_temp in range(N_reps):
    temp_cols_list.append('Sim_meth_frac_'+str(i_temp))

for i_ in range(len(AllRepsMethFracs_Dor10_CoopStrength_df_list)):
    temp_df = AllRepsMethFracs_Dor10_CoopStrength_df_list[i_]
    hist_temp, bins_temp = np.histogram(temp_df[temp_cols_list].values*values_scaler, bins=n_bin, range=(x_start,x_end))
    hist_temp = hist_temp/N_reps # normalise
    # ax.axvline(np.nanmean(temp_df[temp_cols_list].values*values_scaler), linestyle='--', color=colors_list_Dor10_CoopStrength[i_],linewidth = linewidth_2_Sim)
    ax.plot(x_hist, hist_temp,  linewidth = linewidth_2_Sim, color=colors_list_Dor10_CoopStrength[i_],
                    label='%s: Mean = %.1f %%' % (Labels_Dor10_CoopStrength[i_],np.nanmean(temp_df[temp_cols_list].values*values_scaler)) )

ax.legend(framealpha=1)

ax.set_xlim(-10,110)
ax.set_ylim(0,1000)

ax.xaxis.set_tick_params(labelsize=12)
ax.yaxis.set_tick_params(labelsize=12)

ax.set_xlabel("Methylation percentage", fontsize=12)
ax.set_ylabel("Number of loci", fontsize=12)

#fig.tight_layout()
# fig.subplots_adjust(top=0.72)

fig.savefig( os.path.join(GraphFolder,filename_start +"ExtemeAccns_Dor10_CoopStrength.png"), bbox_inches = 'tight')

plt.show()
# End Graph



##############################

# Dor10_Delta
# Start Graph
fig, ax = plt.subplots(1,1,figsize=(4,4))

for spine in ['left','right','top','bottom']:
    ax.spines[spine].set_color('k')
    ax.spines[spine].set_linewidth(0.8)
ax.set_facecolor('white')

#ax.grid(False)
ax.grid(b=True, which='major', color='lightgrey', linestyle=':',linewidth=1)

x_start = 0.
x_end = 100. 
x_hist = np.arange(x_start+x_end/(2.*n_bin),x_end,x_end/n_bin)

# Dor10
hist_Dor10, bins_Dor10 = np.histogram(LocusProperties_Dor10_df[['Col0_meth_frac']].values*values_scaler, bins=n_bin, range=(x_start,x_end))
# ax.axvline(np.nanmean(LocusProperties_Dor10_df[['Col0_meth_frac']].values*values_scaler), linestyle='--', color='k',  linewidth = linewidth_2_Dat)
ax.plot(x_hist, hist_Dor10,  linewidth = linewidth_2_Dat, color='k',
                  label='Dor10: Mean = %.1f %%' % np.nanmean(LocusProperties_Dor10_df[['Col0_meth_frac']].values*values_scaler))

# Single Sim
temp_cols_list = []
for i_temp in range(N_reps):
    temp_cols_list.append('Sim_meth_frac_'+str(i_temp))

for i_ in range(len(AllRepsMethFracs_Dor10_Delta_df_list)):
    temp_df = AllRepsMethFracs_Dor10_Delta_df_list[i_]
    hist_temp, bins_temp = np.histogram(temp_df[temp_cols_list].values*values_scaler, bins=n_bin, range=(x_start,x_end))
    hist_temp = hist_temp/N_reps # normalise
    # ax.axvline(np.nanmean(temp_df[temp_cols_list].values*values_scaler), linestyle='--', color=colors_list_Dor10_Delta[i_],linewidth = linewidth_2_Sim)
    ax.plot(x_hist, hist_temp,  linewidth = linewidth_2_Sim, color=colors_list_Dor10_Delta[i_],
                    label='%s: Mean = %.1f %%' % (Labels_Dor10_Delta[i_],np.nanmean(temp_df[temp_cols_list].values*values_scaler)) )

ax.legend(framealpha=1)

ax.set_xlim(-10,110)
ax.set_ylim(0,1000)

ax.xaxis.set_tick_params(labelsize=12)
ax.yaxis.set_tick_params(labelsize=12)

ax.set_xlabel("Methylation percentage", fontsize=12)
ax.set_ylabel("Number of loci", fontsize=12)

#fig.tight_layout()
# fig.subplots_adjust(top=0.72)

fig.savefig( os.path.join(GraphFolder,filename_start +"ExtemeAccns_Dor10_Delta.png"), bbox_inches = 'tight')

plt.show()
# End Graph


##############################

# Dor10_Both
# Start Graph
fig, ax = plt.subplots(1,1,figsize=(4,4))

for spine in ['left','right','top','bottom']:
    ax.spines[spine].set_color('k')
    ax.spines[spine].set_linewidth(0.8)
ax.set_facecolor('white')

#ax.grid(False)
ax.grid(b=True, which='major', color='lightgrey', linestyle=':',linewidth=1)

x_start = 0.
x_end = 100. 
x_hist = np.arange(x_start+x_end/(2.*n_bin),x_end,x_end/n_bin)

# Col0
hist_Col0, bins_Col0 = np.histogram(LocusProperties_Col0_df[['Col0_meth_frac']].values*values_scaler, bins=n_bin, range=(x_start,x_end))
# ax.axvline(np.nanmean(LocusProperties_Col0_df[['Col0_X_frac']].values*values_scaler), linestyle='--', color='k',  linewidth = linewidth_2_Dat)
ax.plot(x_hist, hist_Col0,  linewidth = linewidth_2_Dat, color='green',
                  label='Col0' )

# Dor10
hist_Dor10, bins_Dor10 = np.histogram(LocusProperties_Dor10_df[['Col0_meth_frac']].values*values_scaler, bins=n_bin, range=(x_start,x_end))
# ax.axvline(np.nanmean(LocusProperties_Dor10_df[['Col0_meth_frac']].values*values_scaler), linestyle='--', color='k',  linewidth = linewidth_2_Dat)
ax.plot(x_hist, hist_Dor10,  linewidth = linewidth_2_Dat, color='k',
                  label='Dor10' )

# Single Sim
temp_cols_list = []
for i_temp in range(N_reps):
    temp_cols_list.append('Sim_meth_frac_'+str(i_temp))

for i_ in range(len(AllRepsMethFracs_Dor10_Both_df_list)):
    temp_df = AllRepsMethFracs_Dor10_Both_df_list[i_]
    hist_temp, bins_temp = np.histogram(temp_df[temp_cols_list].values*values_scaler, bins=n_bin, range=(x_start,x_end))
    hist_temp = hist_temp/N_reps # normalise
    # ax.axvline(np.nanmean(temp_df[temp_cols_list].values*values_scaler), linestyle='--', color=colors_list_Dor10_Both[i_],linewidth = linewidth_2_Sim)
    ax.plot(x_hist, hist_temp,  linewidth = linewidth_2_Sim, color=colors_list_Dor10_Both[i_],
                    label='%s' % (Labels_Dor10_Both[i_]) )

ax.legend(framealpha=1)

ax.set_xlim(-10,110)
ax.set_ylim(0,1000)

ax.xaxis.set_tick_params(labelsize=12)
ax.yaxis.set_tick_params(labelsize=12)

ax.set_xlabel("Methylation percentage", fontsize=12)
ax.set_ylabel("Number of loci", fontsize=12)

#fig.tight_layout()
# fig.subplots_adjust(top=0.72)

fig.savefig( os.path.join(GraphFolder,filename_start +"ExtemeAccns_Dor10_Both_v1.png"), bbox_inches = 'tight')
fig.savefig( os.path.join(GraphFolder,filename_start +"ExtemeAccns_Dor10_Both_v1.pdf"), bbox_inches = 'tight')

plt.show()
# End Graph


##############################

# NS_CoopStrength
# Start Graph
fig, ax = plt.subplots(1,1,figsize=(4,4))

for spine in ['left','right','top','bottom']:
    ax.spines[spine].set_color('k')
    ax.spines[spine].set_linewidth(0.8)
ax.set_facecolor('white')

#ax.grid(False)
ax.grid(b=True, which='major', color='lightgrey', linestyle=':',linewidth=1)

x_start = 0.
x_end = 100. 
x_hist = np.arange(x_start+x_end/(2.*n_bin),x_end,x_end/n_bin)

# NSFilt
temp_cols_list = []
for i_temp in range(N_ExptState_NSFilt):
    temp_cols_list.append('Dset1_meth_frac_'+str(i_temp))
hist_NSFilt, bins_NSFilt = np.histogram(AllRepsMethFracs_NSFilt_df[temp_cols_list].values*values_scaler, bins=n_bin, range=(x_start,x_end))
hist_NSFilt = hist_NSFilt/N_ExptState_NSFilt # normalise
# ax.axvline(np.nanmean(AllRepsMethFracs_NSFilt_df[temp_cols_list].values*values_scaler), linestyle='--', color='k',  linewidth = linewidth_2_Dat)
ax.plot(x_hist, hist_NSFilt,  linewidth = linewidth_2_Dat, color='k',
                  label='NS: Mean = %.1f %%' % np.nanmean(AllRepsMethFracs_NSFilt_df[temp_cols_list].values*values_scaler))

# Single Sim
temp_cols_list = []
for i_temp in range(N_reps):
    temp_cols_list.append('Sim_meth_frac_'+str(i_temp))

for i_ in range(len(AllRepsMethFracs_NS_CoopStrength_df_list)):
    temp_df = AllRepsMethFracs_NS_CoopStrength_df_list[i_]
    hist_temp, bins_temp = np.histogram(temp_df[temp_cols_list].values*values_scaler, bins=n_bin, range=(x_start,x_end))
    hist_temp = hist_temp/N_reps # normalise
    # ax.axvline(np.nanmean(temp_df[temp_cols_list].values*values_scaler), linestyle='--', color=colors_list_NS_CoopStrength[i_],linewidth = linewidth_2_Sim)
    ax.plot(x_hist, hist_temp,  linewidth = linewidth_2_Sim, color=colors_list_NS_CoopStrength[i_],
                    label='%s: Mean = %.1f %%' % (Labels_NS_CoopStrength[i_],np.nanmean(temp_df[temp_cols_list].values*values_scaler)) )

ax.legend(framealpha=1)

ax.set_xlim(-10,110)
ax.set_ylim(0,1000)

ax.xaxis.set_tick_params(labelsize=12)
ax.yaxis.set_tick_params(labelsize=12)

ax.set_xlabel("Methylation percentage", fontsize=12)
ax.set_ylabel("Number of loci", fontsize=12)

#fig.tight_layout()
# fig.subplots_adjust(top=0.72)

fig.savefig( os.path.join(GraphFolder,filename_start +"ExtemeAccns_NS_CoopStrength.png"), bbox_inches = 'tight')

plt.show()
# End Graph




##############################

# NS_Delta
# Start Graph
fig, ax = plt.subplots(1,1,figsize=(4,4))

for spine in ['left','right','top','bottom']:
    ax.spines[spine].set_color('k')
    ax.spines[spine].set_linewidth(0.8)
ax.set_facecolor('white')

#ax.grid(False)
ax.grid(b=True, which='major', color='lightgrey', linestyle=':',linewidth=1)

x_start = 0.
x_end = 100. 
x_hist = np.arange(x_start+x_end/(2.*n_bin),x_end,x_end/n_bin)

# NSFilt
temp_cols_list = []
for i_temp in range(N_ExptState_NSFilt):
    temp_cols_list.append('Dset1_meth_frac_'+str(i_temp))
hist_NSFilt, bins_NSFilt = np.histogram(AllRepsMethFracs_NSFilt_df[temp_cols_list].values*values_scaler, bins=n_bin, range=(x_start,x_end))
hist_NSFilt = hist_NSFilt/N_ExptState_NSFilt # normalise
# ax.axvline(np.nanmean(AllRepsMethFracs_NSFilt_df[temp_cols_list].values*values_scaler), linestyle='--', color='k',  linewidth = linewidth_2_Dat)
ax.plot(x_hist, hist_NSFilt,  linewidth = linewidth_2_Dat, color='k',
                  label='NS: Mean = %.1f %%' % np.nanmean(AllRepsMethFracs_NSFilt_df[temp_cols_list].values*values_scaler))

# Single Sim
temp_cols_list = []
for i_temp in range(N_reps):
    temp_cols_list.append('Sim_meth_frac_'+str(i_temp))

for i_ in range(len(AllRepsMethFracs_NS_Delta_df_list)):
    temp_df = AllRepsMethFracs_NS_Delta_df_list[i_]
    hist_temp, bins_temp = np.histogram(temp_df[temp_cols_list].values*values_scaler, bins=n_bin, range=(x_start,x_end))
    hist_temp = hist_temp/N_reps # normalise
    # ax.axvline(np.nanmean(temp_df[temp_cols_list].values*values_scaler), linestyle='--', color=colors_list_NS_Delta[i_],linewidth = linewidth_2_Sim)
    ax.plot(x_hist, hist_temp,  linewidth = linewidth_2_Sim, color=colors_list_NS_Delta[i_],
                    label='%s: Mean = %.1f %%' % (Labels_NS_Delta[i_],np.nanmean(temp_df[temp_cols_list].values*values_scaler)) )

ax.legend(framealpha=1)

ax.set_xlim(-10,110)
ax.set_ylim(0,1000)

ax.xaxis.set_tick_params(labelsize=12)
ax.yaxis.set_tick_params(labelsize=12)

ax.set_xlabel("Methylation percentage", fontsize=12)
ax.set_ylabel("Number of loci", fontsize=12)

#fig.tight_layout()
# fig.subplots_adjust(top=0.72)

fig.savefig( os.path.join(GraphFolder,filename_start +"ExtemeAccns_NS_Delta.png"), bbox_inches = 'tight')

plt.show()
# End Graph


##############################

# NS_Both
# Start Graph
fig, ax = plt.subplots(1,1,figsize=(4,4))

for spine in ['left','right','top','bottom']:
    ax.spines[spine].set_color('k')
    ax.spines[spine].set_linewidth(0.8)
ax.set_facecolor('white')

#ax.grid(False)
ax.grid(b=True, which='major', color='lightgrey', linestyle=':',linewidth=1)

x_start = 0.
x_end = 100. 
x_hist = np.arange(x_start+x_end/(2.*n_bin),x_end,x_end/n_bin)

# Col0
hist_Col0, bins_Col0 = np.histogram(LocusProperties_Col0_df[['Col0_meth_frac']].values*values_scaler, bins=n_bin, range=(x_start,x_end))
# ax.axvline(np.nanmean(LocusProperties_Col0_df[['Col0_X_frac']].values*values_scaler), linestyle='--', color='k',  linewidth = linewidth_2_Dat)
ax.plot(x_hist, hist_Col0,  linewidth = linewidth_2_Dat, color='green',
                  label='Col0' )

# NSFilt
temp_cols_list = []
for i_temp in range(N_ExptState_NSFilt):
    temp_cols_list.append('Dset1_meth_frac_'+str(i_temp))
hist_NSFilt, bins_NSFilt = np.histogram(AllRepsMethFracs_NSFilt_df[temp_cols_list].values*values_scaler, bins=n_bin, range=(x_start,x_end))
hist_NSFilt = hist_NSFilt/N_ExptState_NSFilt # normalise
# ax.axvline(np.nanmean(AllRepsMethFracs_NSFilt_df[temp_cols_list].values*values_scaler), linestyle='--', color='k',  linewidth = linewidth_2_Dat)
ax.plot(x_hist, hist_NSFilt,  linewidth = linewidth_2_Dat, color='k',
                  label='NS' )

# Single Sim
temp_cols_list = []
for i_temp in range(N_reps):
    temp_cols_list.append('Sim_meth_frac_'+str(i_temp))

for i_ in range(len(AllRepsMethFracs_NS_Both_df_list)):
    temp_df = AllRepsMethFracs_NS_Both_df_list[i_]
    hist_temp, bins_temp = np.histogram(temp_df[temp_cols_list].values*values_scaler, bins=n_bin, range=(x_start,x_end))
    hist_temp = hist_temp/N_reps # normalise
    # ax.axvline(np.nanmean(temp_df[temp_cols_list].values*values_scaler), linestyle='--', color=colors_list_NS_Both[i_],linewidth = linewidth_2_Sim)
    ax.plot(x_hist, hist_temp,  linewidth = linewidth_2_Sim, color=colors_list_NS_Both[i_],
                    label='%s' % (Labels_NS_Both[i_]) )

ax.legend(framealpha=1)

ax.set_xlim(-10,110)
ax.set_ylim(0,1000)

ax.xaxis.set_tick_params(labelsize=12)
ax.yaxis.set_tick_params(labelsize=12)

ax.set_xlabel("Methylation percentage", fontsize=12)
ax.set_ylabel("Number of loci", fontsize=12)

#fig.tight_layout()
# fig.subplots_adjust(top=0.72)

fig.savefig( os.path.join(GraphFolder,filename_start +"ExtemeAccns_NS_Both_v1.png"), bbox_inches = 'tight')
fig.savefig( os.path.join(GraphFolder,filename_start +"ExtemeAccns_NS_Both_v1.pdf"), bbox_inches = 'tight')

plt.show()
# End Graph

#####################################

# RL_CoopStrength
# Start Graph
fig, ax = plt.subplots(1,1,figsize=(4,4))

for spine in ['left','right','top','bottom']:
    ax.spines[spine].set_color('k')
    ax.spines[spine].set_linewidth(0.8)
ax.set_facecolor('white')

#ax.grid(False)
ax.grid(b=True, which='major', color='lightgrey', linestyle=':',linewidth=1)

x_start = 0.
x_end = 100. 
x_hist = np.arange(x_start+x_end/(2.*n_bin),x_end,x_end/n_bin)

# RLFilt
temp_cols_list = []
for i_temp in range(N_ExptState_RLFilt):
    temp_cols_list.append('Dset1_meth_frac_'+str(i_temp))
hist_RLFilt, bins_RLFilt = np.histogram(AllRepsMethFracs_RLFilt_df[temp_cols_list].values*values_scaler, bins=n_bin, range=(x_start,x_end))
hist_RLFilt = hist_RLFilt/N_ExptState_RLFilt # normalise
# ax.axvline(np.nanmean(AllRepsMethFracs_RLFilt_df[temp_cols_list].values*values_scaler), linestyle='--', color='k',  linewidth = linewidth_2_Dat)
ax.plot(x_hist, hist_RLFilt,  linewidth = linewidth_2_Dat, color='k',
                  label='RL: Mean = %.1f %%' % np.nanmean(AllRepsMethFracs_RLFilt_df[temp_cols_list].values*values_scaler))

# Single Sim
temp_cols_list = []
for i_temp in range(N_reps):
    temp_cols_list.append('Sim_meth_frac_'+str(i_temp))

for i_ in range(len(AllRepsMethFracs_RL_CoopStrength_df_list)):
    temp_df = AllRepsMethFracs_RL_CoopStrength_df_list[i_]
    hist_temp, bins_temp = np.histogram(temp_df[temp_cols_list].values*values_scaler, bins=n_bin, range=(x_start,x_end))
    hist_temp = hist_temp/N_reps # normalise
    # ax.axvline(np.nanmean(temp_df[temp_cols_list].values*values_scaler), linestyle='--', color=colors_list_RL_CoopStrength[i_],linewidth = linewidth_2_Sim)
    ax.plot(x_hist, hist_temp,  linewidth = linewidth_2_Sim, color=colors_list_RL_CoopStrength[i_],
                    label='%s: Mean = %.1f %%' % (Labels_RL_CoopStrength[i_],np.nanmean(temp_df[temp_cols_list].values*values_scaler)) )

ax.legend(framealpha=1)

ax.set_xlim(-10,110)
ax.set_ylim(0,2000)

ax.xaxis.set_tick_params(labelsize=12)
ax.yaxis.set_tick_params(labelsize=12)

ax.set_xlabel("Methylation percentage", fontsize=12)
ax.set_ylabel("Number of loci", fontsize=12)

#fig.tight_layout()
# fig.subplots_adjust(top=0.72)

fig.savefig( os.path.join(GraphFolder,filename_start +"ExtemeAccns_RL_CoopStrength.png"), bbox_inches = 'tight')

plt.show()
# End Graph



##############################

# RL_Delta
# Start Graph
fig, ax = plt.subplots(1,1,figsize=(4,4))

for spine in ['left','right','top','bottom']:
    ax.spines[spine].set_color('k')
    ax.spines[spine].set_linewidth(0.8)
ax.set_facecolor('white')

#ax.grid(False)
ax.grid(b=True, which='major', color='lightgrey', linestyle=':',linewidth=1)

x_start = 0.
x_end = 100. 
x_hist = np.arange(x_start+x_end/(2.*n_bin),x_end,x_end/n_bin)

# RLFilt
temp_cols_list = []
for i_temp in range(N_ExptState_RLFilt):
    temp_cols_list.append('Dset1_meth_frac_'+str(i_temp))
hist_RLFilt, bins_RLFilt = np.histogram(AllRepsMethFracs_RLFilt_df[temp_cols_list].values*values_scaler, bins=n_bin, range=(x_start,x_end))
hist_RLFilt = hist_RLFilt/N_ExptState_RLFilt # normalise
# ax.axvline(np.nanmean(AllRepsMethFracs_RLFilt_df[temp_cols_list].values*values_scaler), linestyle='--', color='k',  linewidth = linewidth_2_Dat)
ax.plot(x_hist, hist_RLFilt,  linewidth = linewidth_2_Dat, color='k',
                  label='RL: Mean = %.1f %%' % np.nanmean(AllRepsMethFracs_RLFilt_df[temp_cols_list].values*values_scaler))

# Single Sim
temp_cols_list = []
for i_temp in range(N_reps):
    temp_cols_list.append('Sim_meth_frac_'+str(i_temp))

for i_ in range(len(AllRepsMethFracs_RL_Delta_df_list)):
    temp_df = AllRepsMethFracs_RL_Delta_df_list[i_]
    hist_temp, bins_temp = np.histogram(temp_df[temp_cols_list].values*values_scaler, bins=n_bin, range=(x_start,x_end))
    hist_temp = hist_temp/N_reps # normalise
    # ax.axvline(np.nanmean(temp_df[temp_cols_list].values*values_scaler), linestyle='--', color=colors_list_RL_Delta[i_],linewidth = linewidth_2_Sim)
    ax.plot(x_hist, hist_temp,  linewidth = linewidth_2_Sim, color=colors_list_RL_Delta[i_],
                    label='%s: Mean = %.1f %%' % (Labels_RL_Delta[i_],np.nanmean(temp_df[temp_cols_list].values*values_scaler)) )

ax.legend(framealpha=1)

ax.set_xlim(-10,110)
ax.set_ylim(0,2000)

ax.xaxis.set_tick_params(labelsize=12)
ax.yaxis.set_tick_params(labelsize=12)

ax.set_xlabel("Methylation percentage", fontsize=12)
ax.set_ylabel("Number of loci", fontsize=12)

#fig.tight_layout()
# fig.subplots_adjust(top=0.72)

fig.savefig( os.path.join(GraphFolder,filename_start +"ExtemeAccns_RL_Delta.png"), bbox_inches = 'tight')

plt.show()
# End Graph


##############################

# RL_Both
# Start Graph
fig, ax = plt.subplots(1,1,figsize=(4,4))

for spine in ['left','right','top','bottom']:
    ax.spines[spine].set_color('k')
    ax.spines[spine].set_linewidth(0.8)
ax.set_facecolor('white')

#ax.grid(False)
ax.grid(b=True, which='major', color='lightgrey', linestyle=':',linewidth=1)

x_start = 0.
x_end = 100. 
x_hist = np.arange(x_start+x_end/(2.*n_bin),x_end,x_end/n_bin)


# Col0
hist_Col0, bins_Col0 = np.histogram(LocusProperties_Col0_df[['Col0_meth_frac']].values*values_scaler, bins=n_bin, range=(x_start,x_end))
# ax.axvline(np.nanmean(LocusProperties_Col0_df[['Col0_X_frac']].values*values_scaler), linestyle='--', color='k',  linewidth = linewidth_2_Dat)
ax.plot(x_hist, hist_Col0,  linewidth = linewidth_2_Dat, color='green',
                  label='Col0' )

# RLFilt
temp_cols_list = []
for i_temp in range(N_ExptState_RLFilt):
    temp_cols_list.append('Dset1_meth_frac_'+str(i_temp))
hist_RLFilt, bins_RLFilt = np.histogram(AllRepsMethFracs_RLFilt_df[temp_cols_list].values*values_scaler, bins=n_bin, range=(x_start,x_end))
hist_RLFilt = hist_RLFilt/N_ExptState_RLFilt # normalise
# ax.axvline(np.nanmean(AllRepsMethFracs_RLFilt_df[temp_cols_list].values*values_scaler), linestyle='--', color='k',  linewidth = linewidth_2_Dat)
ax.plot(x_hist, hist_RLFilt,  linewidth = linewidth_2_Dat, color='k',
                  label='RL' )

# Single Sim
temp_cols_list = []
for i_temp in range(N_reps):
    temp_cols_list.append('Sim_meth_frac_'+str(i_temp))

for i_ in range(len(AllRepsMethFracs_RL_Both_df_list)):
    temp_df = AllRepsMethFracs_RL_Both_df_list[i_]
    hist_temp, bins_temp = np.histogram(temp_df[temp_cols_list].values*values_scaler, bins=n_bin, range=(x_start,x_end))
    hist_temp = hist_temp/N_reps # normalise
    # ax.axvline(np.nanmean(temp_df[temp_cols_list].values*values_scaler), linestyle='--', color=colors_list_RL_Both[i_],linewidth = linewidth_2_Sim)
    ax.plot(x_hist, hist_temp,  linewidth = linewidth_2_Sim, color=colors_list_RL_Both[i_],
                    label='%s' % (Labels_RL_Both[i_]) )

ax.legend(framealpha=1)

ax.set_xlim(-10,110)
ax.set_ylim(0,2000)

ax.xaxis.set_tick_params(labelsize=12)
ax.yaxis.set_tick_params(labelsize=12)

ax.set_xlabel("Methylation percentage", fontsize=12)
ax.set_ylabel("Number of loci", fontsize=12)

#fig.tight_layout()
# fig.subplots_adjust(top=0.72)

fig.savefig( os.path.join(GraphFolder,filename_start +"ExtemeAccns_RL_Both_v1.png"), bbox_inches = 'tight')
fig.savefig( os.path.join(GraphFolder,filename_start +"ExtemeAccns_RL_Both_v1.pdf"), bbox_inches = 'tight')

plt.show()
# End Graph

###############
# insets

fig, ax = plt.subplots(2,2,figsize=(4,4))

for axi_ in range(0,2):
    for axj_ in range(0,2):
        for spine in ['left','right','top','bottom']:
            ax[axi_,axj_].spines[spine].set_color('k')
            ax[axi_,axj_].spines[spine].set_linewidth(0.8)
        ax[axi_,axj_].set_facecolor('white')

        #ax.grid(False)
        ax[axi_,axj_].grid(b=True, which='major', color='lightgrey', linestyle=':',linewidth=1)

x_start = 0.
x_end = 100. 
x_hist = np.arange(x_start+x_end/(2.*n_bin),x_end,x_end/n_bin)

# Dor10

# Dor10
hist_Dor10, bins_Dor10 = np.histogram(LocusProperties_Dor10_df[['Col0_meth_frac']].values*values_scaler, bins=n_bin, range=(x_start,x_end))
ax[0,0].plot(x_hist, hist_Dor10,  linewidth = linewidth_2_Dat, color='k',
                  label='Dor10' )

# Single Sim
temp_cols_list = []
for i_temp in range(N_reps):
    temp_cols_list.append('Sim_meth_frac_'+str(i_temp))

for i_ in range(len(AllRepsMethFracs_Dor10_Both_df_list)):
    temp_df = AllRepsMethFracs_Dor10_Both_df_list[i_]
    hist_temp, bins_temp = np.histogram(temp_df[temp_cols_list].values*values_scaler, bins=n_bin, range=(x_start,x_end))
    hist_temp = hist_temp/N_reps # normalise
    ax[0,0].plot(x_hist, hist_temp,  linewidth = linewidth_2_Sim, color=colors_list_Dor10_Both[i_],
                    label='%s' % (Labels_Dor10_Both[i_]) )

ax[0,0].set_xlim(-2,10)
ax[0,0].set_ylim(0,400)

ax[0,0].xaxis.set_tick_params(labelsize=10)
ax[0,0].yaxis.set_tick_params(labelsize=10)

# NS

# NSFilt
temp_cols_list = []
for i_temp in range(N_ExptState_NSFilt):
    temp_cols_list.append('Dset1_meth_frac_'+str(i_temp))
hist_NSFilt, bins_NSFilt = np.histogram(AllRepsMethFracs_NSFilt_df[temp_cols_list].values*values_scaler, bins=n_bin, range=(x_start,x_end))
hist_NSFilt = hist_NSFilt/N_ExptState_NSFilt # normalise
ax[0,1].plot(x_hist, hist_NSFilt,  linewidth = linewidth_2_Dat, color='k',
                  label='NS' )

# Single Sim
temp_cols_list = []
for i_temp in range(N_reps):
    temp_cols_list.append('Sim_meth_frac_'+str(i_temp))

for i_ in range(len(AllRepsMethFracs_NS_Both_df_list)):
    temp_df = AllRepsMethFracs_NS_Both_df_list[i_]
    hist_temp, bins_temp = np.histogram(temp_df[temp_cols_list].values*values_scaler, bins=n_bin, range=(x_start,x_end))
    hist_temp = hist_temp/N_reps # normalise
    ax[0,1].plot(x_hist, hist_temp,  linewidth = linewidth_2_Sim, color=colors_list_NS_Both[i_],
                    label='%s' % (Labels_NS_Both[i_]) )

ax[0,1].set_xlim(-2,10)
ax[0,1].set_ylim(100,600)

ax[0,1].xaxis.set_tick_params(labelsize=10)
ax[0,1].yaxis.set_tick_params(labelsize=10)


# RL

# RLFilt
temp_cols_list = []
for i_temp in range(N_ExptState_RLFilt):
    temp_cols_list.append('Dset1_meth_frac_'+str(i_temp))
hist_RLFilt, bins_RLFilt = np.histogram(AllRepsMethFracs_RLFilt_df[temp_cols_list].values*values_scaler, bins=n_bin, range=(x_start,x_end))
hist_RLFilt = hist_RLFilt/N_ExptState_RLFilt # normalise
ax[1,0].plot(x_hist, hist_RLFilt,  linewidth = linewidth_2_Dat, color='k',
                  label='RL' )

# Single Sim
temp_cols_list = []
for i_temp in range(N_reps):
    temp_cols_list.append('Sim_meth_frac_'+str(i_temp))

for i_ in range(len(AllRepsMethFracs_RL_Both_df_list)):
    temp_df = AllRepsMethFracs_RL_Both_df_list[i_]
    hist_temp, bins_temp = np.histogram(temp_df[temp_cols_list].values*values_scaler, bins=n_bin, range=(x_start,x_end))
    hist_temp = hist_temp/N_reps # normalise
    ax[1,0].plot(x_hist, hist_temp,  linewidth = linewidth_2_Sim, color=colors_list_RL_Both[i_],
                    label='%s' % (Labels_RL_Both[i_]) )

ax[1,0].set_xlim(-2,10)
ax[1,0].set_ylim(600,1200)

ax[1,0].xaxis.set_tick_params(labelsize=10)
ax[1,0].yaxis.set_tick_params(labelsize=10)


# Cvi etc

# Cvi0
hist_Cvi0, bins_Cvi0 = np.histogram(LocusProperties_Cvi0_df[['Col0_meth_frac']].values*values_scaler, bins=n_bin, range=(x_start,x_end))
ax[1,1].plot(x_hist, hist_Cvi0,  linewidth = linewidth_2_Dat, color='k', linestyle='-',
                  label='Cvi0')
print('Cvi0',hist_Cvi0[0:3])




# # Can0
# hist_Can0, bins_Can0 = np.histogram(LocusProperties_Can0_df[['Col0_meth_frac']].values*values_scaler, bins=n_bin, range=(x_start,x_end))
# ax[1,1].plot(x_hist, hist_Can0,  linewidth = linewidth_2_Dat, color='k', linestyle=':',
#                   label='Can0' )
print('Can0', hist_Can0[0:3])

# Single Sim
temp_cols_list = []
for i_temp in range(N_reps):
    temp_cols_list.append('Sim_meth_frac_'+str(i_temp))

for i_ in range(len(AllRepsMethFracs_Can0UKID116Cvi0_Both_df_list)):
    temp_df = AllRepsMethFracs_Can0UKID116Cvi0_Both_df_list[i_]
    hist_temp, bins_temp = np.histogram(temp_df[temp_cols_list].values*values_scaler, bins=n_bin, range=(x_start,x_end))
    hist_temp = hist_temp/N_reps # normalise
    ax[1,1].plot(x_hist, hist_temp,  linewidth = linewidth_2_Sim, color=colors_list_Can0UKID116Cvi0_Both[i_],
                    label='%s' % (Labels_Can0UKID116Cvi0_Both[i_]) )

# UKID116
# hist_UKID116, bins_UKID116 = np.histogram(LocusProperties_UKID116_df[['Col0_meth_frac']].values*values_scaler, bins=n_bin, range=(x_start,x_end))
# ax[1,1].plot(x_hist, hist_UKID116,  linewidth = linewidth_2_Dat, color='k', linestyle='--',
#                   label='UKID116' )
print('UKID116', hist_UKID116[0:3])


ax[1,1].set_xlim(-2,10)
# ax[1,1].set_ylim(0,2500)

ax[1,1].xaxis.set_tick_params(labelsize=10)
ax[1,1].yaxis.set_tick_params(labelsize=10)



fig.text(0.5, -0.00, 'Methylation percentage', ha='center', va='center', fontsize=12)
fig.text(-0.00, 0.5, 'Number of loci', ha='center', va='center', rotation='vertical', fontsize=12)



fig.tight_layout()
plt.subplots_adjust(wspace=0.45)
# fig.subplots_adjust(top=0.72)

# fig.savefig( os.path.join(GraphFolder,filename_start +"ExtemeAccns_Dor10_Both_v1.png"), bbox_inches = 'tight')
fig.savefig( os.path.join(GraphFolder,filename_start +"ExtemeAccns_Insets.png"), bbox_inches = 'tight')

###############################

# Start Graph
fig, ax = plt.subplots(1,1,figsize=(4,4))

for spine in ['left','right','top','bottom']:
    ax.spines[spine].set_color('k')
    ax.spines[spine].set_linewidth(0.8)
ax.set_facecolor('white')

#ax.grid(False)
ax.grid(b=True, which='major', color='lightgrey', linestyle=':',linewidth=1)

x_start = 0.
x_end = 100. 
x_hist = np.arange(x_start+x_end/(2.*n_bin),x_end,x_end/n_bin)

# temp_cols_list = []
# for i_temp in range(N_reps):
#     temp_cols_list.append('Sim_meth_frac_'+str(i_temp))
# hist_Sim, bins_Sim = np.histogram(AllRepsMethFracs_df[temp_cols_list].values, bins=n_bin, range=(0,1))
# hist_Sim = hist_Sim/N_reps # normalise
# ax[0].axvline(np.nanmean(AllRepsMethFracs_df[temp_cols_list].values), linestyle='--', color='r')
# ax[1].axvline(np.nanmean(AllRepsMethFracs_df[temp_cols_list].values), linestyle='--', color='r')
# ax[0].plot(x_hist, hist_Sim,  linewidth = 2, color='r',
#                   label='N_bins = %d\nSim. $\mu$ = %.0f' % (n_bin,np.nanmean(AllRepsMethFracs_df[temp_cols_list].values)))
# ax[1].plot(x_hist, hist_Sim,  linewidth = 2, color='r',
#                   label='Normalised\nSim. $\mu$ = %.0f' % np.nanmean(AllRepsMethFracs_df[temp_cols_list].values))




# Cvi0
hist_Cvi0, bins_Cvi0 = np.histogram(LocusProperties_Cvi0_df[['Col0_meth_frac']].values*values_scaler, bins=n_bin, range=(x_start,x_end))
# ax.axvline(np.nanmean(LocusProperties_Cvi0_df[['Col0_meth_frac']].values*values_scaler), linestyle='--', color=colors_list_sweep[0],  linewidth = linewidth_1)
ax.plot(x_hist, hist_Cvi0,  linewidth = linewidth_1, color=colors_list_sweep[0],
                  label='Cvi0: Mean = %.0f %%' % np.nanmean(LocusProperties_Cvi0_df[['Col0_meth_frac']].values*values_scaler))


# UKID116
hist_UKID116, bins_UKID116 = np.histogram(LocusProperties_UKID116_df[['Col0_meth_frac']].values*values_scaler, bins=n_bin, range=(x_start,x_end))
# ax.axvline(np.nanmean(LocusProperties_UKID116_df[['Col0_meth_frac']].values*values_scaler), linestyle='--', color=colors_list_sweep[1],  linewidth = linewidth_1)
ax.plot(x_hist, hist_UKID116,  linewidth = linewidth_1, color=colors_list_sweep[1],
                  label='UKID116: Mean = %.0f %%' % np.nanmean(LocusProperties_UKID116_df[['Col0_meth_frac']].values*values_scaler))


# Can0
hist_Can0, bins_Can0 = np.histogram(LocusProperties_Can0_df[['Col0_meth_frac']].values*values_scaler, bins=n_bin, range=(x_start,x_end))
# ax.axvline(np.nanmean(LocusProperties_Can0_df[['Col0_meth_frac']].values*values_scaler), linestyle='--', color=colors_list_sweep[2],  linewidth = linewidth_1)
ax.plot(x_hist, hist_Can0,  linewidth = linewidth_1, color=colors_list_sweep[2],
                  label='Can0: Mean = %.0f %%' % np.nanmean(LocusProperties_Can0_df[['Col0_meth_frac']].values*values_scaler))




ax.legend(framealpha=1)

ax.set_xlim(-10,110)

ax.set_ylim(0,2500)

ax.xaxis.set_tick_params(labelsize=12)
ax.yaxis.set_tick_params(labelsize=12)


ax.set_xlabel("Methylation percentage", fontsize=12)
ax.set_ylabel("Number of loci", fontsize=12)

#fig.tight_layout()
# fig.subplots_adjust(top=0.72)

fig.savefig( os.path.join(GraphFolder,filename_start +"Cvi0Can0UKID116.png"), bbox_inches = 'tight')

plt.show()
# End Graph

############################

