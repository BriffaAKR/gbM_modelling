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

# locus_type = params.locus_type
# TE_filt = params.TE_filt

N_corrns_bins = params.N_corrns_bins

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

# unlikely to edit
N_total_IDs = params.N_total_IDs
n_cc = params.n_cc


n_bin = 50
GraphFolder = "Graphs"

filename_start_Sim = params.filename_start_Sim
# filename_start_Dat = params.filename_start_Dat
filename_start_Graph = params.filename_start_Graph
filename_ending = '.tsv'

# file_in_LocusProperties_Data = os.path.join('Output_files', filename_start_Dat + 'LocusProperties'+filename_ending)
file_in_LocusProperties_Sim = os.path.join('Output_files', filename_start_Sim + 'LocusProperties'+filename_ending)

# file_in_Mcorrn_Data = os.path.join('Output_files', filename_start_Dat + 'Mcorrn_Dat'+filename_ending)
# file_in_NonX_CGcorrn_Data = os.path.join('Output_files', filename_start_Dat + 'NonX_CGcorrn_Dat'+filename_ending)
file_in_All_CGcorrn_Sim = os.path.join('Output_files', filename_start_Sim + 'All_CGcorrn_Sim'+filename_ending)
file_in_Mcorrn_Sim = os.path.join('Output_files', filename_start_Sim + 'Mcorrn_Sim'+filename_ending)

file_in_CodeProgress_Sim = os.path.join('Output_files', filename_start_Sim + 'CodeProgress'+filename_ending)

file_in_RealCGsites_Sim = os.path.join('Output_files','MAF5_H2AZWT_EqbrmOutput_FinalFit_Flucts_Col0likeAccnsNonRedCov50_Sim_Pchoice_Excl_InitialState_100U_Mcorrn_Sim.tsv')

# LocusProperties_Data_df = pd.read_csv(file_in_LocusProperties_Data, sep="\t{1}",engine='python')
# LocusProperties_Data_df.set_index("gene_ID", inplace = True)
# N_gene_Data = len(LocusProperties_Data_df)

LocusProperties_Sim_df = pd.read_csv(file_in_LocusProperties_Sim, sep="\t{1}",engine='python')
LocusProperties_Sim_df.set_index("gene_ID", inplace = True)
N_gene_Sim = len(LocusProperties_Sim_df)
# IDs_Sim_list = LocusProperties_Sim_df.index.values.tolist()

# print('N_gene',N_gene_Data, N_gene_Sim)
print('N_gene',N_gene_Sim)

print()

# add up total number of Non-missing data pairs in data for normalisation
CodeProgress_Sim_df = pd.read_csv(file_in_CodeProgress_Sim, sep="\t{1}",engine='python')
# NonX_pairs_total_Data = CodeProgress_Data_df['NonX_pairs_total'].sum()
All_pairs_total_Sim = CodeProgress_Sim_df['All_pairs_total'].sum()
# print('NonX_pairs_total_Data', NonX_pairs_total_Data)
# print('All_pairs_total_Data', All_pairs_total_Data)
print()
Pairs_total_Sim = (LocusProperties_Sim_df['N_CG']**2).sum()*N_reps
print('Norm_Sim',Pairs_total_Sim, 'All_pairs_total_Sim', All_pairs_total_Sim)
print()

# Mcorrn_Data_df = pd.read_csv(file_in_Mcorrn_Data, sep="\t{1}",engine='python',header=None)
# Mcorrn_Data_array = Mcorrn_Data_df.sum(axis=0).values
# NonX_CGcorrn_Data_df = pd.read_csv(file_in_NonX_CGcorrn_Data, sep="\t{1}",engine='python',header=None)
# NonX_CGcorrn_Data_array = NonX_CGcorrn_Data_df.sum(axis=0).values
All_CGcorrn_Sim_df = pd.read_csv(file_in_All_CGcorrn_Sim, sep="\t{1}",engine='python',header=None)
All_CGcorrn_Sim_array = All_CGcorrn_Sim_df.sum(axis=0).values

Mcorrn_Sim_df = pd.read_csv(file_in_Mcorrn_Sim, sep="\t{1}",engine='python',header=None)
Mcorrn_Sim_array = Mcorrn_Sim_df.sum(axis=0).values

Mcorrn_RealCGsites_Sim_df = pd.read_csv(file_in_RealCGsites_Sim, sep="\t{1}",engine='python',header=None)
Mcorrn_RealCGsites_Sim_array = Mcorrn_RealCGsites_Sim_df.sum(axis=0).values

############

# # Normalise Pairs_Data and Pairs_Sim individually
# PairSep_MM_Data_array = PairSep_MM_Data_array/Total_Pairs_Data
# PairSep_MU_Data_array = PairSep_MU_Data_array/Total_Pairs_Data
# PairSep_UU_Data_array = PairSep_UU_Data_array/Total_Pairs_Data
# PairSep_XX_Data_array = PairSep_XX_Data_array/Total_Pairs_Data

# PairSep_MM_Sim_array = PairSep_MM_Sim_array/Total_Pairs_Sim
# PairSep_MU_Sim_array = PairSep_MU_Sim_array/Total_Pairs_Sim
# PairSep_UU_Sim_array = PairSep_UU_Sim_array/Total_Pairs_Sim

x_vals_corrns = np.arange(N_corrns_bins)



############

# Start Graph
fig, ax = plt.subplots(1,1,figsize=(6,3))

for spine in ['left','right','top','bottom']:
    ax.spines[spine].set_color('k')
    ax.spines[spine].set_linewidth(0.8)
ax.set_facecolor('white')

#ax.grid(False)
ax.grid(b=True, which='major', color='lightgrey', linestyle=':',linewidth=1)

y_scale = 1.0e-3
y_scale_label = ' ($10^{-3})$'

ax.plot(x_vals_corrns[2:], (Mcorrn_RealCGsites_Sim_array[2:]/(All_pairs_total_Sim*740/30))/y_scale,  linewidth = 1, color='b',label='M-M pair; 740 reps. simulated')
ax.plot(x_vals_corrns[2:], (Mcorrn_Sim_array[2:]/All_pairs_total_Sim)/y_scale,  linewidth = 1, color='k',label='M-M pair; 30 reps. simulated')
# ax.plot(x_vals_corrns[2:], (Mcorrn_Data_array[2:]/NonX_pairs_total_Data)/y_scale,  linewidth = 1, color='green',label='M-M pair; 30 Col-like data')
# ax.plot(x_vals_corrns[2:], (NonX_CGcorrn_Data_array[2:]/NonX_pairs_total_Data)/y_scale,  linewidth = 1, color='orange',label='NonX CG-site correlation')
ax.plot(x_vals_corrns[2:], (All_CGcorrn_Sim_array[2:]/All_pairs_total_Sim)/y_scale,  linewidth = 1, color='darkorange',label='CG-CG pair; 30 reps. simulated')



# from mpl_toolkits.axes_grid1.inset_locator import inset_axes
# # gains inset
# axins_0 = inset_axes(ax, width="60%", height="50%", borderpad=1)
# axins_0.tick_params(axis='both', which='major', labelsize=14)
# axins_0.set_xlim(0,20)
# axins_0.set_ylim(0,0.02)
# axins_0.plot(x_vals_corrns, Mcorrn_Sim_array,  linewidth = 1, color='k',label='$MM$-pairs Sim.')
# axins_0.fill_between(x_vals_corrns, Mcorrn_Data_array,Mcorrn_Data_array+PairSep_XX_Data_array,  
#                     alpha = 0.5, color='green',label='$MM$-pairs Data')


ax.legend(fontsize=12,ncol=1,facecolor='white', framealpha=1.0)

ax.set_xlim(0,500)
ax.set_ylim(0,1.)


ax.xaxis.set_tick_params(labelsize=12)
ax.yaxis.set_tick_params(labelsize=12)

ax.set_xlabel("CG-site pair separation (bp)", fontsize=12)
ax.set_ylabel("Density"+y_scale_label, fontsize=12)

fig.tight_layout()
#fig.subplots_adjust(top=0.5)

fig.savefig( os.path.join(GraphFolder,filename_start_Graph +"MM_CGcorrn_rhoCGcontrol.png"))


plt.show()
# End Graph