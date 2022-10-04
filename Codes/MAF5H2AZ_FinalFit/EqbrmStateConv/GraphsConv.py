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

N_corrns_bins = params.N_corrns_bins

N_reps = params.N_reps
initial_seed = params.initial_seed
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

GraphFolder = "Graphs"

filename_start = params.filename_start
filename_ending = '.tsv'


# define colours
col_data_Col0 = 'green'
col_Sim_Expt = 'k'
col_Sim_100M = 'r'
col_Sim_100U = 'b'

fig_label_Converged = 'F4F'

# title_sting_single_time = "No. sims. per gene = %d      Cell cycle duration = %.1f       No. of genes simulated = %d        Sim. time = %d gens.\n\
# 'P' choice = %s      $\epsilon$ = variable      $\gamma_0$ = varialbe      $r_{div\\,\gamma} =$ %d      $r_{\gamma}$ = %d      $\lambda_{\gamma} =$ %.3f\n\
# $\delta$ = %.3e      $u_{M PL}^+ =$ varialbe      $r_{div.} =$ %d      $r_{plat.} =$ %d      $\lambda_{M PL} =$ %.3f      $u_{Scale\\,LR1}^+ =$ %.3f\n \
# $u_{M PL\\,LR1}^+ =$ varialbe      $r_{div\\,LR1} =$ %d      $r_{plat.\\,LR1} =$ %d      $\lambda_{M\\,In\\,PL\\,LR1} =$ %.3f      $\lambda_{M\\,Out\\,PL\\,LR1} =$ %.3f      $r_{LR} =$ %d\n \
# $m_u$ = %.2f     $c_u$ = %.2f     $m_{\epsilon}$ = %.2f     $c_{\epsilon}$ = %.2f    $m_{\gamma}$ = %.2f     $c_{\gamma}$ = %.2f     $\\rho_{cap}$ = 1./%.2f\n\
# Initial seed = %d      Sim. initial state: %s      %d $\leq N_{CG\\,min}$ < %d      %.3f $\leq \\rho_{CG\\,min}$ < %.3f\
# \nLocus type: %s      Annotation file: %s\n\
# Filter 1: %s      Filter 2: %s      ExcludeFilter 1: %s\n\
# $N_{bin}$ = %d"\
#              % (N_reps, rep_time, N_gene, N_gen_burn_in,  P_choice,
#                 r_div_gamma, r_gamma, lambda_gamma, delta, 
#                r_div, r_plat,lambda_coop, u_scale_val,
#                 r_div_LR1, r_plat_LR1, lambda_coopIn_LR1, lambda_coopOut_LR1,
#                 r_LR1, u_m_val,u_c_val,e_m_val,e_c_val,g_m_val,g_c_val,spacing_cap,
#                 initial_seed, initial_state_choice,N_CG_min,N_CG_max, N_CG_density_min,N_CG_density_max, locus_type, 
#                 file_in_annotation, file_in_anno_filt_1, file_in_anno_filt_2, file_in_exclude_filt_1,n_bin)

# example filename
# MAF5_H2AZWT_EqbrmOutput_FinalFit_Conv_Pchoice_Excl_InitialState_100U_Ngen_1E2_LocusProperties.tsv



# define list of time-values
#time_scaler = 1.e+4
#time_scaler_label = '($\\times10^{4}$ gens.)'
time_scaler = 1.
time_scaler_label = '(generations)'

TimeValues_list_100U = [1.e+0/time_scaler, 4.e+0/time_scaler, 1.e+1/time_scaler, 
                        2.e+1/time_scaler, 4.e+1/time_scaler, 5.e+1/time_scaler, 
                        6.e+1/time_scaler, 1.e+2/time_scaler, 2.e+2/time_scaler, 
                        3.e+2/time_scaler, 5.e+2/time_scaler, 8.e+2/time_scaler, 
                        1.e+3/time_scaler, 2.e+3/time_scaler, 3.e+3/time_scaler, 
                        4.e+3/time_scaler, 5.e+3/time_scaler, 6.e+3/time_scaler, 
                        7.e+3/time_scaler, 8.e+3/time_scaler, 9.e+3/time_scaler, 
                        1.e+4/time_scaler, 1.5e+4/time_scaler, 2.e+4/time_scaler, 
                        4.e+4/time_scaler, 8.e+4/time_scaler, 1.e+5/time_scaler,
                        1.e+6/time_scaler]
TimeValues_labels_100U = ['1E0_', '4E0_', '1E1_', 
                            '2E1_', '4E1_', '5E1_', 
                            '6E1_', '1E2_', '2E2_', 
                            '3E2_', '5E2_', '8E2_', 
                            '1E3_', '2E3_', '3E3_', 
                            '4E3_', '5E3_', '6E3_', 
                            '7E3_', '8E3_', '9E3_', 
                            '1E4_', '15E3_', '2E4_', 
                            '4E4_', '8E4_', '1E5_', '1E6_']

TimeValues_list_100M = TimeValues_list_100U
TimeValues_labels_100M = TimeValues_labels_100U

# end filename for LocusProperties
filename_end_LocusProperties = 'LocusProperties.tsv'

# define current initial-state
current_initial_state = '100U_'
# define start of filename
filename_begining = 'MAF5_H2AZWT_EqbrmOutput_FinalFit_Conv_Pchoice_Excl_InitialState_'+current_initial_state+'Ngen_'
MeanMethFrac_vals_list_100U = []

# loop over simulation lengths
for i_sim in range(len(TimeValues_list_100U)):
    # LocusProperties used for Col0 data - only need to load in once! 
    file_in_LocusProperties_temp = os.path.join('Output_files', filename_begining+TimeValues_labels_100U[i_sim]+filename_end_LocusProperties)
    print(file_in_LocusProperties_temp)
    LocusProperties_df_temp = pd.read_csv(file_in_LocusProperties_temp, sep="\t{1}",engine='python')
    LocusProperties_df_temp.set_index("gene_ID", inplace = True)
    MeanMethFrac_val_temp = np.nanmean(LocusProperties_df_temp['Sim_mu_meth_frac'])*100
    MeanMethFrac_vals_list_100U.append(MeanMethFrac_val_temp)
    print(len(LocusProperties_df_temp))
    print(TimeValues_list_100U[i_sim],MeanMethFrac_val_temp)

print()

# define current initial-state
current_initial_state = '100M_'
# define start of filename
filename_begining = 'MAF5_H2AZWT_EqbrmOutput_FinalFit_Conv_Pchoice_Excl_InitialState_'+current_initial_state+'Ngen_'
MeanMethFrac_vals_list_100M = []

# loop over simulation lengths
for i_sim in range(len(TimeValues_list_100M)):
    # LocusProperties used for Col0 data - only need to load in once! 
    file_in_LocusProperties_temp = os.path.join('Output_files', filename_begining+TimeValues_labels_100M[i_sim]+filename_end_LocusProperties)
    print(file_in_LocusProperties_temp)
    LocusProperties_df_temp = pd.read_csv(file_in_LocusProperties_temp, sep="\t{1}",engine='python')
    LocusProperties_df_temp.set_index("gene_ID", inplace = True)
    MeanMethFrac_val_temp = np.nanmean(LocusProperties_df_temp['Sim_mu_meth_frac'])*100
    MeanMethFrac_vals_list_100M.append(MeanMethFrac_val_temp)
    print(len(LocusProperties_df_temp))
    print(TimeValues_list_100U[i_sim],MeanMethFrac_val_temp)





# Start Graph
fig, ax = plt.subplots(1,1,figsize=(5,5))

for spine in ['left','right','top','bottom']:
    ax.spines[spine].set_color('k')
    ax.spines[spine].set_linewidth(0.8)
ax.set_facecolor('white')

#ax.grid(False)
ax.grid(b=True, which='major', color='lightgrey', linestyle=':',linewidth=1)

ax.plot(TimeValues_list_100M, MeanMethFrac_vals_list_100M,  linewidth = 1.5, color=col_Sim_100M,
                  label='$\mu_{{%s}}$ = %.1f' % 
                  (TimeValues_labels_100M[-5][:-1],MeanMethFrac_vals_list_100M[-5]))

ax.plot(TimeValues_list_100U, MeanMethFrac_vals_list_100U,  linewidth = 1.5, color=col_Sim_100U,
                  label='100U $\mu^*$ = %.1f\n0.95$\mu^*$ = %.1f\n$\mu_{{%s}}$ = %.1f' % 
                  (MeanMethFrac_vals_list_100U[-1],0.95*MeanMethFrac_vals_list_100U[-1],
                  TimeValues_labels_100U[-11][:-1],MeanMethFrac_vals_list_100U[-11]))
ax.set_xscale('log')
# ax.axhline(y=MeanMethFrac_vals_list_100U[-1], color='k', linestyle=':')
# ax.vlines(x=TimeValues_list_100U[-7], ymin=0.0, ymax=MeanMethFrac_vals_list_100U[-1], color='k', linestyle=':')
# ax.hlines(y=0.95*MeanMethFrac_vals_list_100U[-1], xmin=0.0, xmax=TimeValues_list_100U[-11], color='k', linestyle=':')
# ax.vlines(x=TimeValues_list_100U[-11], ymin=0.0, ymax=0.95*MeanMethFrac_vals_list_100U[-1], color='k', linestyle=':')



#fig.suptitle(title_sting_single_time, fontsize=12)


# ax.legend(fontsize=8,ncol=1,facecolor='white', framealpha=1.0,handlelength=1.0)


#ax.set_xlim(0,1)
ax.set_ylim(0,100)

ax.xaxis.set_tick_params(labelsize=12)
ax.yaxis.set_tick_params(labelsize=12)

ax.set_xlabel("Simulation time "+time_scaler_label, fontsize=12)
ax.set_ylabel("Mean methylation\npercentage", fontsize=12)


#fig.tight_layout()
fig.subplots_adjust(top=0.5)

fig.savefig( os.path.join(GraphFolder,'MAF5_H2AZWT_EqbrmOutput_FinalFit_Conv_Pchoice_Excl_' +"Mfrac_Conv_"+fig_label_Converged+".png"), bbox_inches = 'tight')

plt.show()
# End Graph


##############



