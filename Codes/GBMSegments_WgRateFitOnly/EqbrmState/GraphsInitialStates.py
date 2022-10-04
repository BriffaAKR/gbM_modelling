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

n_bin = 50
GraphFolder = "Graphs"

filename_start = params.filename_start
filename_ending = '.tsv'

filename_start_Graph = 'GBMSegments_EqbrmOutput_WgRatesOnlyFit_Pchoice_Excl_InitialState_Expt_Ngen_1E5_Graph_'

# LocusProperties used for Col0 data - only need to load in once! 
file_in_LocusProperties = os.path.join('Output_files', 'GBMSegments_EqbrmOutput_WgRatesOnlyFit_Pchoice_Excl_InitialState_Expt_Ngen_1E5_LocusProperties.tsv')

LocusProperties_df = pd.read_csv(file_in_LocusProperties, sep="\t{1}",engine='python')
LocusProperties_df.set_index("gene_ID", inplace = True)
N_gene = len(LocusProperties_df)
IDs_list = LocusProperties_df.index.values.tolist()
print('N_gene',N_gene)


# Manually load in the Sims from all three initial states. 
file_in_AllRepsMethFracs_Expt = os.path.join('Output_files', 
'GBMSegments_EqbrmOutput_WgRatesOnlyFit_Pchoice_Excl_InitialState_Expt_Ngen_1E5_AllRepsMethFracs.tsv')
AllRepsMethFracs_Expt_df = pd.read_csv(file_in_AllRepsMethFracs_Expt, sep="\t{1}",engine='python')
AllRepsMethFracs_Expt_df.set_index("gene_ID", inplace = True)
print(len(AllRepsMethFracs_Expt_df))

# file_in_AllRepsMethFracs_100M = os.path.join('Output_files', 
# 'MAF5_H2AZWT_EqbrmOutput_FinalFit_Pchoice_Excl_InitialState_100M_Ngen_1E5_AllRepsMethFracs.tsv')
# AllRepsMethFracs_100M_df = pd.read_csv(file_in_AllRepsMethFracs_100M, sep="\t{1}",engine='python')
# AllRepsMethFracs_100M_df.set_index("gene_ID", inplace = True)

# file_in_AllRepsMethFracs_100U = os.path.join('Output_files', 
# 'MAF5_H2AZWT_EqbrmOutput_FinalFit_Pchoice_Excl_InitialState_100U_Ngen_1E5_AllRepsMethFracs.tsv')
# AllRepsMethFracs_100U_df = pd.read_csv(file_in_AllRepsMethFracs_100U, sep="\t{1}",engine='python')
# AllRepsMethFracs_100U_df.set_index("gene_ID", inplace = True)

# define colours
col_data_Col0 = 'green'
col_Sim_Expt = 'k'
# col_Sim_100M = 'r'
# col_Sim_100U = 'b'

fig_label_Expt = 'FS6A_right'

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


# convert to percentage
values_scaler = 100

# Start Graph
fig, ax = plt.subplots(1,1,figsize=(3,3))

for spine in ['left','right','top','bottom']:
    ax.spines[spine].set_color('k')
    ax.spines[spine].set_linewidth(0.8)
ax.set_facecolor('white')

#ax.grid(False)
ax.grid(b=True, which='major', color='lightgrey', linestyle=':',linewidth=1)

x_start = 0.
x_end = 100. 
x_hist = np.arange(x_start+x_end/(2.*n_bin),x_end,x_end/n_bin)


hist_Data_label = 'Col-0 Data'


hist_Col0, bins_Col0 = np.histogram(LocusProperties_df[['Col0_meth_frac']].values*values_scaler, bins=n_bin, range=(x_start,x_end))
# ax.axvline(np.nanmean(LocusProperties_df[['Col0_meth_frac']].values*values_scaler), linestyle='-', color=col_data_Col0)
ax.plot(x_hist, hist_Col0,  linewidth = 2, color=col_data_Col0,
                  label='%s\nMean = %.1f %%' % (hist_Data_label,np.nanmean(LocusProperties_df[['Col0_meth_frac']].values*values_scaler)) )

hist_Sim_label = 'Simulated'


temp_cols_list = []
for i_temp in range(N_reps):
    temp_cols_list.append('Sim_meth_frac_'+str(i_temp))
hist_Sims_Expt, bins_Sim = np.histogram(AllRepsMethFracs_Expt_df[temp_cols_list].values*values_scaler, bins=n_bin, range=(x_start,x_end))
hist_Sims_Expt = hist_Sims_Expt/N_reps # normalise
# ax.axvline(np.nanmean(AllRepsMethFracs_Expt_df[temp_cols_list].values*values_scaler), linestyle='-', color=col_Sim_Expt)
ax.plot(x_hist, hist_Sims_Expt,  linewidth = 2, color=col_Sim_Expt,linestyle='-',
                  label='%s\nMean = %.1f %%' % (hist_Sim_label,np.nanmean(AllRepsMethFracs_Expt_df[temp_cols_list].values*values_scaler)) )



ax.legend(fontsize=10,ncol=1,facecolor='white', framealpha=1.0)

ax.set_xlim(-10,110)

ax.set_ylim(0,1000)

ax.xaxis.set_tick_params(labelsize=12)
ax.yaxis.set_tick_params(labelsize=12)

ax.set_xlabel("Methylation percentage", fontsize=12)
ax.set_ylabel("Number of loci", fontsize=12)

#fig.tight_layout()
#fig.subplots_adjust(top=0.5)
hspace = 0.
fig.savefig( os.path.join(GraphFolder,filename_start_Graph +"Mfrac_dist_"+fig_label_Expt+".png"), bbox_inches = 'tight')
fig.savefig( os.path.join(GraphFolder,filename_start_Graph +"Mfrac_dist_"+fig_label_Expt+".pdf"), bbox_inches = 'tight')

plt.show()
# End Graph
