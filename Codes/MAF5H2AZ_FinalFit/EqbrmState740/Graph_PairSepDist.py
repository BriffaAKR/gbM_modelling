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

filename_start_Sim = params.filename_start_Sim
filename_start_Dat = params.filename_start_Dat
filename_start_Graph = params.filename_start_Graph
filename_ending = '.tsv'

file_in_LocusProperties_Data = os.path.join('Output_files', filename_start_Dat + 'LocusProperties'+filename_ending)

file_in_PairSep_MM_Data = os.path.join('Output_files', filename_start_Dat + 'PairSeparationsMM_Data'+filename_ending)
file_in_PairSep_MU_Data = os.path.join('Output_files', filename_start_Dat + 'PairSeparationsMU_Data'+filename_ending)
file_in_PairSep_UU_Data = os.path.join('Output_files', filename_start_Dat + 'PairSeparationsUU_Data'+filename_ending)
file_in_PairSep_XX_Data = os.path.join('Output_files', filename_start_Dat + 'PairSeparationsXX_Data'+filename_ending)
file_in_PairSep_MM_Sim = os.path.join('Output_files', filename_start_Sim + 'PairSeparationsMM_Sim'+filename_ending)
file_in_PairSep_MU_Sim = os.path.join('Output_files', filename_start_Sim + 'PairSeparationsMU_Sim'+filename_ending)
file_in_PairSep_UU_Sim = os.path.join('Output_files', filename_start_Sim + 'PairSeparationsUU_Sim'+filename_ending)

LocusProperties_df = pd.read_csv(file_in_LocusProperties_Data, sep="\t{1}",engine='python')
LocusProperties_df.set_index("gene_ID", inplace = True)
N_gene = len(LocusProperties_df)
IDs_list = LocusProperties_df.index.values.tolist()

# add up total number of pairs
Total_Pairs_Data = 0
Total_Pairs_Sim = 0

PairSep_MM_Data_df = pd.read_csv(file_in_PairSep_MM_Data, sep="\t{1}",engine='python',header=None)
print(PairSep_MM_Data_df.head())
print()
print(PairSep_MM_Data_df.sum(axis=0))
PairSep_MM_Data_array = PairSep_MM_Data_df.sum(axis=0).values
print()
print(PairSep_MM_Data_array[:20])
print(len(PairSep_MM_Data_array),np.sum(PairSep_MM_Data_array))
Total_Pairs_Data += np.sum(PairSep_MM_Data_array)
PairSep_MU_Data_df = pd.read_csv(file_in_PairSep_MU_Data, sep="\t{1}",engine='python',header=None)
PairSep_MU_Data_array = PairSep_MU_Data_df.sum(axis=0).values
Total_Pairs_Data += np.sum(PairSep_MU_Data_array)
PairSep_UU_Data_df = pd.read_csv(file_in_PairSep_UU_Data, sep="\t{1}",engine='python',header=None)
PairSep_UU_Data_array = PairSep_UU_Data_df.sum(axis=0).values
Total_Pairs_Data += np.sum(PairSep_UU_Data_array)
PairSep_XX_Data_df = pd.read_csv(file_in_PairSep_XX_Data, sep="\t{1}",engine='python',header=None)
PairSep_XX_Data_array = PairSep_XX_Data_df.sum(axis=0).values
Total_Pairs_Data += np.sum(PairSep_XX_Data_array)

PairSep_MM_Sim_df = pd.read_csv(file_in_PairSep_MM_Sim, sep="\t{1}",engine='python',header=None)
PairSep_MM_Sim_array = PairSep_MM_Sim_df.sum(axis=0).values
Total_Pairs_Sim += np.sum(PairSep_MM_Sim_array)
PairSep_MU_Sim_df = pd.read_csv(file_in_PairSep_MU_Sim, sep="\t{1}",engine='python',header=None)
PairSep_MU_Sim_array = PairSep_MU_Sim_df.sum(axis=0).values
Total_Pairs_Sim += np.sum(PairSep_MU_Sim_array)
PairSep_UU_Sim_df = pd.read_csv(file_in_PairSep_UU_Sim, sep="\t{1}",engine='python',header=None)
PairSep_UU_Sim_array = PairSep_UU_Sim_df.sum(axis=0).values
Total_Pairs_Sim += np.sum(PairSep_UU_Sim_array)
print()
print(Total_Pairs_Data, Total_Pairs_Sim, Total_Pairs_Data - Total_Pairs_Sim, 'len(AT4G23000) = 81', (Total_Pairs_Data - Total_Pairs_Sim)/N_reps)

# Normalise Pairs_Data and Pairs_Sim individually
PairSep_MM_Data_array = PairSep_MM_Data_array/Total_Pairs_Data
PairSep_MU_Data_array = PairSep_MU_Data_array/Total_Pairs_Data
PairSep_UU_Data_array = PairSep_UU_Data_array/Total_Pairs_Data
PairSep_XX_Data_array = PairSep_XX_Data_array/Total_Pairs_Data

PairSep_MM_Sim_array = PairSep_MM_Sim_array/Total_Pairs_Sim
PairSep_MU_Sim_array = PairSep_MU_Sim_array/Total_Pairs_Sim
PairSep_UU_Sim_array = PairSep_UU_Sim_array/Total_Pairs_Sim

x_vals_corrns = np.arange(N_corrns_bins)


title_sting_single_time = "No. sims. per gene = %d      Cell cycle duration = %.1f       No. of genes simulated = %d        Sim. time = %d gens.\n\
'P' choice = %s      $\epsilon$ = variable      $\gamma_0$ = varialbe      $r_{div\\,\gamma} =$ %d      $r_{\gamma}$ = %d      $\lambda_{\gamma} =$ %.3f\n\
$\delta$ = %.3e      $u_{M PL}^+ =$ varialbe      $r_{div.} =$ %d      $r_{plat.} =$ %d      $\lambda_{M PL} =$ %.3f      $u_{Scale\\,LR1}^+ =$ %.3f\n \
$u_{M PL\\,LR1}^+ =$ varialbe      $r_{div\\,LR1} =$ %d      $r_{plat.\\,LR1} =$ %d      $\lambda_{M\\,In\\,PL\\,LR1} =$ %.3f      $\lambda_{M\\,Out\\,PL\\,LR1} =$ %.3f      $r_{LR} =$ %d\n \
$m_u$ = %.2f     $c_u$ = %.2f     $m_{\epsilon}$ = %.2f     $c_{\epsilon}$ = %.2f    $m_{\gamma}$ = %.2f     $c_{\gamma}$ = %.2f     $\\rho_{cap}$ = 1./%.2f\n\
Initial seed = %d      Sim. initial state: %s      %d $\leq N_{CG\\,min}$ < %d      %.3f $\leq \\rho_{CG\\,min}$ < %.3f\
\nLocus type: %s      Annotation file: %s\n\
Filter 1: %s      Filter 2: %s      ExcludeFilter 1: %s\n\
$N_{bin}$ = %d"\
             % (N_reps, rep_time, N_gene, N_gen_burn_in,  P_choice,
                r_div_gamma, r_gamma, lambda_gamma, delta, 
               r_div, r_plat,lambda_coop, u_scale_val,
                r_div_LR1, r_plat_LR1, lambda_coopIn_LR1, lambda_coopOut_LR1,
                r_LR1, u_m_val,u_c_val,e_m_val,e_c_val,g_m_val,g_c_val,spacing_cap,
                initial_seed, initial_state_choice,N_CG_min,N_CG_max, N_CG_density_min,N_CG_density_max, locus_type, 
                file_in_annotation, file_in_anno_filt_1, file_in_anno_filt_2, file_in_exclude_filt_1,n_bin)


# Start Graph
fig, ax = plt.subplots(1,3,figsize=(12,3))

for spine in ['left','right','top','bottom']:
    ax[0].spines[spine].set_color('k')
    ax[0].spines[spine].set_linewidth(0.8)
    ax[1].spines[spine].set_color('k')
    ax[1].spines[spine].set_linewidth(0.8)
    ax[2].spines[spine].set_color('k')
    ax[2].spines[spine].set_linewidth(0.8)

ax[0].set_facecolor('white')
ax[1].set_facecolor('white')
ax[2].set_facecolor('white')



#ax[0].grid(False)
ax[0].grid(b=True, which='major', color='lightgrey', linestyle=':',linewidth=1)
ax[1].grid(b=True, which='major', color='lightgrey', linestyle=':',linewidth=1)
ax[2].grid(b=True, which='major', color='lightgrey', linestyle=':',linewidth=1)


ax[0].plot(x_vals_corrns, PairSep_MM_Sim_array,  linewidth = 1, color='k',label='$MM$-pairs Sim.')
ax[0].fill_between(x_vals_corrns, PairSep_MM_Data_array,PairSep_MM_Data_array+PairSep_XX_Data_array,  
                    alpha = 0.5, color='green',label='$MM$-pairs Data')

from mpl_toolkits.axes_grid1.inset_locator import inset_axes
# gains inset
axins_0 = inset_axes(ax[0], width="60%", height="50%", borderpad=1)
axins_0.tick_params(axis='both', which='major', labelsize=14)
axins_0.set_xlim(0,20)
axins_0.set_ylim(0,0.02)
axins_0.plot(x_vals_corrns, PairSep_MM_Sim_array,  linewidth = 1, color='k',label='$MM$-pairs Sim.')
axins_0.fill_between(x_vals_corrns, PairSep_MM_Data_array,PairSep_MM_Data_array+PairSep_XX_Data_array,  
                    alpha = 0.5, color='green',label='$MM$-pairs Data')


ax[1].plot(x_vals_corrns, PairSep_MU_Sim_array,  linewidth = 1, color='k',label='$MU$-pairs Sim.')
ax[1].fill_between(x_vals_corrns, PairSep_MU_Data_array,PairSep_MU_Data_array+PairSep_XX_Data_array,  
                    alpha = 0.5, color='green',label='$MU$-pairs Data')

ax[2].plot(x_vals_corrns, PairSep_UU_Sim_array,  linewidth = 1, color='k',label='Simulated')
ax[2].fill_between(x_vals_corrns, PairSep_UU_Data_array,PairSep_UU_Data_array+PairSep_XX_Data_array,  
                    alpha = 0.5, color='green',label='740 Col-like Data')

#fig.suptitle(title_sting_single_time, fontsize=12)
#ax[0].legend(fontsize=12,ncol=1,facecolor='white', framealpha=1.0)
#ax[1].legend(fontsize=12,ncol=1,facecolor='white', framealpha=1.0)
ax[2].legend(fontsize=12,ncol=1,facecolor='white', framealpha=1.0)

ax[0].set_xlim(0,150)
ax[0].set_ylim(0,0.02)

ax[1].set_xlim(0,150)
ax[1].set_ylim(0,0.01*2)

ax[2].set_xlim(0,150)
ax[2].set_ylim(0,0.01*2)

ax[0].xaxis.set_tick_params(labelsize=12)
ax[0].yaxis.set_tick_params(labelsize=12)
ax[1].xaxis.set_tick_params(labelsize=12)
ax[1].yaxis.set_tick_params(labelsize=12)
ax[2].xaxis.set_tick_params(labelsize=12)
ax[2].yaxis.set_tick_params(labelsize=12)

ax[0].set_xlabel("$MM$-Pair separation (bp)", fontsize=12)
ax[0].set_ylabel("Density", fontsize=12)
ax[1].set_xlabel("$MU$-Pair separation (bp)", fontsize=12)
ax[1].set_ylabel("Density", fontsize=12)
ax[2].set_xlabel("$UU$-Pair separation (bp)", fontsize=12)
ax[2].set_ylabel("Density", fontsize=12)

fig.tight_layout()
#fig.subplots_adjust(top=0.5)

#fig.savefig( os.path.join(GraphFolder,filename_start_Graph +"PairSepDist_SF9E.png"))
fig.savefig( os.path.join(GraphFolder,filename_start_Graph +"PairSepDist_SF9E.png"), bbox_inches = 'tight')
fig.savefig( os.path.join(GraphFolder,filename_start_Graph +"PairSepDist_SF9E.pdf"), bbox_inches = 'tight')


plt.show()
# End Graph


