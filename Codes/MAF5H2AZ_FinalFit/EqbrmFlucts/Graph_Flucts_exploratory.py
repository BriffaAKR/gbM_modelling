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

n_cc = params.n_cc

N_gen_total = params.N_gen_total
# convert to cell cycles
N_gen_total_cc = N_gen_total*n_cc

TimeAverage_StartGen = params.TimeAverage_StartGen

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

N_ExptState_Decile3 = 10
N_ExptState_Decile4 = 10
N_ExptState_Decile5 = 10

n_bin = 50
GraphFolder = "Graphs"

filename_start_Sim = 'MAF5_H2AZWT_EqbrmOutput_FinalFit_Flucts_Col0likeAccnsNonRedCov50_Sim_Pchoice_Excl_InitialState_100U_'
filename_start_Dat = 'MAF5_H2AZWT_EqbrmOutput_FinalFit_Flucts_Col0likeAccnsNonRedCov50_Dat_Pchoice_Excl_InitialState_100U_'
filename_start_Graph = 'MAF5_H2AZWT_EqbrmOutput_FinalFit_Flucts_Col0likeAccnsNonRedCov50_Pchoice_Excl_InitialState_100U_'

# filename_start_Sim = params.filename_start_Sim
# filename_start_Dat = params.filename_start_Dat
# filename_start_Graph = params.filename_start_Graph
filename_ending = '.tsv'

# file_in_LocusProperties_Data = os.path.join('Output_files', filename_start_Dat + 'LocusProperties'+filename_ending)

# LocusProperties_Data_df = pd.read_csv(file_in_LocusProperties_Data, sep="\t{1}",engine='python')
# LocusProperties_Data_df.set_index("gene_ID", inplace = True)
# N_gene = len(LocusProperties_Data_df)
# IDs_list = LocusProperties_Data_df.index.values.tolist()


file_in_LocusProperties_Sim = os.path.join('Output_files', filename_start_Sim + 'LocusProperties'+filename_ending)

LocusProperties_Sim_df = pd.read_csv(file_in_LocusProperties_Sim, sep="\t{1}",engine='python')
LocusProperties_Sim_df.set_index("gene_ID", inplace = True)
N_gene = len(LocusProperties_Sim_df)
IDs_list = LocusProperties_Sim_df.index.values.tolist()
print('N_gene', N_gene)







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


print('N_CG', 'Counts')
length_val = []
loci_number = []

for i_ in range(5,350+1):
    length_val.append(i_)
    loci_number_temp = LocusProperties_Sim_df.loc[ LocusProperties_Sim_df['N_CG'] == i_ ].shape[0]
    loci_number.append(loci_number_temp)
    print(length_val[-1], loci_number[-1])

TimeScaleFactor = 1E-3
TimeScaleFactorLabel = '(1000 gens.)'

###########

# Start Graph
fig, ax = plt.subplots(1,2,figsize=(6,3))

for spine in ['left','right','top','bottom']:
    ax[0].spines[spine].set_color('k')
    ax[0].spines[spine].set_linewidth(0.8)
    ax[1].spines[spine].set_color('k')
    ax[1].spines[spine].set_linewidth(0.8)

ax[0].set_facecolor('white')
ax[1].set_facecolor('white')


#ax.grid(False)
ax[0].grid(b=True, which='major', color='lightgrey', linestyle=':',linewidth=1)
ax[1].grid(b=True, which='major', color='lightgrey', linestyle=':',linewidth=1)

ax[0].scatter(LocusProperties_Sim_df['N_CG'].values,
            LocusProperties_Sim_df['Sim_Mfrac_FirstTimeSteadyState_AccnMu'].values*TimeScaleFactor/n_cc, s=1,color='k')
ax[1].scatter(LocusProperties_Sim_df['N_CG'].values,
            LocusProperties_Sim_df['Sim_Mfrac_FirstTimeSteadyState_AccnSigma'].values*TimeScaleFactor/n_cc, s=1,color='k')

# ax.plot(x_vals_corrns, PairSep_MM_Sim_array,  linewidth = 1, color='k',label='Sim.')
# ax.fill_between(x_vals_corrns, PairSep_MM_Data_array,PairSep_MM_Data_array+PairSep_XX_Data_array,
#                     alpha = 0.5, color='k',)

# from mpl_toolkits.axes_grid1.inset_locator import inset_axes
# # gains inset
# axins_0 = inset_axes(ax, width="60%", height="50%", borderpad=1)
# axins_0.tick_params(axis='both', which='major', labelsize=14)
# axins_0.set_xlim(0,20)
# axins_0.set_ylim(0,0.02)
# axins_0.plot(x_vals_corrns, PairSep_MM_Sim_array,  linewidth = 1, color='k',label='$MM$-pairs Sim.')
# axins_0.fill_between(x_vals_corrns, PairSep_MM_Data_array,PairSep_MM_Data_array+PairSep_XX_Data_array,
#                     alpha = 0.5, color='green',label='$MM$-pairs Data')



#fig.suptitle(title_sting_single_time, fontsize=12)
#ax.legend(fontsize=12,ncol=1,facecolor='white', framealpha=1.0)

# ax.set_xlim(0,150)
# ax.set_ylim(0,0.02)

ax[0].xaxis.set_tick_params(labelsize=10)
ax[0].yaxis.set_tick_params(labelsize=10)
ax[0].set_xlabel("Locus length (CG sites)", fontsize=10)
ax[0].set_ylabel("Mean time to first reach steady-state\nfrom all-$U$ initial state "+TimeScaleFactorLabel, fontsize=10)
ax[1].xaxis.set_tick_params(labelsize=10)
ax[1].yaxis.set_tick_params(labelsize=10)
ax[1].set_xlabel("Locus length (CG sites)", fontsize=10)
ax[1].set_ylabel("sd of time (over locus replicates)\nto first reach steady-state\nfrom all-$U$ initial state "+TimeScaleFactorLabel, fontsize=10)

ax[0].set_title("Mean over all loci: %.2f" % np.mean(LocusProperties_Sim_df['Sim_Mfrac_FirstTimeSteadyState_AccnMu'].values*TimeScaleFactor/n_cc),fontsize=10)
ax[1].set_title("Mean over all loci: %.2f" % np.mean(LocusProperties_Sim_df['Sim_Mfrac_FirstTimeSteadyState_AccnSigma'].values*TimeScaleFactor/n_cc),fontsize=10)

fig.tight_layout()
#fig.subplots_adjust(top=0.5)

fig.savefig( os.path.join(GraphFolder,filename_start_Graph +"FirstSteadyStateTime_scatter.png"), bbox_inches = 'tight')


plt.show()
# End Graph
###########

# Start Graph
fig, ax = plt.subplots(1,2,figsize=(6,3))

for spine in ['left','right','top','bottom']:
    ax[0].spines[spine].set_color('k')
    ax[0].spines[spine].set_linewidth(0.8)
    ax[1].spines[spine].set_color('k')
    ax[1].spines[spine].set_linewidth(0.8)

ax[0].set_facecolor('white')
ax[1].set_facecolor('white')


#ax.grid(False)
ax[0].grid(b=True, which='major', color='lightgrey', linestyle=':',linewidth=1)
ax[1].grid(b=True, which='major', color='lightgrey', linestyle=':',linewidth=1)

ax[0].scatter(LocusProperties_Sim_df['N_CG'].values,
            LocusProperties_Sim_df['Sim_Mfrac_FirstTimeSteadyState_AccnMu'].values*TimeScaleFactor/n_cc, s=1,color='k')
ax[1].scatter(LocusProperties_Sim_df['N_CG'].values,
            LocusProperties_Sim_df['Sim_Mfrac_FirstTimeSteadyState_AccnSigma'].values*TimeScaleFactor/n_cc, s=1,color='k')

# ax.plot(x_vals_corrns, PairSep_MM_Sim_array,  linewidth = 1, color='k',label='Sim.')
# ax.fill_between(x_vals_corrns, PairSep_MM_Data_array,PairSep_MM_Data_array+PairSep_XX_Data_array,
#                     alpha = 0.5, color='k',)

# from mpl_toolkits.axes_grid1.inset_locator import inset_axes
# # gains inset
# axins_0 = inset_axes(ax, width="60%", height="50%", borderpad=1)
# axins_0.tick_params(axis='both', which='major', labelsize=14)
# axins_0.set_xlim(0,20)
# axins_0.set_ylim(0,0.02)
# axins_0.plot(x_vals_corrns, PairSep_MM_Sim_array,  linewidth = 1, color='k',label='$MM$-pairs Sim.')
# axins_0.fill_between(x_vals_corrns, PairSep_MM_Data_array,PairSep_MM_Data_array+PairSep_XX_Data_array,
#                     alpha = 0.5, color='green',label='$MM$-pairs Data')



#fig.suptitle(title_sting_single_time, fontsize=12)
#ax.legend(fontsize=12,ncol=1,facecolor='white', framealpha=1.0)

ax[0].set_xlim(0,40)
ax[1].set_xlim(0,40)
# ax.set_ylim(0,0.02)

ax[0].xaxis.set_tick_params(labelsize=10)
ax[0].yaxis.set_tick_params(labelsize=10)
ax[0].set_xlabel("Locus length (CG sites)", fontsize=10)
ax[0].set_ylabel("Mean time (over locus replicates)\nto first reach steady-state\nfrom all-$U$ initial state "+TimeScaleFactorLabel, fontsize=10)
ax[1].xaxis.set_tick_params(labelsize=10)
ax[1].yaxis.set_tick_params(labelsize=10)
ax[1].set_xlabel("Locus length (CG sites)", fontsize=10)
ax[1].set_ylabel("Std of time (over locus replicates)\nto first reach steady-state\nfrom all-$U$ initial state "+TimeScaleFactorLabel, fontsize=10)

ax[0].set_title("Mean over all loci: %.2f" % np.mean(LocusProperties_Sim_df['Sim_Mfrac_FirstTimeSteadyState_AccnMu'].values*TimeScaleFactor/n_cc),fontsize=10)
ax[1].set_title("Mean over all loci: %.2f" % np.mean(LocusProperties_Sim_df['Sim_Mfrac_FirstTimeSteadyState_AccnSigma'].values*TimeScaleFactor/n_cc),fontsize=10)


fig.tight_layout()
#fig.subplots_adjust(top=0.5)

fig.savefig( os.path.join(GraphFolder,filename_start_Graph +"FirstSteadyStateTime_scatterA.png"), bbox_inches = 'tight')


plt.show()
# End Graph

###########

# Start Graph
fig, ax = plt.subplots(1,1,figsize=(3,3))

for spine in ['left','right','top','bottom']:
    ax.spines[spine].set_color('k')
    ax.spines[spine].set_linewidth(0.8)

ax.set_facecolor('white')


#ax.grid(False)
ax.grid(b=True, which='major', color='lightgrey', linestyle=':',linewidth=1)

ax.scatter(LocusProperties_Sim_df['N_CG'].values,
            LocusProperties_Sim_df['Sim_Mfrac_star_TimeSigma_AccnMu'].values, s=1,color='k')

# ax.plot(x_vals_corrns, PairSep_MM_Sim_array,  linewidth = 1, color='k',label='Sim.')
# ax.fill_between(x_vals_corrns, PairSep_MM_Data_array,PairSep_MM_Data_array+PairSep_XX_Data_array,
#                     alpha = 0.5, color='k',)

from mpl_toolkits.axes_grid1.inset_locator import inset_axes
# gains inset
axins_0 = inset_axes(ax, width="60%", height="50%", borderpad=0.6)
axins_0.tick_params(axis='both', which='major', labelsize=10)
axins_0.set_xlim(0,20)
#axins_0.set_ylim(0,0.02)
# axins_0.plot(x_vals_corrns, PairSep_MM_Sim_array,  linewidth = 1, color='k',label='$MM$-pairs Sim.')
# axins_0.fill_between(x_vals_corrns, PairSep_MM_Data_array,PairSep_MM_Data_array+PairSep_XX_Data_array,
#                     alpha = 0.5, color='green',label='$MM$-pairs Data')
axins_0.scatter(LocusProperties_Sim_df['N_CG'].values,
            LocusProperties_Sim_df['Sim_Mfrac_star_TimeSigma_AccnMu'].values, s=1,color='k')



#fig.suptitle(title_sting_single_time, fontsize=12)
#ax.legend(fontsize=12,ncol=1,facecolor='white', framealpha=1.0)

# ax.set_xlim(0,40)
# ax.set_ylim(0,0.02)

ax.xaxis.set_tick_params(labelsize=10)
ax.yaxis.set_tick_params(labelsize=10)
ax.set_xlabel("Locus length (CG sites)", fontsize=10)
ax.set_ylabel("Typical fluctuation magnitude OR\nsd of mCG fluctuations about\nthe mean steady-state level\n(methylation percentage)", fontsize=10)

ax.set_title("Mean over all loci: %.2f" % np.mean(LocusProperties_Sim_df['Sim_Mfrac_star_TimeSigma_AccnMu'].values),fontsize=10)


fig.tight_layout()
#fig.subplots_adjust(top=0.5)

fig.savefig( os.path.join(GraphFolder,filename_start_Graph +"Mfrac_star_TimeSigma_scatter.png"), bbox_inches = 'tight')


plt.show()
# End Graph


###########

# Start Graph
fig, ax = plt.subplots(1,2,figsize=(6,3))

for spine in ['left','right','top','bottom']:
    ax[0].spines[spine].set_color('k')
    ax[0].spines[spine].set_linewidth(0.8)
    ax[1].spines[spine].set_color('k')
    ax[1].spines[spine].set_linewidth(0.8)

ax[0].set_facecolor('white')
ax[1].set_facecolor('white')



#ax.grid(False)
ax[0].grid(b=True, which='major', color='lightgrey', linestyle=':',linewidth=1)
ax[1].grid(b=True, which='major', color='lightgrey', linestyle=':',linewidth=1)

ax[0].scatter(LocusProperties_Sim_df['N_CG'].values,
            LocusProperties_Sim_df['Sim_Mfrac_MaxFluct_Magnitude_AccnMu'].values, s=1,color='k')
ax[1].scatter(LocusProperties_Sim_df['N_CG'].values,
            LocusProperties_Sim_df['Sim_Mfrac_MaxFluct_Magnitude_AccnSigma'].values, s=1,color='k')

# ax.plot(x_vals_corrns, PairSep_MM_Sim_array,  linewidth = 1, color='k',label='Sim.')
# ax.fill_between(x_vals_corrns, PairSep_MM_Data_array,PairSep_MM_Data_array+PairSep_XX_Data_array,
#                     alpha = 0.5, color='k',)

from mpl_toolkits.axes_grid1.inset_locator import inset_axes
# gains inset
axins_0 = inset_axes(ax[0], width="60%", height="50%", borderpad=0.6)
axins_0.tick_params(axis='both', which='major', labelsize=10)
axins_0.set_xlim(0,20)
#axins_0.set_ylim(0,0.02)
# axins_0.plot(x_vals_corrns, PairSep_MM_Sim_array,  linewidth = 1, color='k',label='$MM$-pairs Sim.')
# axins_0.fill_between(x_vals_corrns, PairSep_MM_Data_array,PairSep_MM_Data_array+PairSep_XX_Data_array,
#                     alpha = 0.5, color='green',label='$MM$-pairs Data')
axins_0.scatter(LocusProperties_Sim_df['N_CG'].values,
            LocusProperties_Sim_df['Sim_Mfrac_MaxFluct_Magnitude_AccnMu'].values, s=1,color='k')

axins_1 = inset_axes(ax[1], width="60%", height="50%", borderpad=0.6)
axins_1.tick_params(axis='both', which='major', labelsize=10)
axins_1.set_xlim(0,20)
#axins_0.set_ylim(0,0.02)
# axins_1.plot(x_vals_corrns, PairSep_MM_Sim_array,  linewidth = 1, color='k',label='$MM$-pairs Sim.')
# axins_1.fill_between(x_vals_corrns, PairSep_MM_Data_array,PairSep_MM_Data_array+PairSep_XX_Data_array,
#                     alpha = 0.5, color='green',label='$MM$-pairs Data')
axins_1.scatter(LocusProperties_Sim_df['N_CG'].values,
            LocusProperties_Sim_df['Sim_Mfrac_MaxFluct_Magnitude_AccnSigma'].values, s=1,color='k')



#fig.suptitle(title_sting_single_time, fontsize=12)
#ax.legend(fontsize=12,ncol=1,facecolor='white', framealpha=1.0)

# ax.set_xlim(0,150)
# ax.set_ylim(0,0.02)

ax[0].xaxis.set_tick_params(labelsize=10)
ax[0].yaxis.set_tick_params(labelsize=10)
ax[0].set_xlabel("Locus length (CG sites)", fontsize=10)
ax[0].set_ylabel("Mean magnitude (over locus replicates)\nof greatest mCG fluctuation about\nthe mean steady-state level\n(methylation percentage)", fontsize=10)
ax[1].xaxis.set_tick_params(labelsize=10)
ax[1].yaxis.set_tick_params(labelsize=10)
ax[1].set_xlabel("Locus length (CG sites)", fontsize=10)
ax[1].set_ylabel("sd (over locus replicates) of greatest\nmagnitude mCG fluctuation about\nthe mean steady-state level\n(methylation percentage)", fontsize=10)

ax[0].set_title("Mean over all loci: %.2f" % np.mean(LocusProperties_Sim_df['Sim_Mfrac_MaxFluct_Magnitude_AccnMu'].values),fontsize=10)
ax[1].set_title("Mean over all loci: %.2f" % np.mean(LocusProperties_Sim_df['Sim_Mfrac_MaxFluct_Magnitude_AccnSigma'].values),fontsize=10)


fig.tight_layout()
#fig.subplots_adjust(top=0.5)

fig.savefig( os.path.join(GraphFolder,filename_start_Graph +"MaxFluct_Magnitude_scatter.png"), bbox_inches = 'tight')


plt.show()
# End Graph
###########

# Start Graph
fig, ax = plt.subplots(1,2,figsize=(6,3))

for spine in ['left','right','top','bottom']:
    ax[0].spines[spine].set_color('k')
    ax[0].spines[spine].set_linewidth(0.8)
    ax[1].spines[spine].set_color('k')
    ax[1].spines[spine].set_linewidth(0.8)

ax[0].set_facecolor('white')
ax[1].set_facecolor('white')


#ax.grid(False)
ax[0].grid(b=True, which='major', color='lightgrey', linestyle=':',linewidth=1)
ax[1].grid(b=True, which='major', color='lightgrey', linestyle=':',linewidth=1)

ax[0].scatter(LocusProperties_Sim_df['N_CG'].values,
            LocusProperties_Sim_df['Sim_Mfrac_MaxFluct_Duration_AccnMu'].values*TimeScaleFactor/n_cc, s=1,color='k')
ax[1].scatter(LocusProperties_Sim_df['N_CG'].values,
            LocusProperties_Sim_df['Sim_Mfrac_MaxFluct_Duration_AccnSigma'].values*TimeScaleFactor/n_cc, s=1,color='k')

# ax.plot(x_vals_corrns, PairSep_MM_Sim_array,  linewidth = 1, color='k',label='Sim.')
# ax.fill_between(x_vals_corrns, PairSep_MM_Data_array,PairSep_MM_Data_array+PairSep_XX_Data_array,
#                     alpha = 0.5, color='k',)

# from mpl_toolkits.axes_grid1.inset_locator import inset_axes
# # gains inset
# axins_0 = inset_axes(ax, width="60%", height="50%", borderpad=1)
# axins_0.tick_params(axis='both', which='major', labelsize=14)
# axins_0.set_xlim(0,20)
# axins_0.set_ylim(0,0.02)
# axins_0.plot(x_vals_corrns, PairSep_MM_Sim_array,  linewidth = 1, color='k',label='$MM$-pairs Sim.')
# axins_0.fill_between(x_vals_corrns, PairSep_MM_Data_array,PairSep_MM_Data_array+PairSep_XX_Data_array,
#                     alpha = 0.5, color='green',label='$MM$-pairs Data')



#fig.suptitle(title_sting_single_time, fontsize=12)
#ax.legend(fontsize=12,ncol=1,facecolor='white', framealpha=1.0)

# ax[0].set_xlim(0,40)
# ax[1].set_xlim(0,40)
# ax.set_ylim(0,0.02)

ax[0].xaxis.set_tick_params(labelsize=10)
ax[0].yaxis.set_tick_params(labelsize=10)
ax[0].set_xlabel("Locus length (CG sites)", fontsize=10)
ax[0].set_ylabel('Mean durration (over locus replicates)\nof greatest mCG fluctuation about\nthe mean steady-state level '+TimeScaleFactorLabel, fontsize=10)
ax[1].xaxis.set_tick_params(labelsize=10)
ax[1].yaxis.set_tick_params(labelsize=10)
ax[1].set_xlabel("Locus length (CG sites)", fontsize=10)
ax[1].set_ylabel("sd (over locus replicates) of durration of greatest\nmagnitude mCG fluctuation about\nthe mean steady-state level "+TimeScaleFactorLabel, fontsize=10)

ax[0].set_title("Mean over all loci: %.2f" % np.mean(LocusProperties_Sim_df['Sim_Mfrac_MaxFluct_Duration_AccnMu'].values*TimeScaleFactor/n_cc),fontsize=10)
ax[1].set_title("Mean over all loci: %.2f" % np.mean(LocusProperties_Sim_df['Sim_Mfrac_MaxFluct_Duration_AccnSigma'].values*TimeScaleFactor/n_cc),fontsize=10)


fig.tight_layout()
#fig.subplots_adjust(top=0.5)

fig.savefig( os.path.join(GraphFolder,filename_start_Graph +"MaxFluct_Duration_scatter.png"), bbox_inches = 'tight')


plt.show()
# End Graph

######

# Start Graph
fig, ax = plt.subplots(1,2,figsize=(6,3))

for spine in ['left','right','top','bottom']:
    ax[0].spines[spine].set_color('k')
    ax[0].spines[spine].set_linewidth(0.8)
    ax[1].spines[spine].set_color('k')
    ax[1].spines[spine].set_linewidth(0.8)

ax[0].set_facecolor('white')
ax[1].set_facecolor('white')


#ax.grid(False)
ax[0].grid(b=True, which='major', color='lightgrey', linestyle=':',linewidth=1)
ax[1].grid(b=True, which='major', color='lightgrey', linestyle=':',linewidth=1)

ax[0].scatter(LocusProperties_Sim_df['N_CG'].values,
            LocusProperties_Sim_df['Sim_Mfrac_MaxFluct_Duration_AccnMu'].values*TimeScaleFactor/n_cc, s=1,color='k')
ax[1].scatter(LocusProperties_Sim_df['N_CG'].values,
            LocusProperties_Sim_df['Sim_Mfrac_MaxFluct_Duration_AccnSigma'].values*TimeScaleFactor/n_cc, s=1,color='k')

# ax.plot(x_vals_corrns, PairSep_MM_Sim_array,  linewidth = 1, color='k',label='Sim.')
# ax.fill_between(x_vals_corrns, PairSep_MM_Data_array,PairSep_MM_Data_array+PairSep_XX_Data_array,
#                     alpha = 0.5, color='k',)

# from mpl_toolkits.axes_grid1.inset_locator import inset_axes
# # gains inset
# axins_0 = inset_axes(ax, width="60%", height="50%", borderpad=1)
# axins_0.tick_params(axis='both', which='major', labelsize=14)
# axins_0.set_xlim(0,20)
# axins_0.set_ylim(0,0.02)
# axins_0.plot(x_vals_corrns, PairSep_MM_Sim_array,  linewidth = 1, color='k',label='$MM$-pairs Sim.')
# axins_0.fill_between(x_vals_corrns, PairSep_MM_Data_array,PairSep_MM_Data_array+PairSep_XX_Data_array,
#                     alpha = 0.5, color='green',label='$MM$-pairs Data')



#fig.suptitle(title_sting_single_time, fontsize=12)
#ax.legend(fontsize=12,ncol=1,facecolor='white', framealpha=1.0)

ax[0].set_xlim(0,40)
ax[1].set_xlim(0,40)
# ax.set_ylim(0,0.02)

ax[0].xaxis.set_tick_params(labelsize=10)
ax[0].yaxis.set_tick_params(labelsize=10)
ax[0].set_xlabel("Locus length (CG sites)", fontsize=10)
ax[0].set_ylabel('Mean durration (over locus replicates)\nof greatest mCG fluctuation about\nthe mean steady-state level '+TimeScaleFactorLabel, fontsize=10)
ax[1].xaxis.set_tick_params(labelsize=10)
ax[1].yaxis.set_tick_params(labelsize=10)
ax[1].set_xlabel("Locus length (CG sites)", fontsize=10)
ax[1].set_ylabel("sd (over locus replicates) of durration of greatest\nmagnitude mCG fluctuation about\nthe mean steady-state level "+TimeScaleFactorLabel, fontsize=10)

ax[0].set_title("Mean over all loci: %.2f" % np.mean(LocusProperties_Sim_df['Sim_Mfrac_MaxFluct_Duration_AccnMu'].values*TimeScaleFactor/n_cc),fontsize=10)
ax[1].set_title("Mean over all loci: %.2f" % np.mean(LocusProperties_Sim_df['Sim_Mfrac_MaxFluct_Duration_AccnSigma'].values*TimeScaleFactor/n_cc),fontsize=10)


fig.tight_layout()
#fig.subplots_adjust(top=0.5)

fig.savefig( os.path.join(GraphFolder,filename_start_Graph +"MaxFluct_Duration_scatterA.png"), bbox_inches = 'tight')


plt.show()
# End Graph

######

N_CG_vals = np.linspace(5,500,500-5+1)
FirstMTimes_vals = []
Gillespie_FirstMTimes_vals = []
for i_ in range(len(N_CG_vals)):
    FirstMTimes_vals.append( 1./(delta*N_CG_vals[i_]) )
    Gillespie_FirstMTimes_vals.append( (-1./(delta))*np.log( N_CG_vals[i_]/(N_CG_vals[i_]+1.) ) )
FirstMTimes_vals = np.array(FirstMTimes_vals)
Gillespie_FirstMTimes_vals = np.array(Gillespie_FirstMTimes_vals)



# Start Graph
fig, ax = plt.subplots(1,2,figsize=(6,3))

for spine in ['left','right','top','bottom']:
    ax[0].spines[spine].set_color('k')
    ax[0].spines[spine].set_linewidth(0.8)
    ax[1].spines[spine].set_color('k')
    ax[1].spines[spine].set_linewidth(0.8)

ax[0].set_facecolor('white')
ax[1].set_facecolor('white')

#ax.grid(False)
ax[0].grid(b=True, which='major', color='lightgrey', linestyle=':',linewidth=1)
ax[1].grid(b=True, which='major', color='lightgrey', linestyle=':',linewidth=1)

ax[0].scatter(LocusProperties_Sim_df['N_CG'].values,
             LocusProperties_Sim_df['Sim_FirstTimeChangedSite_AccnMu'].values/n_cc, s=1,color='k')
ax[1].scatter(LocusProperties_Sim_df['N_CG'].values,
            LocusProperties_Sim_df['Sim_FirstTimeChangedSite_AccnSigma'].values/n_cc, s=1,color='k')

# ax[0].errorbar(LocusProperties_Sim_df['N_CG'].values,
#             LocusProperties_Sim_df['Sim_FirstTimeChangedSite_AccnMu'].values/n_cc,
            #  yerr=LocusProperties_Sim_df['Sim_FirstTimeChangedSite_AccnSigma'].values/n_cc, fmt="o",markersize=1)

ax[0].plot(N_CG_vals, FirstMTimes_vals/n_cc,  linewidth = 1, color='r')
ax[1].plot(N_CG_vals, FirstMTimes_vals/n_cc,  linewidth = 1, color='r')

# ax[0].plot(N_CG_vals, Gillespie_FirstMTimes_vals/n_cc,  linewidth = 1, color='b')


# ax.plot(x_vals_corrns, PairSep_MM_Sim_array,  linewidth = 1, color='k',label='Sim.')
# ax.fill_between(x_vals_corrns, PairSep_MM_Data_array,PairSep_MM_Data_array+PairSep_XX_Data_array,
#                     alpha = 0.5, color='k',)

# from mpl_toolkits.axes_grid1.inset_locator import inset_axes
# # gains inset
# axins_0 = inset_axes(ax, width="60%", height="50%", borderpad=1)
# axins_0.tick_params(axis='both', which='major', labelsize=14)
# axins_0.set_xlim(0,20)
# axins_0.set_ylim(0,0.02)
# axins_0.plot(x_vals_corrns, PairSep_MM_Sim_array,  linewidth = 1, color='k',label='$MM$-pairs Sim.')
# axins_0.fill_between(x_vals_corrns, PairSep_MM_Data_array,PairSep_MM_Data_array+PairSep_XX_Data_array,
#                     alpha = 0.5, color='green',label='$MM$-pairs Data')



#fig.suptitle(title_sting_single_time, fontsize=12)
#ax.legend(fontsize=12,ncol=1,facecolor='white', framealpha=1.0)

# ax.set_xlim(0,150)
ax[0].set_ylim(0,500)
ax[1].set_ylim(0,500)


ax[0].xaxis.set_tick_params(labelsize=10)
ax[0].yaxis.set_tick_params(labelsize=10)
ax[0].set_xlabel("Locus length (CG sites)", fontsize=10)
ax[0].set_ylabel("Mean time (over locus replicates)\nto first $M$ gain (gens.)", fontsize=10)
ax[1].xaxis.set_tick_params(labelsize=10)
ax[1].yaxis.set_tick_params(labelsize=10)
ax[1].set_xlabel("Locus length (CG sites)", fontsize=10)
ax[1].set_ylabel("sd of time (over locus replicates)\nto first $M$ gain (gens.)", fontsize=10)

ax[0].set_title("Mean over all loci: %.2f" % np.mean(LocusProperties_Sim_df['Sim_FirstTimeChangedSite_AccnMu'].values/n_cc),fontsize=10)
ax[1].set_title("Mean over all loci: %.2f" % np.mean(LocusProperties_Sim_df['Sim_FirstTimeChangedSite_AccnSigma'].values/n_cc),fontsize=10)

fig.tight_layout()
#fig.subplots_adjust(top=0.5)

fig.savefig( os.path.join(GraphFolder,filename_start_Graph +"FirstTimeChangedSite_scatter.png"), bbox_inches = 'tight')


plt.show()
# End Graph


###############
# Start Graph
fig, ax = plt.subplots(1,2,figsize=(6,3))

for spine in ['left','right','top','bottom']:
    ax[0].spines[spine].set_color('k')
    ax[0].spines[spine].set_linewidth(0.8)
    ax[1].spines[spine].set_color('k')
    ax[1].spines[spine].set_linewidth(0.8)

ax[0].set_facecolor('white')
ax[1].set_facecolor('white')

#ax.grid(False)
ax[0].grid(b=True, which='major', color='lightgrey', linestyle=':',linewidth=1)
ax[1].grid(b=True, which='major', color='lightgrey', linestyle=':',linewidth=1)

ax[0].scatter(LocusProperties_Sim_df['N_CG'].values,
             LocusProperties_Sim_df['Sim_FirstTimeChangedSite_AccnMu'].values/n_cc, s=1,color='k')
ax[1].scatter(LocusProperties_Sim_df['N_CG'].values,
            LocusProperties_Sim_df['Sim_FirstTimeChangedSite_AccnSigma'].values/n_cc, s=1,color='k')

# ax[0].errorbar(LocusProperties_Sim_df['N_CG'].values,
#             LocusProperties_Sim_df['Sim_FirstTimeChangedSite_AccnMu'].values/n_cc,
            #  yerr=LocusProperties_Sim_df['Sim_FirstTimeChangedSite_AccnSigma'].values/n_cc, fmt="o",markersize=1)

ax[0].plot(N_CG_vals, FirstMTimes_vals/n_cc,  linewidth = 1, color='r')
ax[1].plot(N_CG_vals, FirstMTimes_vals/n_cc,  linewidth = 1, color='r')

# ax[0].plot(N_CG_vals, Gillespie_FirstMTimes_vals/n_cc,  linewidth = 1, color='b')


# ax.plot(x_vals_corrns, PairSep_MM_Sim_array,  linewidth = 1, color='k',label='Sim.')
# ax.fill_between(x_vals_corrns, PairSep_MM_Data_array,PairSep_MM_Data_array+PairSep_XX_Data_array,
#                     alpha = 0.5, color='k',)

# from mpl_toolkits.axes_grid1.inset_locator import inset_axes
# # gains inset
# axins_0 = inset_axes(ax, width="60%", height="50%", borderpad=1)
# axins_0.tick_params(axis='both', which='major', labelsize=14)
# axins_0.set_xlim(0,20)
# axins_0.set_ylim(0,0.02)
# axins_0.plot(x_vals_corrns, PairSep_MM_Sim_array,  linewidth = 1, color='k',label='$MM$-pairs Sim.')
# axins_0.fill_between(x_vals_corrns, PairSep_MM_Data_array,PairSep_MM_Data_array+PairSep_XX_Data_array,
#                     alpha = 0.5, color='green',label='$MM$-pairs Data')



#fig.suptitle(title_sting_single_time, fontsize=12)
#ax.legend(fontsize=12,ncol=1,facecolor='white', framealpha=1.0)

# ax[0].set_xlim(0,200)
# ax[0].set_ylim(0,500)
# ax[1].set_ylim(0,500)


ax[0].xaxis.set_tick_params(labelsize=10)
ax[0].yaxis.set_tick_params(labelsize=10)
ax[0].set_xlabel("Locus length (CG sites)", fontsize=10)
ax[0].set_ylabel("Mean time (over locus replicates)\nto first $M$ gain (gens.)", fontsize=10)
ax[1].xaxis.set_tick_params(labelsize=10)
ax[1].yaxis.set_tick_params(labelsize=10)
ax[1].set_xlabel("Locus length (CG sites)", fontsize=10)
ax[1].set_ylabel("sd of time (over locus replicates)\nto first $M$ gain (gens.)", fontsize=10)

ax[0].set_title("Mean over all loci: %.2f" % np.mean(LocusProperties_Sim_df['Sim_FirstTimeChangedSite_AccnMu'].values/n_cc),fontsize=10)
ax[1].set_title("Mean over all loci: %.2f" % np.mean(LocusProperties_Sim_df['Sim_FirstTimeChangedSite_AccnSigma'].values/n_cc),fontsize=10)

fig.tight_layout()
#fig.subplots_adjust(top=0.5)

fig.savefig( os.path.join(GraphFolder,filename_start_Graph +"FirstTimeChangedSite_scatterA.png"), bbox_inches = 'tight')


plt.show()
# End Graph