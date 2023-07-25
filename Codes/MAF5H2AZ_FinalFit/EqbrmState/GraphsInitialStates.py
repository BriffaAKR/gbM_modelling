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

filename_start = "MAF5_H2AZWT_EqbrmOutput_FinalFit_Pchoice_Excl_InitialState_Expt_Ngen_1E5_"
filename_ending = '.tsv'

# LocusProperties used for Col0 data - only need to load in once! 
file_in_LocusProperties = os.path.join('Output_files', 'MAF5_H2AZWT_EqbrmOutput_FinalFit_Pchoice_Excl_InitialState_Expt_Ngen_1E5_LocusProperties.tsv')

LocusProperties_df = pd.read_csv(file_in_LocusProperties, sep="\t{1}",engine='python')
LocusProperties_df.set_index("gene_ID", inplace = True)
N_gene = len(LocusProperties_df)
IDs_list = LocusProperties_df.index.values.tolist()
print('N_gene_Expt',N_gene)

# Manually load in the Sims from all three initial states. 
file_in_AllRepsMethFracs_Expt = os.path.join('Output_files', 
'MAF5_H2AZWT_EqbrmOutput_FinalFit_Pchoice_Excl_InitialState_Expt_Ngen_1E5_AllRepsMethFracs.tsv')
AllRepsMethFracs_Expt_df = pd.read_csv(file_in_AllRepsMethFracs_Expt, sep="\t{1}",engine='python')
AllRepsMethFracs_Expt_df.set_index("gene_ID", inplace = True)
print('N_gene_Expt',len(AllRepsMethFracs_Expt_df))

file_in_AllRepsMethFracs_100M = os.path.join('Output_files', 
'MAF5_H2AZWT_EqbrmOutput_FinalFit_Pchoice_Excl_InitialState_100M_Ngen_1E5_AllRepsMethFracs.tsv')
AllRepsMethFracs_100M_df = pd.read_csv(file_in_AllRepsMethFracs_100M, sep="\t{1}",engine='python')
AllRepsMethFracs_100M_df.set_index("gene_ID", inplace = True)
print('N_gene_100M',len(AllRepsMethFracs_100M_df))

file_in_AllRepsMethFracs_100U = os.path.join('Output_files', 
'MAF5_H2AZWT_EqbrmOutput_FinalFit_Pchoice_Excl_InitialState_100U_Ngen_1E5_AllRepsMethFracs.tsv')
AllRepsMethFracs_100U_df = pd.read_csv(file_in_AllRepsMethFracs_100U, sep="\t{1}",engine='python')
AllRepsMethFracs_100U_df.set_index("gene_ID", inplace = True)
print('N_gene_100U',len(AllRepsMethFracs_100U_df))

# Load in RatesOnlyFit
# MAF5_H2AZWT_EqbrmOutput_RatesOnlyFit_Pchoice_Excl_InitialState_Expt_Ngen_1E5_AllRepsMethFracs.tsv
file_in_AllRepsMethFracs_Expt_RatesOnlyFit = os.path.join('Output_files', 
'MAF5_H2AZWT_EqbrmOutput_RatesOnlyFit_Pchoice_Excl_InitialState_Expt_Ngen_1E5_AllRepsMethFracs.tsv')
AllRepsMethFracs_Expt_RatesOnlyFit_df = pd.read_csv(file_in_AllRepsMethFracs_Expt_RatesOnlyFit, sep="\t{1}",engine='python')
AllRepsMethFracs_Expt_RatesOnlyFit_df.set_index("gene_ID", inplace = True)
print('N_gene_Expt',len(AllRepsMethFracs_Expt_RatesOnlyFit_df))


# define colours
col_data_Col0 = 'green'
col_Sim_Expt = 'k'
col_Sim_100M = 'r'
col_Sim_100U = 'b'

fig_label_Expt = 'F4A'
fig_label_Converged = 'F4B'

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

# Make lists of loci underand over methylated by model. 
Meth_threshold = 0.2

Loci_OverMeth_list = LocusProperties_df.loc[ LocusProperties_df['D3D4D5_mu_meth_diff'] >= Meth_threshold ].index.values
Loci_CorrectMeth_list = LocusProperties_df.loc[ (LocusProperties_df['D3D4D5_mu_meth_diff'] < Meth_threshold) & 
                                                (LocusProperties_df['D3D4D5_mu_meth_diff'] >= -Meth_threshold) ].index.values
Loci_UnderMeth_list = LocusProperties_df.loc[ LocusProperties_df['D3D4D5_mu_meth_diff'] < -Meth_threshold ].index.values
print('Over','Correct','Under')
print(len(Loci_OverMeth_list),len(Loci_CorrectMeth_list),len(Loci_UnderMeth_list))
print(len(Loci_OverMeth_list)+len(Loci_CorrectMeth_list)+len(Loci_UnderMeth_list))

file_name_output = os.path.join(GraphFolder,filename_start +'OverMethLoci'+'.tsv')
open(file_name_output,'w').close
file_name_output = os.path.join(GraphFolder,filename_start +'CorrectMethLoci'+'.tsv')
open(file_name_output,'w').close
file_name_output = os.path.join(GraphFolder,filename_start +'UnderMethLoci'+'.tsv')
open(file_name_output,'w').close



for i_ in range(len(Loci_OverMeth_list)):
    # Write out to file
    file_name_output = os.path.join(GraphFolder,filename_start +'OverMethLoci'+'.tsv')
    with open(file_name_output,'a',newline='') as output_file:
        line_writer = csv.writer(output_file,delimiter='\t')
        line_writer.writerow( [Loci_OverMeth_list[i_]] )

for i_ in range(len(Loci_CorrectMeth_list)):
    # Write out to file
    file_name_output = os.path.join(GraphFolder,filename_start +'CorrectMethLoci'+'.tsv')
    with open(file_name_output,'a',newline='') as output_file:
        line_writer = csv.writer(output_file,delimiter='\t')
        line_writer.writerow( [Loci_CorrectMeth_list[i_]] )

for i_ in range(len(Loci_UnderMeth_list)):
    # Write out to file
    file_name_output = os.path.join(GraphFolder,filename_start +'UnderMethLoci'+'.tsv')
    with open(file_name_output,'a',newline='') as output_file:
        line_writer = csv.writer(output_file,delimiter='\t')
        line_writer.writerow( [Loci_UnderMeth_list[i_]] )



# Start Graph
fig, ax = plt.subplots(1,1,figsize=(3,5))

for spine in ['left','right','top','bottom']:
    ax.spines[spine].set_color('k')
    ax.spines[spine].set_linewidth(0.8)
ax.set_facecolor('white')

#ax.grid(False)
ax.grid(b=True, which='major', color='lightgrey', linestyle=':',linewidth=1)


# convert to percentage
values_scaler = 100
x_start = 0.
x_end = 100. 

x_hist = np.arange(x_start+x_end/(2.*n_bin),x_end,x_end/n_bin)


# # hist_Col0_label = 'Col-0 Data'

# # hist_Col0, bins_Col0 = np.histogram(LocusProperties_df[['Col0_meth_frac']].values*values_scaler, bins=n_bin, range=(x_start,x_end))
# # # ax.axvline(np.nanmean(LocusProperties_df[['Col0_meth_frac']].values*values_scaler), linestyle='-', color=col_data_Col0)
# # ax.plot(x_hist, hist_Col0,  linewidth = 2, color=col_data_Col0,
# #                   label='%s\nMean = %.1f %%' % (hist_Col0_label,np.nanmean(LocusProperties_df[['Col0_meth_frac']].values*values_scaler)) )


# decided we should show 30 Col-like here, rather than just Col0

hist_Data_label = '30 Col-like Data'

temp_cols_list = []
for i_temp in range(10):
    temp_cols_list.append('D3_meth_frac_'+str(i_temp))
    temp_cols_list.append('D4_meth_frac_'+str(i_temp))
    temp_cols_list.append('D5_meth_frac_'+str(i_temp))
hist_Data, bins_Data = np.histogram(AllRepsMethFracs_Expt_df[temp_cols_list].values*values_scaler, bins=n_bin, range=(x_start,x_end))
hist_Data = hist_Data/N_reps # normalise
# ax.axvline(np.nanmean(AllRepsMethFracs_Expt_df[temp_cols_list].values*values_scaler), linestyle='--', color='g')
ax.plot(x_hist, hist_Data,  linewidth = 2, color=col_data_Col0, linestyle='-',
                  label='%s\nMean = %.1f %%' % (hist_Data_label,np.nanmean(AllRepsMethFracs_Expt_df[temp_cols_list].values*values_scaler)) )



##

hist_Sims_Expt_label = 'Simulated'

temp_cols_list = []
for i_temp in range(N_reps):
    temp_cols_list.append('Sim_meth_frac_'+str(i_temp))
hist_Sims_Expt, bins_Sim = np.histogram(AllRepsMethFracs_Expt_df[temp_cols_list].values*values_scaler, bins=n_bin, range=(x_start,x_end))
hist_Sims_Expt = hist_Sims_Expt/N_reps # normalise
# ax.axvline(np.nanmean(AllRepsMethFracs_Expt_df[temp_cols_list].values*values_scaler), linestyle='-', color=col_Sim_Expt)
ax.plot(x_hist, hist_Sims_Expt,  linewidth = 2, color=col_Sim_Expt,linestyle='-',
                  label='%s\nMean = %.1f %%' % (hist_Sims_Expt_label,np.nanmean(AllRepsMethFracs_Expt_df[temp_cols_list].values*values_scaler)) )



# fig.suptitle(title_sting_single_time, fontsize=12)


ax.legend(fontsize=10,ncol=1,facecolor='white', framealpha=1.0)


ax.set_xlim(-10,110)



ax.set_ylim(0,1000)


ax.xaxis.set_tick_params(labelsize=12)
ax.yaxis.set_tick_params(labelsize=12)


ax.set_xlabel("Methylation percentage", fontsize=12)
ax.set_ylabel("Number of loci", fontsize=12)


#fig.tight_layout()
fig.subplots_adjust(top=0.5)

fig.savefig( os.path.join(GraphFolder,filename_start +"Mfrac_InitialState_Expt_"+fig_label_Expt+".png"), bbox_inches = 'tight')
fig.savefig( os.path.join(GraphFolder,filename_start +"Mfrac_InitialState_Expt_"+fig_label_Expt+".pdf"), bbox_inches = 'tight')


plt.show()
# End Graph


################################################


# Start Graph
fig, ax = plt.subplots(1,1,figsize=(3,5))

for spine in ['left','right','top','bottom']:
    ax.spines[spine].set_color('k')
    ax.spines[spine].set_linewidth(0.8)
ax.set_facecolor('white')

#ax.grid(False)
ax.grid(b=True, which='major', color='lightgrey', linestyle=':',linewidth=1)


# convert to percentage
values_scaler = 100
x_start = 0.
x_end = 100. 

x_hist = np.arange(x_start+x_end/(2.*n_bin),x_end,x_end/n_bin)

# hist_Col0, bins_Col0 = np.histogram(LocusProperties_df[['Col0_meth_frac']].values*values_scaler, bins=n_bin, range=(x_start,x_end))
# ax.axvline(np.nanmean(LocusProperties_df[['Col0_meth_frac']].values*values_scaler), linestyle='-', color='g')
# ax.plot(x_hist, hist_Col0,  linewidth = 2, color='g',
#                   label='Data Col0\n$\mu$ = %.4f' % np.nanmean(LocusProperties_df[['Col0_meth_frac']].values*values_scaler))

hist_Sims_Expt_label = 'Col-0 Simulated'

temp_cols_list = []
for i_temp in range(N_reps):
    temp_cols_list.append('Sim_meth_frac_'+str(i_temp))
hist_Sims_Expt, bins_Sim = np.histogram(AllRepsMethFracs_Expt_df[temp_cols_list].values*values_scaler, bins=n_bin, range=(x_start,x_end))
hist_Sims_Expt = hist_Sims_Expt/N_reps # normalise
# ax.axvline(np.nanmean(AllRepsMethFracs_Expt_df[temp_cols_list].values*values_scaler), linestyle='-', color=col_Sim_Expt)
ax.plot(x_hist, hist_Sims_Expt,  linewidth = 2, color=col_Sim_Expt,linestyle='-',
                  label='%s\nMean = %.1f %%' % (hist_Sims_Expt_label,np.nanmean(AllRepsMethFracs_Expt_df[temp_cols_list].values*values_scaler)) )

hist_Sims_100M_label = 'All $M$ Simulated'

temp_cols_list = []
for i_temp in range(N_reps):
    temp_cols_list.append('Sim_meth_frac_'+str(i_temp))
hist_Sims_100M, bins_Sim = np.histogram(AllRepsMethFracs_100M_df[temp_cols_list].values*values_scaler, bins=n_bin, range=(x_start,x_end))
hist_Sims_100M = hist_Sims_100M/N_reps # normalise
# ax.axvline(np.nanmean(AllRepsMethFracs_100M_df[temp_cols_list].values*values_scaler), linestyle='--', color=col_Sim_100M)
ax.plot(x_hist, hist_Sims_100M,  linewidth = 2, color=col_Sim_100M,linestyle='--',
                  label='%s\nMean = %.1f %%' % (hist_Sims_100M_label,np.nanmean(AllRepsMethFracs_100M_df[temp_cols_list].values*values_scaler)) )

hist_Sims_100U_label = 'All $U$ Simulated'

temp_cols_list = []
for i_temp in range(N_reps):
    temp_cols_list.append('Sim_meth_frac_'+str(i_temp))
hist_Sims_100U, bins_Sim = np.histogram(AllRepsMethFracs_100U_df[temp_cols_list].values*values_scaler, bins=n_bin, range=(x_start,x_end))
hist_Sims_100U = hist_Sims_100U/N_reps # normalise
# ax.axvline(np.nanmean(AllRepsMethFracs_100U_df[temp_cols_list].values*values_scaler), linestyle=':', color=col_Sim_100U)
ax.plot(x_hist, hist_Sims_100U,  linewidth = 2, color=col_Sim_100U,linestyle=':',
                  label='%s\nMean = %.1f %%' % (hist_Sims_100U_label,np.nanmean(AllRepsMethFracs_100U_df[temp_cols_list].values*values_scaler)) )

# fig.suptitle(title_sting_single_time, fontsize=12)


ax.legend(fontsize=9,ncol=1,facecolor='white', framealpha=1.0,handlelength=2.0)


ax.set_xlim(-10,110)



ax.set_ylim(0,1000)


ax.xaxis.set_tick_params(labelsize=12)
ax.yaxis.set_tick_params(labelsize=12)


ax.set_xlabel("Methylation percentage", fontsize=12)
ax.set_ylabel("Number of loci", fontsize=12)


#fig.tight_layout()
fig.subplots_adjust(top=0.5)

fig.savefig( os.path.join(GraphFolder,filename_start +"Mfrac_InitialState_Convgd_"+fig_label_Converged+".png"), bbox_inches = 'tight')
fig.savefig( os.path.join(GraphFolder,filename_start +"Mfrac_InitialState_Convgd_"+fig_label_Converged+".pdf"), bbox_inches = 'tight')


plt.show()
# End Graph


###########################################

# Start Graph
fig, ax = plt.subplots(1,1,figsize=(3,3))

for spine in ['left','right','top','bottom']:
    ax.spines[spine].set_color('k')
    ax.spines[spine].set_linewidth(0.8)
ax.set_facecolor('white')

#ax.grid(False)
ax.grid(b=True, which='major', color='lightgrey', linestyle=':',linewidth=1)


n_bin_2 = 100

# convert to percentage
values_scaler = 100
x_start = 0.
x_end = 100. 

x_hist = np.arange(x_start+x_end/(2.*n_bin_2),x_end,x_end/n_bin_2)

hist_D3D4D5_label = 'Col-0 Data'

hist_D3D4D5, bins_D3D4D5 = np.histogram(LocusProperties_df[['D3D4D5_mu_meth_frac']].values*values_scaler, bins=n_bin_2, range=(x_start,x_end))
# ax.axvline(np.nanmean(LocusProperties_df[['D3D4D5_mu_meth_frac']].values*values_scaler), linestyle='-', color=col_data_Col0)
ax.plot(x_hist, hist_D3D4D5,  linewidth = 2, color=col_data_Col0,
                  label='%s\nMean = %.1f %%' % (hist_D3D4D5_label,np.nanmean(LocusProperties_df[['D3D4D5_mu_meth_frac']].values*values_scaler)) )

hist_Sims_Expt_label = 'Simulated'

hist_Sims_Expt, bins_Sim = np.histogram(LocusProperties_df[['Sim_mu_meth_frac']].values*values_scaler, bins=n_bin_2, range=(x_start,x_end))
# ax.axvline(np.nanmean(LocusProperties_df[['Sim_mu_meth_frac']].values*values_scaler), linestyle='-', color=col_Sim_Expt)
ax.plot(x_hist, hist_Sims_Expt,  linewidth = 2, color=col_Sim_Expt,linestyle='-',
                  label='%s\nMean = %.1f %%' % (hist_Sims_Expt_label,np.nanmean(LocusProperties_df[['Sim_mu_meth_frac']].values*values_scaler)) )



# fig.suptitle(title_sting_single_time, fontsize=12)


# ax.legend(fontsize=10,ncol=1,facecolor='white', framealpha=1.0,handlelength=2.0)


ax.set_xlim(-10,110)


#ax.set_ylim(0,200)


ax.xaxis.set_tick_params(labelsize=12)
ax.yaxis.set_tick_params(labelsize=12)


ax.set_xlabel("Methylation percentage", fontsize=12)
ax.set_ylabel("Number of loci", fontsize=12)


#fig.tight_layout()
# fig.subplots_adjust(top=0.5)

fig.savefig( os.path.join(GraphFolder,filename_start +"Mfrac_InitialState_Expt_MeanOverReps"+".png"), bbox_inches = 'tight')

plt.show()
# End Graph


################################################


#############

N_plots_axi = 1
N_plots_axj = 4

n_bin = 50


CG_density_maxs = [55.,47.,40.,10.]
CG_density_mins = [1000.,55.,47.,40.]

# convert to percentage
values_scaler = 100


fig, ax = plt.subplots(N_plots_axi,N_plots_axj,figsize=(N_plots_axj*3,N_plots_axi*3.8))

x_start = 0.
x_end = 100. 
x_hist = np.arange(x_start+x_end/(2.*n_bin),x_end,x_end/n_bin)

for axi in range(N_plots_axi):
    for axj in range(N_plots_axj):

        for spine in ['left','right','top','bottom']:
            ax[axj].spines[spine].set_color('k')
            ax[axj].spines[spine].set_linewidth(0.8)
        ax[axj].set_facecolor('white')
        #ax.grid(False)
        ax[axj].grid(b=True, which='major', color='lightgrey', linestyle=':',linewidth=1)

        #filter dataframes
        # Simulation including CG density correction
        Temp_Sim_df1 = AllRepsMethFracs_Expt_df.loc[ ((LocusProperties_df['CG_density'] >= (1./CG_density_mins[axj])) &
                                                    (LocusProperties_df['CG_density'] < (1./CG_density_maxs[axj]))) ]
        Temp_Sim1_label = 'Sim-2\nMean = %.1f %%' % (np.nanmean(Temp_Sim_df1[temp_cols_list].values*values_scaler))

        # Simulation with no CG density correction
        Temp_Sim_df0 = AllRepsMethFracs_Expt_RatesOnlyFit_df.loc[ ((LocusProperties_df['CG_density'] >= (1./CG_density_mins[axj])) &
                                                    (LocusProperties_df['CG_density'] < (1./CG_density_maxs[axj]))) ]
        Temp_Sim0_label = 'Sim-1\nMean = %.1f %%' % (np.nanmean(Temp_Sim_df0[temp_cols_list].values*values_scaler))

        Temp_Data_df = LocusProperties_df.loc[ ((LocusProperties_df['CG_density'] >= 1./CG_density_mins[axj]) &
                                                    (LocusProperties_df['CG_density'] < 1./CG_density_maxs[axj])) ]
        Temp_Data_label = 'Col-0\nMean = %.1f %%' % (np.nanmean(Temp_Data_df['Col0_meth_frac'].values*values_scaler))
        Title_label = "$N_{loci}$ = %d" % ( len(Temp_Data_df) )
        #CG_density_label = "$\\frac{{1}}{{%d}} \leq \\rho_{CG} < \\frac{{1}}{{%d}}$" % (CG_density_mins[axj], CG_density_maxs[axj])
        CG_density_label = "%.3f $\leq \\rho_{CG} <$ %.3f" % (1./CG_density_mins[axj], 1./CG_density_maxs[axj])

        # dummy data
        #ax[axj].plot(x_hist,x_hist+1E5,linewidth=0,label=Title_label )

        hist_Col0, bins_Col0 = np.histogram(Temp_Data_df[['Col0_meth_frac']].values*values_scaler, bins=n_bin, range=(x_start,x_end))
        # ax[axj].axvline(np.nanmean(Temp_Data_df[['Col0_meth_frac']].values*values_scaler), linestyle=':', color='g')
        ax[axj].plot(x_hist, hist_Col0,  linewidth = 2, color='g',label=Temp_Data_label)


        temp_cols_list = []
        for i_temp in range(N_reps):
            temp_cols_list.append('Sim_meth_frac_'+str(i_temp))
        hist_Sim0, bins_Sim = np.histogram(Temp_Sim_df0[temp_cols_list].values*values_scaler, bins=n_bin, range=(x_start,x_end))
        hist_Sim0 = hist_Sim0/N_reps # normalise
        # ax[axj].axvline(np.nanmean(Temp_Sim_df0[temp_cols_list].values*values_scaler), linestyle=':', color='k')
        ax[axj].plot(x_hist, hist_Sim0,  linewidth = 2, color='k',label=Temp_Sim0_label)


        temp_cols_list = []
        for i_temp in range(N_reps):
            temp_cols_list.append('Sim_meth_frac_'+str(i_temp))
        hist_Sim1, bins_Sim = np.histogram(Temp_Sim_df1[temp_cols_list].values*values_scaler, bins=n_bin, range=(x_start,x_end))
        hist_Sim1 = hist_Sim1/N_reps # normalise
        # ax[axj].axvline(np.nanmean(Temp_Sim_df1[temp_cols_list].values*values_scaler), linestyle=':', color='m')
        ax[axj].plot(x_hist, hist_Sim1,  linewidth = 2, color='m',label=Temp_Sim1_label, linestyle='--')


        # temp_cols_list = []
        # for i_temp in range(N_reps):
        #     temp_cols_list.append('Dset1_meth_frac_'+str(i_temp))
        # hist_Data, bins_Data = np.histogram(Temp_Data_df[temp_cols_list].values*values_scaler, bins=n_bin, range=(x_start,x_end))fontsize=12
        # hist_Data = hist_Data/N_reps # normalise
        # ax[axj].axvline(np.nanmean(Temp_Data_df[temp_cols_list].values*values_scaler), linestyle=':', color='green')
        # ax[axj].plot(x_hist, hist_Data,  linewidth = 2, color='green', linestyle='--')



        #fig.suptitle(title_sting_single_time, fontsize=12)
        ax[axj].legend(fontsize=10,ncol=1,facecolor='white', framealpha=1.0,handlelength=1.5)

        ax[axj].set_xlim(-10,110)
        ax[axj].set_ylim(0,320)

        ax[axj].xaxis.set_tick_params(labelsize=12)
        ax[axj].yaxis.set_tick_params(labelsize=12)

        ax[axj].set_title(Title_label, fontsize=12)


        ax[axj].set_title(CG_density_label+'\n\n\n'+Title_label, fontsize=12)

        ax[axj].set_xlabel("Methylation percentage", fontsize=12)

        ax[axj].set_ylabel("Number of loci", fontsize=12)

ax[0].arrow(0, 395, 580, 0, fc='k', ec='k', clip_on=False, width=0.8, head_width=8, head_length=12)
ax[-1].annotate("Increasing CG-site density", (-80, 370), fontsize=12, annotation_clip=False)

fig.tight_layout()
fig.subplots_adjust(hspace = 0.4)

fig.savefig( os.path.join(GraphFolder,filename_start +'Mfrac_dist_groups_Rho'+str(N_plots_axj)+'_SF6G.png'), bbox_inches = 'tight')
fig.savefig( os.path.join(GraphFolder,filename_start +'Mfrac_dist_groups_Rho'+str(N_plots_axj)+'_SF6G.pdf'), bbox_inches = 'tight')

plt.show()
#######################################################################

#############

N_plots_axi = 1
N_plots_axj = 4

n_bin = 50


CG_density_maxs = [55.,47.,40.,10.]
CG_density_mins = [1000.,55.,47.,40.]

# convert to percentage
values_scaler = 100


fig, ax = plt.subplots(N_plots_axi,N_plots_axj,figsize=(N_plots_axj*3,N_plots_axi*3.8))

x_start = 0.
x_end = 100. 
x_hist = np.arange(x_start+x_end/(2.*n_bin),x_end,x_end/n_bin)

for axi in range(N_plots_axi):
    for axj in range(N_plots_axj):

        for spine in ['left','right','top','bottom']:
            ax[axj].spines[spine].set_color('k')
            ax[axj].spines[spine].set_linewidth(0.8)
        ax[axj].set_facecolor('white')
        #ax.grid(False)
        ax[axj].grid(b=True, which='major', color='lightgrey', linestyle=':',linewidth=1)

        #filter dataframes
        # Simulation including CG density correction
        Temp_Sim_df1 = AllRepsMethFracs_Expt_df.loc[ ((LocusProperties_df['CG_density'] >= (1./CG_density_mins[axj])) &
                                                    (LocusProperties_df['CG_density'] < (1./CG_density_maxs[axj]))) ]
        Temp_Sim1_label = 'Sim-2\nMean = %.1f %%' % (np.nanmean(Temp_Sim_df1[temp_cols_list].values*values_scaler))

        # Simulation with no CG density correction
        Temp_Sim_df0 = AllRepsMethFracs_Expt_RatesOnlyFit_df.loc[ ((LocusProperties_df['CG_density'] >= (1./CG_density_mins[axj])) &
                                                    (LocusProperties_df['CG_density'] < (1./CG_density_maxs[axj]))) ]
        Temp_Sim0_label = 'Sim-1\nMean = %.1f %%' % (np.nanmean(Temp_Sim_df0[temp_cols_list].values*values_scaler))

        # 30 Col-like accns 
        Temp_Data_df = AllRepsMethFracs_Expt_df.loc[ ((LocusProperties_df['CG_density'] >= 1./CG_density_mins[axj]) &
                                                    (LocusProperties_df['CG_density'] < 1./CG_density_maxs[axj])) ]                                           

        Title_label = "$N_{loci}$ = %d" % ( len(Temp_Data_df) )
        #CG_density_label = "$\\frac{{1}}{{%d}} \leq \\rho_{CG} < \\frac{{1}}{{%d}}$" % (CG_density_mins[axj], CG_density_maxs[axj])
        CG_density_label = "%.3f $\leq \\rho_{CG} <$ %.3f" % (1./CG_density_mins[axj], 1./CG_density_maxs[axj])

        # dummy data
        #ax[axj].plot(x_hist,x_hist+1E5,linewidth=0,label=Title_label )

        temp_cols_list = []
        for i_temp in range(10):
            temp_cols_list.append('D3_meth_frac_'+str(i_temp))
            temp_cols_list.append('D4_meth_frac_'+str(i_temp))
            temp_cols_list.append('D5_meth_frac_'+str(i_temp))
        Temp_Data_label = '30 Col-like'
        hist_Data, bins_Data = np.histogram(Temp_Data_df[temp_cols_list].values*values_scaler, bins=n_bin, range=(x_start,x_end))
        hist_Data = hist_Data/len(temp_cols_list) # normalise
        # ax[axj].axvline(np.nanmean(Temp_Data_df[temp_cols_list].values*values_scaler), linestyle='--', color='g')
        ax[axj].plot(x_hist, hist_Data,  linewidth = 2, color=col_data_Col0, linestyle='-',
                        label='%s\nMean = %.1f %%' % (Temp_Data_label,np.nanmean(Temp_Data_df[temp_cols_list].values*values_scaler)) )

        # hist_Col0, bins_Col0 = np.histogram(Temp_Data_df[['Col0_meth_frac']].values*values_scaler, bins=n_bin, range=(x_start,x_end))
        # # ax[axj].axvline(np.nanmean(Temp_Data_df[['Col0_meth_frac']].values*values_scaler), linestyle=':', color='g')
        # ax[axj].plot(x_hist, hist_Col0,  linewidth = 2, color='g',label=Temp_Data_label)


        temp_cols_list = []
        for i_temp in range(N_reps):
            temp_cols_list.append('Sim_meth_frac_'+str(i_temp))
        hist_Sim0, bins_Sim = np.histogram(Temp_Sim_df0[temp_cols_list].values*values_scaler, bins=n_bin, range=(x_start,x_end))
        hist_Sim0 = hist_Sim0/N_reps # normalise
        # ax[axj].axvline(np.nanmean(Temp_Sim_df0[temp_cols_list].values*values_scaler), linestyle=':', color='k')
        ax[axj].plot(x_hist, hist_Sim0,  linewidth = 2, color='k',label=Temp_Sim0_label)


        temp_cols_list = []
        for i_temp in range(N_reps):
            temp_cols_list.append('Sim_meth_frac_'+str(i_temp))
        hist_Sim1, bins_Sim = np.histogram(Temp_Sim_df1[temp_cols_list].values*values_scaler, bins=n_bin, range=(x_start,x_end))
        hist_Sim1 = hist_Sim1/N_reps # normalise
        # ax[axj].axvline(np.nanmean(Temp_Sim_df1[temp_cols_list].values*values_scaler), linestyle=':', color='m')
        ax[axj].plot(x_hist, hist_Sim1,  linewidth = 2, color='m',label=Temp_Sim1_label, linestyle='--')


        # temp_cols_list = []
        # for i_temp in range(N_reps):
        #     temp_cols_list.append('Dset1_meth_frac_'+str(i_temp))
        # hist_Data, bins_Data = np.histogram(Temp_Data_df[temp_cols_list].values*values_scaler, bins=n_bin, range=(x_start,x_end))fontsize=12
        # hist_Data = hist_Data/N_reps # normalise
        # ax[axj].axvline(np.nanmean(Temp_Data_df[temp_cols_list].values*values_scaler), linestyle=':', color='green')
        # ax[axj].plot(x_hist, hist_Data,  linewidth = 2, color='green', linestyle='--')



        #fig.suptitle(title_sting_single_time, fontsize=12)
        ax[axj].legend(fontsize=10,ncol=1,facecolor='white', framealpha=1.0,handlelength=1.5)

        ax[axj].set_xlim(-10,110)
        ax[axj].set_ylim(0,320)

        ax[axj].xaxis.set_tick_params(labelsize=12)
        ax[axj].yaxis.set_tick_params(labelsize=12)

        ax[axj].set_title(Title_label, fontsize=12)


        ax[axj].set_title(CG_density_label+'\n\n\n'+Title_label, fontsize=12)

        ax[axj].set_xlabel("Methylation percentage", fontsize=12)

        ax[axj].set_ylabel("Number of loci", fontsize=12)

ax[0].arrow(0, 395, 580, 0, fc='k', ec='k', clip_on=False, width=0.8, head_width=8, head_length=12)
ax[-1].annotate("Increasing CG-site density", (-80, 370), fontsize=12, annotation_clip=False)

fig.tight_layout()
fig.subplots_adjust(hspace = 0.4)

fig.savefig( os.path.join(GraphFolder,filename_start +'Mfrac_dist_groups_Rho'+str(N_plots_axj)+'_SF6G_NEW.png'), bbox_inches = 'tight')
fig.savefig( os.path.join(GraphFolder,filename_start +'Mfrac_dist_groups_Rho'+str(N_plots_axj)+'_SF6G_NEW.pdf'), bbox_inches = 'tight')

plt.show()


###########################################

# Simulation with no CG density correction
# # Load in RatesOnlyFit
# MAF5_H2AZWT_EqbrmOutput_RatesOnlyFit_Pchoice_Excl_InitialState_Expt_Ngen_1E5_AllRepsMethFracs.tsv


# Start Graph
fig, ax = plt.subplots(1,1,figsize=(3,5))

for spine in ['left','right','top','bottom']:
    ax.spines[spine].set_color('k')
    ax.spines[spine].set_linewidth(0.8)
ax.set_facecolor('white')

#ax.grid(False)
ax.grid(b=True, which='major', color='lightgrey', linestyle=':',linewidth=1)


# convert to percentage
values_scaler = 100
x_start = 0.
x_end = 100. 

x_hist = np.arange(x_start+x_end/(2.*n_bin),x_end,x_end/n_bin)


hist_Col0_label = 'Col-0 Data'

hist_Col0, bins_Col0 = np.histogram(LocusProperties_df[['Col0_meth_frac']].values*values_scaler, bins=n_bin, range=(x_start,x_end))
# ax.axvline(np.nanmean(LocusProperties_df[['Col0_meth_frac']].values*values_scaler), linestyle='-', color=col_data_Col0)
ax.plot(x_hist, hist_Col0,  linewidth = 2, color=col_data_Col0,
                  label='%s\nMean = %.1f %%' % (hist_Col0_label,np.nanmean(LocusProperties_df[['Col0_meth_frac']].values*values_scaler)) )

hist_Sims_Expt_label = 'Simulated'

temp_cols_list = []
for i_temp in range(N_reps):
    temp_cols_list.append('Sim_meth_frac_'+str(i_temp))
hist_Sims_Expt, bins_Sim = np.histogram(AllRepsMethFracs_Expt_RatesOnlyFit_df[temp_cols_list].values*values_scaler, bins=n_bin, range=(x_start,x_end))
hist_Sims_Expt = hist_Sims_Expt/N_reps # normalise
# ax.axvline(np.nanmean(AllRepsMethFracs_Expt_RatesOnlyFit_df[temp_cols_list].values*values_scaler), linestyle='-', color=col_Sim_Expt)
ax.plot(x_hist, hist_Sims_Expt,  linewidth = 2, color=col_Sim_Expt,linestyle='-',
                  label='%s\nMean = %.1f %%' % (hist_Sims_Expt_label,np.nanmean(AllRepsMethFracs_Expt_RatesOnlyFit_df[temp_cols_list].values*values_scaler)) )



# fig.suptitle(title_sting_single_time, fontsize=12)


ax.legend(fontsize=10,ncol=1,facecolor='white', framealpha=1.0)


ax.set_xlim(-10,110)



ax.set_ylim(0,1000)


ax.xaxis.set_tick_params(labelsize=12)
ax.yaxis.set_tick_params(labelsize=12)


ax.set_xlabel("Methylation percentage", fontsize=12)
ax.set_ylabel("Number of loci", fontsize=12)


#fig.tight_layout()
fig.subplots_adjust(top=0.5)

fig.savefig( os.path.join(GraphFolder,"MAF5_H2AZWT_EqbrmOutput_RatesOnlyFit_Pchoice_Excl_InitialState_Expt_Ngen_1E5_Mfrac_InitialState_Expt_SF6F"+".png"), bbox_inches = 'tight')
fig.savefig( os.path.join(GraphFolder,"MAF5_H2AZWT_EqbrmOutput_RatesOnlyFit_Pchoice_Excl_InitialState_Expt_Ngen_1E5_Mfrac_InitialState_Expt_SF6F"+".pdf"), bbox_inches = 'tight')


plt.show()
# End Graph


################################################
# Simulation with no CG density correction
# # Load in RatesOnlyFit
# MAF5_H2AZWT_EqbrmOutput_RatesOnlyFit_Pchoice_Excl_InitialState_Expt_Ngen_1E5_AllRepsMethFracs.tsv


# Start Graph
fig, ax = plt.subplots(1,1,figsize=(3,5))

for spine in ['left','right','top','bottom']:
    ax.spines[spine].set_color('k')
    ax.spines[spine].set_linewidth(0.8)
ax.set_facecolor('white')

#ax.grid(False)
ax.grid(b=True, which='major', color='lightgrey', linestyle=':',linewidth=1)


# convert to percentage
values_scaler = 100
x_start = 0.
x_end = 100. 

x_hist = np.arange(x_start+x_end/(2.*n_bin),x_end,x_end/n_bin)




hist_Data_label = '30 Col-like Data'

temp_cols_list = []
for i_temp in range(10):
    temp_cols_list.append('D3_meth_frac_'+str(i_temp))
    temp_cols_list.append('D4_meth_frac_'+str(i_temp))
    temp_cols_list.append('D5_meth_frac_'+str(i_temp))
hist_Data, bins_Data = np.histogram(AllRepsMethFracs_Expt_df[temp_cols_list].values*values_scaler, bins=n_bin, range=(x_start,x_end))
hist_Data = hist_Data/N_reps # normalise
# ax.axvline(np.nanmean(AllRepsMethFracs_Expt_df[temp_cols_list].values*values_scaler), linestyle='--', color='g')
ax.plot(x_hist, hist_Data,  linewidth = 2, color=col_data_Col0, linestyle='-',
                  label='%s\nMean = %.1f %%' % (hist_Data_label,np.nanmean(AllRepsMethFracs_Expt_df[temp_cols_list].values*values_scaler)) )



# hist_Col0_label = 'Col-0 Data'

# hist_Col0, bins_Col0 = np.histogram(LocusProperties_df[['Col0_meth_frac']].values*values_scaler, bins=n_bin, range=(x_start,x_end))
# # ax.axvline(np.nanmean(LocusProperties_df[['Col0_meth_frac']].values*values_scaler), linestyle='-', color=col_data_Col0)
# ax.plot(x_hist, hist_Col0,  linewidth = 2, color=col_data_Col0,
#                   label='%s\nMean = %.1f %%' % (hist_Col0_label,np.nanmean(LocusProperties_df[['Col0_meth_frac']].values*values_scaler)) )

hist_Sims_Expt_label = 'Simulated'

temp_cols_list = []
for i_temp in range(N_reps):
    temp_cols_list.append('Sim_meth_frac_'+str(i_temp))
hist_Sims_Expt, bins_Sim = np.histogram(AllRepsMethFracs_Expt_RatesOnlyFit_df[temp_cols_list].values*values_scaler, bins=n_bin, range=(x_start,x_end))
hist_Sims_Expt = hist_Sims_Expt/N_reps # normalise
# ax.axvline(np.nanmean(AllRepsMethFracs_Expt_RatesOnlyFit_df[temp_cols_list].values*values_scaler), linestyle='-', color=col_Sim_Expt)
ax.plot(x_hist, hist_Sims_Expt,  linewidth = 2, color=col_Sim_Expt,linestyle='-',
                  label='%s\nMean = %.1f %%' % (hist_Sims_Expt_label,np.nanmean(AllRepsMethFracs_Expt_RatesOnlyFit_df[temp_cols_list].values*values_scaler)) )



# fig.suptitle(title_sting_single_time, fontsize=12)


ax.legend(fontsize=10,ncol=1,facecolor='white', framealpha=1.0)


ax.set_xlim(-10,110)



ax.set_ylim(0,1000)


ax.xaxis.set_tick_params(labelsize=12)
ax.yaxis.set_tick_params(labelsize=12)


ax.set_xlabel("Methylation percentage", fontsize=12)
ax.set_ylabel("Number of loci", fontsize=12)


#fig.tight_layout()
fig.subplots_adjust(top=0.5)

fig.savefig( os.path.join(GraphFolder,"MAF5_H2AZWT_EqbrmOutput_RatesOnlyFit_Pchoice_Excl_InitialState_Expt_Ngen_1E5_Mfrac_InitialState_Expt_SF6F_NEW"+".png"), bbox_inches = 'tight')
fig.savefig( os.path.join(GraphFolder,"MAF5_H2AZWT_EqbrmOutput_RatesOnlyFit_Pchoice_Excl_InitialState_Expt_Ngen_1E5_Mfrac_InitialState_Expt_SF6F_NEW"+".pdf"), bbox_inches = 'tight')


plt.show()
# End Graph


################################################

# mCG level distribution
# three simulated initial states 
# 30 Col-like accessions


# Start Graph
fig, ax = plt.subplots(1,1,figsize=(3,5))

for spine in ['left','right','top','bottom']:
    ax.spines[spine].set_color('k')
    ax.spines[spine].set_linewidth(0.8)
ax.set_facecolor('white')

#ax.grid(False)
ax.grid(b=True, which='major', color='lightgrey', linestyle=':',linewidth=1)


# convert to percentage
values_scaler = 100
x_start = 0.
x_end = 100. 

x_hist = np.arange(x_start+x_end/(2.*n_bin),x_end,x_end/n_bin)

hist_Data_label = '30 Col-like Data'

temp_cols_list = []
for i_temp in range(10):
    temp_cols_list.append('D3_meth_frac_'+str(i_temp))
    temp_cols_list.append('D4_meth_frac_'+str(i_temp))
    temp_cols_list.append('D5_meth_frac_'+str(i_temp))
hist_Data, bins_Data = np.histogram(AllRepsMethFracs_Expt_df[temp_cols_list].values*values_scaler, bins=n_bin, range=(x_start,x_end))
hist_Data = hist_Data/N_reps # normalise
# ax.axvline(np.nanmean(AllRepsMethFracs_Expt_df[temp_cols_list].values*values_scaler), linestyle='--', color='g')
# ax.plot(x_hist, hist_Data,  linewidth = 2, color=col_data_Col0, linestyle='-',
#                   label='%s\nMean = %.1f %%' % (hist_Data_label,np.nanmean(AllRepsMethFracs_Expt_df[temp_cols_list].values*values_scaler)) )
ax.plot(x_hist, hist_Data,  linewidth = 2, color=col_data_Col0, linestyle='-',
                  label='%s' % (hist_Data_label) )




hist_Sims_Expt_label = 'Simulated (Col-0)'

temp_cols_list = []
for i_temp in range(N_reps):
    temp_cols_list.append('Sim_meth_frac_'+str(i_temp))
hist_Sims_Expt, bins_Sim = np.histogram(AllRepsMethFracs_Expt_df[temp_cols_list].values*values_scaler, bins=n_bin, range=(x_start,x_end))
hist_Sims_Expt = hist_Sims_Expt/N_reps # normalise
# ax.axvline(np.nanmean(AllRepsMethFracs_Expt_df[temp_cols_list].values*values_scaler), linestyle='-', color=col_Sim_Expt)
# ax.plot(x_hist, hist_Sims_Expt,  linewidth = 2, color=col_Sim_Expt,linestyle='-',
#                   label='%s\nMean = %.1f %%' % (hist_Sims_Expt_label,np.nanmean(AllRepsMethFracs_Expt_df[temp_cols_list].values*values_scaler)) )
ax.plot(x_hist, hist_Sims_Expt,  linewidth = 2, color=col_Sim_Expt,linestyle='-',
                  label='%s' % (hist_Sims_Expt_label) )


hist_Sims_100M_label = 'Simulated (All $M$)'

temp_cols_list = []
for i_temp in range(N_reps):
    temp_cols_list.append('Sim_meth_frac_'+str(i_temp))
hist_Sims_100M, bins_Sim = np.histogram(AllRepsMethFracs_100M_df[temp_cols_list].values*values_scaler, bins=n_bin, range=(x_start,x_end))
hist_Sims_100M = hist_Sims_100M/N_reps # normalise
# ax.axvline(np.nanmean(AllRepsMethFracs_100M_df[temp_cols_list].values*values_scaler), linestyle='--', color=col_Sim_100M)
# ax.plot(x_hist, hist_Sims_100M,  linewidth = 2, color=col_Sim_100M,linestyle='--',
#                   label='%s\nMean = %.1f %%' % (hist_Sims_100M_label,np.nanmean(AllRepsMethFracs_100M_df[temp_cols_list].values*values_scaler)) )
ax.plot(x_hist, hist_Sims_100M,  linewidth = 2, color=col_Sim_100M,linestyle='--',
                  label='%s' % (hist_Sims_100M_label) )


hist_Sims_100U_label = 'Simulated (All $U$)'

temp_cols_list = []
for i_temp in range(N_reps):
    temp_cols_list.append('Sim_meth_frac_'+str(i_temp))
hist_Sims_100U, bins_Sim = np.histogram(AllRepsMethFracs_100U_df[temp_cols_list].values*values_scaler, bins=n_bin, range=(x_start,x_end))
hist_Sims_100U = hist_Sims_100U/N_reps # normalise
# ax.axvline(np.nanmean(AllRepsMethFracs_100U_df[temp_cols_list].values*values_scaler), linestyle=':', color=col_Sim_100U)
# ax.plot(x_hist, hist_Sims_100U,  linewidth = 2, color=col_Sim_100U,linestyle=':',
#                   label='%s\nMean = %.1f %%' % (hist_Sims_100U_label,np.nanmean(AllRepsMethFracs_100U_df[temp_cols_list].values*values_scaler)) )
ax.plot(x_hist, hist_Sims_100U,  linewidth = 2, color=col_Sim_100U,linestyle=':',
                  label='%s' % (hist_Sims_100U_label) )


# fig.suptitle(title_sting_single_time, fontsize=12)


ax.legend(fontsize=10,ncol=1,facecolor='white', framealpha=1.0,handlelength=2.0)


ax.set_xlim(-10,110)



ax.set_ylim(0,1000)


ax.xaxis.set_tick_params(labelsize=12)
ax.yaxis.set_tick_params(labelsize=12)


ax.set_xlabel("Methylation percentage", fontsize=12)
ax.set_ylabel("Number of loci", fontsize=12)


#fig.tight_layout()
fig.subplots_adjust(top=0.5)

fig.savefig( os.path.join(GraphFolder,filename_start +"Mfrac_InitialStates_Data"+".png"), bbox_inches = 'tight')


plt.show()
# End Graph


###########################################