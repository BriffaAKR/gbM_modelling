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
CG_density_min = params.N_CG_density_min
CG_density_max = params.N_CG_density_max

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

# Drop this from simulatation dataframes
missing_gene_in_data_ID = 'AT4G23000'

file_in_LocusProperties_Data = os.path.join('Output_files', filename_start_Dat + 'LocusProperties'+filename_ending)
file_in_AllRepsMethFracs_Data = os.path.join('Output_files', filename_start_Dat + 'AllRepsMethFracs'+filename_ending)
file_in_AllRepsXFracs_Data = os.path.join('Output_files', filename_start_Dat + 'AllRepsXFracs'+filename_ending)
file_in_LocusProperties_Sim = os.path.join('Output_files', filename_start_Sim + 'LocusProperties'+filename_ending)
file_in_AllRepsMethFracs_Sim = os.path.join('Output_files', filename_start_Sim + 'AllRepsMethFracs'+filename_ending)

LocusProperties_data_df = pd.read_csv(file_in_LocusProperties_Data, sep="\t{1}",engine='python')
LocusProperties_data_df.set_index("gene_ID", inplace = True)
N_gene = len(LocusProperties_data_df)
IDs_list = LocusProperties_data_df.index.values.tolist()
print(file_in_LocusProperties_Data)
print('N_gene_data',N_gene)

LocusProperties_sim_df = pd.read_csv(file_in_LocusProperties_Sim, sep="\t{1}",engine='python')
LocusProperties_sim_df.set_index("gene_ID", inplace = True)
N_gene_sim = len(LocusProperties_sim_df)
IDs_list_sim = LocusProperties_sim_df.index.values.tolist()
print(file_in_LocusProperties_Sim)
print('N_gene_sim',N_gene)


AllRepsMethFracs_Data_df = pd.read_csv(file_in_AllRepsMethFracs_Data, sep="\t{1}",engine='python')
AllRepsMethFracs_Data_df.set_index("gene_ID", inplace = True)
AllRepsXFracs_Data_df = pd.read_csv(file_in_AllRepsXFracs_Data, sep="\t{1}",engine='python')
AllRepsXFracs_Data_df.set_index("gene_ID", inplace = True)
AllRepsMethFracs_Sim_df = pd.read_csv(file_in_AllRepsMethFracs_Sim, sep="\t{1}",engine='python')
AllRepsMethFracs_Sim_df.set_index("gene_ID", inplace = True)
AllRepsMethFracs_Sim_df = AllRepsMethFracs_Sim_df.drop([missing_gene_in_data_ID])
print(file_in_AllRepsMethFracs_Sim)
print('N_gene_Sim',len(AllRepsMethFracs_Sim_df))


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
#                 initial_seed, initial_state_choice,N_CG_min,N_CG_max, CG_density_min,CG_density_max, locus_type, 
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


hist_Data_label = '740 Col-like Data'

temp_cols_list = []
for i_temp in range(N_reps):
    temp_cols_list.append('Dset1_meth_frac_'+str(i_temp))
hist_Data, bins_Data = np.histogram(AllRepsMethFracs_Data_df[temp_cols_list].values*values_scaler, bins=n_bin, range=(x_start,x_end))
hist_Data = hist_Data/N_reps # normalise
# ax.axvline(np.nanmean(AllRepsMethFracs_Data_df[temp_cols_list].values*values_scaler), linestyle='--', color='g')
ax.plot(x_hist, hist_Data,  linewidth = 2, color='g', linestyle='-',
                  label='%s\nMean = %.1f %%' % (hist_Data_label,np.nanmean(AllRepsMethFracs_Data_df[temp_cols_list].values*values_scaler)) )


hist_Sim_label = 'Simulated'

temp_cols_list = []
for i_temp in range(N_reps):
    temp_cols_list.append('Dset1_meth_frac_'+str(i_temp))
hist_Sim, bins_Sim = np.histogram(AllRepsMethFracs_Sim_df[temp_cols_list].values*values_scaler, bins=n_bin, range=(x_start,x_end))
hist_Sim = hist_Sim/N_reps # normalise
# ax.axvline(np.nanmean(AllRepsMethFracs_Sim_df[temp_cols_list].values*values_scaler), linestyle='--', color='k')
ax.plot(x_hist, hist_Sim,  linewidth = 2, color='k',
                  label='%s\nMean = %.1f %%' % (hist_Sim_label,np.nanmean(AllRepsMethFracs_Sim_df[temp_cols_list].values*values_scaler)) )




#fig.suptitle(title_sting_single_time, fontsize=12)
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
fig.savefig( os.path.join(GraphFolder,filename_start_Graph +"Mfrac_dist_SF7B.png"), bbox_inches = 'tight')
fig.savefig( os.path.join(GraphFolder,filename_start_Graph +"Mfrac_dist_SF7B.pdf"), bbox_inches = 'tight')

plt.show()
# End Graph

###############################################

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

# temp_cols_list = []
# for i_temp in range(N_reps):
#     temp_cols_list.append('Dset1_meth_frac_'+str(i_temp))
# hist_Sim, bins_Sim = np.histogram(AllRepsXFracs_Sim_df[temp_cols_list].values*values_scaler, bins=n_bin, range=(x_start,x_end))
# hist_Sim = hist_Sim/N_reps # normalise
# ax.axvline(np.nanmean(AllRepsMethFracs_Sim_df[temp_cols_list].values*values_scaler), linestyle='--', color='k')
# ax.plot(x_hist, hist_Sim,  linewidth = 1, color='k',
#                   label='Sim.\n$\mu$ = %.2f' % (np.nanmean(AllRepsXFracs_Sim_df[temp_cols_list].values*values_scaler)))

hist_X_label = 'Col-like Data'

temp_cols_list = []
for i_temp in range(N_reps):
    temp_cols_list.append('Dset1_meth_frac_'+str(i_temp))
hist_X, bins_Data = np.histogram(AllRepsXFracs_Data_df[temp_cols_list].values*values_scaler, bins=n_bin, range=(x_start,x_end))
hist_X = hist_X/N_reps # normalise
# ax.axvline(np.nanmean(AllRepsXFracs_Data_df[temp_cols_list].values*values_scaler), linestyle='--', color='g')
ax.plot(x_hist, hist_X,  linewidth = 2, color='g', linestyle='--',
                  label='%s\nMean = %.1f %%' % (hist_X_label,np.nanmean(AllRepsXFracs_Data_df[temp_cols_list].values*values_scaler)) )

#fig.suptitle(title_sting_single_time, fontsize=12)
ax.legend(fontsize=10,ncol=1,facecolor='white', framealpha=1.0)

ax.set_xlim(-10,110)

ax.set_ylim(0,1000)

ax.xaxis.set_tick_params(labelsize=12)
ax.yaxis.set_tick_params(labelsize=12)

ax.set_xlabel("Missing data percentage", fontsize=12)
ax.set_ylabel("Number of loci", fontsize=12)

#fig.tight_layout()
#fig.subplots_adjust(top=0.5)

fig.savefig( os.path.join(GraphFolder,filename_start_Graph +"Xfrac_dist.png"), bbox_inches = 'tight')

plt.show()
# End Graph


#############

N_plots_axi = 4
N_plots_axj = 4

n_bin = 50

#N_CG_mins = [5,10,25,60,100,150]
#N_CG_maxs = [10,25,60,100,150,500]

N_CG_mins = [5,10,25,60]
N_CG_maxs = [10,25,60,500]

#CG_density_maxs = [60.,50.,40.,10.]
#CG_density_mins = [1000.,60.,50.,40.]

#CG_density_maxs = [57.,50.,43.,10.]
#CG_density_mins = [1000.,57.,50.,43.]

CG_density_maxs = [55.,47.,40.,10.]
CG_density_mins = [1000.,55.,47.,40.]

# convert to percentage
values_scaler = 100


# # Start Graph A
# fig, ax = plt.subplots(N_plots_axi,N_plots_axj,figsize=(N_plots_axj*3,N_plots_axi*3))

# x_start = 0.
# x_end = 100. 
# x_hist = np.arange(x_start+x_end/(2.*n_bin),x_end,x_end/n_bin)

# for axi in range(N_plots_axi):
#     for axj in range(N_plots_axj):

#         for spine in ['left','right','top','bottom']:
#             ax[axi,axj].spines[spine].set_color('k')
#             ax[axi,axj].spines[spine].set_linewidth(0.8)
#         ax[axi,axj].set_facecolor('white')
#         #ax.grid(False)
#         ax[axi,axj].grid(b=True, which='major', color='lightgrey', linestyle=':',linewidth=1)

#         #filter dataframes
#         Temp_Sim_df = AllRepsMethFracs_Sim_df.loc[ ((LocusProperties_sim_df['N_CG'] >= N_CG_mins[axi]) &
#                                                     (LocusProperties_sim_df['N_CG'] < N_CG_maxs[axi]) & 
#                                                     (LocusProperties_sim_df['CG_density'] >= (1./CG_density_mins[axj])) &
#                                                     (LocusProperties_sim_df['CG_density'] < (1./CG_density_maxs[axj]))) ]
#         Temp_Sim_label = '$\mu_{{Model}}$ = %.1f' % (np.nanmean(Temp_Sim_df[temp_cols_list].values*values_scaler))

#         Temp_Data_df = AllRepsMethFracs_Data_df.loc[ ((LocusProperties_data_df['N_CG'] >= N_CG_mins[axi]) &
#                                                     (LocusProperties_data_df['N_CG'] < N_CG_maxs[axi]) & 
#                                                     (LocusProperties_data_df['CG_density'] >= 1./CG_density_mins[axj]) &
#                                                     (LocusProperties_data_df['CG_density'] < 1./CG_density_maxs[axj])) ]
#         Temp_Data_label = '$\mu_{{Data}}$ = %.1f' % (np.nanmean(Temp_Data_df[temp_cols_list].values*values_scaler))
#         Title_label = "$N_{loci}$ = %d" % ( len(Temp_Data_df) )
#         N_CG_label = "%d $\leq N_{CG} <$ %d" % (N_CG_mins[axi], N_CG_maxs[axi])
#         CG_density_label = "$\\frac{{1}}{{%d}} \leq \\rho_{CG} < \\frac{{1}}{{%d}}$" % (CG_density_mins[axj], CG_density_maxs[axj])

#         # dummy data
#         ax[axi,axj].plot(x_hist,x_hist+1E5,linewidth=0,label=Title_label+'\nMeans:' )

#         hist_Sim_label = 'Sim.'

#         temp_cols_list = []
#         for i_temp in range(N_reps):
#             temp_cols_list.append('Dset1_meth_frac_'+str(i_temp))
#         hist_Sim, bins_Sim = np.histogram(Temp_Sim_df[temp_cols_list].values*values_scaler, bins=n_bin, range=(x_start,x_end))
#         hist_Sim = hist_Sim/N_reps # normalise
#         # ax[axi,axj].axvline(np.nanmean(Temp_Sim_df[temp_cols_list].values*values_scaler), linestyle=':', color='b')
#         ax[axi,axj].plot(x_hist, hist_Sim,  linewidth = 2, color='b', label='%s: %.1f %%' % (hist_Sim_label,np.nanmean(Temp_Sim_df[temp_cols_list].values*values_scaler)) )

#         hist_Data_label = 'Data'

#         temp_cols_list = []
#         for i_temp in range(N_reps):
#             temp_cols_list.append('Dset1_meth_frac_'+str(i_temp))
#         hist_Data, bins_Data = np.histogram(Temp_Data_df[temp_cols_list].values*values_scaler, bins=n_bin, range=(x_start,x_end))
#         hist_Data = hist_Data/N_reps # normalise
#         # ax[axi,axj].axvline(np.nanmean(Temp_Data_df[temp_cols_list].values*values_scaler), linestyle=':', color='k')
#         ax[axi,axj].plot(x_hist, hist_Data,  linewidth = 2, color='k', linestyle='--', label='%s: %.1f %%' % (hist_Data_label,np.nanmean(Temp_Data_df[temp_cols_list].values*values_scaler)) )



#         #fig.suptitle(title_sting_single_time, fontsize=12)
#         ax[axi,axj].legend(fontsize=14,ncol=1,facecolor='white', framealpha=1.0,handlelength=0)

#         ax[axi,axj].set_xlim(-10,110)
#         ax[axi,axj].set_ylim(0,200)

#         ax[axi,axj].xaxis.set_tick_params(labelsize=12)
#         ax[axi,axj].yaxis.set_tick_params(labelsize=12)


#         if (axi == 0):
#             ax[axi,axj].set_title(CG_density_label+'\n\n\n\n', fontsize=14)
#         if (axi == N_plots_axi -1):
#             ax[axi,axj].set_xlabel("Methylation\npercentage", fontsize=14)
#         if (axj == 0):
#             ax[axi,axj].set_ylabel(N_CG_label+"\n\n\n\n\nNumber of loci", fontsize=14)

# fig.tight_layout()
# #fig.subplots_adjust(top=0.5)

# fig.savefig( os.path.join(GraphFolder,filename_start_Graph +'Mfrac_dist_groups_Rho'+str(N_plots_axj)+'_NCG'+str(N_plots_axi)+'_A.png'), bbox_inches = 'tight')

# plt.show()
# # End Graph A


# Start Graph B
fig, ax = plt.subplots(N_plots_axi,N_plots_axj,figsize=(N_plots_axj*3,N_plots_axi*2.8))

x_start = 0.
x_end = 100. 
x_hist = np.arange(x_start+x_end/(2.*n_bin),x_end,x_end/n_bin)

for axi in range(N_plots_axi):
    for axj in range(N_plots_axj):

        for spine in ['left','right','top','bottom']:
            ax[axi,axj].spines[spine].set_color('k')
            ax[axi,axj].spines[spine].set_linewidth(0.8)
        ax[axi,axj].set_facecolor('white')
        #ax.grid(False)
        ax[axi,axj].grid(b=True, which='major', color='lightgrey', linestyle=':',linewidth=1)

        #filter dataframes
        Temp_Sim_df = AllRepsMethFracs_Sim_df.loc[ ((LocusProperties_sim_df['N_CG'] >= N_CG_mins[axi]) &
                                                    (LocusProperties_sim_df['N_CG'] < N_CG_maxs[axi]) & 
                                                    (LocusProperties_sim_df['CG_density'] >= (1./CG_density_mins[axj])) &
                                                    (LocusProperties_sim_df['CG_density'] < (1./CG_density_maxs[axj]))) ]
        Temp_Sim_label = '$\mu_{{Model}}$ = %.1f' % (np.nanmean(Temp_Sim_df[temp_cols_list].values*values_scaler))

        Temp_Data_df = AllRepsMethFracs_Data_df.loc[ ((LocusProperties_data_df['N_CG'] >= N_CG_mins[axi]) &
                                                    (LocusProperties_data_df['N_CG'] < N_CG_maxs[axi]) & 
                                                    (LocusProperties_data_df['CG_density'] >= 1./CG_density_mins[axj]) &
                                                    (LocusProperties_data_df['CG_density'] < 1./CG_density_maxs[axj])) ]
        Temp_Data_label = '$\mu_{{Data}}$ = %.1f' % (np.nanmean(Temp_Data_df[temp_cols_list].values*values_scaler))
        Title_label = "$N_{loci}$ = %d" % ( len(Temp_Data_df) )
        N_CG_label = "%d $\leq N_{CG} <$ %d" % (N_CG_mins[axi], N_CG_maxs[axi])
        CG_density_label = "$\\frac{{1}}{{%d}} \leq \\rho_{CG} < \\frac{{1}}{{%d}}$" % (CG_density_mins[axj], CG_density_maxs[axj])

        # dummy data
        #ax[axi,axj].plot(x_hist,x_hist+1E5,linewidth=0,label=Title_label )

        temp_cols_list = []
        for i_temp in range(N_reps):
            temp_cols_list.append('Dset1_meth_frac_'+str(i_temp))
        hist_Sim, bins_Sim = np.histogram(Temp_Sim_df[temp_cols_list].values*values_scaler, bins=n_bin, range=(x_start,x_end))
        hist_Sim = hist_Sim/N_reps # normalise
        ax[axi,axj].axvline(np.nanmean(Temp_Sim_df[temp_cols_list].values*values_scaler), linestyle=':', color='k')
        ax[axi,axj].plot(x_hist, hist_Sim,  linewidth = 2, color='k')

        temp_cols_list = []
        for i_temp in range(N_reps):
            temp_cols_list.append('Dset1_meth_frac_'+str(i_temp))
        hist_Data, bins_Data = np.histogram(Temp_Data_df[temp_cols_list].values*values_scaler, bins=n_bin, range=(x_start,x_end))
        hist_Data = hist_Data/N_reps # normalise
        ax[axi,axj].axvline(np.nanmean(Temp_Data_df[temp_cols_list].values*values_scaler), linestyle=':', color='g')
        ax[axi,axj].plot(x_hist, hist_Data,  linewidth = 2, color='g', linestyle='--')



        #fig.suptitle(title_sting_single_time, fontsize=12)
        #ax[axi,axj].legend(fontsize=14,ncol=1,facecolor='white', framealpha=1.0,handlelength=0)

        ax[axi,axj].set_xlim(-10,110)
        #ax[axi,axj].set_ylim(0,200)

        ax[axi,axj].xaxis.set_tick_params(labelsize=12)
        ax[axi,axj].yaxis.set_tick_params(labelsize=12)

        ax[axi,axj].set_title(Title_label, fontsize=14)

        if (axi == 0):
            ax[axi,axj].set_title(CG_density_label+'\n\n\n\n'+Title_label, fontsize=14)
        if (axi == N_plots_axi -1):
            ax[axi,axj].set_xlabel("Methylation\npercentage", fontsize=14)
        if (axj == 0):
            ax[axi,axj].set_ylabel(N_CG_label+"\n\n\n\n\nNumber of loci", fontsize=14)

ax[0,0].arrow(0, 152, 560, 0, fc='k', ec='k', clip_on=False, width=0.8, head_width=8, head_length=12)
ax[0,0].arrow(-90, 100, 0, -540, fc='k', ec='k', clip_on=False, width=0.8, head_width=8, head_length=12)
ax[0,-1].annotate("Increasing CG-site density", (-80, 250), fontsize=14, annotation_clip=False)
ax[-1,0].annotate("Increasing locus size", (-80, 50), fontsize=14, annotation_clip=False,rotation=90)

fig.tight_layout()
#fig.subplots_adjust(top=0.5)

fig.savefig( os.path.join(GraphFolder,filename_start_Graph +'Mfrac_dist_groups_Rho'+str(N_plots_axj)+'_NCG'+str(N_plots_axi)+'_SF7C.png'), bbox_inches = 'tight')
fig.savefig( os.path.join(GraphFolder,filename_start_Graph +'Mfrac_dist_groups_Rho'+str(N_plots_axj)+'_NCG'+str(N_plots_axi)+'_SF7C.pdf'), bbox_inches = 'tight')

plt.show()
# End Graph B


###########################################

# Start Graph
fig, ax = plt.subplots(1,1,figsize=(3,3))

for spine in ['left','right','top','bottom']:
    ax.spines[spine].set_color('k')
    ax.spines[spine].set_linewidth(0.8)
ax.set_facecolor('white')

#ax.grid(False)
ax.grid(b=True, which='major', color='lightgrey', linestyle=':',linewidth=1)


n_bin_2 = 200

# convert to percentage
values_scaler = 100
x_start = 0.
x_end = 100. 

x_hist = np.arange(x_start+x_end/(2.*n_bin_2),x_end,x_end/n_bin_2)


hist_D3D4D5, bins_D3D4D5 = np.histogram(LocusProperties_data_df[['Dset1_mu_meth_frac']].values*values_scaler, bins=n_bin_2, range=(x_start,x_end))
# ax.axvline(np.nanmean(LocusProperties_data_df[['Dset1_mu_meth_frac']].values*values_scaler), linestyle='-', color='green')
ax.plot(x_hist, hist_D3D4D5,  linewidth = 2, color='green',
                  label='Data Col-like\n$\mu$ = %.4f' % (np.nanmean(LocusProperties_data_df[['Dset1_mu_meth_frac']].values*values_scaler)) )


hist_Sims_Expt, bins_Sim = np.histogram(LocusProperties_sim_df[['Sim_mu_meth_frac']].values*values_scaler, bins=n_bin_2, range=(x_start,x_end))
# ax.axvline(np.nanmean(LocusProperties_sim_df[['Sim_mu_meth_frac']].values*values_scaler), linestyle='-', color='k')
ax.plot(x_hist, hist_Sims_Expt,  linewidth = 2, color='k',linestyle='-',
                  label='Mean Over Reps\nSim. Expt.\n$\mu$ = %.4f' % (np.nanmean(LocusProperties_sim_df[['Sim_mu_meth_frac']].values*values_scaler)) )



# fig.suptitle(title_sting_single_time, fontsize=12)


ax.legend(fontsize=10,ncol=1,facecolor='white', framealpha=1.0)


ax.set_xlim(-10,110)


#ax.set_ylim(0,200)


ax.xaxis.set_tick_params(labelsize=12)
ax.yaxis.set_tick_params(labelsize=12)


ax.set_xlabel("Methylation percentage", fontsize=12)
ax.set_ylabel("Number of loci", fontsize=12)


#fig.tight_layout()
# fig.subplots_adjust(top=0.5)

fig.savefig( os.path.join(GraphFolder,filename_start_Graph +"Mfrac_InitialState_Expt_MeanOverReps"+".png"), bbox_inches = 'tight')

plt.show()
# End Graph


################################################