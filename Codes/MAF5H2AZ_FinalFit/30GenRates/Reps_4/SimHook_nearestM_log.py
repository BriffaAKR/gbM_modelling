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

# fig_label = "_FigS4A"
# sim_type = params.sim_type

# fetch parameters

# define annoation files
file_in_annotation = params.file_in_annotation
file_in_anno_filt_1 = params.file_in_anno_filt_1
file_in_anno_filt_2 = params.file_in_anno_filt_2
file_in_exclude_filt_1 = params.file_in_exclude_filt_1

file_in_Parental_Col0_state = params.file_in_Parental_Col0_state

filename_params_code = params.filename_params_code

locus_type = params.locus_type
TE_filt = params.TE_filt

file_in_GainRateData = params.file_in_GainRateData
file_in_LossRateData = params.file_in_LossRateData



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


GraphFolder = "Graphs"

InputFiles_folder = params.InputFiles_folder
filename_start = params.filename_start
filename_ending = '.tsv'

# build up input path string
path_string = os.path.join( os.getcwd(), '..')
path_string = os.path.join( path_string, '..')
path_string = os.path.join( path_string, '..')
path_string = os.path.join( path_string, '..')
path_string = os.path.join( path_string, InputFiles_folder)

# load in gain/loss data
# NearestM
GainRateData_file = os.path.join( path_string, file_in_GainRateData )
LossRateData_file = os.path.join( path_string, file_in_LossRateData ) 

GainRate_data_df = pd.read_csv(GainRateData_file, sep="\t{1}", engine='python')
LossRate_data_df = pd.read_csv(LossRateData_file, sep="\t{1}", engine='python')



# moving average using convolution
def calc_moving_average(N_window_, list_in_):


    list_out_ = np.convolve( list_in_, np.ones((N_window_,))/N_window_, mode='same')


    # If N_window_ = odd
    # start and end edges are both of length (N_window_-1.)/2.
    # range (0, (N_window_-1.)/2.)
    # range (-(N_window_-1.)/2.,0)

    #If N_window_ = even 
    # start edge of length = N_window_/2.
    # end window_ of length = (N_window_/2.-1.)
    # range (0, N_window_/2.)
    # range (-(N_window_/2.-1.),0)

    # normalisation corrections to ends: 
    # If N_window_ = odd (eg. 5)
    # (3,4) and (4,3)
    # ((N_window_+1.)/2. , (N_window_+1.)/2. +1) and ((N_window_+1.)/2. +1 , (N_window_+1.)/2)
    # If N_window_ = even (eg. 6)
    # (3,4,5) and (5,4)
    # (N_window_/2. , N_window_/2.+1 , N_window_/2.+2) and (N_window_/2.+2 , N_window_/2.+1)


    if (N_window_%2) == 0:
        #print('N_window_ even')
        for i_temp in range(0, int(N_window_/2.)):
            #print(i_temp, i_temp + N_window_/2.)
            list_out_[i_temp] = list_out_[i_temp]*N_window_/(i_temp + N_window_/2.)
        for i_temp in range(int(-(N_window_/2.-1.)),0):
            #print(i_temp, N_window_/2. - i_temp)
            list_out_[i_temp] = list_out_[i_temp]*N_window_/(N_window_/2. - i_temp)
    if (N_window_%2) ==1:
        #print('N_window_ odd')
        for i_temp in range(0, int((N_window_-1.)/2.)):
            #print(i_temp, i_temp + (N_window_+1.)/2.)
            list_out_[i_temp] = list_out_[i_temp]*N_window_/(i_temp + (N_window_+1.)/2.)
        for i_temp in range(int(-(N_window_-1.)/2.),0):
            #print(i_temp,  (N_window_-1.)/2. -i_temp)
            list_out_[i_temp] = list_out_[i_temp]*N_window_/((N_window_-1.)/2. -i_temp)
            
    return list_out_


# define linear function to fit to. 
def linear_func(x_, m_lin_, c_lin_):
    return m_lin_*x_ + c_lin_

def fit_trendline_lin(x_,y_, N_):

    # N is number of data points to calc. fit line at
    # fit linear trendlines to delta values

    # find variance
    y_var_ = np.var(y_)

    ## linear
    coef_lin_ , stats_lin_ = np.polynomial.polynomial.polyfit(x_,y_,1,full=True)
    m_lin_ = coef_lin_[-1]; c_lin_ = coef_lin_[-2]
    resids_lin_ = stats_lin_[0][0]
    R_sq_lin_ = 1. - resids_lin_/(len(y_)*y_var_)

    x_vals_ = np.arange(0,N_ + 1, 10)

    lin_vals_ = [ linear_func(i_, m_lin_, c_lin_) for i_ in x_vals_]

    return (x_vals_,lin_vals_,R_sq_lin_, m_lin_, c_lin_)

####
Log10_x_vals = np.log10(GainRate_data_df['x'].tolist()[1:])
x_vals = GainRate_data_df['x'].tolist()[1:]

# Gain
Log10_gains_Lizzie_conservative = np.log10(GainRate_data_df['Lizzie_Conservative_GainRate'].tolist()[1:])
Log10_gains_All_mean = np.log10(GainRate_data_df['Average_GainRate'].tolist()[1:])
Log10_gains_Schmitz_mean = np.log10(GainRate_data_df['Schmitz_GainRate'].tolist()[1:])


# Loss
Log10_losses_Lizzie_conservative = np.log10(LossRate_data_df['Lizzie_Conservative_LossRate'].tolist()[1:])
Log10_losses_All_mean = np.log10(LossRate_data_df['Average_LossRate'].tolist()[1:])
Log10_losses_Schmitz_mean = np.log10(LossRate_data_df['Schmitz_LossRate'].tolist()[1:])

noise_start = 500

gains_start = 25
gains_end = 47
losses_start = 20
losses_end = 200

# remove missing data from gains mean
x_fit_input_gains = []
y_fit_input_gains = []

for i_ in range(gains_start,gains_end):
    if Log10_gains_All_mean[i_] > -1.e+50:
        x_fit_input_gains.append(Log10_x_vals[i_])
        y_fit_input_gains.append(Log10_gains_All_mean[i_])


x_vals_linfit_gains ,lin_vals_linfit_gains ,R_sq_lin_gains , m_lin_gains , c_lin_gains = fit_trendline_lin(
    x_fit_input_gains,y_fit_input_gains, 10)


# remove missing data from loss mean
x_fit_input_losses = []
y_fit_input_losses = []

for i_ in range(losses_start,losses_end):
    if Log10_losses_All_mean[i_] > -1.e+50:
        x_fit_input_losses.append(Log10_x_vals[i_])
        y_fit_input_losses.append(Log10_losses_All_mean[i_])
        
        
x_vals_linfit_losses ,lin_vals_linfit_losses ,R_sq_lin_losses , m_lin_losses , c_lin_losses = fit_trendline_lin(
    x_fit_input_losses,y_fit_input_losses, 10)




# Start Graph

fig0, ax0 = plt.subplots(2,1,figsize=(12,12))

for spine in ['left','right','top','bottom']:
    ax0[0].spines[spine].set_color('k')
    ax0[1].spines[spine].set_color('k')
ax0[0].set_facecolor('white')
ax0[1].set_facecolor('white')

ax0[0].grid(False)
ax0[1].grid(False)

linewidth=linewidth_choice = 2

# Plots NearestM
ax0[0].plot(Log10_x_vals,Log10_gains_All_mean,label='Expt. mean\n gain rate',color='grey',linestyle='-',linewidth=linewidth_choice)
ax0[0].plot(x_vals_linfit_gains, lin_vals_linfit_gains, label='Fit to range:\n%d $< x \leq$ %d' 
         % (gains_start, gains_end), color='r',linewidth=linewidth_choice)

ax0[1].plot(Log10_x_vals,Log10_losses_All_mean,label='Expt. mean\n loss rate', color='grey',linestyle='-',linewidth=linewidth_choice)
ax0[1].plot(x_vals_linfit_losses, lin_vals_linfit_losses, label='Fit to range:\n%d $< x \leq$ %d' 
         % (losses_start, losses_end), color='b',linewidth=linewidth_choice)


ax0[0].axvline(x=np.log10(noise_start), linestyle='--', color='k',label='x = %d bp' % noise_start,linewidth=linewidth_choice)
# ax0[1].axvline(x=np.log10(noise_start), linestyle='--', color='k',label='x = %d bp' % noise_start,linewidth=linewidth_choice)

ax0[0].axvline(x=np.log10(gains_start), linestyle='--', color='r',linewidth=linewidth_choice)
ax0[0].axvline(x=np.log10(gains_end), linestyle='--', color='r',linewidth=linewidth_choice)

ax0[1].axvline(x=np.log10(losses_start), linestyle='--', color='b',linewidth=linewidth_choice)
ax0[1].axvline(x=np.log10(losses_end), linestyle='--', color='b',linewidth=linewidth_choice)



ax0[0].set_xlabel("Log$_{10}$(bp distance to nearest $M$-site)", fontsize=18)
ax0[0].set_ylabel("Log$_{10}$(gain rate per c.c.)", fontsize=18)

ax0[0].legend(fontsize=16,ncol=1,facecolor='white', framealpha=1)
ax0[0].tick_params(axis='both', which='major', labelsize=18)

ax0[1].set_xlabel("Log$_{10}$(bp distance to nearest $M$-site)", fontsize=18)
ax0[1].set_ylabel("Log$_{10}$(loss rate per c.c.)", fontsize=18)

ax0[1].legend(fontsize=16,ncol=1,facecolor='white', framealpha=1)
ax0[1].tick_params(axis='both', which='major', labelsize=18)



ax0[0].set_xlim(0,3.5)
ax0[0].set_ylim(-7,-2)

ax0[1].set_xlim(0,3.5)
ax0[1].set_ylim(-7,-2)

# ax0[0].yaxis.set_major_locator(plt.MaxNLocator(7))
# ax0[1].yaxis.set_major_locator(plt.MaxNLocator(7))



#fig0.tight_layout()
fig0.subplots_adjust(hspace=0.3)

fig0.savefig( os.path.join(GraphFolder,filename_start +'DataRates_NearestM_LogLog_large.png'), bbox_inches = 'tight')

# End Graph


###################
# Start Graph

fig0, ax0 = plt.subplots(2,1,figsize=(6,6))

for spine in ['left','right','top','bottom']:
    ax0[0].spines[spine].set_color('k')
    ax0[1].spines[spine].set_color('k')
ax0[0].set_facecolor('white')
ax0[1].set_facecolor('white')

ax0[0].grid(False)
ax0[1].grid(False)

linewidth=linewidth_choice = 1

# Plots NearestM
ax0[0].plot(Log10_x_vals,Log10_gains_All_mean,label='Expt. mean',color='grey',linestyle='-',linewidth=linewidth_choice)
ax0[0].plot(x_vals_linfit_gains, lin_vals_linfit_gains, label='Fit to range:\n%d $< x \leq$ %d' 
         % (gains_start, gains_end), color='r',linewidth=linewidth_choice)

ax0[1].plot(Log10_x_vals,Log10_losses_All_mean,label='Expt. mean', color='grey',linestyle='-',linewidth=linewidth_choice)
ax0[1].plot(x_vals_linfit_losses, lin_vals_linfit_losses, label='Fit to range:\n%d $< x \leq$ %d' 
         % (losses_start, losses_end), color='b',linewidth=linewidth_choice)


# ax0[0].axvline(x=np.log10(noise_start), linestyle='--', color='k',label='x = %d bp' % noise_start,linewidth=linewidth_choice)
# ax0[1].axvline(x=np.log10(noise_start), linestyle='--', color='k',label='x = %d bp' % noise_start,linewidth=linewidth_choice)

ax0[0].axvline(x=np.log10(noise_start), linestyle='--', color='k', linewidth=linewidth_choice)
# ax0[1].axvline(x=np.log10(noise_start), linestyle='--', color='k', linewidth=linewidth_choice)

ax0[0].axvline(x=np.log10(gains_start), linestyle='--', color='r',linewidth=linewidth_choice)
ax0[0].axvline(x=np.log10(gains_end), linestyle='--', color='r',linewidth=linewidth_choice)

ax0[1].axvline(x=np.log10(losses_start), linestyle='--', color='b',linewidth=linewidth_choice)
ax0[1].axvline(x=np.log10(losses_end), linestyle='--', color='b',linewidth=linewidth_choice)



ax0[0].set_xlabel("Log$_{10}$(bp distance to nearest $M$-site)", fontsize=12)
ax0[0].set_ylabel("Log$_{10}$(gain rate per c.c.)", fontsize=12)

ax0[0].legend(fontsize=10,ncol=1,facecolor='white', framealpha=1)
ax0[0].tick_params(axis='both', which='major', labelsize=12)

ax0[1].set_xlabel("Log$_{10}$(bp distance to nearest $M$-site)", fontsize=12)
ax0[1].set_ylabel("Log$_{10}$(loss rate per c.c.)", fontsize=12)

ax0[1].legend(fontsize=10,ncol=1,facecolor='white', framealpha=1)
ax0[1].tick_params(axis='both', which='major', labelsize=12)



ax0[0].set_xlim(0,3.5)
ax0[0].set_ylim(-7,-2)

ax0[1].set_xlim(0,3.5)
ax0[1].set_ylim(-7,-2)

# ax0[0].yaxis.set_major_locator(plt.MaxNLocator(7))
# ax0[1].yaxis.set_major_locator(plt.MaxNLocator(7))



#fig0.tight_layout()
fig0.subplots_adjust(hspace=0.3)

fig0.savefig( os.path.join(GraphFolder,filename_start +'DataRates_NearestM_LogLog_Pair.png'), bbox_inches = 'tight')

# End Graph

###################
###################

# Start Graph

fig0, ax0 = plt.subplots(2,1,figsize=(12,12))

for spine in ['left','right','top','bottom']:
    ax0[0].spines[spine].set_color('k')
    ax0[1].spines[spine].set_color('k')
ax0[0].set_facecolor('white')
ax0[1].set_facecolor('white')

ax0[0].grid(False)
ax0[1].grid(False)

linewidth=linewidth_choice = 2

# Plots NearestM
ax0[0].plot(x_vals,Log10_gains_All_mean,label='Expt. mean\n gain rate',color='grey',linestyle='-',linewidth=linewidth_choice)
# ax0[0].plot(x_vals_linfit_gains, lin_vals_linfit_gains, label='Fit to range:\n%d $< x \leq$ %d' 
#          % (gains_start, gains_end), color='r',linewidth=linewidth_choice)

ax0[1].plot(x_vals,Log10_losses_All_mean,label='Expt. mean\n loss rate', color='grey',linestyle='-',linewidth=linewidth_choice)
# ax0[1].plot(x_vals_linfit_losses, lin_vals_linfit_losses, label='Fit to range:\n%d $< x \leq$ %d' 
#          % (losses_start, losses_end), color='b',linewidth=linewidth_choice)


# ax0[0].axvline(x=np.log10(noise_start), linestyle='--', color='k',label='x = %d bp' % noise_start,linewidth=linewidth_choice)
# ax0[1].axvline(x=np.log10(noise_start), linestyle='--', color='k',label='x = %d bp' % noise_start,linewidth=linewidth_choice)

ax0[0].axvline(x=gains_start, linestyle='--', color='r',linewidth=linewidth_choice)
ax0[0].axvline(x=gains_end, linestyle='--', color='r',linewidth=linewidth_choice)

ax0[1].axvline(x=losses_start, linestyle='--', color='b',linewidth=linewidth_choice)
ax0[1].axvline(x=losses_end, linestyle='--', color='b',linewidth=linewidth_choice)



ax0[0].set_xlabel("distance to nearest $M$-site (bp)", fontsize=18)
ax0[0].set_ylabel("Log$_{10}$(gain rate per c.c.)", fontsize=18)

ax0[0].legend(fontsize=16,ncol=1,facecolor='white', framealpha=1)
ax0[0].tick_params(axis='both', which='major', labelsize=18)

ax0[1].set_xlabel("distance to nearest $M$-ite (bp)", fontsize=18)
ax0[1].set_ylabel("Log$_{10}$(loss rate per c.c.)", fontsize=18)

ax0[1].legend(fontsize=16,ncol=1,facecolor='white', framealpha=1)
ax0[1].tick_params(axis='both', which='major', labelsize=18)



ax0[0].set_xlim(0,500)
ax0[0].set_ylim(-7,-2)

ax0[1].set_xlim(0,500)
ax0[1].set_ylim(-7,-2)

# ax0[0].yaxis.set_major_locator(plt.MaxNLocator(7))
# ax0[1].yaxis.set_major_locator(plt.MaxNLocator(7))



#fig0.tight_layout()
fig0.subplots_adjust(hspace=0.3)

fig0.savefig( os.path.join(GraphFolder,filename_start +'DataRates_NearestM_LogLin_large.png'), bbox_inches = 'tight')

# End Graph


# # Start Graph

# fig0, ax0 = plt.subplots(2,1,figsize=(6,6))

# for spine in ['left','right','top','bottom']:
#     ax0[0].spines[spine].set_color('k')
#     ax0[1].spines[spine].set_color('k')
#     ax0[0].spines[spine].set_linewidth(0.8)
#     ax0[1].spines[spine].set_linewidth(0.8)
# ax0[0].set_facecolor('white')
# ax0[1].set_facecolor('white')

# #ax0[0].grid(False)
# #ax0[1].grid(False)
# ax0[0].grid(b=True, which='major', color='lightgrey', linestyle=':',linewidth=1)
# ax0[1].grid(b=True, which='major', color='lightgrey', linestyle=':',linewidth=1)

# # Plots NearestM
# #Data
# # gains
# print(len(x_vals[2:]),len(MovAv_gains_All_Min),len(MovAv_gains_All_Max))
# ax0[0].fill_between(x_vals[2:],y_scale*MovAv_gains_All_Min, y_scale*MovAv_gains_All_Max, color='lightgrey', alpha = 0.8, label='All range $\pm 1\sigma$')
# ax0[0].fill_between(x_vals[2:],y_scale*MovAv_gains_Schmitz_Min, y_scale*MovAv_gains_Schmitz_Max, color='plum', alpha = 0.6, label='Schmitz range $\pm 1\sigma$')
# ax0[0].plot(x_vals[2:], y_scale*MovAv_gains_All_mean, linewidth=2, color='grey',alpha=0.8, label = 'All mean')
# ax0[0].plot(x_vals[2:], y_scale*MovAv_gains_Schmitz_mean, linewidth=2, color='k', label = 'Schmitz mean')
# ax0[0].plot(x_vals[2:], y_scale*MovAv_gains_Lizzie_conservative, linewidth=1, color='k',linestyle='--', label = 'Hollwey')


# # losses 
# ax0[1].fill_between(x_vals[2:],y_scale*MovAv_losses_All_Min, y_scale*MovAv_losses_All_Max, color='lightgrey', alpha = 0.8, label='All range $\pm 1\sigma$')
# ax0[1].fill_between(x_vals[2:],y_scale*MovAv_losses_Schmitz_Min, y_scale*MovAv_losses_Schmitz_Max, color='lightseagreen', alpha = 0.6, label='Schmitz range $\pm 1\sigma$')
# ax0[1].plot(x_vals[2:], y_scale*MovAv_losses_All_mean, linewidth=2, color='grey',alpha=0.8, label = 'All mean')
# ax0[1].plot(x_vals[2:], y_scale*MovAv_losses_Schmitz_mean, linewidth=2, color='k', label = 'Schmitz mean')
# ax0[1].plot(x_vals[2:], y_scale*MovAv_losses_Lizzie_conservative, linewidth=1, color='k',linestyle='--', label = 'Hollwey')

# #Model
# ax0[0].plot(x_vals[2:], y_scale*MovAv_SimGainRate_nn_dist_M,label='Simulated\nSmoothed n = %d bp' % N_window, linewidth=2, color='r')
# ax0[1].plot(x_vals[2:], y_scale*MovAv_SimLossRate_nn_dist_M,label='Simulated\nSmoothed n = %d bp' % N_window, linewidth=2, color='b')




# ax0[0].set_ylim(0,3)

# #ax0[0].set_xlabel("Distance to nearest M site (bp)", fontsize=18)
# ax0[0].set_ylabel("Gain rate\n per cell cycle"+y_scale_label, fontsize=18)
# # ,loc=(0.18,0.40)
# #ax0[0].legend(fontsize=18,ncol=2,facecolor='white', framealpha=0.6)
# ax0[0].tick_params(axis='both', which='major', labelsize=16)

# ax0[1].set_ylim(0,3)
# ax0[1].set_xlabel("Distance to nearest M site (bp)", fontsize=18)
# ax0[1].set_ylabel("Loss rate\n per cell cylce"+y_scale_label, fontsize=18)
# #ax0[1].legend(fontsize=18,ncol=2,loc=2,facecolor='white', framealpha=0.6)
# ax0[1].tick_params(axis='both', which='major', labelsize=16)

# # fig0.suptitle(title_sting_single_time, fontsize=20)

# #ax0[0].set_xlim(0,800)
# #ax0[1].set_xlim(0,250)
# ax0[0].set_xlim(0,500)
# ax0[1].set_xlim(0,500)



# #ax0[0].yaxis.set_major_locator(plt.MaxNLocator(7))
# #ax0[1].yaxis.set_major_locator(plt.MaxNLocator(7))



# #fig0.tight_layout()
# #fig0.subplots_adjust(hspace=0.2)
# # fig0.subplots_adjust(top=0.68)

# fig0.savefig( os.path.join(GraphFolder,filename_start +'SimRates_NearestM_Pair'+'.png'), bbox_inches = 'tight')

# # End Graph


# # Start Graph

# fig0, ax0 = plt.subplots(1,1,figsize=(6,3))

# for spine in ['left','right','top','bottom']:
#     ax0.spines[spine].set_color('k')
#     ax0.spines[spine].set_linewidth(0.8)
# ax0.set_facecolor('white')

# #ax0.grid(False)
# ax0.grid(b=True, which='major', color='lightgrey', linestyle=':',linewidth=1)

# # Plots NearestM
# #Data
# # gains
# print(len(x_vals[2:]),len(MovAv_gains_All_Min),len(MovAv_gains_All_Max))
# ax0.fill_between(x_vals[2:],y_scale*MovAv_gains_All_Min, y_scale*MovAv_gains_All_Max, color='lightgrey', alpha = 0.8, label='All range $\pm 1\sigma$')
# ax0.fill_between(x_vals[2:],y_scale*MovAv_gains_Schmitz_Min, y_scale*MovAv_gains_Schmitz_Max, color='plum', alpha = 0.6, label='Schmitz range $\pm 1\sigma$')
# ax0.plot(x_vals[2:], y_scale*MovAv_gains_All_mean, linewidth=2, color='grey',alpha=0.8, label = 'All mean')
# ax0.plot(x_vals[2:], y_scale*MovAv_gains_Schmitz_mean, linewidth=2, color='k', label = 'Schmitz mean')
# ax0.plot(x_vals[2:], y_scale*MovAv_gains_Lizzie_conservative, linewidth=1, color='k',linestyle='--', label = 'Hollwey')


# #Model
# ax0.plot(x_vals[2:], y_scale*MovAv_SimGainRate_nn_dist_M,label='Simulated\nSmoothed n = %d bp' % N_window, linewidth=2, color='r')



# ax0.set_ylim(0,3)

# ax0.set_xlabel("Distance to nearest M site (bp)", fontsize=18)
# ax0.set_ylabel("Gain rate\n per cell cycle"+y_scale_label, fontsize=18)
# # ,loc=(0.18,0.40)
# #ax0.legend(fontsize=18,ncol=2,facecolor='white', framealpha=0.6)
# ax0.tick_params(axis='both', which='major', labelsize=16)


# ax0.set_xlim(0,500)



# #ax0.yaxis.set_major_locator(plt.MaxNLocator(7))

# #fig0.tight_layout()
# #fig0.subplots_adjust(hspace=0.2)

# fig0.savefig( os.path.join(GraphFolder,filename_start +'SimRates_NearestM_Gains'+'.png'), bbox_inches = 'tight')

# # End Graph

# #######################

# # Start Graph

# fig0, ax0 = plt.subplots(1,1,figsize=(6,3))

# for spine in ['left','right','top','bottom']:
#     ax0.spines[spine].set_color('k')
#     ax0.spines[spine].set_linewidth(0.8)
# ax0.set_facecolor('white')

# #ax0.grid(False)
# ax0.grid(b=True, which='major', color='lightgrey', linestyle=':',linewidth=1)

# # Plots NearestM
# #Data

# # losses 
# ax0.fill_between(x_vals[2:],y_scale*MovAv_losses_All_Min, y_scale*MovAv_losses_All_Max, color='lightgrey', alpha = 0.8, label='All range $\pm 1\sigma$')
# ax0.fill_between(x_vals[2:],y_scale*MovAv_losses_Schmitz_Min, y_scale*MovAv_losses_Schmitz_Max, color='lightseagreen', alpha = 0.6, label='Schmitz range $\pm 1\sigma$')
# ax0.plot(x_vals[2:], y_scale*MovAv_losses_All_mean, linewidth=2, color='grey',alpha=0.8, label = 'All mean')
# ax0.plot(x_vals[2:], y_scale*MovAv_losses_Schmitz_mean, linewidth=2, color='k', label = 'Schmitz mean')
# ax0.plot(x_vals[2:], y_scale*MovAv_losses_Lizzie_conservative, linewidth=1, color='k',linestyle='--', label = 'Hollwey')

# #Model
# ax0.plot(x_vals[2:], y_scale*MovAv_SimLossRate_nn_dist_M,label='Simulated\nSmoothed n = %d bp' % N_window, linewidth=2, color='b')




# ax0.set_xlim(0,500)

# #ax0.yaxis.set_major_locator(plt.MaxNLocator(7))


# #fig0.tight_layout()
# #fig0.subplots_adjust(hspace=0.2)
# # fig0.subplots_adjust(top=0.68)

# fig0.savefig( os.path.join(GraphFolder,filename_start +'SimRates_NearestM_Losses'+'.png'), bbox_inches = 'tight')

# # End Graph