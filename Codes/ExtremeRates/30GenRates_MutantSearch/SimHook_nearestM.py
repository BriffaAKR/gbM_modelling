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
coop_strength = params.coop_strength

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
# path_string = os.path.join( path_string, '..')
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

N_window = 10

N_window_gain = N_window
N_window_loss = N_window

N_window_gain_Jay = N_window_gain
N_window_loss_Jay = N_window_loss
N_window_gain_Lizzie = N_window_gain
N_window_loss_Lizzie = N_window_loss
N_window_gain_Sim = N_window_gain
N_window_loss_Sim = N_window_loss

x_vals = np.arange(0,N_corrns_bins,1)

# NearestM
# Gain smoothing data
MovAv_gains_Lizzie_conservative = calc_moving_average(N_window_gain_Lizzie, GainRate_data_df['Lizzie_Conservative_GainRate'].tolist()[1:])
MovAv_gains_All_mean = calc_moving_average(N_window_gain_Jay, GainRate_data_df['Average_GainRate'].tolist()[1:])
MovAv_gains_All_Max = calc_moving_average(N_window_gain_Jay, GainRate_data_df['All_Max_GainRate'].tolist()[1:])
MovAv_gains_All_Min = calc_moving_average(N_window_gain_Jay, GainRate_data_df['All_Min_GainRate'].tolist()[1:])
MovAv_gains_Schmitz_mean = calc_moving_average(N_window_gain_Jay, GainRate_data_df['Schmitz_GainRate'].tolist()[1:])
MovAv_gains_Schmitz_Max = calc_moving_average(N_window_gain_Jay, GainRate_data_df['Schmitz_Max_GainRate'].tolist()[1:])
MovAv_gains_Schmitz_Min = calc_moving_average(N_window_gain_Jay, GainRate_data_df['Schmitz_Min_GainRate'].tolist()[1:])

# Loss smoothing data
MovAv_losses_Lizzie_conservative = calc_moving_average(N_window_loss_Lizzie, LossRate_data_df['Lizzie_Conservative_LossRate'].tolist()[1:])
MovAv_losses_All_mean = calc_moving_average(N_window_loss_Jay, LossRate_data_df['Average_LossRate'].tolist()[1:])
MovAv_losses_All_Max = calc_moving_average(N_window_loss_Jay, LossRate_data_df['All_Max_LossRate'].tolist()[1:])
MovAv_losses_All_Min = calc_moving_average(N_window_loss_Jay, LossRate_data_df['All_Min_LossRate'].tolist()[1:])
MovAv_losses_Schmitz_mean = calc_moving_average(N_window_loss_Jay, LossRate_data_df['Schmitz_LossRate'].tolist()[1:])
MovAv_losses_Schmitz_Max = calc_moving_average(N_window_loss_Jay, LossRate_data_df['Schmitz_Max_LossRate'].tolist()[1:])
MovAv_losses_Schmitz_Min = calc_moving_average(N_window_loss_Jay, LossRate_data_df['Schmitz_Min_LossRate'].tolist()[1:])



# Load in CodeProgress to get total number of simulated IDs
file_in_CodeProgress = os.path.join('Output_files', filename_start + 'CodeProgress'+filename_ending)
CodeProgress_df = pd.read_csv(file_in_CodeProgress, sep="\t{1}",engine='python')
# sum N_gene column
N_gene = CodeProgress_df['N_gene'].sum(axis=0)
file_in_Skipped_IDs = os.path.join('Output_files', filename_start + 'Skipped_IDs'+filename_ending)
Skipped_IDs_df = pd.read_csv(file_in_Skipped_IDs, sep="\t{1}",engine='python')
N_gene = N_gene - len(Skipped_IDs_df) - 1
print('N_gene',N_gene)
print()

# Calc. simulated overall gains and losses. 
Total_Gains = CodeProgress_df['Sim_total_gains'].sum()
Total_Gains_per_rep = Total_Gains/N_reps
Total_Losses = CodeProgress_df['Sim_total_losses'].sum()
Total_Losses_per_rep = Total_Losses/N_reps
Total_StayM = CodeProgress_df['Sim_total_stayM'].sum()
Total_StayU = CodeProgress_df['Sim_total_statyU'].sum()
Total_GainRate = Total_Gains/(Total_Gains+Total_StayU)
Total_LossRate = Total_Losses/(Total_Losses+Total_StayM)
Total_Gains_per_cc = Total_Gains_per_rep/(N_gen_burn_in*n_cc)
Total_Losses_per_cc = Total_Losses_per_rep/(N_gen_burn_in*n_cc)
Total_GainRate_per_cc = Total_GainRate/(N_gen_burn_in*n_cc)
Total_LossRate_per_cc = Total_LossRate/(N_gen_burn_in*n_cc)
# Fraction of M-sites in initial state
Total_M_0 = (Total_Losses + Total_StayM)/(Total_Losses + Total_StayM + Total_Gains + Total_StayU)
# Write out to file
file_name_output_Totals = os.path.join(GraphFolder,filename_start +'RatesTotals'+'.tsv')
with open(file_name_output_Totals,'w',newline='') as output_file:
    line_writer = csv.writer(output_file,delimiter='\t')
    line_writer.writerow( ['delta','CoopStrength','epsilon','Total_GainRate_per_cc','Total_LossRate_per_cc','Total_M_0'] )
    line_writer.writerow( [delta,coop_strength,e_c_val,Total_GainRate_per_cc,Total_LossRate_per_cc,Total_M_0] )




# Load in Simulated Changes
file_in_SimGains_nn_dist_M = os.path.join('Output_files', filename_start + 'SimGains_nn_dist_M'+filename_ending)
file_in_SimLosses_nn_dist_M = os.path.join('Output_files', filename_start + 'SimLosses_nn_dist_M'+filename_ending)
file_in_SimStayM_nn_dist_M = os.path.join('Output_files', filename_start + 'SimStayM_nn_dist_M'+filename_ending)
file_in_SimStayU_nn_dist_M = os.path.join('Output_files', filename_start + 'SimStayU_nn_dist_M'+filename_ending)



SimGains_nn_dist_M_df = pd.read_csv(file_in_SimGains_nn_dist_M, sep="\t{1}",engine='python')
SimLosses_nn_dist_M_df = pd.read_csv(file_in_SimLosses_nn_dist_M, sep="\t{1}",engine='python')
SimStayM_nn_dist_M_df = pd.read_csv(file_in_SimStayM_nn_dist_M, sep="\t{1}",engine='python')
SimStayU_nn_dist_M_df = pd.read_csv(file_in_SimStayU_nn_dist_M, sep="\t{1}",engine='python')


# sum columns of SimGains/Losses and convert to array
SimGains_nn_dist_M_array = np.array(SimGains_nn_dist_M_df.sum(axis=0).tolist())
SimLosses_nn_dist_M_array = np.array(SimLosses_nn_dist_M_df.sum(axis=0).tolist())
SimStayM_nn_dist_M_array = np.array(SimStayM_nn_dist_M_df.sum(axis=0).tolist())
SimStayU_nn_dist_M_array = np.array(SimStayU_nn_dist_M_df.sum(axis=0).tolist())


# NearestM
MovAv_SimGains_nn_dist_M = calc_moving_average(N_window, SimGains_nn_dist_M_array[2:])
MovAv_SimLosses_nn_dist_M = calc_moving_average(N_window, SimLosses_nn_dist_M_array[2:])
MovAv_SimInitiallyM_nn_dist_M = calc_moving_average(N_window, SimStayM_nn_dist_M_array[2:]+SimLosses_nn_dist_M_array[2:])
MovAv_SimInitiallyU_nn_dist_M = calc_moving_average(N_window, SimStayU_nn_dist_M_array[2:]+SimGains_nn_dist_M_array[2:])

MovAv_SimGainRate_nn_dist_M = MovAv_SimGains_nn_dist_M/MovAv_SimInitiallyU_nn_dist_M
MovAv_SimLossRate_nn_dist_M = MovAv_SimLosses_nn_dist_M/MovAv_SimInitiallyM_nn_dist_M
# rescale to per cell cycle
MovAv_SimGainRate_nn_dist_M = MovAv_SimGainRate_nn_dist_M/(n_cc*N_gen_burn_in)
MovAv_SimLossRate_nn_dist_M = MovAv_SimLossRate_nn_dist_M/(n_cc*N_gen_burn_in)
print(len(MovAv_SimGainRate_nn_dist_M),len(MovAv_SimLossRate_nn_dist_M))


# # Calc. gains and losses totaled from hook graphs
# # Sim 
# Total_Hook_Sim_Gains = np.sum(SimGains_nn_dist_M_array[2:])
# Total_Hook_Sim_Gains_per_rep = Total_Hook_Sim_Gains/N_reps
# Total_Hook_Sim_Losses = np.sum(SimLosses_nn_dist_M_array[2:])
# Total_Hook_Sim_Losses_per_rep = Total_Hook_Sim_Losses/N_reps
# Total_Hook_Sim_StayM = np.sum(SimStayM_nn_dist_M_array)
# Total_Hook_Sim_StayU = np.sum(SimStayU_nn_dist_M_array)
# Total_Hook_Sim_GainRate = Total_Hook_Sim_Gains/(Total_Hook_Sim_Gains+Total_Hook_Sim_StayU)
# Total_Hook_Sim_LossRate = Total_Hook_Sim_Losses/(Total_Hook_Sim_Losses+Total_Hook_Sim_StayM)
# Total_Hook_Sim_Gains_per_cc = Total_Hook_Sim_Gains_per_rep/(N_gen_burn_in*n_cc)
# Total_Hook_Sim_Losses_per_cc = Total_Hook_Sim_Losses_per_rep/(N_gen_burn_in*n_cc)
# Total_Hook_Sim_GainRate_per_cc = Total_Hook_Sim_GainRate/(N_gen_burn_in*n_cc)
# Total_Hook_Sim_LossRate_per_cc = Total_Hook_Sim_LossRate/(N_gen_burn_in*n_cc)
# # Write out to file
# file_name_output_Totals = os.path.join(GraphFolder,filename_start +'Hook_Sim_RatesTotals'+'.tsv')
# with open(file_name_output_Totals,'w',newline='') as output_file:
#     line_writer = csv.writer(output_file,delimiter='\t')
#     line_writer.writerow( ['Index', 'Total_Hook_Sim_Gains_per_cc','Total_Hook_Sim_Losses_per_cc',
#     'Total_Hook_Sim_GainRate_per_cc','Total_Hook_Sim_LossRate_per_cc','Gain_Loss_Ratio_per_cc'] )
#     line_writer.writerow( [sim_type,Total_Hook_Sim_Gains_per_cc,Total_Hook_Sim_Losses_per_cc,
#     Total_Hook_Sim_GainRate_per_cc,Total_Hook_Sim_LossRate_per_cc,Total_Hook_Sim_GainRate_per_cc/Total_Hook_Sim_LossRate_per_cc] )

# # Data_Lizzie
# Total_Hook_Data_Lizzie_Gains_per_cc = np.nansum(np.multiply(GainRate_data_df['Lizzie_Conservative_GainRate'].values[1:],GainRate_data_df['InitialStateU_Lizzie'].values[1:]))
# Total_Hook_Data_Lizzie_Losses_per_cc = np.nansum(np.multiply(LossRate_data_df['Lizzie_Conservative_LossRate'].values[1:],GainRate_data_df['InitialStateM_Lizzie'].values[1:]))
# Total_Hook_Data_Lizzie_InitialStateM = np.nansum(GainRate_data_df['InitialStateM_Lizzie'].values[1:])
# Total_Hook_Data_Lizzie_InitialStateU = np.nansum(GainRate_data_df['InitialStateU_Lizzie'].values[1:])
# Total_Hook_Data_Lizzie_GainRate_per_cc = Total_Hook_Data_Lizzie_Gains_per_cc/(Total_Hook_Data_Lizzie_InitialStateU)
# Total_Hook_Data_Lizzie_LossRate_per_cc = Total_Hook_Data_Lizzie_Losses_per_cc/(Total_Hook_Data_Lizzie_InitialStateM)
# # Write out to file
# file_name_output_Totals = os.path.join(GraphFolder,filename_start +'Hook_Data_Lizzie_RatesTotals'+'.tsv')
# with open(file_name_output_Totals,'w',newline='') as output_file:
#     line_writer = csv.writer(output_file,delimiter='\t')
#     line_writer.writerow( ['Index', 'Total_Hook_Data_Lizzie_Gains_per_cc','Total_Hook_Data_Lizzie_Losses_per_cc'
#     'Total_Hook_Data_Lizzie_GainRate_per_cc','Total_Hook_Data_Lizzie_LossRate_per_cc','Gain_Loss_Ratio_per_cc'] )
#     line_writer.writerow( ['Data_Lizzie',Total_Hook_Data_Lizzie_Gains_per_cc,Total_Hook_Data_Lizzie_Losses_per_cc,
#     Total_Hook_Data_Lizzie_GainRate_per_cc,Total_Hook_Data_Lizzie_LossRate_per_cc,Total_Hook_Data_Lizzie_GainRate_per_cc/Total_Hook_Data_Lizzie_LossRate_per_cc] )

# # Data_Schmitz
# Total_Hook_Data_Schmitz_Gains_per_cc = np.nansum(np.multiply(GainRate_data_df['Schmitz_GainRate'].values[1:],GainRate_data_df['InitialStateU_Schmitz'].values[1:]))
# Total_Hook_Data_Schmitz_Losses_per_cc = np.nansum(np.multiply(LossRate_data_df['Schmitz_LossRate'].values[1:],GainRate_data_df['InitialStateM_Schmitz'].values[1:]))
# Total_Hook_Data_Schmitz_InitialStateM = np.nansum(GainRate_data_df['InitialStateM_Schmitz'].values[1:])
# Total_Hook_Data_Schmitz_InitialStateU = np.nansum(GainRate_data_df['InitialStateU_Schmitz'].values[1:])
# Total_Hook_Data_Schmitz_GainRate_per_cc = Total_Hook_Data_Schmitz_Gains_per_cc/(Total_Hook_Data_Schmitz_InitialStateU)
# Total_Hook_Data_Schmitz_LossRate_per_cc = Total_Hook_Data_Schmitz_Losses_per_cc/(Total_Hook_Data_Schmitz_InitialStateM)
# # Write out to file
# file_name_output_Totals = os.path.join(GraphFolder,filename_start +'Hook_Data_Schmitz_RatesTotals'+'.tsv')
# with open(file_name_output_Totals,'w',newline='') as output_file:
#     line_writer = csv.writer(output_file,delimiter='\t')
#     line_writer.writerow( ['Index', 'Total_Hook_Data_Schmitz_Gains_per_cc','Total_Hook_Data_Schmitz_Losses_per_cc',
#     'Total_Hook_Data_Schmitz_GainRate_per_cc','Total_Hook_Data_Schmitz_LossRate_per_cc','Gain_Loss_Ratio_per_cc'] )
#     line_writer.writerow( ['Data_Schmitz',Total_Hook_Data_Schmitz_Gains_per_cc,Total_Hook_Data_Schmitz_Losses_per_cc,
#     Total_Hook_Data_Schmitz_GainRate_per_cc,Total_Hook_Data_Schmitz_LossRate_per_cc,Total_Hook_Data_Schmitz_GainRate_per_cc/Total_Hook_Data_Schmitz_LossRate_per_cc] )

# title_sting_single_time = "No. sims. per gene = %d      Cell cycle duration = %.1f       No. of genes simulated = %d        Sim. time = %d gens.\n\
# 'P' choice = %s      $\epsilon$ = variable      $\gamma_0$ = varialbe      $r_{div\\,\gamma} =$ %d      $r_{\gamma}$ = %d      $\lambda_{\gamma} =$ %.3f\n\
# $\delta$ = %.3e      $u_{M PL}^+ =$ varialbe      $r_{div.} =$ %d      $r_{plat.} =$ %d      $\lambda_{M PL} =$ %.3f      $u_{Scale\\,LR1}^+ =$ %.3f\n \
# $u_{M PL\\,LR1}^+ =$ varialbe      $r_{div\\,LR1} =$ %d      $r_{plat.\\,LR1} =$ %d      $\lambda_{M\\,In\\,PL\\,LR1} =$ %.3f      $\lambda_{M\\,Out\\,PL\\,LR1} =$ %.3f      $r_{LR} =$ %d\n \
# $m_u$ = %.2f     $c_u$ = %.2f     $m_{\epsilon}$ = %.2f     $c_{\epsilon}$ = %.2f    $m_{\gamma}$ = %.2f     $c_{\gamma}$ = %.2f     $\\rho_{cap}$ = 1./%.2f     $N_{smooth}$ = %d bp\n\
# Initial seed = %d      Sim. initial state: %s      %d $\leq N_{CG\\,min}$ < %d      %.3f $\leq \\rho_{CG\\,min}$ < %.3f\
# \nLocus type: %s      Annotation file: %s\n\
# Filter 1: %s      Filter 2: %s      ExcludeFilter 1: %s\nExperimental Initial State: %s"\
#              % (N_reps, rep_time, N_gene, N_gen_burn_in,  P_choice,
#                 r_div_gamma, r_gamma, lambda_gamma, delta, 
#                r_div, r_plat,lambda_coop, u_scale_val,
#                 r_div_LR1, r_plat_LR1, lambda_coopIn_LR1, lambda_coopOut_LR1,
#                 r_LR1, u_m_val,u_c_val,e_m_val,e_c_val,g_m_val,g_c_val,spacing_cap,N_window,
#                 initial_seed, initial_state_choice,N_CG_min,N_CG_max, N_CG_density_min,N_CG_density_max, locus_type, 
#                 file_in_annotation, file_in_anno_filt_1, file_in_anno_filt_2, file_in_exclude_filt_1,file_in_Parental_Col0_state)

# y-scale factor for plots
y_scale = 1.0e+4
y_scale_label = " ($\\times 10^{-4}$)"


# # Start Graph

# fig0, ax0 = plt.subplots(2,1,figsize=(12,12))

# for spine in ['left','right','top','bottom']:
#     ax0[0].spines[spine].set_color('k')
#     ax0[1].spines[spine].set_color('k')
# ax0[0].set_facecolor('white')
# ax0[1].set_facecolor('white')

# ax0[0].grid(False)
# ax0[1].grid(False)


# # Plots NearestM
# #Data
# # gains
# print(len(x_vals[2:]),len(MovAv_gains_All_Min),len(MovAv_gains_All_Max))
# ax0[0].fill_between(x_vals[2:],y_scale*MovAv_gains_All_Min, y_scale*MovAv_gains_All_Max, color='lightgrey', alpha = 0.8, label='All range $\pm 1\sigma$')
# ax0[0].fill_between(x_vals[2:],y_scale*MovAv_gains_Schmitz_Min, y_scale*MovAv_gains_Schmitz_Max, color='plum', alpha = 0.6, label='Schmitz range $\pm 1\sigma$')
# ax0[0].plot(x_vals[2:], y_scale*MovAv_gains_All_mean, linewidth=2, color='grey',alpha=0.8, label = 'All mean\nSmoothed n = %d bp' % N_window_gain_Jay)
# ax0[0].plot(x_vals[2:], y_scale*MovAv_gains_Schmitz_mean, linewidth=2, color='k', label = 'Schmitz mean\nSmoothed n = %d bp' % N_window_gain_Jay)
# ax0[0].plot(x_vals[2:], y_scale*MovAv_gains_Lizzie_conservative, linewidth=1, color='k',linestyle='--', label = 'Hollwey\nSmoothed n = %d bp' % N_window_gain_Lizzie)


# # losses 
# ax0[1].fill_between(x_vals[2:],y_scale*MovAv_losses_All_Min, y_scale*MovAv_losses_All_Max, color='lightgrey', alpha = 0.8, label='All range $\pm 1\sigma$')
# ax0[1].fill_between(x_vals[2:],y_scale*MovAv_losses_Schmitz_Min, y_scale*MovAv_losses_Schmitz_Max, color='lightseagreen', alpha = 0.6, label='Schmitz range $\pm 1\sigma$')
# ax0[1].plot(x_vals[2:], y_scale*MovAv_losses_All_mean, linewidth=2, color='grey',alpha=0.8, label = 'All mean\nSmoothed n = %d bp' % N_window_loss_Jay)
# ax0[1].plot(x_vals[2:], y_scale*MovAv_losses_Schmitz_mean, linewidth=2, color='k', label = 'Schmitz mean\nSmoothed n = %d bp' % N_window_loss_Jay)
# ax0[1].plot(x_vals[2:], y_scale*MovAv_losses_Lizzie_conservative, linewidth=1, color='k',linestyle='--', label = 'Hollwey\nSmoothed n = %d bp' % N_window_loss_Lizzie)

# #Model
# ax0[0].plot(x_vals[2:], y_scale*MovAv_SimGainRate_nn_dist_M,label='Simulated\nSmoothed n = %d bp' % N_window, linewidth=2, color='r')
# ax0[1].plot(x_vals[2:], y_scale*MovAv_SimLossRate_nn_dist_M,label='Simulated\nSmoothed n = %d bp' % N_window, linewidth=2, color='b')

# from mpl_toolkits.axes_grid1.inset_locator import inset_axes
# # gains inset
# axins_0 = inset_axes(ax0[0], width="85%", height="65%", borderpad=1)
# axins_0.tick_params(axis='both', which='major', labelsize=20)
# axins_0.set_xlim(64,489)
# axins_0.set_ylim(0,0.4)

# axins_0.fill_between(x_vals[2:],y_scale*MovAv_gains_All_Min, y_scale*MovAv_gains_All_Max, color='lightgrey', alpha = 0.8, label='All range $\pm 1\sigma$')
# axins_0.fill_between(x_vals[2:],y_scale*MovAv_gains_Schmitz_Min, y_scale*MovAv_gains_Schmitz_Max, color='plum', alpha = 0.6, label='Schmitz range $\pm 1\sigma$')
# axins_0.plot(x_vals[2:], y_scale*MovAv_gains_All_mean, linewidth=2, color='grey',alpha=0.8, label = 'All mean')
# axins_0.plot(x_vals[2:], y_scale*MovAv_gains_Schmitz_mean, linewidth=2, color='k', label = 'Schmitz mean')
# axins_0.plot(x_vals[2:], y_scale*MovAv_gains_Lizzie_conservative, linewidth=1, color='k',linestyle='--', label = 'Hollwey')
# axins_0.plot(x_vals[2:], y_scale*MovAv_SimGainRate_nn_dist_M, color='r',linestyle='-',linewidth=2,label='Simulated\nSmoothed n = %d bp' % N_window)

# # losses inset
# axins_1 = inset_axes(ax0[1], width="50%", height="50%", borderpad=1, loc=2)
# axins_1.tick_params(axis='both', which='major', labelsize=20)
# axins_1.yaxis.tick_right()
# axins_1.set_xlim(0,100)
# axins_1.set_ylim(0,0.75)

# axins_1.fill_between(x_vals[2:],y_scale*MovAv_losses_All_Min, y_scale*MovAv_losses_All_Max, color='lightgrey', alpha = 0.8, label='All range $\pm 1\sigma$')
# axins_1.fill_between(x_vals[2:],y_scale*MovAv_losses_Schmitz_Min, y_scale*MovAv_losses_Schmitz_Max, color='lightseagreen', alpha = 0.6, label='Schmitz range $\pm 1\sigma$')
# axins_1.plot(x_vals[2:], y_scale*MovAv_losses_All_mean, linewidth=2, color='grey',alpha=0.8, label = 'All mean')
# axins_1.plot(x_vals[2:], y_scale*MovAv_losses_Schmitz_mean, linewidth=2, color='k', label = 'Schmitz mean')
# axins_1.plot(x_vals[2:], y_scale*MovAv_losses_Lizzie_conservative, linewidth=1, color='k',linestyle='--', label = 'Hollwey')
# axins_1.plot(x_vals[2:], y_scale*MovAv_SimLossRate_nn_dist_M,label='Simulated\nSmoothed n = %d bp' % N_window, linewidth=2, color='b')





# ax0[0].set_ylim(0,3)
# ax0[1].set_ylim(0,3)

# ax0[0].set_xlabel("Distance to nearest M site (bps)", fontsize=18)
# ax0[0].set_ylabel("Gain rate per cell cycle"+y_scale_label, fontsize=20)
# # ,loc=(0.18,0.40)
# # ax0[0].legend(fontsize=18,ncol=4,facecolor='white', framealpha=0.6)
# ax0[0].tick_params(axis='both', which='major', labelsize=18)

# ax0[1].set_xlabel("Distance to nearest M site (bps)", fontsize=18)
# ax0[1].set_ylabel("Loss rate per cell cylce"+y_scale_label, fontsize=18)
# # ax0[1].legend(fontsize=18,ncol=4,loc=2,facecolor='white', framealpha=0.6)
# ax0[1].tick_params(axis='both', which='major', labelsize=18)

# # fig0.suptitle(title_sting_single_time, fontsize=20)

# #ax0[0].set_xlim(0,800)
# #ax0[1].set_xlim(0,250)
# ax0[0].set_xlim(0,500)
# ax0[1].set_xlim(0,500)

# # ## Try to format y-axis in scientific notation ##

# # class MathTextSciFormatter(mtick.Formatter):
# #     def __init__(self, fmt="%1.2e"):
# #         self.fmt = fmt
# #     def __call__(self, x, pos=None):
# #         s = self.fmt % x
# #         decimal_point = '.'
# #         positive_sign = '+'
# #         tup = s.split('e')
# #         significand = tup[0].rstrip(decimal_point)
# #         sign = tup[1][0].replace(positive_sign, '')
# #         exponent = tup[1][1:].lstrip('0')
# #         if exponent:
# #             exponent = '10^{%s%s}' % (sign, exponent)
# #         if significand and exponent:
# #             s =  r'%s{\times}%s' % (significand, exponent)
# #         else:
# #             s =  r'%s%s' % (significand, exponent)
# #         return "${}$".format(s)

# # # Format with 2 decimal places
# # ax0[0].yaxis.set_major_formatter(MathTextSciFormatter("%1.1e"))
# # ax0[1].yaxis.set_major_formatter(MathTextSciFormatter("%1.1e"))

# # ax0[0].yaxis.set_major_locator(plt.MaxNLocator(7))
# # ax0[1].yaxis.set_major_locator(plt.MaxNLocator(7))

# # #################################################


# #fig0.tight_layout()
# fig0.subplots_adjust(hspace=0.3)
# # fig0.subplots_adjust(top=0.8)

# fig0.savefig( os.path.join(GraphFolder,filename_start +'SimRates_NearestM_large.png'), bbox_inches = 'tight')

# # End Graph


# ###################


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



# from mpl_toolkits.axes_grid1.inset_locator import inset_axes
# # gains inset
# axins_0 = inset_axes(ax0[0], width="85%", height="65%", borderpad=1)
# axins_0.tick_params(axis='both', which='major', labelsize=14)
# axins_0.set_xlim(64,489)
# axins_0.set_ylim(0,0.4)

# axins_0.fill_between(x_vals[2:],y_scale*MovAv_gains_All_Min, y_scale*MovAv_gains_All_Max, color='lightgrey', alpha = 0.8, label='All range $\pm 1\sigma$')
# axins_0.fill_between(x_vals[2:],y_scale*MovAv_gains_Schmitz_Min, y_scale*MovAv_gains_Schmitz_Max, color='plum', alpha = 0.6, label='Schmitz range $\pm 1\sigma$')
# axins_0.plot(x_vals[2:], y_scale*MovAv_gains_All_mean, linewidth=2, color='grey',alpha=0.8, label = 'All mean')
# axins_0.plot(x_vals[2:], y_scale*MovAv_gains_Schmitz_mean, linewidth=2, color='k', label = 'Schmitz mean')
# axins_0.plot(x_vals[2:], y_scale*MovAv_gains_Lizzie_conservative, linewidth=1, color='k',linestyle='--', label = 'Hollwey')
# axins_0.plot(x_vals[2:], y_scale*MovAv_SimGainRate_nn_dist_M, color='r',linestyle='-',linewidth=2,label='Simulated\nSmoothed n = %d bp' % N_window)

# # losses inset
# axins_1 = inset_axes(ax0[1], width="50%", height="50%", borderpad=1, loc=2)
# axins_1.tick_params(axis='both', which='major', labelsize=14)
# axins_1.yaxis.tick_right()
# axins_1.set_xlim(0,100)
# axins_1.set_ylim(0,0.75)

# axins_1.fill_between(x_vals[2:],y_scale*MovAv_losses_All_Min, y_scale*MovAv_losses_All_Max, color='lightgrey', alpha = 0.8, label='All range $\pm 1\sigma$')
# axins_1.fill_between(x_vals[2:],y_scale*MovAv_losses_Schmitz_Min, y_scale*MovAv_losses_Schmitz_Max, color='lightseagreen', alpha = 0.6, label='Schmitz range $\pm 1\sigma$')
# axins_1.plot(x_vals[2:], y_scale*MovAv_losses_All_mean, linewidth=2, color='grey',alpha=0.8, label = 'All mean')
# axins_1.plot(x_vals[2:], y_scale*MovAv_losses_Schmitz_mean, linewidth=2, color='k', label = 'Schmitz mean')
# axins_1.plot(x_vals[2:], y_scale*MovAv_losses_Lizzie_conservative, linewidth=1, color='k',linestyle='--', label = 'Hollwey')
# axins_1.plot(x_vals[2:], y_scale*MovAv_SimLossRate_nn_dist_M,label='Simulated\nSmoothed n = %d bp' % N_window, linewidth=2, color='b')


# ax0[0].set_ylim(0,3)

# ax0[0].set_xlabel("Distance to nearest M site (bps)", fontsize=12)
# ax0[0].set_ylabel("Gain rate\n per cell cycle"+y_scale_label, fontsize=12)
# # ,loc=(0.18,0.40)
# #ax0[0].legend(fontsize=16,ncol=2,facecolor='white', framealpha=0.6)
# ax0[0].tick_params(axis='both', which='major', labelsize=12)

# ax0[1].set_ylim(0,3)
# ax0[1].set_xlabel("Distance to nearest M site (bps)", fontsize=12)
# ax0[1].set_ylabel("Loss rate\n per cell cylce"+y_scale_label, fontsize=12)
# #ax0[1].legend(fontsize=16,ncol=2,loc=2,facecolor='white', framealpha=0.6)
# ax0[1].tick_params(axis='both', which='major', labelsize=12)

# # fig0.suptitle(title_sting_single_time, fontsize=20)

# #ax0[0].set_xlim(0,800)
# #ax0[1].set_xlim(0,250)
# ax0[0].set_xlim(0,500)
# ax0[1].set_xlim(0,500)



# #ax0[0].yaxis.set_major_locator(plt.MaxNLocator(7))
# #ax0[1].yaxis.set_major_locator(plt.MaxNLocator(7))



# #fig0.tight_layout()
# fig0.subplots_adjust(hspace=0.3)
# # fig0.subplots_adjust(top=0.68)

# fig0.savefig( os.path.join(GraphFolder,filename_start +'SimRates_NearestM_Pair'+'.png'), bbox_inches = 'tight')

# # End Graph

# ##########################################

# # Start Graph

# fig0, ax0 = plt.subplots(1,1,figsize=(6,3))

# for spine in ['left','right','top','bottom']:
#     ax0.spines[spine].set_color('k')
#     ax0.spines[spine].set_linewidth(0.8)
# ax0.set_facecolor('white')

# #ax0[0].grid(False)

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



# from mpl_toolkits.axes_grid1.inset_locator import inset_axes
# # gains inset
# axins_0 = inset_axes(ax0, width="75%", height="65%", borderpad=1)
# axins_0.tick_params(axis='both', which='major', labelsize=12)
# axins_0.set_xlim(200,1400)
# axins_0.set_ylim(0,0.2)

# axins_0.fill_between(x_vals[2:],y_scale*MovAv_gains_All_Min, y_scale*MovAv_gains_All_Max, color='lightgrey', alpha = 0.8, label='All range $\pm 1\sigma$')
# axins_0.fill_between(x_vals[2:],y_scale*MovAv_gains_Schmitz_Min, y_scale*MovAv_gains_Schmitz_Max, color='plum', alpha = 0.6, label='Schmitz range $\pm 1\sigma$')
# axins_0.plot(x_vals[2:], y_scale*MovAv_gains_All_mean, linewidth=2, color='grey',alpha=0.8, label = 'All mean')
# axins_0.plot(x_vals[2:], y_scale*MovAv_gains_Schmitz_mean, linewidth=2, color='k', label = 'Schmitz mean')
# axins_0.plot(x_vals[2:], y_scale*MovAv_gains_Lizzie_conservative, linewidth=1, color='k',linestyle='--', label = 'Hollwey')
# axins_0.plot(x_vals[2:], y_scale*MovAv_SimGainRate_nn_dist_M, color='r',linestyle='-',linewidth=2,label='Simulated\nSmoothed n = %d bp' % N_window)



# ax0.set_ylim(0,3)

# ax0.set_xlabel("Distance to nearest M site (bps)", fontsize=12)
# ax0.set_ylabel("Gain rate\n per cell cycle"+y_scale_label, fontsize=12)
# # ,loc=(0.18,0.40)
# #ax0[0].legend(fontsize=12,ncol=2,facecolor='white', framealpha=0.6)
# ax0.tick_params(axis='both', which='major', labelsize=12)


# #fig0.suptitle(title_sting_single_time, fontsize=20)

# #ax0[0].set_xlim(0,800)
# #ax0[1].set_xlim(0,250)
# ax0.set_xlim(0,1400)



# #ax0[0].yaxis.set_major_locator(plt.MaxNLocator(7))
# #ax0[1].yaxis.set_major_locator(plt.MaxNLocator(7))



# #fig0.tight_layout()
# #fig0.subplots_adjust(hspace=0.3)
# #fig0.subplots_adjust(top=0.68)

# fig0.savefig( os.path.join(GraphFolder,filename_start +'SimGainsTail_NearestM'+'.png'), bbox_inches = 'tight')

# # End Graph


# ######################

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


# from mpl_toolkits.axes_grid1.inset_locator import inset_axes
# # gains inset
# axins_0 = inset_axes(ax0, width="85%", height="65%", borderpad=1)
# axins_0.tick_params(axis='both', which='major', labelsize=14)
# axins_0.set_xlim(64,489)
# axins_0.set_ylim(0,0.4)

# axins_0.fill_between(x_vals[2:],y_scale*MovAv_gains_All_Min, y_scale*MovAv_gains_All_Max, color='lightgrey', alpha = 0.8, label='All range $\pm 1\sigma$')
# axins_0.fill_between(x_vals[2:],y_scale*MovAv_gains_Schmitz_Min, y_scale*MovAv_gains_Schmitz_Max, color='plum', alpha = 0.6, label='Schmitz range $\pm 1\sigma$')
# axins_0.plot(x_vals[2:], y_scale*MovAv_gains_All_mean, linewidth=2, color='grey',alpha=0.8, label = 'All mean')
# axins_0.plot(x_vals[2:], y_scale*MovAv_gains_Schmitz_mean, linewidth=2, color='k', label = 'Schmitz mean')
# axins_0.plot(x_vals[2:], y_scale*MovAv_gains_Lizzie_conservative, linewidth=1, color='k',linestyle='--', label = 'Hollwey')
# axins_0.plot(x_vals[2:], y_scale*MovAv_SimGainRate_nn_dist_M, color='r',linestyle='-',linewidth=2,label='Simulated\nSmoothed n = %d bp' % N_window)


# ax0.set_ylim(0,3)

# ax0.set_xlabel("Distance to nearest M site (bps)", fontsize=12)
# ax0.set_ylabel("Gain rate\n per cell cycle"+y_scale_label, fontsize=12)
# # ,loc=(0.18,0.40)
# #ax0.legend(fontsize=12,ncol=2,facecolor='white', framealpha=0.6)
# ax0.tick_params(axis='both', which='major', labelsize=12)


# ax0.set_xlim(0,500)



# #ax0.yaxis.set_major_locator(plt.MaxNLocator(7))

# #fig0.tight_layout()
# #fig0.subplots_adjust(hspace=0.3)

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



# from mpl_toolkits.axes_grid1.inset_locator import inset_axes
# # losses inset
# axins_1 = inset_axes(ax0, width="50%", height="50%", borderpad=1, loc=2)
# axins_1.tick_params(axis='both', which='major', labelsize=14)
# axins_1.yaxis.tick_right()
# axins_1.set_xlim(0,100)
# axins_1.set_ylim(0,0.75)

# axins_1.fill_between(x_vals[2:],y_scale*MovAv_losses_All_Min, y_scale*MovAv_losses_All_Max, color='lightgrey', alpha = 0.8, label='All range $\pm 1\sigma$')
# axins_1.fill_between(x_vals[2:],y_scale*MovAv_losses_Schmitz_Min, y_scale*MovAv_losses_Schmitz_Max, color='lightseagreen', alpha = 0.6, label='Schmitz range $\pm 1\sigma$')
# axins_1.plot(x_vals[2:], y_scale*MovAv_losses_All_mean, linewidth=2, color='grey',alpha=0.8, label = 'All mean')
# axins_1.plot(x_vals[2:], y_scale*MovAv_losses_Schmitz_mean, linewidth=2, color='k', label = 'Schmitz mean')
# axins_1.plot(x_vals[2:], y_scale*MovAv_losses_Lizzie_conservative, linewidth=1, color='k',linestyle='--', label = 'Hollwey')
# axins_1.plot(x_vals[2:], y_scale*MovAv_SimLossRate_nn_dist_M,label='Simulated\nSmoothed n = %d bp' % N_window, linewidth=2, color='b')

# ax0.set_ylim(0,3)
# ax0.set_xlabel("Distance to nearest M site (bps)", fontsize=12)
# ax0.set_ylabel("Loss rate\n per cell cylce"+y_scale_label, fontsize=12)
# #ax0.legend(fontsize=12,ncol=2,loc=2,facecolor='white', framealpha=0.6)
# ax0.tick_params(axis='both', which='major', labelsize=12)

# ax0.set_xlim(0,500)

# #ax0.yaxis.set_major_locator(plt.MaxNLocator(7))


# #fig0.tight_layout()
# #fig0.subplots_adjust(hspace=0.3)
# # fig0.subplots_adjust(top=0.68)

# fig0.savefig( os.path.join(GraphFolder,filename_start +'SimRates_NearestM_Losses'+'.png'), bbox_inches = 'tight')

# # End Graph