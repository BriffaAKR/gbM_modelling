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
N_gen_total = params.N_gen_total
N_gen_burn_in = N_gen_total

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

coop_strength = params.coop_strength

# unlikely to edit
N_total_IDs = params.N_total_IDs
n_cc = params.n_cc


n_bin = 50
GraphFolder = "Graphs"

mutant_label = params.mutant_label

filename_start = params.filename_start
filename_ending = '.tsv'

# LocusProperties used for Col0 data - only need to load in once! 
file_in_LocusProperties = os.path.join('Output_files', 'MAF5_H2AZWT_EqbrmOutput_FinalFit_Pchoice_Excl_InitialState_100U_Ngen_1E5_LocusProperties.tsv')

LocusProperties_df = pd.read_csv(file_in_LocusProperties, sep="\t{1}",engine='python')
LocusProperties_df.set_index("gene_ID", inplace = True)
N_gene = len(LocusProperties_df)
IDs_list = LocusProperties_df.index.values.tolist()
print('N_gene_Expt',N_gene)

# Manually load in the Sims for 'mutant rates'
file_in_AllRepsMethFracs_mutant = os.path.join('Output_files', filename_start+'AllRepsMethFracs.tsv')
AllRepsMethFracs_mutant_df = pd.read_csv(file_in_AllRepsMethFracs_mutant, sep="\t{1}",engine='python')
AllRepsMethFracs_mutant_df.set_index("gene_ID", inplace = True)
print('N_gene_mutant',len(AllRepsMethFracs_mutant_df))


# define colours
col_data_Col0 = 'green'
col_Sim_mutant = 'k'


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
                  label='%s\nMean = %.0f %%' % (hist_Col0_label,np.nanmean(LocusProperties_df[['Col0_meth_frac']].values*values_scaler)) )

hist_Sims_mutant_label =  '$\it{%s}$ Sim.' % mutant_label[:-6]

temp_cols_list = []
for i_temp in range(N_reps):
    temp_cols_list.append('Dset1_meth_frac_'+str(i_temp))
hist_Sims_mutant, bins_Sim = np.histogram(AllRepsMethFracs_mutant_df[temp_cols_list].values*values_scaler, bins=n_bin, range=(x_start,x_end))
hist_Sims_mutant = hist_Sims_mutant/N_reps # normalise
# ax.axvline(np.nanmean(AllRepsMethFracs_mutant_df[temp_cols_list].values*values_scaler), linestyle='-', color=col_Sim_mutant)
ax.plot(x_hist, hist_Sims_mutant,  linewidth = 2, color=col_Sim_mutant,linestyle='-',
                  label='%s\nMean = %.0f %%' % (hist_Sims_mutant_label,np.nanmean(AllRepsMethFracs_mutant_df[temp_cols_list].values*values_scaler)) )

ax.set_title('$r^*$ = %.2f,  $\epsilon$ = %.1f$\\times10^{-4}$' % (coop_strength, e_c_val), y=1.05 )

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

fig.savefig( os.path.join(GraphFolder,filename_start +"Mfrac.png"), bbox_inches = 'tight')
fig.savefig( os.path.join(GraphFolder,filename_start +"Mfrac.pdf"), bbox_inches = 'tight')


plt.show()
# End Graph


