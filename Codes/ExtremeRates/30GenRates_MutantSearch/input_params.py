# define annoation files
file_in_annotation = "MAF5_H2AZWT_End1p2_Mean1p2.tsv"
N_total_IDs = 10900 # This must be more than the total number of IDs in the annotation!!!

# 'UseU' 'Excl'
P_choice = "Excl"
# "100U"  "100M"  'Expt'
initial_state_choice = "Expt"

N_gen_burn_in = 30
#N_gen_code = '1E5'

N_reps = 100 # Number of simulation reps

InputFiles_folder = 'Sim_Input_Files_Final'
sim_type = 'ExtremeRates'

# Sim options
coop_strength_int = 1.0
coop_strength_dec = 0.5
coop_strength = coop_strength_int + coop_strength_dec

epsilon_int = 4.0
epsilon_dec = 0.1

delta_int = 4.0
delta_dec = 0.0



file_in_anno_filt_1 = 'non_te_genes_0p01.txt' # USE 'NONE' if no filtering required
file_in_anno_filt_2 = 'NONE' # USE 'NONE' if no filtering required
file_in_exclude_filt_1 = 'SegmentOverlap_gt20perc_smallestgone.tsv' # IDs to remove from annotation

file_in_Parental_Col0_state = 'SeqRun9All_CG_2021-01-18_WTconsensus_meth_status_allgenes_withID.tsv'


filename_params_code = 'InitialState_'+initial_state_choice+'_'+'Ngen_'+str(N_gen_burn_in)+'_'+'Nreps_'+str(N_reps)+'_'+'output_Sim_Coop'+str(int(coop_strength_int))+'p'+str(coop_strength_dec)[2:]+'_delta'+str(int(delta_int))+'p'+str(delta_dec)[2:]+'_'+'_epsilon'+str(int(epsilon_int))+'p'+str(epsilon_dec)[2:]+'_'
locus_type = 'ExclSegmentOverlap_gt20perc_SmallestGone'
TE_filt = 'TEfilt0p01'

#filename_start = file_in_annotation[:-4] + '_' + locus_type + '_' + TE_filt + '_Col0' + filename_params_code
filename_start = "MAF5_H2AZWT_RatesOutput_"+filename_params_code

# files for graphs
file_in_GainRateData = 'MAF5Perc_TransH2az1p2_MeanH2az1p2_ExclTE0p01_ExclOlapGT20_GainRateRanges_ExclOutlier.tsv'
file_in_LossRateData = 'MAF5Perc_TransH2az1p2_MeanH2az1p2_ExclTE0p01_ExclOlapGT20_LossRateRanges_ExclOutlier.tsv'

file_in_GainRateData_NearestU = 'MAF5Perc_TransH2az1p2_MeanH2az1p2_ExclTE0p01_ExclOlapGT20_GainRateRanges_ExclOutlier_NearestU.tsv'
file_in_LossRateData_NearestU = 'MAF5Perc_TransH2az1p2_MeanH2az1p2_ExclTE0p01_ExclOlapGT20_LossRateRanges_ExclOutlier_NearestU.tsv'


N_corrns_bins = 2000 + 1 # No. of pb bins to calculate gaps correlations for. 




initial_seed = 0



rep_time = 1




off_rate = 1.0e-50

N_CG_min = 5
N_CG_max = 5000
#N_CG_max = 500
N_CG_density_min = 0.0
N_CG_density_max = 1.0
#N_CG_density_max = 0.1

# define linear relationships:
u_scale_val = 0.1
spacing_cap = 26.

u_m_val = -33.0
u_c_val = coop_strength*1.95 # 1.34 - (u_m_val/53.)
u_cap_val = (u_m_val/spacing_cap) + u_c_val

e_m_val = 0.0 # 0.0
# e_c_val = 2.0 # e_m_val0*(1./83.) + e_c_val0
e_cap_val = 1.0E+50

g_m_val = -15.0
g_c_val = coop_strength*0.90 # 0.59 - (g_m_val/53.)
g_cap_val = (g_m_val/spacing_cap) + g_c_val

delta = (delta_int + delta_dec)*1.0e-6
e_c_val = epsilon_int + epsilon_dec

lambda_gamma = 0.6
r_div_gamma = 11
r_gamma = 13

# SR (short-range) parameters
lambda_coop = 1.6 
r_div = 15
r_plat = 18 

# LR (long-range) parameters
lambda_coopIn_LR1 = lambda_coop
lambda_coopOut_LR1 = lambda_coop
r_div_LR1 = 15
r_plat_LR1 = 18

r_LR1 = 167


n_cc = 34


