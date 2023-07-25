# define annoation files
file_in_annotation = "Whole_GBM_Genes.tsv"
N_total_IDs = 14600 # This must be more than the total number of IDs in the annotation!!!

# 'UseU' 'Excl'
P_choice = "Excl"
# "100U"  "100M"  'Expt'
initial_state_choice = "100U"

N_gen_burn_in = 100*1000
N_gen_code = '1E5'

InputFiles_folder = 'Sim_Input_Files_Final'

file_in_anno_filt_1 = 'non_te_genes_0p01.txt' # USE 'NONE' if no filtering required
file_in_anno_filt_2 = 'NONE' # USE 'NONE' if no filtering required
file_in_exclude_filt_1 = 'SegmentOverlap_gt20perc_smallestgone.tsv' # IDs to remove from annotation

file_in_Parental_Col0_state = 'SeqRun9All_CG_2021-01-18_WTconsensus_meth_status_allgenes_withID.tsv'


#filename_params_code = '_output_Try1_Um33_Gm15_Floor26_Cu1p95_Cg1p00_'+'Pchoice_'+P_choice+'_'
filename_params_code = 'Pchoice_'+P_choice+'_'+'InitialState_'+initial_state_choice+'_'+'Ngen_'+N_gen_code+'_'

locus_type = 'ExclSegmentOverlap_gt20perc_SmallestGone'
TE_filt = 'TEfilt0p01'

#filename_start = file_in_annotation[:-4] + '_' + locus_type + '_' + TE_filt + '_Col0' + filename_params_code
filename_start = "WholeGBMgenes_EqbrmOutput_RatesOnlyFit_"+filename_params_code

# files for graphs
file_in_GainRateData = 'WholeGbmGenes_GainRateRanges_ExclOutlier.tsv'
file_in_LossRateData = 'WholeGbmGenes_LossRateRanges_ExclOutlier.tsv'




N_corrns_bins = 2000 + 1 # No. of pb bins to calculate gaps correlations for. 


N_reps = 30 # Number of simulation reps

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

u_m_val = 0.0 # -33.0
u_c_val = 1.3 # 1.95 # 1.34 - (u_m_val/53.)
u_cap_val = (u_m_val/spacing_cap) + u_c_val

e_m_val = 0.0 # 0.0
e_c_val = 2.0 # e_m_val0*(1./83.) + e_c_val0
e_cap_val = 1.0E+50

g_m_val = 0.0 # -15.0
g_c_val = 0.5 # 0.90 # 0.59 - (g_m_val/53.)
g_cap_val = (g_m_val/spacing_cap) + g_c_val

delta =  5.0e-7

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


