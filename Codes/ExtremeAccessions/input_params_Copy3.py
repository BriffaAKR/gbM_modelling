import os

# define annoation files
file_in_annotation = "MAF5_H2AZWT_End1p2_Mean1p2.tsv"
N_total_IDs = 10900 # This must be more than the total number of IDs in the annotation!!!

# 'UseU' 'Excl'
P_choice = "Excl"
# "100U"  "100M"  'Expt'
initial_state_choice = "100U"

N_gen_burn_in = 100*1000
N_gen_code = '1E5'

file_in_anno_filt_1 = 'non_te_genes_0p01.txt' # USE 'NONE' if no filtering required
file_in_anno_filt_2 = 'NONE' # USE 'NONE' if no filtering required
file_in_exclude_filt_1 = 'SegmentOverlap_gt20perc_smallestgone.tsv' # IDs to remove from annotation

# build up input path string
InputFiles_folder = 'Sim_Input_Files_Final'
path_string = os.path.join( os.getcwd(), '..')
path_string = os.path.join( path_string, '..')
#path_string = os.path.join( path_string, '..')
#path_string = os.path.join( path_string, '..')
path_string = os.path.join( path_string, InputFiles_folder)

####

# Sim options
coop_strength_int = 0.0
coop_strength_dec = 0.98
coop_strength = coop_strength_int + coop_strength_dec
delta_int = 4.0
delta_dec = 0.0

filename_params_code = 'output_Sim_Coop'+str(int(coop_strength_int))+'p'+str(coop_strength_dec)[2:]+'_delta'+str(int(delta_int))+'p'+str(delta_dec)[2:]+'_'


# RLFilt options
# input_accessions_folder = 'm1001_states_Relict'
# file_in_accessions_codes = "Accessions_RelictsFilt_MinCov_0p7.tsv" # list of accessions to load in
# filename_params_code = 'output_RLFilt_'

# NSFilt options
# input_accessions_folder = 'm1001_states_NS'
# file_in_accessions_codes = "Accessions_NorthenSwedishFilt_MinCov_0p7.tsv" # list of accessions to load in
# filename_params_code = 'output_NSFilt_'


chosen_accession_file_end = '3.tsv'
chosen_accession_file_start = 'm1001_samples_meth_status_mCG_density_'
path_string_single_accession = os.path.join( path_string, 'm1001_states_Extreme')


#for single accession code
# SRX248646 Cvi_0
# SRX2190740 UKID116
# SRX2190758 Can-0
# SRX445897 Dör-10

SRXcode_Cvi0 = '0_SRX248646_'
SRXcode_UKID116 = '0_SRX2190740_'
SRXcode_Can0 = '0_SRX2190758_'
SRXcode_Dor10 = '11_SRX445897_'


# chosen_accession_file = os.path.join( path_string_single_accession, chosen_accession_file_start+SRXcode_Cvi0+chosen_accession_file_end)
# filename_params_code = 'output_Cvi0_'

# chosen_accession_file = os.path.join( path_string_single_accession, chosen_accession_file_start+SRXcode_UKID116+chosen_accession_file_end)
# filename_params_code = 'output_UKID116_'

# chosen_accession_file = os.path.join( path_string_single_accession, chosen_accession_file_start+SRXcode_Can0+chosen_accession_file_end)
# filename_params_code = 'output_Can0_'

# chosen_accession_file = os.path.join( path_string_single_accession, chosen_accession_file_start+SRXcode_Dor10+chosen_accession_file_end)
# filename_params_code = 'output_Dor10_'



locus_type = 'ExclSegmentOverlap_gt20perc_SmallestGone'
TE_filt = 'TEfilt0p01'

filename_start = file_in_annotation[:-4] + '_' 'ExtremeAccns_' + filename_params_code


N_corrns_bins = 1000 + 1 # No. of pb bins to calculate gaps correlations for. 


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

u_m_val = -33.0
u_c_val = coop_strength*1.95 # 1.34 - (u_m_val/53.)
u_cap_val = (u_m_val/spacing_cap) + u_c_val

e_m_val = 0.0 # 0.0
e_c_val = 2.0 # e_m_val0*(1./83.) + e_c_val0
e_cap_val = 1.0E+50

g_m_val = -15.0
g_c_val = coop_strength*0.90 # 0.59 - (g_m_val/53.)
g_cap_val = (g_m_val/spacing_cap) + g_c_val

delta = (delta_int + delta_dec)*1.0e-6

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


