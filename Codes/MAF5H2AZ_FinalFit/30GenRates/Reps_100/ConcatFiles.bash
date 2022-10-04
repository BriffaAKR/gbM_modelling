#!/bin/bash -e
#SBATCH --mail-type=END,FAIL # notifications for job done & fail
#SBATCH --mail-user=amy.briffa@jic.ac.uk # send-to address
#SBATCH --mem 4000 # memory pool for all cores
# #SBATCH -p rg-mh
#SBATCH -e slurm_logs/slurm-%j.err 
#SBATCH -o slurm_logs/slurm-%j.out

# script to concatinate some files

# declare filename_start
filename_start="MAF5_H2AZWT_RatesOutput_FinalFit_Pchoice_Excl_InitialState_Expt_Nreps_100_"

folder_1="Output_batch_files/CodeProgress/"

#folder_2a="Output_batch_files/Correlations_Col0/NearestPairSeparationsMM_Col0/"
#folder_2b="Output_batch_files/Correlations_Col0/NearestPairSeparationsMU_Col0/"
#folder_2c="Output_batch_files/Correlations_Col0/NearestPairSeparationsUU_Col0/"
#folder_2d="Output_batch_files/Correlations_Col0/NearestPairSeparationsXX_Col0/"
#folder_2e="Output_batch_files/Correlations_Col0/PairSeparationsMM_Col0/"
#folder_2f="Output_batch_files/Correlations_Col0/PairSeparationsMU_Col0/"
#folder_2g="Output_batch_files/Correlations_Col0/PairSeparationsUU_Col0/"
#folder_2h="Output_batch_files/Correlations_Col0/PairSeparationsXX_Col0/"

#folder_3a="Output_batch_files/Correlations_D3/NearestPairSeparationsMM_D3/"
#folder_3b="Output_batch_files/Correlations_D3/NearestPairSeparationsMU_D3/"
#folder_3c="Output_batch_files/Correlations_D3/NearestPairSeparationsUU_D3/"
#folder_3d="Output_batch_files/Correlations_D3/NearestPairSeparationsXX_D3/"
#folder_3e="Output_batch_files/Correlations_D3/PairSeparationsMM_D3/"
#folder_3f="Output_batch_files/Correlations_D3/PairSeparationsMU_D3/"
#folder_3g="Output_batch_files/Correlations_D3/PairSeparationsUU_D3/"
#folder_3h="Output_batch_files/Correlations_D3/PairSeparationsXX_D3/"

#folder_4a="Output_batch_files/Correlations_D4/NearestPairSeparationsMM_D4/"
#folder_4b="Output_batch_files/Correlations_D4/NearestPairSeparationsMU_D4/"
#folder_4c="Output_batch_files/Correlations_D4/NearestPairSeparationsUU_D4/"
#folder_4d="Output_batch_files/Correlations_D4/NearestPairSeparationsXX_D4/"
#folder_4e="Output_batch_files/Correlations_D4/PairSeparationsMM_D4/"
#folder_4f="Output_batch_files/Correlations_D4/PairSeparationsMU_D4/"
#folder_4g="Output_batch_files/Correlations_D4/PairSeparationsUU_D4/"
#folder_4h="Output_batch_files/Correlations_D4/PairSeparationsXX_D4/"

#folder_5a="Output_batch_files/Correlations_D5/NearestPairSeparationsMM_D5/"
#folder_5b="Output_batch_files/Correlations_D5/NearestPairSeparationsMU_D5/"
#folder_5c="Output_batch_files/Correlations_D5/NearestPairSeparationsUU_D5/"
#folder_5d="Output_batch_files/Correlations_D5/NearestPairSeparationsXX_D5/"
#folder_5e="Output_batch_files/Correlations_D5/PairSeparationsMM_D5/"
#folder_5f="Output_batch_files/Correlations_D5/PairSeparationsMU_D5/"
#folder_5g="Output_batch_files/Correlations_D5/PairSeparationsUU_D5/"
#folder_5h="Output_batch_files/Correlations_D5/PairSeparationsXX_D5/"

#folder_6a="Output_batch_files/Correlations_Sim/NearestPairSeparationsMM_Sim/"
#folder_6b="Output_batch_files/Correlations_Sim/NearestPairSeparationsMU_Sim/"
#folder_6c="Output_batch_files/Correlations_Sim/NearestPairSeparationsUU_Sim/"
#folder_6d="Output_batch_files/Correlations_Sim/NearestPairSeparationsXX_Sim/"
#folder_6e="Output_batch_files/Correlations_Sim/PairSeparationsMM_Sim/"
#folder_6f="Output_batch_files/Correlations_Sim/PairSeparationsMU_Sim/"
#folder_6g="Output_batch_files/Correlations_Sim/PairSeparationsUU_Sim/"
#folder_6h="Output_batch_files/Correlations_Sim/PairSeparationsXX_Sim/"

#folder_7="Output_batch_files/LocusProperties/"
folder_8="Output_batch_files/Output_TEST/"
#folder_9="Output_batch_files/SimulatedStates/"
folder_10="Output_batch_files/Skipped_IDs/"
#folder_11="Output_batch_files/AllRepsMethFracs/"

folder_12a="Output_batch_files/ChangedSites/gains_nn_dist_M/"
folder_12b="Output_batch_files/ChangedSites/gains_nn_dist_U/"
folder_12c="Output_batch_files/ChangedSites/losses_nn_dist_M/"
folder_12d="Output_batch_files/ChangedSites/losses_nn_dist_U/"
folder_12e="Output_batch_files/ChangedSites/stayM_nn_dist_M/"
folder_12f="Output_batch_files/ChangedSites/stayM_nn_dist_U/"
folder_12g="Output_batch_files/ChangedSites/stayU_nn_dist_M/"
folder_12h="Output_batch_files/ChangedSites/stayU_nn_dist_U/"

filename_end_1="CodeProgress_batch_*"

#filename_end_2a="NearestPairSeparationsMM_Col0_batch_*"
#filename_end_2b="NearestPairSeparationsMU_Col0_batch_*"
#filename_end_2c="NearestPairSeparationsUU_Col0_batch_*"
#filename_end_2d="NearestPairSeparationsXX_Col0_batch_*"
#filename_end_2e="PairSeparationsMM_Col0_batch_*"
#filename_end_2f="PairSeparationsMU_Col0_batch_*"
#filename_end_2g="PairSeparationsUU_Col0_batch_*"
#filename_end_2h="PairSeparationsXX_Col0_batch_*"

#filename_end_3a="NearestPairSeparationsMM_D3_batch_*"
#filename_end_3b="NearestPairSeparationsMU_D3_batch_*"
#filename_end_3c="NearestPairSeparationsUU_D3_batch_*"
#filename_end_3d="NearestPairSeparationsXX_D3_batch_*"
#filename_end_3e="PairSeparationsMM_D3_batch_*"
#filename_end_3f="PairSeparationsMU_D3_batch_*"
#filename_end_3g="PairSeparationsUU_D3_batch_*"
#filename_end_3h="PairSeparationsXX_D3_batch_*"

#filename_end_4a="NearestPairSeparationsMM_D4_batch_*"
#filename_end_4b="NearestPairSeparationsMU_D4_batch_*"
#filename_end_4c="NearestPairSeparationsUU_D4_batch_*"
#filename_end_4d="NearestPairSeparationsXX_D4_batch_*"
#filename_end_4e="PairSeparationsMM_D4_batch_*"
#filename_end_4f="PairSeparationsMU_D4_batch_*"
#filename_end_4g="PairSeparationsUU_D4_batch_*"
#filename_end_4h="PairSeparationsXX_D4_batch_*"

#filename_end_5a="NearestPairSeparationsMM_D5_batch_*"
#filename_end_5b="NearestPairSeparationsMU_D5_batch_*"
#filename_end_5c="NearestPairSeparationsUU_D5_batch_*"
#filename_end_5d="NearestPairSeparationsXX_D5_batch_*"
#filename_end_5e="PairSeparationsMM_D5_batch_*"
#filename_end_5f="PairSeparationsMU_D5_batch_*"
#filename_end_5g="PairSeparationsUU_D5_batch_*"
#filename_end_5h="PairSeparationsXX_D5_batch_*"

#filename_end_6a="NearestPairSeparationsMM_Sim_batch_*"
#filename_end_6b="NearestPairSeparationsMU_Sim_batch_*"
#filename_end_6c="NearestPairSeparationsUU_Sim_batch_*"
#filename_end_6d="NearestPairSeparationsXX_Sim_batch_*"
#filename_end_6e="PairSeparationsMM_Sim_batch_*"
#filename_end_6f="PairSeparationsMU_Sim_batch_*"
#filename_end_6g="PairSeparationsUU_Sim_batch_*"
#filename_end_6h="PairSeparationsXX_Sim_batch_*"

#filename_end_7="LocusProperties_batch_*"
filename_end_8="TEST_batch_*"
filename_end_10="Skipped_IDs_batch_*"
#filename_end_11="AllRepsMethFracs_batch_*"

filename_end_12a="SimGains_nn_dist_M_batch_*"
filename_end_12b="SimGains_nn_dist_U_batch_*"
filename_end_12c="SimLosses_nn_dist_M_batch_*"
filename_end_12d="SimLosses_nn_dist_U_batch_*"
filename_end_12e="SimStayM_nn_dist_M_batch_*"
filename_end_12f="SimStayM_nn_dist_U_batch_*"
filename_end_12g="SimStayU_nn_dist_M_batch_*"
filename_end_12h="SimStayU_nn_dist_U_batch_*"


output_folder="Output_files/"

output_filename_end_1="CodeProgress.tsv"

#output_filename_end_2a="NearestPairSeparationsMM_Col0.tsv"
#output_filename_end_2b="NearestPairSeparationsMU_Col0.tsv"
#output_filename_end_2c="NearestPairSeparationsUU_Col0.tsv"
#output_filename_end_2d="NearestPairSeparationsXX_Col0.tsv"
#output_filename_end_2e="PairSeparationsMM_Col0.tsv"
#output_filename_end_2f="PairSeparationsMU_Col0.tsv"
#output_filename_end_2g="PairSeparationsUU_Col0.tsv"
#output_filename_end_2h="PairSeparationsXX_Col0.tsv"

#output_filename_end_3a="NearestPairSeparationsMM_D3.tsv"
#output_filename_end_3b="NearestPairSeparationsMU_D3.tsv"
#output_filename_end_3c="NearestPairSeparationsUU_D3.tsv"
#output_filename_end_3d="NearestPairSeparationsXX_D3.tsv"
#output_filename_end_3e="PairSeparationsMM_D3.tsv"
#output_filename_end_3f="PairSeparationsMU_D3.tsv"
#output_filename_end_3g="PairSeparationsUU_D3.tsv"
#output_filename_end_3h="PairSeparationsXX_D3.tsv"

#output_filename_end_4a="NearestPairSeparationsMM_D4.tsv"
#output_filename_end_4b="NearestPairSeparationsMU_D4.tsv"
#output_filename_end_4c="NearestPairSeparationsUU_D4.tsv"
#output_filename_end_4d="NearestPairSeparationsXX_D4.tsv"
#output_filename_end_4e="PairSeparationsMM_D4.tsv"
#output_filename_end_4f="PairSeparationsMU_D4.tsv"
#output_filename_end_4g="PairSeparationsUU_D4.tsv"
#output_filename_end_4h="PairSeparationsXX_D4.tsv"

#output_filename_end_5a="NearestPairSeparationsMM_D5.tsv"
#output_filename_end_5b="NearestPairSeparationsMU_D5.tsv"
#output_filename_end_5c="NearestPairSeparationsUU_D5.tsv"
#output_filename_end_5d="NearestPairSeparationsXX_D5.tsv"
#output_filename_end_5e="PairSeparationsMM_D5.tsv"
#output_filename_end_5f="PairSeparationsMU_D5.tsv"
#output_filename_end_5g="PairSeparationsUU_D5.tsv"
#output_filename_end_5h="PairSeparationsXX_D5.tsv"

#output_filename_end_6a="NearestPairSeparationsMM_Sim.tsv"
#output_filename_end_6b="NearestPairSeparationsMU_Sim.tsv"
#output_filename_end_6c="NearestPairSeparationsUU_Sim.tsv"
#output_filename_end_6d="NearestPairSeparationsXX_Sim.tsv"
#output_filename_end_6e="PairSeparationsMM_Sim.tsv"
#output_filename_end_6f="PairSeparationsMU_Sim.tsv"
#output_filename_end_6g="PairSeparationsUU_Sim.tsv"
#output_filename_end_6h="PairSeparationsXX_Sim.tsv"

#output_filename_end_7="LocusProperties.tsv"
output_filename_end_8="TEST.tsv"
output_filename_end_10="Skipped_IDs.tsv"
#output_filename_end_11="AllRepsMethFracs.tsv"

output_filename_end_12a="SimGains_nn_dist_M.tsv"
output_filename_end_12b="SimGains_nn_dist_U.tsv"
output_filename_end_12c="SimLosses_nn_dist_M.tsv"
output_filename_end_12d="SimLosses_nn_dist_U.tsv"
output_filename_end_12e="SimStayM_nn_dist_M.tsv"
output_filename_end_12f="SimStayM_nn_dist_U.tsv"
output_filename_end_12g="SimStayU_nn_dist_M.tsv"
output_filename_end_12h="SimStayU_nn_dist_U.tsv"


input_filename_1="$folder_1$filename_start$filename_end_1"

#input_filename_2a="$folder_2a$filename_start$filename_end_2a"
#input_filename_2b="$folder_2b$filename_start$filename_end_2b"
#input_filename_2c="$folder_2c$filename_start$filename_end_2c"
#input_filename_2d="$folder_2d$filename_start$filename_end_2d"
#input_filename_2e="$folder_2e$filename_start$filename_end_2e"
#input_filename_2f="$folder_2f$filename_start$filename_end_2f"
#input_filename_2g="$folder_2g$filename_start$filename_end_2g"
#input_filename_2h="$folder_2h$filename_start$filename_end_2h"

#input_filename_3a="$folder_3a$filename_start$filename_end_3a"
#input_filename_3b="$folder_3b$filename_start$filename_end_3b"
#input_filename_3c="$folder_3c$filename_start$filename_end_3c"
#input_filename_3d="$folder_3d$filename_start$filename_end_3d"
#input_filename_3e="$folder_3e$filename_start$filename_end_3e"
#input_filename_3f="$folder_3f$filename_start$filename_end_3f"
#input_filename_3g="$folder_3g$filename_start$filename_end_3g"
#input_filename_3h="$folder_3h$filename_start$filename_end_3h"

#input_filename_4a="$folder_4a$filename_start$filename_end_4a"
#input_filename_4b="$folder_4b$filename_start$filename_end_4b"
#input_filename_4c="$folder_4c$filename_start$filename_end_4c"
#input_filename_4d="$folder_4d$filename_start$filename_end_4d"
#input_filename_4e="$folder_4e$filename_start$filename_end_4e"
#input_filename_4f="$folder_4f$filename_start$filename_end_4f"
#input_filename_4g="$folder_4g$filename_start$filename_end_4g"
#input_filename_4h="$folder_4h$filename_start$filename_end_4h"

#input_filename_5a="$folder_5a$filename_start$filename_end_5a"
#input_filename_5b="$folder_5b$filename_start$filename_end_5b"
#input_filename_5c="$folder_5c$filename_start$filename_end_5c"
#input_filename_5d="$folder_5d$filename_start$filename_end_5d"
#input_filename_5e="$folder_5e$filename_start$filename_end_5e"
#input_filename_5f="$folder_5f$filename_start$filename_end_5f"
#input_filename_5g="$folder_5g$filename_start$filename_end_5g"
#input_filename_5h="$folder_5h$filename_start$filename_end_5h"

#input_filename_6a="$folder_6a$filename_start$filename_end_6a"
#input_filename_6b="$folder_6b$filename_start$filename_end_6b"
#input_filename_6c="$folder_6c$filename_start$filename_end_6c"
#input_filename_6d="$folder_6d$filename_start$filename_end_6d"
#input_filename_6e="$folder_6e$filename_start$filename_end_6e"
#input_filename_6f="$folder_6f$filename_start$filename_end_6f"
#input_filename_6g="$folder_6g$filename_start$filename_end_6g"
#input_filename_6h="$folder_6h$filename_start$filename_end_6h"

#input_filename_7="$folder_7$filename_start$filename_end_7"
input_filename_8="$folder_8$filename_start$filename_end_8"
input_filename_10="$folder_10$filename_start$filename_end_10"
#input_filename_11="$folder_11$filename_start$filename_end_11"

input_filename_12a="$folder_12a$filename_start$filename_end_12a"
input_filename_12b="$folder_12b$filename_start$filename_end_12b"
input_filename_12c="$folder_12c$filename_start$filename_end_12c"
input_filename_12d="$folder_12d$filename_start$filename_end_12d"
input_filename_12e="$folder_12e$filename_start$filename_end_12e"
input_filename_12f="$folder_12f$filename_start$filename_end_12f"
input_filename_12g="$folder_12g$filename_start$filename_end_12g"
input_filename_12h="$folder_12h$filename_start$filename_end_12h"

output_filename_1="$output_folder$filename_start$output_filename_end_1"

#output_filename_2a="$output_folder$filename_start$output_filename_end_2a"
#output_filename_2b="$output_folder$filename_start$output_filename_end_2b"
#output_filename_2c="$output_folder$filename_start$output_filename_end_2c"
#output_filename_2d="$output_folder$filename_start$output_filename_end_2d"
#output_filename_2e="$output_folder$filename_start$output_filename_end_2e"
#output_filename_2f="$output_folder$filename_start$output_filename_end_2f"
#output_filename_2g="$output_folder$filename_start$output_filename_end_2g"
#output_filename_2h="$output_folder$filename_start$output_filename_end_2h"

#output_filename_3a="$output_folder$filename_start$output_filename_end_3a"
#output_filename_3b="$output_folder$filename_start$output_filename_end_3b"
#output_filename_3c="$output_folder$filename_start$output_filename_end_3c"
#output_filename_3d="$output_folder$filename_start$output_filename_end_3d"
#output_filename_3e="$output_folder$filename_start$output_filename_end_3e"
#output_filename_3f="$output_folder$filename_start$output_filename_end_3f"
#output_filename_3g="$output_folder$filename_start$output_filename_end_3g"
#output_filename_3h="$output_folder$filename_start$output_filename_end_3h"

#output_filename_4a="$output_folder$filename_start$output_filename_end_4a"
#output_filename_4b="$output_folder$filename_start$output_filename_end_4b"
#output_filename_4c="$output_folder$filename_start$output_filename_end_4c"
#output_filename_4d="$output_folder$filename_start$output_filename_end_4d"
#output_filename_4e="$output_folder$filename_start$output_filename_end_4e"
#output_filename_4f="$output_folder$filename_start$output_filename_end_4f"
#output_filename_4g="$output_folder$filename_start$output_filename_end_4g"
#output_filename_4h="$output_folder$filename_start$output_filename_end_4h"

#output_filename_5a="$output_folder$filename_start$output_filename_end_5a"
#output_filename_5b="$output_folder$filename_start$output_filename_end_5b"
#output_filename_5c="$output_folder$filename_start$output_filename_end_5c"
#output_filename_5d="$output_folder$filename_start$output_filename_end_5d"
#output_filename_5e="$output_folder$filename_start$output_filename_end_5e"
#output_filename_5f="$output_folder$filename_start$output_filename_end_5f"
#output_filename_5g="$output_folder$filename_start$output_filename_end_5g"
#output_filename_5h="$output_folder$filename_start$output_filename_end_5h"

#output_filename_6a="$output_folder$filename_start$output_filename_end_6a"
#output_filename_6b="$output_folder$filename_start$output_filename_end_6b"
#output_filename_6c="$output_folder$filename_start$output_filename_end_6c"
#output_filename_6d="$output_folder$filename_start$output_filename_end_6d"
#output_filename_6e="$output_folder$filename_start$output_filename_end_6e"
#output_filename_6f="$output_folder$filename_start$output_filename_end_6f"
#output_filename_6g="$output_folder$filename_start$output_filename_end_6g"
#output_filename_6h="$output_folder$filename_start$output_filename_end_6h"

#output_filename_7="$output_folder$filename_start$output_filename_end_7"
output_filename_8="$output_folder$filename_start$output_filename_end_8"
output_filename_10="$output_folder$filename_start$output_filename_end_10"
#output_filename_11="$output_folder$filename_start$output_filename_end_11"

output_filename_12a="$output_folder$filename_start$output_filename_end_12a"
output_filename_12b="$output_folder$filename_start$output_filename_end_12b"
output_filename_12c="$output_folder$filename_start$output_filename_end_12c"
output_filename_12d="$output_folder$filename_start$output_filename_end_12d"
output_filename_12e="$output_folder$filename_start$output_filename_end_12e"
output_filename_12f="$output_folder$filename_start$output_filename_end_12f"
output_filename_12g="$output_folder$filename_start$output_filename_end_12g"
output_filename_12h="$output_folder$filename_start$output_filename_end_12h"


cat $input_filename_1 > $output_filename_1

#cat $input_filename_2a > $output_filename_2a
#cat $input_filename_2b > $output_filename_2b
#cat $input_filename_2c > $output_filename_2c
#cat $input_filename_2d > $output_filename_2d
#cat $input_filename_2e > $output_filename_2e
#cat $input_filename_2f > $output_filename_2f
#cat $input_filename_2g > $output_filename_2g
#cat $input_filename_2h > $output_filename_2g

#cat $input_filename_3a > $output_filename_3a
#cat $input_filename_3b > $output_filename_3b
#cat $input_filename_3c > $output_filename_3c
#cat $input_filename_3d > $output_filename_3d
#cat $input_filename_3e > $output_filename_3e
#cat $input_filename_3f > $output_filename_3f
#cat $input_filename_3g > $output_filename_3g
#cat $input_filename_3h > $output_filename_3g

#cat $input_filename_4a > $output_filename_4a
#cat $input_filename_4b > $output_filename_4b
#cat $input_filename_4c > $output_filename_4c
#cat $input_filename_4d > $output_filename_4d
#cat $input_filename_4e > $output_filename_4e
#cat $input_filename_4f > $output_filename_4f
#cat $input_filename_4g > $output_filename_4g
#cat $input_filename_4h > $output_filename_4g

#cat $input_filename_5a > $output_filename_5a
#cat $input_filename_5b > $output_filename_5b
#cat $input_filename_5c > $output_filename_5c
#cat $input_filename_5d > $output_filename_5d
#cat $input_filename_5e > $output_filename_5e
#cat $input_filename_5f > $output_filename_5f
#cat $input_filename_5g > $output_filename_5g
#cat $input_filename_5h > $output_filename_5g

#cat $input_filename_6a > $output_filename_6a
#cat $input_filename_6b > $output_filename_6b
#cat $input_filename_6c > $output_filename_6c
#cat $input_filename_6d > $output_filename_6d
#cat $input_filename_6e > $output_filename_6e
#cat $input_filename_6f > $output_filename_6f
#cat $input_filename_6g > $output_filename_6g
#cat $input_filename_6h > $output_filename_6g

#cat $input_filename_7 > $output_filename_7
cat $input_filename_8 > $output_filename_8

#for i in {0..29}
#do
#    filename_end_9="SimState_$i_batch_*"
#    output_filename_end_9="SimState_$i.tsv"
#    input_filename_9="$folder_9$filename_start$filename_end_9"
#    output_filename_9="$output_folder$filename_start$output_filename_end_9"
#    cat $input_filename_9 > $output_filename_9
#done

cat $input_filename_10 > $output_filename_10
#cat $input_filename_11 > $output_filename_11


cat $input_filename_12a > $output_filename_12a
cat $input_filename_12b > $output_filename_12b
cat $input_filename_12c > $output_filename_12c
cat $input_filename_12d > $output_filename_12d
cat $input_filename_12e > $output_filename_12e
cat $input_filename_12f > $output_filename_12f
cat $input_filename_12g > $output_filename_12g
cat $input_filename_12h > $output_filename_12h