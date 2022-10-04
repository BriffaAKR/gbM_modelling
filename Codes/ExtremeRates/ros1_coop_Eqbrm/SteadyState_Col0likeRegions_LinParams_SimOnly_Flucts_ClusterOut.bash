#!/bin/bash -e
#SBATCH --mail-type=END,FAIL # notifications for job done & fail
#SBATCH --mail-user=amy.briffa@jic.ac.uk # send-to address
#SBATCH --mem 4000 # memory pool for all cores
##SBATCH -p rg-mh
#SBATCH -p jic-medium
##SBATCH -t 2-00:00
#SBATCH -e slurm_logs/slurm-%j.err 
#SBATCH -o slurm_logs/slurm-%j.out
##SBATCH --array=0-2 # end of array (N).  2nd argument of python code should be: N . start of array must by 0
#SBATCH --array=0-1000 # end of array (N).  2nd argument of python code should be: N . start of array must by 0
# plan to make output_file_0 to contain just headers, and then have N output files

source anaconda3-5.2.0

python3 -u SteadyState_Col0likeRegions_LinParams_SimOnly_Flucts_ClusterOut.py ${SLURM_ARRAY_TASK_ID} 1000 # check last integer compatible with array end
##python3 -u SteadyState_Col0likeRegions_LinParams_SimOnly_Flucts_ClusterOut.py ${SLURM_ARRAY_TASK_ID} 10000 # check last integer compatible with array end