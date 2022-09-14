#!/bin/bash -e
#SBATCH --mail-type=END,FAIL # notifications for job done & fail
#SBATCH --mail-user=amy.briffa@jic.ac.uk # send-to address
#SBATCH --mem 4000 # memory pool for all cores
##SBATCH -p rg-mh
#SBATCH -e slurm_logs/slurm-%j.err 
#SBATCH -o slurm_logs/slurm-%j.out
## No array job needed for this code! Designed to run on single core. 
##SBATCH --array=0-1 # end of array (N).  2nd argument of python code should be: N . start of array must by 0. 
# plan to make output_file_0 to contain just headers, and then have N output files

source anaconda3-5.2.0

python3 -u HighM_accessions_screen.py