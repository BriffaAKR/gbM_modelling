#!/bin/bash -e
#SBATCH --mail-type=END,FAIL # notifications for job done & fail
#SBATCH --mail-user=amy.briffa@jic.ac.uk # send-to address
#SBATCH --mem 4000 # memory pool for all cores
# #SBATCH -p rg-mh
#SBATCH -e slurm_logs/slurm-%j.err 
#SBATCH -o slurm_logs/slurm-%j.out

# script to run graphs code

source anaconda3-5.2.0

python3 GraphsInitialStates.py