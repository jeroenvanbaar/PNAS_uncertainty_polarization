#!/bin/bash

# Job Name
#SBATCH -J dyad_reg_run3

# Core and memory
#SBATCH -c 1
#SBATCH --mem=32G
#SBATCH --account=carney-ofeldman-condo

# Walltime requested
#SBATCH -t 24:00:00

# Provide index values (voxel set indices used to split the jobs)
#SBATCH --array=0-72
# SBATCH --array=0,
# Use '%A' for array-job ID, '%J' for job ID and '%a' for task ID
#SBATCH -e Logs/dyad_reg_run-3_voxel-set-%a_err.txt
#SBATCH -o Logs/dyad_reg_run-3_voxel-set-%a_out.txt

# Messages to
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jvb@brown.edu

# Use the $SLURM_ARRAY_TASK_ID variable to provide different inputs for each job

run=3
filter_TR=0 # Option to run the analysis on only a segment of the video (indicated by volume (TR) number)
TR_start=1
TR_end=711

module load R/3.5.2

Rscript --vanilla /users/jvanbaar/data/jvanbaar/polarization/NeuroPolitics/Analyses/voxelwise_analyses/dyadic_regression/2.Run_dyadic_regression_cluster.R $SLURM_ARRAY_TASK_ID 1 $run 1000 $filter_TR $TR_start $TR_end