#!/bin/bash

# Job Name
#SBATCH -J niftify_ISRSA

# Core and memory
#SBATCH -c 4
#SBATCH --mem=32G
#SBATCH --account=carney-ofeldman-condo

# Walltime requested
#SBATCH -t 24:00:00

# Provide index values (runs)
#SBATCH --array=3
# Use '%A' for array-job ID, '%J' for job ID and '%a' for task ID
#SBATCH -e Logs/niftify_ISRSA_run-%a_err.txt
#SBATCH -o Logs/niftify_ISRSA_run-%a_out.txt

# Messages to
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jvb@brown.edu

# Use the $SLURM_ARRAY_TASK_ID variable to provide different inputs for each job

filter_TR=0
TR_start=1
TR_end=711

echo "Collecting ISRSA results for run " $SLURM_ARRAY_TASK_ID

module load R/3.5.2

Rscript --vanilla /users/jvanbaar/data/jvanbaar/polarization/NeuroPolitics/Analyses/voxelwise_analyses/dyadic_regression/3.Collect_results_in_nifti.R $SLURM_ARRAY_TASK_ID $filter_TR $TR_start $TR_end