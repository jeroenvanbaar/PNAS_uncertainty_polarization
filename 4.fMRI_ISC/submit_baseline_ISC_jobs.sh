#!/bin/bash

# Job Name
#SBATCH -J baseline_ISC_run1

# Core and memory
#SBATCH -c 4
#SBATCH --mem=16G
#SBATCH --account=carney-ofeldman-condo

# Walltime requested
#SBATCH -t 1:00:00

# Provide index values (subject IDs)
#SBATCH --array=2,4-8,10-24,28,30-38,43-55
# SBATCH --array=2,
# Use '%A' for array-job ID, '%J' for job ID and '%a' for task ID
#SBATCH -e Logs/baseline_ISC_run1_sub-%a_err.txt
#SBATCH -o Logs/baseline_ISC_run1_sub-%a_out.txt

# Messages to
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jvb@brown.edu

# Use the $SLURM_ARRAY_TASK_ID variable to provide different inputs for each job

run=1

module load anaconda/3-5.2.0
source activate pp_fmri_new

python /users/jvanbaar/data/jvanbaar/polarization/NeuroPolitics/Analyses/Inter-subject_correlation/baseline_ISC_subject.py $SLURM_ARRAY_TASK_ID $run