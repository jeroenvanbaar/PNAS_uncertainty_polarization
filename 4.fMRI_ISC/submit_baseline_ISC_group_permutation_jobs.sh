#!/bin/bash

# Job Name
#SBATCH -J baseline_ISC_groupperm_run1

# Core and memory
#SBATCH -c 2
#SBATCH --mem=16G
#SBATCH --account=carney-ofeldman-condo

# Walltime requested
#SBATCH -t 4:00:00

# Provide index values (voxel group IDs)
# SBATCH --array=0,
#SBATCH --array=0-72
# Use '%A' for array-job ID, '%J' for job ID and '%a' for task ID
#SBATCH -e Logs/baseline_ISC_groupperm_run1_voxgr-%a_err.txt
#SBATCH -o Logs/baseline_ISC_groupperm_run1_voxgr-%a_out.txt

# Messages to
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jvb@brown.edu

# Use the $SLURM_ARRAY_TASK_ID variable to provide different inputs for each job

run=1
n_vox_per_group=1000

module load anaconda/3-5.2.0
source activate pp_fmri_new

python ~/data/jvanbaar/polarization/NeuroPolitics/Analyses/Inter-subject_correlation/baseline_ISC_group_permutation.py $run $SLURM_ARRAY_TASK_ID $n_vox_per_group