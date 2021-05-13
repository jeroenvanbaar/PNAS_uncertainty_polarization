import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

import nltools, scipy, glob, os, sys
from nltools import Brain_Data, Design_Matrix, expand_mask, collapse_mask
from nltools.plotting import plot_brain, plot_glass_brain
from nltools.stats import one_sample_permutation

# from helpers import *

# Import variables
if len(sys.argv) < 2:
    exit("Needs run, voxel group, n voxels per group as input argument")
else:
    run = int(sys.argv[1])
    voxgr = int(sys.argv[2])
    n_vox_per_group = int(sys.argv[3])

# Settings
base_dir = '/gpfs_home/jvanbaar/data/jvanbaar/polarization/NeuroPolitics/'
mask_dir = '/gpfs_home/jvanbaar/data/jvanbaar/polarization/NeuroPolitics/Analyses/voxelwise_analyses/brain_map_consistency'
mask_fname = (mask_dir + '/80pct_brain_mask.nii.gz')
all_subs = pd.read_csv(base_dir + 'Data/Subjects_and_exclusions/all_subjects.csv')['sub'].values.tolist()
print(len(all_subs))

# Exclusions
exclude_motion = pd.read_csv(base_dir + '/Data/Subjects_and_exclusions/exclude_video-watching_motion.csv'
                            ).query('run == @run')['sub'].values.tolist()
exclude_attention = pd.read_csv(base_dir + '/Data/Subjects_and_exclusions/exclude_video-watching_attention.csv'
                            ).query('run == @run')['sub'].values.tolist()
exclude = exclude_motion + exclude_attention
if run == 1:
    last_TR = 390
elif run == 2:
    last_TR = 307
elif run == 3:
    last_TR = 720
print('Exclusions: %s'%exclude)
sub_list = [i for i in all_subs if i not in exclude]
print('%i subjects remain: %s'%(len(sub_list),sub_list))

# Load r maps for each subject
print('Loading subject-wise correlation maps...')
r_maps = Brain_Data(mask = mask_fname)
for si,sub in enumerate(sub_list):
    print(sub, end = ', ')
    r_maps = r_maps.append(Brain_Data(
        base_dir + 'Results/Inter-subject_correlations/baseline_ISC/run-%i_sub-%03d_corr_mean_others.nii.gz'%(run,sub),
        mask = mask_fname))
print('...done')

# Find voxels
voxels = np.arange(voxgr*n_vox_per_group, (voxgr+1)*n_vox_per_group)
voxels = voxels[voxels<len(r_maps[0].data)] # Only relevant for the last group, which gets truncated
print('Running %i voxels, from %s to %s'%(len(voxels), voxels[:3], voxels[-3:]))

# Take raw data and compute mean r per voxel
# We leave brain space here to enter np arrays
raw_data = r_maps.data[:,voxels]
mean_r = np.mean(raw_data,0) # One value per voxel, mean across rows (subjects)

# Compute p-value per voxel
ps = []
print('Running permutation test per voxel...')
for vi in range(raw_data.shape[1]):
    if np.mod(vi,10) == 0:
        print(vi, end = ',')
    ps.extend([one_sample_permutation(raw_data[:,vi], n_permute=5001)['p']])
out = pd.DataFrame({'r':mean_r,'p':ps,'vox':voxels})

# Store
print('Done. Storing data...')
out.to_csv(base_dir + 'Results/Inter-subject_correlations/baseline_ISC/run-%i_r_permutation-test_voxgr-%i.csv'%(run,voxgr))
print('Done.')