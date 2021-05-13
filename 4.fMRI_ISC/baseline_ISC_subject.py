# This script computes the 'baseline' inter-subject correlations of neural activity during video watching.
# 'Baseline' is meant in the sense of 'not associated with predictors such as ideology', simply shared activity across all subjects regardless of creed etc.
# This is the first step in many ISC papers.

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os, sys
import nltools, scipy, glob
from nltools import Brain_Data, Design_Matrix, expand_mask, collapse_mask
from nltools.plotting import plot_brain, plot_glass_brain

from helpers import load_cleaned_funx_in_80pct_mask

# Import variables
if len(sys.argv) < 3:
    exit("Needs subject number, run as input argument")
else:
    sub = int(sys.argv[1])
    run = int(sys.argv[2])

# Settings
base_dir = '/gpfs_home/jvanbaar/data/jvanbaar/polarization/NeuroPolitics/'
ses = 1
task = 'videoWatching'
space = 'MNI152NLin2009cAsym'
TR = 1.5
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
print(sub_list)
if (sub not in sub_list) or (sub in exclude):
    print('Subject %i must be excluded for this run... aborting...'%sub)
    sys.exit()
else:
    other_list = [i for i in sub_list if i is not sub]
    print('Comparing sub %i to %i others: %s'%(sub,len(other_list),other_list))

# Load subject's own data
n_TRs = []
print('Loading subject %is own data...'%sub, end = '')
own_image = load_cleaned_funx_in_80pct_mask(sub, ses, task, run, space)
print(' with %i TRs'%own_image.shape()[0])
n_TRs.extend([own_image.shape()[0]])
min_n_TR = np.min(n_TRs)
print('\nBrain mask has the following affine matrix:')
print(own_image.mask.affine)
print('\n')

# We are trying to correlate sub against mean of others. Since corr to mean == corr to sum,
# we simply sum all signals from the other subjects. We thus dont stack the 'others' data but aggregate
# instantly using the sum to avoid overloading working memory.
# Start with the first
print('Loading other sub 1 out of %i: sub-%03d...'%(len(other_list),other_list[0]), end = '')
others_image = load_cleaned_funx_in_80pct_mask(other_list[0],ses,task,run,space) # Start with the first
print(' with %i TRs'%others_image.shape()[0])
n_TRs.extend([others_image.shape()[0]])
if np.min(n_TRs) < min_n_TR:
    min_n_TR = np.min(n_TRs)
    print('Truncating images to %i TRs.'%min_n_TR)
    others_image = others_image[:min_n_TR]
# Add the rest
for otheri,other in enumerate(other_list[1:]):
    print('Loading other sub %i out of %i: sub-%03d...'%(otheri+2,len(other_list),other), end = '')
    yet_another_image = load_cleaned_funx_in_80pct_mask(other,ses,task,run,space)
    print(' with %i TRs'%yet_another_image.shape()[0])
    n_TRs.extend([yet_another_image.shape()[0]])
    if np.min(n_TRs) < min_n_TR:
        min_n_TR = np.min(n_TRs)
        print('Truncating images to %i TRs.'%min_n_TR)
    others_image = sum([others_image[:min_n_TR],yet_another_image[:min_n_TR]])
# Take off the first 4 TRs (countdown screen):
print('Removing first 4 volumes...')
others_image = others_image[4:]
print('Shape of others image is now:')
print(others_image.shape())

# Truncate own image TRs if needed
own_image = own_image[4:min_n_TR]
print('Shape of own image is now:')
print(own_image.shape())

# Correlate own image to others summed image (correlation is invariant to scale)
print('Correlating voxel-wise time series of sub versus others...')
rs = []
for vi in range(len(own_image[0])):
    if np.mod(vi,1000) == 0:
        print(vi, end = ',')
    r = np.corrcoef(own_image.data[:,vi],others_image.data[:,vi])[0,1]
    rs.extend([r])
print('Done correlating.')

# Store r map
print('Creating correlation map...')
r_map = own_image[0].copy()
r_map.data = np.array(rs)
r_map.write(base_dir + '/Results/Inter-subject_correlations/baseline_ISC/run-%i_sub-%03d_corr_mean_others.nii.gz'%(run,sub))
print('Done.')