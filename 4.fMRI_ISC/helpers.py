import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

import nltools, scipy, glob, os, sys
from nltools import Brain_Data, Design_Matrix, Adjacency, expand_mask, collapse_mask
from nltools.plotting import plot_brain, plot_glass_brain
from nltools.prefs import MNI_Template, resolve_mni_path
from nilearn.plotting import plot_img, plot_roi, plot_glass_brain, plot_stat_map

def load_confounds(sub, ses, task, run,
                   prep_dir = ('/gpfs_home/jvanbaar/data/jvanbaar/' + 
                               'polarization/derivatives/fmriprep')
                  ):
    func_dir = prep_dir + '/sub-%03d/ses-%i/func/'%(sub,ses)
    confounds = pd.read_csv(func_dir + 
       'sub-%03d_ses-%i_task-%s_run-%i_desc-confounds_regressors.tsv'%(
       sub, ses, task, run), sep='\t')
    return confounds

def load_cleaned_funx(sub, ses, task, run, space,
                    prep_dir = ('/gpfs_home/jvanbaar/data/jvanbaar/' + 
                               'polarization/derivatives/fmriprep'),
                   cleaned_dir = ('/gpfs_home/jvanbaar/data/jvanbaar/' + 
                               'polarization/derivatives/cleaning')
                  ):
    cleaned_bold_fname = (cleaned_dir + '/sub-%03d/ses-%i/func/'%(sub,ses) +
       'sub-%03d_ses-%i_task-%s_run-%i_space-%s_desc-cleaned_bold.nii.gz'%(
       sub, ses, task, run, space))
    mask_fname = (prep_dir + '/sub-%03d/ses-%i/func/'%(sub,ses) +
                   '/sub-%03d_ses-%i_task-%s_run-%i_space-%s_desc-brain_mask.nii.gz'%(
                   sub, ses, task, run, space))
    dat = Brain_Data(cleaned_bold_fname, mask = mask_fname)
    return dat

def load_cleaned_funx_in_80pct_mask(sub, ses, task, run, space,
                    mask_dir = ('/gpfs_home/jvanbaar/data/jvanbaar/' + 
                               'polarization/NeuroPolitics/Analyses/voxelwise_analyses/brain_map_consistency'),
                   cleaned_dir = ('/gpfs_home/jvanbaar/data/jvanbaar/' + 
                               'polarization/derivatives/cleaning')
                  ):
    cleaned_bold_fname = (cleaned_dir + '/sub-%03d/ses-%i/func/'%(sub,ses) +
       'sub-%03d_ses-%i_task-%s_run-%i_space-%s_desc-cleaned_bold.nii.gz'%(
       sub, ses, task, run, space))
    mask_fname = (mask_dir + '/80pct_brain_mask.nii.gz')
    dat = Brain_Data(cleaned_bold_fname, mask = mask_fname)
    return dat

def plot_sub_timeseries(dat, confounds, sub = None):
    
    nrows = 7
    lw = .5
    fig, axes = plt.subplots(nrows = nrows, ncols = 1, figsize = [6,1.5*nrows])
    
    # Global signal
    axes[0].plot(confounds['global_signal'], lw = lw)
    axes[0].set(title = '%sGS from confounds file'%(
        ('sub-%03d '%sub) if sub is not None else ''))
    
    # Mean voxel signal
    axes[1].plot(dat.data.mean(axis=1),'g', lw = lw)
    axes[1].set(title = 'sub-%03d mean voxel signal'%sub)
    # FD
    axes[2].plot(confounds['framewise_displacement'], lw = lw)
    axes[2].set(title = 'sub-%03d framewise displacement'%sub)
    # DVARS
    axes[3].plot(confounds['std_dvars'], lw = lw)
    axes[3].set(title = 'sub-%03d std_dvars'%sub)
    # CSF
    axes[4].plot(confounds['csf'], lw = lw)
    axes[4].set(title = 'sub-%03d csf'%sub)
    # WM
    axes[5].plot(confounds['white_matter'], lw = lw)
    axes[5].set(title = 'sub-%03d white_matter'%sub)
    # Some voxels
    axes[6].plot(dat.data[:,1000],'g', lw = lw)
    axes[6].set(title = 'sub-%03d voxel 1000'%sub)
    
    plt.tight_layout()
    
    return fig, axes

def make_volume_regressor(amplitude_envelope, rate, TR = 1.5, video_onset = 6,
                          post_pad = 30, visual_check = True):
    
    # Pre-pad with zeros until video onset:
    print('Pre-padding with a %.1f-second pause before video...'%video_onset)
    pre_onset = np.zeros([1,video_onset*rate])[0]
    print('Post-padding with a %.1f-second pause after video...'%post_pad)
    post_audio = np.zeros([1,rate*post_pad])[0] # Add 20 seconds of silence
    volume_signal = np.hstack([pre_onset, amplitude_envelope, post_audio])

    # Convolve
    print('Convolving...')
    dm = Design_Matrix(pd.DataFrame(volume_signal, columns = ['volume']),
                       sampling_freq = rate).convolve()
    volume_regressor = dm.values.flatten()
    len_mins = np.floor((len(volume_regressor)/rate)/60)
    len_secs = len(volume_regressor)/rate - len_mins*60
    print('Regressor duration is %i:%.2f.'%(len_mins, len_secs))

    # Average by TR to align with BOLD sampling frequency
    print('Downsampling...')
    volume_regressor_TRfreq = [np.mean(volume_regressor[int(i*TR*rate):int((i+1)*TR*rate)])
                             for i in range(int(np.floor(len(volume_signal)/(rate*TR))))]
    n_TRs = len(volume_regressor_TRfreq)
    print('Regressor length is %i TRs.'%(n_TRs))

    # Visually check audio envelope
    if visual_check:
        fig,ax = plt.subplots(1,1)
        ax.plot(amplitude_envelope)
        ax.plot(np.arange(0,n_TRs*TR*rate, TR*rate),
                 np.multiply(volume_regressor_TRfreq,1))

    return volume_regressor_TRfreq

def store_brain_plots(obj, out_file, title = None, include_glass = True, include_ortho = True):
    
    if (out_file[-4] == '.') & ((out_file[-3:] == 'png') | (out_file[-3:] == 'pdf')):
        out_file = out_file[:-4]
    
    obj_nii = obj.to_nifti()
    
    cmap = 'RdBu_r'
    
    if include_glass:
        print('Writing glass brain plot...')
        if title is not None:
            plot_glass_brain(obj_nii, display_mode='lzry', colorbar=True, cmap=cmap, plot_abs=False,
                         title = title, output_file = out_file + '_glass.pdf')
        else:
            plot_glass_brain(obj_nii, display_mode='lzry', colorbar=True, cmap=cmap, plot_abs=False,
                         output_file = out_file + '_glass.pdf')
        
        
    if include_ortho:
        print('Writing ortho plots...')
        views = ['x', 'y', 'z']
        cut_coords = [range(-40, 50, 10),
                      [-88, -72, -58, -38, -26, 8, 20, 34, 46],
                      [-34, -22, -10, 0, 16, 34, 46, 56, 66]]
        
        for v, c in zip(views, cut_coords):
            plot_stat_map(obj_nii, cut_coords=c, display_mode=v,
                          cmap=cmap, bg_img=resolve_mni_path(MNI_Template)['brain'],
                          output_file = out_file + '_ortho_%s.pdf'%v)
    