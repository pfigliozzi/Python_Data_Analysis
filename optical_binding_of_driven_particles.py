import string
import scipy
import Tkinter, tkFileDialog
import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist
import matplotlib.pyplot as plt
import os
import sys
import skimage.viewer
from PIL import Image, ImageDraw, ImageColor
import trackpy
import trackpy.predict

def plot_nn_dist_par_perp(data_frame,fit_params,plot_name):
    ang_sep,ang_bins=nn_distance_angle_seperation_ver_2(data_frame,8,1)
    print len(ang_sep)
    ang_sep=np.array([np.concatenate((ang_sep[v],ang_sep[v+4])) for v in range(len(ang_sep)/2)])
    um_conv=6.5/60/1.6/2
    print len(ang_sep)
    ang_sep_conv=ang_sep*um_conv
    #labels=['Parellel','Perpendicular']
    print ang_bins
    #print ang_sep
    labels=[str(v[0])+'-'+str(v[1]) for v in ang_bins]
    plt.figure(figsize=[4,3])
    n,b,p=plt.hist([v for v in ang_sep_conv],bins=100,
                   range=(0,3),normed=True,histtype='step',label=labels)
    plt.legend()
    plt.xlabel('Seperation (um)')
    plt.ylabel('Probability Density')
    plt.title(plot_name)
    #fig_dir="J:\Pat's Projects\Dynamical Phase Transition\Figures for Paper\\nn_histogram_comparison\\"
    #plt.savefig(fig_dir+plot_name[:-4]+".png", dpi=200)
    plt.show()

def filter_df_by_particles_in_frame(data_frame, num_particles, mode='equal'):
    '''Return a DataFrame where just the frames with the requested number
    of particles are present. This only works on DataFrames that have gone
    through find_nn_ver_2
    
    :param data_frame: The input data_frame
    :param num_particles: The number of particles you want in each frame
    :param mode: Can be 'equal' means only frames with num_particles
    are returned. 'less' means only frames less than or equal to num_particles are
    returned. 'greater' means only frames greater than or equal to num_particles 
    are returned.
    :return data_frame:
    '''
    # The function does not compute the number of particles in each frame
    # but instead the number of entries for a given number of particles
    # The number of entries = 2 * (num_particles nCr 2)
    from scipy.special import comb
    if num_particles!=1:
        num_particles = 2 * (comb(num_particles,2))
    data = data_frame.copy()
    part_num_in_frame = data.groupby('frame').apply(len)
    if mode == 'equal':
        return data.set_index('frame')[part_num_in_frame==num_particles].reset_index()
    elif mode == 'less':
        return data.set_index('frame')[part_num_in_frame<=num_particles].reset_index()
    elif mode == 'greater':
        return data.set_index('frame')[part_num_in_frame>=num_particles].reset_index()


def filter_by_n_nn_with_more_d_distance_away(data_frame, nn_num_allowed=1, dist_extra_nn_allowed=2.0, min_nn_dist_allowed=float('inf'), um_conv=6.5/60/1.6/2):
    '''Filter a group of n nearest neighbors where the n+1 neighbor is a minimum distance away from 
    the group.
    
    :params data_frame: A data frame containing all the particle trajectory information
    :params nn_num_allowed: The number of nearest neighbors you want in the group (1 nn means you are 
    looking for pairs of particles)
    :params (um) dist_extra_nn_allowed: The minimum distance a nearest neighbor of nn_num_allowed+1 has to
    be away from the nearest neighbor group
    :params min_nn_dist_allowed: The minimum distance allowed for the nn of interest to be in. Given 
    that you are targeting n nearest neighbors if their interparticle distance is not less than 
    :params um_conv: The conversion factor to change the distances from pixels to microns
    '''
    um_data_frame = data_frame.copy()
    # Convert the NN Dist to um
    um_data_frame['nn_dist'] *= um_conv
    # Pivot the data to look down the nn_num and view the nn_dis
    frm_trkid_by_nnnum = um_data_frame.pivot_table(index=['frame','track id'], columns='nn_num', values='nn_dist')
    # Find valid nn less than the min_nn_distance_allowed
    valid_nn = frm_trkid_by_nnnum.loc[:,:nn_num_allowed] < min_nn_dist_allowed
    # Find valid nn where nn_num+1 is greater than dist_extra_nn_allowed away
    valid_greater_nn_away = frm_trkid_by_nnnum.loc[:,nn_num_allowed+1:] > dist_extra_nn_allowed
    # NaN values give False when they should be True
    valid_greater_nn_away[pd.isnull(frm_trkid_by_nnnum.loc[:,nn_num_allowed+1:])] = True
    valid_frm_trkid = valid_nn.all(axis=1) & valid_greater_nn_away.all(axis=1)
    return data_frame.set_index(['frame','track id']).ix[valid_frm_trkid[valid_frm_trkid].index].reset_index()
    
    