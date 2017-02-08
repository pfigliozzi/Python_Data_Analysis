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
import glob
sys.path.append(os.path.abspath("C:\Users\Scherer Lab E\Documents\GitHub\Python_Data_Analysis"))
import common_functions as cf


def mode_for_rolling(df):
    '''Pandas does not support calculating the mode for a rolling object. This
    function will allow me to calculate the mode of a rolling object.
    '''

    mode, count = scipy.stats.mode(df)
    return mode

def find_frames_particle_number_diff_rolling_mode(df, window=5):
    '''Returns the frame numbers from a DataFrame where the number of 
    particles in the frame is different from the mode around that frame.

    :param df: DataFrame object without duplicate values that has frames
    marked by a 'frame' column
    :param (int) window: The size of the window to calculate the mode around
    each frame. A window of 5 would calculate the mode including the two
    frames before and after the current frame.
    '''

    particles_per_frame = df.groupby('frame').apply(len)
    rolling_mode = particles_per_frame.rolling(window, center=True).apply(mode_for_rolling)

    # Add values of to the rolling mode to the ends with smaller windows
    for frame in range(window/2):
        rolling_mode.iloc[frame] = rolling_mode.iloc[:window].mode().values
        rolling_mode.iloc[-(frame+1)] = rolling_mode.iloc[-window:].mode().values


    different_num_particles = particles_per_frame != rolling_mode
    frames = different_num_particles[different_num_particles==True].index

    return frames

def trackpy_gaussian_rot_motion_linker(data_frame, search_range, rot_velocity=0.0, fwhm=6, memory=0, **kwargs):
    '''A wrapper for trackpy linking that includes a predictor for rotating particles

    :params data_frame: DataFrame containing all the particle position information
    :params search_range: Max distance a particle can move between frames
    :params rot_velocity: The bias (in degrees) that a candidate particle should be
    found at for each frame. This value reflects the maximum bias applied to the 
    particle based on its position.
    :param fwhm: The full-width-half-maximum of the guassian function which applies
    the rotational bias. Should be set to the full-width-half-maximum of the ring 
    trap in r.
    :params memory: The number of frames a particle can disappear for and still be 
    considered the same particle.
    :params kwrgs: Additional keyword arguments passed to trackpy.link_df
    '''

    # Find the particle locations in polar coords
    xf, yf, rf = cf.least_sq_fit_circle(data_frame)
    cf.polar_coor_data_frame(data_frame, xf, yf)
    
    # Setup the gaussian profile of the bias using the full width half
    # max and the radius for the circle fit.
    std = fwhm/2.3548
    mean = rf
    guass = lambda x: np.exp(-(x - mean)**2/(2 * std**2)) 
    
    # Generate the predictor function
    @trackpy.predict.predictor
    def predict(t1, particle):
        theta = cf.calc_angle(particle.pos[0], particle.pos[1], xf, yf)
        r = cf.calc_radius(particle.pos[0], particle.pos[1], xf, yf)
        
        new_theta = theta + rot_velocity * guass(r) * (t1 - particle.t)
        new_theta %= 360.0
        new_x = cf.calc_x_from_polar(r, new_theta, xf)
        new_y = cf.calc_y_from_polar(r, new_theta, yf)
        return np.array((new_x,new_y))
        
    
    # Track the data and restructure the resulting DataFrame
    trackpy.link_df(data_frame, search_range, memory=memory, pos_columns=['x pos', 'y pos'],
                    retain_index=True, link_strategy='numba', predictor=predict, **kwargs)
    data_frame['track id'] = data_frame['particle']
    del data_frame['particle']
