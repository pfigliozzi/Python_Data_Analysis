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
