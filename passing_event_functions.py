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

def add_del_r_del_theta(dataframe):
    """A function that adds delta theta and delta r columns to a data 
    frame for a particles nearest neighbors in theta.

    Using the nearest neighbors in theta it will calculate the delta r 
    and delta theta for each of those nearest neighbors in theta. The theta
    nearest neighbor is not always the same as the cartesian nearest 
    neighbor. The delta r and delta theta will be aligned with the 
    theta nearest neighbor id. The 'del_theta' and 'del_r' columns describe
    the position of the nearest neighbor relative to the current particle
    ('track id').

    :params dataframe: DataFrame that contains columns for ['frame', 'track id',
    'x pos', 'y pos', 'r', 'theta', 'theta_nn_id'] that you want to add the 
    'del_r' and 'del_theta' columns to.
    """
    df_copy = dataframe.copy()
    # Align nearest neighbors with the theta_nn_id col.
    frame_track_indexed = df_copy.set_index(['frame', 'track id'])
    frame_track_indexed = frame_track_indexed.drop_duplicates(['x pos', 'y pos'])

    frame_track_keys = df_copy[['frame', 'theta_nn_id']].values
    zipped_keys = zip(frame_track_keys[:,0], frame_track_keys[:,1])
    frame_track_aligned = frame_track_indexed.ix[zipped_keys].reset_index()
    
    del_theta = frame_track_aligned['theta'] - df_copy['theta'] 
    del_theta = del_theta - 360 * np.round(del_theta/360)
    df_copy['del_theta'] = del_theta
    
    df_copy['del_r'] = frame_track_aligned['r'] - df_copy['r']
    return df_copy

def find_passing_events(df):
    """A function that finds when a passing event occurs by looking
    at when the relative theta coordinate changes between a pair of particles
    
    This function will call 'add_del_r_del_theta' to calculate the relative
    coordinates in 'r' and 'theta' for all particle pairs in the DataFrame.
    This function will look at the sign of 'del_theta' for a given pair and
    determine when 'del_theta' goes from negative to positive (or vice versa)
    to declare it as a passing event. The passing events will be denoted with
    booling column, 'passing_event', appended on the returned DataFrame. Note:
    the relative coordinates 'del_r' and 'del_theta' are also appended to the 
    DataFrame returned from this function.
    
    :param df: DataFrame that contains columns for ['frame', 'track id', 
    'theta', 'r', 'x pos', 'y pos'].
    """
    
    # Add relative coordinates
    df_rel = add_del_r_del_theta(df)
    
    # Create a DataFrame with a 1 frame shift
    df_new_index = df_rel.set_index(['frame', 'track id', 'theta_nn_id'])
    new_index = df_new_index.index
    new_index_array = np.asarray(list(new_index.values))
    new_index_array[:, 0] += 1 # Shift the 'frame' index by 1
    selector = list(map(tuple, new_index_array))
    df_frame_shift = df_new_index.loc[selector]
    
    # Find where 'del_theta' changes sign
    del_theta_change = df_new_index['del_theta'].values * df_frame_shift['del_theta'].values
    pass_transition_state = (del_theta_change < 0)
    within_90_deg = (abs(del_theta_change) < 90**2) # Don't count >180 passing events
    df_rel['passing_event'] = pass_transition_state & within_90_deg
    
    return df_rel

def add_dwell_time_and_travel_distance_passing_events(df, del_theta_separation, min_r=40, dev_from_center_start=2):
    '''Add columns to passing events indicating the travel distance of a pair
    of particles and the dwell time of a pair before they pass.
    
    :param df: DataFrame object created for analyzing the passing event data
    :param del_theta_separation: The separation (in degrees) that a pair of
    particles need to be to be considered a pair and begin calculating travel
    distance and dwell time.
    :param min_r: The minimum absolute delta r particles are allowed to have
    to be considered a legitimate passing event.
    :param dev_from_center_start: The deviation that each particle in the pair
    can have from the mean radius to begin calculating dwell time and travel 
    distance.
    '''
    
    df_copy = df.copy()
    
    df_copy['begin_paired_up_frame'] = np.nan
    df_copy['distance_traveled_as_pair'] = np.nan
    df_copy['dwell_time_to_pass'] = np.nan
    
    # Find all the passing events that meet the criteria
    df_events = df_copy[df_copy.last_pass_event == True]
    df_events = df_events[df_events.del_theta < 0]
    df_events = df_events[abs(df_events.del_r) < min_r]    
    
    for idx, series in df_events.iterrows():
        track_id = series['track id']
        pair_id = series['pair_id']
        last_pass_frame = series['frame']
        df_valid = df_copy.query('frame < @last_pass_frame')
        df_valid = df_valid[df_valid.pair_id == pair_id]
        df_valid = df_valid[df_valid['track id'] == track_id]
        
        #print df_valid.query('del_theta < @del_theta_separation').drop(['frame', 'track id', 'x pos', 'y pos', 'nn_num', 'nn_id', 'nn_dist'], axis=1)
        df_less_than_sep = df_valid.del_theta < -del_theta_separation
        df_less_than_sep = df_less_than_sep > -90
        df_greater_than_sep = df_valid.del_theta > -del_theta_separation
        df_greater_than_sep = df_greater_than_sep < 90
        df_particle_close_trap = abs(df_valid.dist_avg_r) < dev_from_center_start
        df_other_particle_close_trap = abs(df_valid.dist_avg_r_other_part) < dev_from_center_start
        
        entering_range = df_less_than_sep.shift(1) * df_greater_than_sep * df_particle_close_trap * df_other_particle_close_trap
        if len(entering_range[entering_range == 1]) == 0:
            continue
        latest_entering_range_index = entering_range[entering_range == 1].index[-1]
        
        frame_begin_pair = df_valid.loc[latest_entering_range_index, 'frame']
        
        assignment_sel = (df.frame == last_pass_frame) & (df.pair_id == pair_id)
        df_copy.loc[assignment_sel, 'begin_paired_up_frame'] = frame_begin_pair
        
        df_find_travel = df.query('@frame_begin_pair <= frame <= @last_pass_frame')
        df_find_travel = df_find_travel[df_find_travel.pair_id == pair_id]
        df_find_travel = df_find_travel[df_find_travel['track id'] == track_id]
        theta_midpoint = ((df_find_travel.del_theta/2.0) + df_find_travel.theta) % 360
        theta_disp = theta_midpoint.shift(-1) - theta_midpoint
        theta_disp = theta_disp - 360 * np.round(theta_disp/360)
        
        df_copy.loc[assignment_sel, 'distance_traveled_as_pair'] = theta_disp.sum()
        df_copy.loc[assignment_sel, 'dwell_time_to_pass'] = last_pass_frame - frame_begin_pair
    
    return df_copy

def select_frame_from_experiment(frame, key, store, new_base_dir=None):
    '''Selects a particlular image frame from an experiment when you give
    the data base, the key, and the frame from an HDF store that has appended
    paths of image data.
    '''
    image_path = store.index.image_path[store.index['key'] == key].values[0]
    if not new_base_dir == None:
        image_path_components = image_path.split(os.sep)[-2:]
        image_path = os.path.join(new_base_dir, *image_path_components)
    image = plt.imread(image_path + '\\' + str(10000+frame) + '.tif')
    return image

def high_res_points(x, y, factor=10):
    '''Takes a trajectory and makes a high resolution trajectory by
    linear interpolation between each point
    
    :param x (arr): x-values of the trajectory
    :param y (arr): x-values of the trajectory
    :param factor (int): The factor of the number of points to increase
    the trajectory.
    '''
    new_x = []
    new_y = []
    for i in range(len(x)-1):
        m = (y[i+1] - y[i]) / (x[i+1] - x[i])
        b = y[i] - m*x[i]
        calc_x = np.linspace(x[i], x[i+1], factor)
        new_x += list(calc_x[:-1])
        new_y += list(calc_x[:-1]*m + b)
        last_x = x[i+1]
        last_y = y[i+1]
    new_x.append(last_x)
    new_y.append(last_y)
    return np.asarray(new_x), np.asarray(new_y)

def two_color_cmap(color1, color2):
    '''Creates a matplotlib color map between two colors. Colors
    entered can be strings such as hex or matplotlib color codes (e.g.
    'r' or 'C0')
    
    :param color1: The first color in the color map.
    :param color2: The end color of the color map
    '''
    rgb_c1 = colors.to_rgb(color1)
    rgb_c2 = colors.to_rgb(color2)
    
    cdict={'red': ((0.0, rgb_c1[0], rgb_c1[0]),
                  (1.0, rgb_c2[0], rgb_c2[0])),
          'green': ((0.0, rgb_c1[1], rgb_c1[1]),
                   (1.0, rgb_c2[1], rgb_c2[1])),
          'blue': ((0.0, rgb_c1[2], rgb_c1[2]),
                  (1.0, rgb_c2[2], rgb_c2[2]))}
    return colors.LinearSegmentedColormap('new_cmap', cdict)


def plot_trajectory(x, y, cmap):
    '''Plots a single trajectory'''
    if isinstance(x, pd.Series):
        x = x.values
    if isinstance(y, pd.Series):
        y = y.values
    x, y = high_res_points(x, y, 2)
    color_array = cmap(np.linspace(0, 1, len(x)))
    colored_plot(x, y, color_array, lw=3)
    plt.xlim(min(x), max(x))
    plt.ylim(min(y), max(y))
    
def plot_df_pass_events(key, store, frame_range=(-5, 10)):
    '''Creates a plot of the trajectory on top of the raw movie data'''
    
    df = store.get(key).copy()
    
    pass_events = df[df['last_pass_event'] == True]
    
    grouped = pass_events.groupby(['frame', 'pair_id'])
    
    for (frame, pair_id), data in grouped:
        image = select_frame_from_experiment(frame + frame_range[0], key, store)
        plt.imshow(image, cmap='gray', origin='upper')

        
        track_behind = data[data.del_theta < 0]['track id'].iloc[0]
        track_ahead = data[data.del_theta > 0]['track id'].iloc[0]
        
        frame_start = frame + frame_range[0]
        frame_end = frame + frame_range[1]
        print key
        print 'Frame '+str(frame_start)+' to '+str(frame_end)
        print 'L = '+str(store.index[store.index.key == key].L.values[0])

        plot_values = df.query('@frame_start <= frame <= @frame_end')
        
        traj_behind = plot_values[plot_values['track id'] == track_behind]
        traj_behind = traj_behind.drop_duplicates(['frame', 'track id'])
        
        traj_ahead = plot_values[plot_values['track id'] == track_ahead]
        traj_ahead = traj_ahead.drop_duplicates(['frame', 'track id'])
        
        cmap1 = two_color_cmap('C0', 'g')
        cmap2 = two_color_cmap('C4', '#CC0000')
        
        y_size = image.shape[1]
        plot_circle(traj_behind.iloc[0], y_size)
        
        plot_trajectory(traj_behind['x pos'], image.shape[1] - traj_behind['y pos'], cmap1)
        plot_trajectory(traj_ahead['x pos'], image.shape[1] - traj_ahead['y pos'], cmap2)
        
        both_trajs = pd.concat([traj_behind, traj_ahead]).copy()
        both_trajs['y pos'] = image.shape[1] - both_trajs['y pos']
        min_max_x = (both_trajs['x pos'].min()-10, both_trajs['x pos'].max()+10)
        min_max_y = (both_trajs['y pos'].max()+10, both_trajs['y pos'].min()-10)
        
        plt.xlim(min_max_x)
        plt.ylim(min_max_y)
        
        plt.show()

def plot_circle(series, y_size):
    '''Plot a circle on an image indicating the average radius of 
    the circle. The values are calculated from a single row from
    a DataFrame representing one experiment. Performs the y-flip as well.
    '''
    x = series['x pos']
    y = series['y pos']
    r = series.r
    theta = series.theta
    
    x_cent, y_cent = cf.calc_cent_from_polar(x, y, r, theta)
    
    angle = np.linspace(0, 360, 1000)
    
    y_cent = y_size - y_cent
    
    x_plot = cf.calc_x_from_polar(series.r_avg, angle, x_cent)
    y_plot = cf.calc_y_from_polar(series.r_avg, angle, y_cent)
        
    plt.plot(x_plot, y_plot, 'y')
