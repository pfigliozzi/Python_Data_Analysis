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

def matlab_gui_to_data_frame_with_intesity(m):
    return pd.DataFrame(
        {'x pos': m[0][:], 'y pos': m[1][:], 'frame': m[4][:], 'particle id': m[3][:], 'track id': m[5][:], 'intensity': m[2][:], 'sigma': m[-2][:]})


def find_k_nn(grp, k=2):
    '''A helper function to find the 1st NN for the function
    filter_close_false_particles_over_plate.'''
    if len(grp) == 1:
        grp['nn_id'] = np.nan
        grp['nn_dist'] = np.nan
        return grp
    xy_data = grp[['x pos','y pos']]
    tree = scipy.spatial.KDTree(xy_data)
    nn_d, nn_i=tree.query(xy_data, k=k)
    particle_ids = grp['particle id'].values
    nn_id = particle_ids[nn_i[:,1]]
    nn_dist = nn_d[:,1]
    grp['nn_id'] = nn_id
    grp['nn_dist'] = nn_dist
    return grp

def sorted_id_string(series):
    '''A helper function for filter_close_false_particles_over_plate
    that will give each pair of NN a unique string (the same string) to
    identify them. '''
    particle_list = [series['particle id'].values, series.nn_id.values]
    particle_list = sorted(particle_list)
    string_list = [str(int(i[0])) for i in particle_list]
    string = ''
    for num, part_id in enumerate(string_list):
        if num == 0:
            string += part_id
        else:
            string += '+'+part_id
    return string

def filter_close_false_particles_over_plate(data_frame, max_nn_dist, angle_bounds=(0,360), nn_num=1):
    '''This function filters out false particles due to halo effect around 
    particles over the nanoplate.
    
    The identification of the Half Nanoplate data produce several false particles
    due to a halo around the particles over the nanoplate. Usually these false
    particles appear less than one particle diameter from the real one. Removing
    them will give more consistent trajectories for the tracking. The real 
    particles are assumed to be the ones closer to the average radius of particles
    over the nanoplate and are also within some NN distance around one particle
    diameter.
    
    :param data_frame: The data_frame that contains the particle xy data and has
    particles over the nanoplate labeled
    :param max_nn_dist: The maximum distance two particles will be flagged with one
    potentially being false. Should be set to be about the particle diameter.
    :param angle_bounds: If you want to limit the search for false particles to
    a particular angle range (beyond just checking over the nanoplate) this can
    be set.
    '''
    over_plate_particles = data_frame[data_frame.over_plate == True].copy()
    lower_angle = angle_bounds[0]
    upper_angle = angle_bounds[1]
    valid_particles = over_plate_particles.query('@lower_angle < theta < @upper_angle')
    # Find first NN of all particles in all frames (over nanoplate)
    df_nn = valid_particles.groupby('frame', group_keys=False).apply(find_k_nn, k=nn_num+1)
    r_mean = df_nn.r.mean()
    # Find particles with less than nn_dist separation
    particles_within_bounds = df_nn[df_nn.nn_dist < max_nn_dist].copy()
    particles_within_bounds['del_r'] = particles_within_bounds['r']-r_mean
    paired_df = particles_within_bounds.groupby(['frame','particle id'], group_keys=False).apply(lambda x: x.assign(pair = sorted_id_string(x)))
    bad_indices = paired_df.groupby(['frame','pair']).apply(lambda x: abs(x.del_r).idxmax())
    return data_frame.drop(bad_indices.values)

def filter_close_false_particles_over_plate_intensity(data_frame, max_nn_dist, angle_bounds=(0,360), nn_num=1):
    '''This function filters out false particles due to halo effect around 
    particles over the nanoplate with a preference of keeping particles with higher
    integrated intensity.
    
    The identification of the Half Nanoplate data produce several false particles
    due to a halo around the particles over the nanoplate. Usually these false
    particles appear less than one particle diameter from the real one. Removing
    them will give more consistent trajectories for the tracking. The real 
    particles are assumed to be the ones with a higher integrated intensity and
    also within some NN distance around one particle diameter.
    
    :param data_frame: The data_frame that contains the particle xy data and has
    particles over the nanoplate labeled
    :param max_nn_dist: The maximum distance two particles will be flagged with one
    potentially being false. Should be set to be about the particle diameter.
    :param angle_bounds: If you want to limit the search for false particles to
    a particular angle range (beyond just checking over the nanoplate) this can
    be set.
    :param nn_num: The number of nearest neighbors to account for when determining 
    which particles to get rid of.
    '''
    over_plate_particles = data_frame[data_frame.over_plate == True].copy()
    lower_angle = angle_bounds[0]
    upper_angle = angle_bounds[1]
    valid_particles = over_plate_particles.query('@lower_angle < theta < @upper_angle')
    # Find first NN of all particles in all frames (over nanoplate)
    df_nn = valid_particles.groupby('frame', group_keys=False).apply(find_k_nn, k=nn_num+1)
    r_mean = df_nn.r.mean()
    # Find particles with less than nn_dist separation
    particles_within_bounds = df_nn[df_nn.nn_dist < max_nn_dist].copy()
    particles_within_bounds['del_r'] = particles_within_bounds['r']-r_mean
    paired_df = particles_within_bounds.groupby(['frame','particle id'], group_keys=False).apply(lambda x: x.assign(pair = sorted_id_string(x)))
    bad_indices = paired_df.groupby(['frame','pair']).apply(lambda x: abs(x.intensity).sort)
    return data_frame.drop(bad_indices.values)

def particle_number_difference_between_data_sets_on_raw_data(df1, df2, raw_image_path=None, theta_subset=(240, 310), image_size=(390, 390), particle_size=6.0):
    '''Draws the positions of two different data frames on top of raw image data for 
    visualizing the positions from frame to frame. Shows only frames where the number
    of particles in df1 differ from df2.

    :params df1: The first DataFrame with particle positions. These positions
    will appear in red.
    :params df2: The second DataFrame of particle positions. These positions 
    will appear in yellow
    :params raw_image_path: The full path to the raw data that corresponds 
    with df1 and df2
    :params theta_subset: Only positions within this theta value are shown
    :params image_size: The size of the raw image data
    :params particle_size: The size to draw the particles on the raw data
    (in pixels)
    '''

    df1_theta = df1.query('@theta_subset[0] < theta < @theta_subset[1]')
    df2_theta = df2.query('@theta_subset[0] < theta < @theta_subset[1]')
    part_in_frame1 = df1_theta.groupby('frame').apply(len)
    part_in_frame2 = df2_theta.groupby('frame').apply(len)
    
    concat = pd.concat([part_in_frame1, part_in_frame2], axis=1)
    frame_mismatch_bool = concat.loc[:,0] != concat.loc[:,1]
    frame_mismatch = frame_mismatch_bool.index[frame_mismatch_bool]
    if raw_image_path != None:
        raw_image_list = glob.glob(raw_image_path+"*.tif")
    image_frames = []
    for frame in frame_mismatch:
        if raw_image_path != None:
            curr_image = Image.open(raw_image_list[int(frame)-1])
            img_arr = np.array(curr_image.getdata()).reshape(curr_image.size)
            min_pixel = np.min(img_arr)
            img_arr = img_arr - min_pixel
            max_pixel = np.max(img_arr)
            scale = 255/float(max_pixel)
            img_arr = np.round(img_arr * scale)
            img_arr = np.flipud(img_arr)
            curr_image = Image.fromarray(img_arr)
            img_frame = curr_image.convert(mode='RGB')
        else:
            img_frame = Image.new('RGB', image_size, 'black')
        img_draw = ImageDraw.Draw(img_frame)
        img_draw.text((3,3), 'Frame '+str(frame),  fill='gray')
        for idx, particle in df1_theta[df1_theta['frame'] == frame].iterrows():
            # Draw a circle for a particle in each frame in df1
            ellipse_bb = [tuple(particle['x pos':'y pos'].values-particle_size),
                          tuple(particle['x pos':'y pos'].values+particle_size)]
            img_draw.ellipse(ellipse_bb, outline='red')
            track_num = particle['track id']
            
        for idx, particle in df2_theta[df2_theta['frame'] == frame].iterrows():
            # Draw a circle for a particle in each frame in df2
            ellipse_bb = [tuple(particle['x pos':'y pos'].values-particle_size),
                          tuple(particle['x pos':'y pos'].values+particle_size)]
            img_draw.ellipse(ellipse_bb, outline='yellow')
            track_num = particle['track id']
        
        df1_frame = np.floor(df1_theta[df1_theta['frame'] == frame].loc[:,'x pos':'y pos'])
        df2_frame = np.floor(df2_theta[df2_theta['frame'] == frame].loc[:,'x pos':'y pos'])
        for num, (idx, row) in enumerate(df1_frame.iterrows()):
            if (row == df2_frame).all(axis=1).any(axis=0):
                ellipse_bb = [tuple(row.values-particle_size),
                          tuple(row.values+particle_size)]
                img_draw.ellipse(ellipse_bb, outline='purple')
        
        image_frame_array = np.array(img_frame)
        img_frame.close()
        image_frames.append(image_frame_array)
    
    skimage.viewer.CollectionViewer(image_frames).show()
    gc.collect()

def draw_particles_on_raw_data(df1, df2, raw_image_path=None, theta_subset=(240, 310), image_size=(390, 390), particle_size=6.0, frame_range=None):
    '''Draws the positions of two different data frames on top of raw image data for 
    visualizing the positions from frame to frame.

    :params df1: The first DataFrame with particle positions. These positions
    will appear in red.
    :params df2: The second DataFrame of particle positions. These positions 
    will appear in yellow
    :params raw_image_path: The full path to the raw data that corresponds 
    with df1 and df2
    :params theta_subset: Only positions within this theta value are shown
    :params image_size: The size of the raw image data
    :params particle_size: The size to draw the particles on the raw data
    (in pixels)
    :params frame_range: The subset of frames to show the positions for
    '''

    df1_theta = df1.query('@theta_subset[0] < theta < @theta_subset[1]')
    df2_theta = df2.query('@theta_subset[0] < theta < @theta_subset[1]')
    
    if raw_image_path != None:
        raw_image_list = glob.glob(raw_image_path+"*.tif")
    if frame_range == None:
        min_frame = pd.concat([df1_theta, df2_theta]).frame.min()
        max_frame = pd.concat([df1_theta, df2_theta]).frame.max()
        frame_range = (min_frame, max_frame)
    frames = xrange(frame_range[0], frame_range[1]+1)
    image_frames = []
    for frame in frames:
        if raw_image_path != None:
            curr_image = Image.open(raw_image_list[int(frame)-1])
            img_arr = np.array(curr_image.getdata()).reshape(curr_image.size)
            min_pixel = np.min(img_arr)
            img_arr = img_arr - min_pixel
            max_pixel = np.max(img_arr)
            scale = 255/float(max_pixel)
            img_arr = np.round(img_arr * scale)
            img_arr = np.flipud(img_arr)
            curr_image = Image.fromarray(img_arr)
            img_frame = curr_image.convert(mode='RGB')
        else:
            img_frame = Image.new('RGB', image_size, 'black')
        img_draw = ImageDraw.Draw(img_frame)
        img_draw.text((3,3), 'Frame '+str(frame),  fill='gray')
        for idx, particle in df1_theta[df1_theta['frame'] == frame].iterrows():
            # Draw a circle for a particle in each frame in df1
            ellipse_bb = [tuple(particle['x pos':'y pos'].values-particle_size),
                          tuple(particle['x pos':'y pos'].values+particle_size)]
            img_draw.ellipse(ellipse_bb, outline='red')
            track_num = particle['track id']
            
        for idx, particle in df2_theta[df2_theta['frame'] == frame].iterrows():
            # Draw a circle for a particle in each frame in df2
            ellipse_bb = [tuple(particle['x pos':'y pos'].values-particle_size),
                          tuple(particle['x pos':'y pos'].values+particle_size)]
            img_draw.ellipse(ellipse_bb, outline='yellow')
            track_num = particle['track id']
        
        df1_frame = np.floor(df1_theta[df1_theta['frame'] == frame].loc[:,'x pos':'y pos'])
        df2_frame = np.floor(df2_theta[df2_theta['frame'] == frame].loc[:,'x pos':'y pos'])
        for num, (idx, row) in enumerate(df1_frame.iterrows()):
            if (row == df2_frame).all(axis=1).any(axis=0):
                ellipse_bb = [tuple(row.values-particle_size),
                          tuple(row.values+particle_size)]
                img_draw.ellipse(ellipse_bb, outline='purple')
        
        image_frame_array = np.array(img_frame)
        img_frame.close()
        image_frames.append(image_frame_array)
    
    skimage.viewer.CollectionViewer(image_frames).show()
    gc.collect()

def plot_positions_on_image(image_path, dfs_positions, save_path=None, particle_size=6.0):
    """Draws particle positions from one or more DataFrames onto raw data
    
    :param image_path: The path where the raw data is. Expect an image sequence of tifs that is ordered when glob
    is called.
    :param dfs_positions: A list of DataFrame objects that have particle positions in "frame", "x pos" and "y pos"
    that will be drawn on the raw data
    :param save_path: The path where the new images will be saved. If the directory does not exist it will be created
    :param particle_size: The size to draw each particle (in pixels) around each position.
    """
    
    try:
        os.chdir(save_path)
    except WindowsError:
        os.mkdir(save_path)
        os.chdir(save_path)
    image_list = glob.glob(image_path+"*.tif")

    trajectory_palette = ['yellow', 'firebrick', 'lime', 'cyan', 'peachpuff', 'mediumaquamarine', 'lavenderblush', 'plum', 'turquoise', 'wheat', 'palevioletred']
    for num, img_name in enumerate(image_list):
        # Load the image
        file_name = os.path.split(img_name)[-1]
        file_name = os.path.splitext(file_name)[0]
        curr_image = Image.open(img_name)
        
        # Rescale image so it shows up in a png
        img_arr = np.array(curr_image.getdata()).reshape(curr_image.size)
        min_pixel = np.min(img_arr)
        img_arr = img_arr - min_pixel
        max_pixel = np.max(img_arr)
        scale = 255/float(max_pixel)
        img_arr = np.round(img_arr * scale)
        curr_image = Image.fromarray(img_arr)
        rgb_image = curr_image.convert(mode='RGB')
        img_draw = ImageDraw.Draw(rgb_image)
        for num2,df in enumerate(dfs_positions):
            frame = df[df.frame == num+1]
            points = list(frame.loc[:,'x pos':'y pos'].values.flatten())
            img_draw.point(points, fill=trajectory_palette[num2])
            for idx, particle in frame.iterrows():
                ellipse_bb = [tuple(particle['x pos':'y pos'].values-particle_size),
                              tuple(particle['x pos':'y pos'].values+particle_size)]
                img_draw.ellipse(ellipse_bb, outline=trajectory_palette[num2])
        rgb_image.save(file_name+'.png')


def hist_bin_optimization(data, data_range=None, max_bins=None):
    '''Uses the Shimazaki histogram binwidth optimization for a given
    data set of integers and returns the optimal number of bins
    
    This function attempts to find the optimal number of bins for a 
    data set you want to histogram. The optimization method was 
    developed by Shimazaki with ref:
    
    Shimazaki and Shinomoto. Neural Comput, 2007, 19(6), 1503-1527
    
    The function is hard coded to work with integer data and constrains
    the bin width such that it cannot go below the sampling frequency
    (in this case is one).
    
    :param data: List or array of your data that you will histogram
    :param (tuple) data_range: The range that you plan on making a
    histogram of the data. You should use this range when you make your
    histogram after finding the optimal number of bins.
    :param (int) max_bins: The maximum number of bins you will allow for
    the given range. Note, if max_bins is greater than the difference of 
    data_range then the difference of data_range is used of max_bins. This
    function is hard coded to work with data that has a sampling frequency 
    equal to 1. Your max bins must be less than your sampling frequency.
    '''
    
    if data_range == None:
        data_range = [min(data), max(data)]
    if max_bins !=None and max_bins < data_range[1]-data_range[0]:
        bin_numbers = range(1,max_bins)
    else:
        bin_numbers = range(1,int(data_range[1]-data_range[0])+1)
    cost_results = []
    for bin_num in bin_numbers:
        bin_width = (data_range[1] - data_range[0])/float(bin_num)
        counts, bins = np.histogram(data, bins=bin_num, range=data_range)
        mean_counts = np.mean(counts)
        variance = np.sum((counts - mean_counts)**2)/float(bin_num)
        c = (2*mean_counts - variance)/float(bin_width)**2
        cost_results.append(c)
    cost_results = np.array(cost_results)
    idx_min_cost = np.argmin(cost_results)
    return bin_numbers[idx_min_cost]

def hist_bin_optimization_continuous(data, data_range=None, max_bins=200):
    '''Uses the Shimazaki histogram binwidth optimization for a given
    data set of continuous data and returns the optimal number of bins
    
    This function attempts to find the optimal number of bins for a 
    data set you want to histogram. The optimization method was 
    developed by Shimazaki with ref:
    
    Shimazaki and Shinomoto. Neural Comput, 2007, 19(6), 1503-1527
    
    The function should work with any continuous data set.
    
    :param data: List or array of your data that you will histogram
    :param (tuple) data_range: The range that you plan on making a
    histogram of the data. You should use this range when you make your
    histogram after finding the optimal number of bins.
    :param (int) max_bins: The maximum number of bins you will allow for
    the given range. Your max bins must be less than your sampling
    frequency. If you data looks over binned try reducing this number
    '''
    
    if data_range == None:
        data_range = [min(data), max(data)]
    bin_numbers = range(1,max_bins)
    cost_results = []
    for bin_num in bin_numbers:
        bin_width = (data_range[1] - data_range[0])/float(bin_num)
        counts, bins = np.histogram(data, bins=bin_num, range=data_range)
        mean_counts = np.mean(counts)
        variance = np.sum((counts - mean_counts)**2)/float(bin_num)
        c = (2*mean_counts - variance)/float(bin_width)**2
        cost_results.append(c)
    cost_results = np.array(cost_results)
    idx_min_cost = np.argmin(cost_results)
    return bin_numbers[idx_min_cost]

def displacement_calc(group):
    '''Calculates the displacement of the group with an x and y position. Use with
    pandas groupby and apply.

    :param group: DataFrame with 'x pos' and 'y pos' columns and contains only one 'track id'
    '''
    group_non_consec_index = group[group.frame - group.shift(1).frame > 1].index
    xy_data = group[['x pos', 'y pos']]
    xy_disp = xy_data - xy_data.shift(1)
    xy_disp.loc[group_non_consec_index] = np.nan
    xy_disp = xy_disp.dropna()
    disp = np.sqrt(np.sum((xy_disp)**2, axis=1))
    return disp

def displacement_in_theta_calc(group):
    '''Calculates the displacement of the group with a theta position. Use with
    pandas groupby and apply.

    :param group: DataFrame with 'theta' column and contains only one 'track id'
    '''
    group_non_consec_index = group[group.frame - group.shift(1).frame > 1].index
    theta_data = group['theta']
    theta_disp = theta_data - theta_data.shift(1)
    theta_disp = theta_disp - 360.0 * np.round(theta_disp/360.0)
    theta_disp.loc[group_non_consec_index] = np.nan
    theta_disp = theta_disp.dropna()
    return theta_disp

default_um_conv = 6.5/60.0/1.6/2.00
default_frame_rate = 90.00
default_radius = 0.075
default_viscosity = 1.002

def calc_velocities_consec_frames(df, frame_rate=default_frame_rate, um_conv=default_um_conv, theta_range=[120,240], only_over_plate=True):
    '''Calculates the velocity of particles in the data frame within a specified part of the ring
    trap and only counts consecutive frames (where the particle does not disappear).

    :param df: DataFrame that contains the trajectory information with keys 
    ['frame','track id','x pos','y pos', 'over_plate']
    :param frame_rate: The frame rate that the experiment was recorded at
    :param um_conv: The conversion factor to convert position from pixels to um.
    :param theta_range: Tuple (2 elements) which describe the lower and upper limits of theta that you want
    to find the velocity of particles over (e.g. the theta range of particles over the nanoplate)
    :param over_plate: If True then only particles over the nanoplate are considered. If False all particles
    velocities are calculated.
    '''
    df_temp = df.query('@theta_range[0] < theta < @theta_range[1]')
    if only_over_plate == True:
        df_temp = df_temp[df_temp['over_plate'] == True]
    df_temp = df_temp.drop_duplicates(subset=['frame', 'track id'], keep='first')
    displacements = df_temp.groupby('track id', group_keys=False).apply(displacement_calc) 
    velocities = displacements * um_conv * frame_rate
    return velocities

def calc_velocities_consec_frames_in_theta(df, frame_rate=default_frame_rate, um_conv=default_um_conv, theta_range=[120,240], r_range=[128,140], only_over_plate=True):
    '''Calculates the velocity of particles in theta from a data frame within a specified part of the ring
    trap and only counts consecutive frames (where the particle does not disappear). The average radius
    of the particles within the designated theta region is used to calculate the arc traveled
    in the ring trap.

    :param df: DataFrame that contains the trajectory information with keys 
    ['frame','track id', 'theta', 'over_plate']
    :param frame_rate: The frame rate that the experiment was recorded at
    :param um_conv: The conversion factor to convert position from pixels to um.
    :param theta_range: Tuple (2 elements) which describe the lower and upper limits of theta that you want
    to find the velocity of particles over (e.g. the theta range of particles over the nanoplate)
    :param r_range: Tuple (2 elements) which describes the upper and lower cutoff in radius for 
    allowed trajectories. This prevents outliers of particles moving great distances
    in theta in the center of the ring (not in the ring trap).
    :param over_plate: If True then only particles over the nanoplate are considered. If False all particles
    velocities are calculated.
    '''
    df_temp = df.query('@theta_range[0] < theta < @theta_range[1]')
    df_temp = df_temp.query('@r_range[0] < r < @r_range [1]')
    if only_over_plate == True:
        df_temp = df_temp[df_temp['over_plate'] == True]
    df_temp = df_temp.drop_duplicates(subset=['frame', 'track id'], keep='first')
    r_avg = df_temp.r.mean()
    displacements_theta = df_temp.groupby('track id', group_keys=False).apply(displacement_in_theta_calc) 
    velocities = displacements_theta* r_avg * (np.pi/180) * um_conv * frame_rate
    return velocities

def calc_force_stokes_drag(velocities, radius=default_radius, viscosity=default_viscosity):
    '''Calculates the drag force on a particle with a measured velocity

    :param velocities: Series of velocities to calculate the force from. Must be in um/s
    :param radius: Radius of particle. Must be in um
    :param viscosity: Viscosity of fluid. Must be in centipoise
    '''
    R = radius * 10**(-6)
    mu = viscosity * 10**(-3)
    si_velocities = velocities * 10**(-6)
    stokes_law = lambda v: 6*np.pi*mu*R*v
    force = stokes_law(si_velocities)
    return force

def check_over_plate_in_past(group):
    '''Returns the entries where particles that are over glass were just over 
    the plate last time they could be identified. This should be used with a 
    groupby on a data frame grouped by 'track id'. 
    
    :pram group: The group from the groupby operation, should be no duplicate
    entries w.r.t. frame and track id.
    '''
    new_group = group.iloc[1:]
    cur_bool = group.over_plate.iloc[1:]
    fut_bool = group.over_plate.shift(1).iloc[1:].astype(np.bool)
    return new_group[(~cur_bool & fut_bool)]

def check_over_glass_future(group):
    '''Returns the entries where particles that are over the plate will be found
    over the glass the next time they could be identified. This should be used 
    with a groupby on a data frame grouped by 'track id'. 
    
    :pram group: The group from the groupby operation, should be no duplicate
    entries w.r.t. frame and track id.
    '''
    new_group = group.iloc[:-1]
    cur_bool = group.over_plate.iloc[:-1]
    fut_bool = group.over_plate.shift(-1).iloc[:-1].astype(np.bool)
    return new_group[(cur_bool & ~fut_bool)]

def add_barrier_crossing_index(df, theta_range=(220,300)):
    '''Returns the data frame with an over_plate designation with two additional columns designating
    a particular particle's last frame over the plate and it's first frame over the glass.

    :param df: DataFrame that contains trajectory information. Must contain keys
    ['frame','track id','x pos','y pos', 'over_plate]
    :param theta_range: Tuple (2 elements) that defines the range where to consider particles leaving the
    plate or entering the glass.
    '''
    df_copy = df.copy()
    df_subset = df_copy.query('@theta_range[0] < theta < @theta_range[1]')
    df_subset = df_subset.drop_duplicates(subset=['frame', 'track id'], keep='first')
    
    # I call a pd.merge in the next two blocks in order to make all the duplicate entries in the original
    # DataFrame also be correctly marked True for last_frame_plate and first_frame_glass

    # Find the a transitioning particles' last frame on plate
    leaving_plate = df_subset.groupby('track id', group_keys=False).apply(check_over_glass_future)
    index_leaving_plate = pd.merge(df_copy, leaving_plate, on=['frame', 'track id'], right_index=True).index
    df_copy['last_frame_plate'] = False
    df_copy.loc[index_leaving_plate, 'last_frame_plate'] = True
    
    # Find the a transitioning particles' first frame on glass
    entering_glass = df_subset.groupby('track id', group_keys=False).apply(check_over_plate_in_past)
    index_entering_glass = pd.merge(df_copy, entering_glass, on=['frame', 'track id'], right_index=True).index
    df_copy['first_frame_glass'] = False
    df_copy.loc[index_entering_glass, 'first_frame_glass'] = True
    
    return df_copy

def trajectory_dwell_time_to_index(df, index, max_frames_absent=1):
    '''Finds the dwell time (in frames) for a track up until a selected 
    index. This will look for the longest consecutive track as long as the
    doesn't leave the field of view for some number of frames.
    
    :param df: The DataFrame that is already filtered for your dwell time
    criteria. That is, if I want to calculated the dwell time within a certain
    area, the DataFrame I input here would only contain positions in that area.
    :param index: The index number for the last position that I want to
    calculated the dwell time up to. This is the index where the end of the
    dwell time is.
    :param max_frames_absent: The number of frames a specific track is allowed
    to be absent and still count towards the dwell time. If the particle
    disappears for longer than max_frames_absent then the point where it
    disappears will be considered the first frame for the dwell time
    calculation.
    '''
    
    final_frame = df.ix[index]['frame']
    track_id = df.ix[index]['track id']
    track_df = df[(df['track id'] == track_id) & (df['frame'] <= final_frame)]
    frame_diff = track_df['frame'] - track_df['frame'].shift(1)
    valid_diffs = frame_diff[frame_diff > max_frames_absent]
    if len(valid_diffs) == 0:
        first_frame = track_df.iloc[0]['frame']
    else:
        first_frame = track_df.ix[valid_diffs.index[-1]]['frame']
    return final_frame-first_frame + 1

def check_exiting_region(group):
    '''Returns the entries where particles that are in the region will be found
    outside the region the next time they could be identified. This should be used 
    with a groupby on a data frame grouped by 'track id'. 
    
    :pram group: The group from the groupby operation, should be no duplicate
    entries w.r.t. frame and track id. Must have a boolean column 'in_region'.
    '''
    new_group = group.iloc[:-1]
    cur_bool = group.in_region.iloc[:-1]
    fut_bool = group.in_region.shift(-1).iloc[:-1].astype(np.bool)
    return new_group[(cur_bool & ~fut_bool)]

def check_entering_region(group):
    '''Returns the entries where particles that are in the region that were just 
    outside the region last time they could be identified. This should be used 
    with a groupby on a data frame grouped by 'track id'. 
    
    :pram group: The group from the groupby operation, should be no duplicate
    entries w.r.t. frame and track id. Must have a boolean column 'in_region'
    '''
    new_group = group.iloc[1:]
    cur_bool = group.in_region.iloc[1:]
    past_bool = group.in_region.shift(1).iloc[1:].astype(np.bool)
    return new_group[(cur_bool & ~past_bool)]

def find_longest_trajs_in_group_in_region(group):
    '''Returns only the entries in the data frame that represent the longest
    trajectory in a region for each 'track id'. This should be used with 
    groupby on a DataFrame grouped by 'track id' with boolean values for
    particles entering and leaving the region.
    
    :param group: A 'track id' group from the groupby operation. Should have
    no duplicate entries w.r.t. frame and 'track id'. Must have boolean 
    columns 'first_frame_region' and 'last_frame_region'
    '''
    # Find trajectories that enter and exit region
    ent_group = group.first_frame_region == True
    exit_group = group.last_frame_region == True
    ent_exit_group = group[ent_group | exit_group]
    
    # Return if particle doesn't have both and entrance and exit
    if len(ent_exit_group) == 1:
        return pd.DataFrame(columns=group.columns)
    
    # Find start index of longest traj
    frames = ent_exit_group.shift(-1).frame - ent_exit_group.frame
    if len(frames) == 0:
        return pd.DataFrame(columns=group.columns)
    start = frames.idxmax()
    if group.loc[start, 'last_frame_region'] == True:
        return pd.DataFrame(columns=group.columns)
    
    # Find end index of longest traj
    frames = ent_exit_group.frame - ent_exit_group.shift(1).frame
    if len(frames) == 0:
        return pd.DataFrame(columns=group.columns)
    end = frames.idxmax()
    if group.loc[end, 'first_frame_region'] == True:
        return pd.DataFrame(columns=group.columns)
    
    return group.ix[start:end]


def find_longest_trajs_in_region(df, theta_range=(120,240)):
    '''Returns a data frame of trajectories of only the longest consecutive 
    trajectories in a defined region.
    
    This function accepts a DataFrame of rotational trajectories and allows you
    to find the longest consecutive trajectories in a selected region. First 
    particle trajectories are cut to just the ones that are in the region. Next,
    for each trajectory (track id) that passes through the region only the longest
    one is kept. The trajectories that are returned must move from outside the
    region into it and also must move from inside the region to outside of it in
    order for the trajectory to be considered. Having these critera prevents 
    counting trajectories of particles that 'appear' in the region (and not seen
    outside first) or particles that 'disappear' out fo the region (and not seen
    again). This should give trajectories that are continuous and move all the 
    way through the region.
    
    Note: This function does not respect periodic boundary conditions. You 
    cannot use a theta range that spans across 0 degrees.
    
    :param df: DataFrame of trajectories to be analyzed. Does not need to have
    unique values w.r.t. 'frame' and 'track id'.
    :param (tuple) theta_range: Tuple (2 elements) wich describes the region you
    will be selecting the longest trajectories from.
    '''
    # Designate if particle positions are in the region with 'in_region" column
    df_unique = df.copy().drop_duplicates(subset=['frame', 'track id'], keep='first')
    df_in_region = (theta_range[0] < df_unique.theta) & (df_unique.theta < theta_range[1])
    df_unique.loc[df_in_region, 'in_region'] = True
    df_unique.loc[~df_in_region, 'in_region'] = False
    
    # Find positions where particles just enter or exit the region
    entering_region = df_unique.groupby('track id', group_keys=False).apply(check_entering_region)
    exiting_region = df_unique.groupby('track id', group_keys=False).apply(check_exiting_region)
    
    # Designate which positions are the first/last frame in the region
    df_unique['first_frame_region'] = False
    df_unique['last_frame_region'] = False
    df_unique.loc[entering_region.index, 'first_frame_region'] = True
    df_unique.loc[exiting_region.index, 'last_frame_region'] = True
    
    # Determine the longest consecutive trajectory in the region
    df_long_traj_in_region = df_unique.groupby('track id', group_keys=False).apply(find_longest_trajs_in_group_in_region)
    return df_long_traj_in_region