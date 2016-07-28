import string
import scipy
import Tkinter, tkFileDialog
import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist
import matplotlib.pyplot as plt
import os
import sys
import re
import skimage.viewer
from PIL import Image, ImageDraw, ImageColor
import trackpy
import trackpy.predict
import gc


def import_matlab_gui(file_path=None):
    '''This function will import data files that are generated from
	Rahgu's tracking GUI. Pick a matlab file and the output will be
	an array'''
    import scipy.io
    if file_path == None:
        import Tkinter, tkFileDialog
        root = Tkinter.Tk()
        root.withdraw()
        file_path = tkFileDialog.askopenfilename(filetypes=[('matlab files', '.mat')])
    print 'Loading ' + file_path
    m = scipy.io.loadmat(str(file_path))['objs_link']
    if len(m) == 0:
        m = scipy.io.loadmat(str(file_path))['objs']
    return m


def matlab_gui_to_data_frame(m):
    import pandas
    return pandas.DataFrame(
        {'x pos': m[0][:], 'y pos': m[1][:], 'frame': m[4][:], 'particle id': m[3][:], 'track id': m[5][:]})


def find_nn(grp):
    '''Used to find the nearest neighbors in a data frame generated
	from matlab_gui_to_data_frame. This adds both the nearest
	neighbor track ids and the distances to each particle in each
	frame. Use this by using df.groupby('frame').apply(find_nn)
	where df is the data frame generated from matlab_gui_to_data_frame.
	Set the output to what you want to call your new data frame (or write
	over your old variable name).
	
	Note: This version of the function works when the input data frame
	lacks track id information. That is, one does not know anything
	about where particles go from frame to frame.'''
    tree = scipy.spatial.KDTree(grp[['x pos', 'y pos']])
    nn = tree.query(grp[['x pos', 'y pos']], k=12)
    reindex_value = grp.index[0]
    grp['nn_part'] = [grp['track id'][v].values for v in nn[1] + reindex_value]
    grp['nn_dist'] = [v for v in nn[0]]
    return grp

def find_nn_ver_2(grp):
    '''This function will find all the nearest neighbors in x-y and add 
    columns to the data frame that include the nn number (1st, 2nd,...), the
    particle id of the nn, and the distance of that nn.
    The way to use this function is as such:
    
    df.groupby('frame',group_keys=False).apply(find_nn).reset_index()
    df=df[['frame','track id','x pos','y pos','nn_num','nn_id','nn_dist','theta','r']]
    
    The group_keys kwarg prevents a redundant frames column. Reseting the index
    will give the data frame a regular integer index like before the function
    is applied. The second line rearranges the columns to the correct order.
    '''
    xy_data=grp[['x pos','y pos']]
    tree=scipy.spatial.KDTree(xy_data)
    nn_d, nn_i=tree.query(xy_data, k=len(xy_data))
    if len(nn_d)==1: # If only one particle return the group
        #grp['nn_num'],grp['nn_id'],grp['nn_dist']=[np.nan,np.nan,np.nan]
        return grp
    particle_ids=grp['track id'].values
    track_ids=np.tile(particle_ids, (len(particle_ids),1))
    track_ids=track_ids.T[:,1:].flatten()
    nn_ids=particle_ids[nn_i]
    # Create nn number column (1st, 2nd, etc)
    nn_num=np.arange(len(particle_ids))
    nn_num=np.tile(nn_num,(len(particle_ids),1))[:,1:]
    nn_num=nn_num.flatten()
    # Create corresponding nn track id
    nn_ids=nn_ids[:,1:].flatten()
    nn_dist=nn_d[:,1:].flatten()
    # Merge with current group
    nn_df=pd.DataFrame(np.vstack((track_ids,nn_num,nn_ids,nn_dist)).T, columns=['track id','nn_num','nn_id','nn_dist'])
    new_df=pd.merge(grp,nn_df, left_on='track id', right_on='track id')
    return new_df

def restructure_nn_data_frame(data_frame):
    '''Restructures the data frame after using the function find_nn such that 
    each nn has its on line in the dataframe. If there are no nearest neighbors
    the first one is shown as NaN particle id and inf distance. This function
    returns the new dataframe with a reset index.
    '''
    index_list = data_frame.index
    nn_part = np.concatenate(data_frame['nn_part'].values)
    nn_dist = np.concatenate(data_frame['nn_dist'].values)
    nn_num = np.tile(np.arange(12), len(nn_part) / len(data_frame['nn_part'].iloc[0]))
    index_new_nn_part = np.repeat(data_frame.index, len(data_frame['nn_part'].iloc[0]))
    new_nn_part = pd.DataFrame(np.vstack([nn_num, nn_part, nn_dist]).T,
                               columns=['nn_num2', 'nn_part2', 'nn_dist2'],
                               index=index_new_nn_part)
    # Filter out nn_num=0
    new_nn_part = new_nn_part.query('nn_num2>0')
    # Filter all NaNs and infs above nn_num=2
    new_nn_part = new_nn_part[(new_nn_part.nn_num2 >= 1) &
                              -((new_nn_part.nn_num2 >= 2) &
                                (new_nn_part.nn_part2.isnull()) &
                                (new_nn_part.nn_dist2 == np.inf))]
    return data_frame.join(new_nn_part).reset_index(drop=True)


def plot_particle_positions(data_frame):
    '''If the particle data is in a data frame (from matlab_gui_to_data_frame)
	then this function will plot the position of all the particles.
	Keep in mind that the y axis is flipped after images are run through
	tracking software'''
    plt.plot(data_frame['x pos'], data_frame['y pos'], linestyle='', marker='o', markerfacecolor='b')
    plt.axis('equal')
    plt.xlabel('Pixels')
    plt.ylabel('Pixels')

def polar_transform_points(data_frame, origin=None):
    '''Polar transforms the particle positions in all frames around
	a defined center (origin).'''
    xy = data_frame[['x pos', 'y pos']].values
    if origin is None:
        origin = np.mean(xy, axis=0)
    output = np.zeros(xy.shape)
    # This finds the radius of each point fromt the center
    output[:, 0] = np.sqrt((origin[0] - xy[:, 0]) ** 2 + (origin[1] - xy[:, 1]) ** 2)
    # This finds the angle of each point starting with 0 at 12 o'clock
    # and going counter clockwise
    output[:, 1] = 180 + 180 * np.arctan2(xy[:, 0] - origin[0], xy[:, 1] - origin[1]) / np.pi
    return output


def y_axis_flip(data_frame, img_height):
    '''The Raghu code spits out particle positions with the a shift
	of the origin from the top left to the bottom left resulting in
	a flip in the y-axis for all the data. This function flips it 
	back to values that match the image. It requires knowledge of
	the height of the original data.
	'''
    temp_df = data_frame.copy()
    temp_df['y pos'] = abs(temp_df['y pos'] - img_height)
    return temp_df


def least_sq_fit_circle(data_frame):
    '''Calculate the best fit circle by least squares for the particle
	locations for all frames in a data frame. This funtion finds
	the center via center of mass'''
    from scipy.spatial.distance import cdist
    import scipy.optimize
    xy = data_frame[['x pos', 'y pos']].values
    # Find center via center of mass
    center = np.mean(xy, axis=0).reshape(1, 2)
    radius = cdist(center, xy).mean()
    # Define the error function to minimize
    def err_function((x_center, y_center, radius)):
        err = np.sqrt((xy[:, 0] - x_center) ** 2 + (xy[:, 1] - y_center) ** 2) - radius
        return (err ** 2).sum()

    xf, yf, rf = scipy.optimize.fmin(err_function, [center[0, 0], center[0, 1], radius], disp=False)
    return xf, yf, rf


def filter_data_frame_by_radius(data_frame, x_cent, y_cent, r_in, r_out):
    '''A simple function that will create a data frame that will
	exclude a set of points not within two set radii'''
    # outer_box = data_frame[['x pos','y pos']].applymap(lambda x: cdist([x_cent,y_cent],[x,x])>r_out)
    outer_box = np.sqrt((data_frame['x pos'] - x_cent) ** 2 + (data_frame['y pos'] - y_cent) ** 2) < r_out
    inner_box = np.sqrt((data_frame['x pos'] - x_cent) ** 2 + (data_frame['y pos'] - y_cent) ** 2) > r_in
    return data_frame[outer_box & inner_box]


def calc_angle(x_pos, y_pos, x_cent, y_cent):
    '''This function calculates the angle of the positions with respect
	to the centers x_cent and y_cent. This function is designed so that
	the origin is at the 12 o'clock position and increases going counter-
	clockwise
	'''
    return 180 + 180 * np.arctan2(x_pos - x_cent, y_cent - y_pos) / np.pi


def calc_radius(x_pos, y_pos, x_cent, y_cent):
    '''This function finds the radius in polar coordinates with respect to
	the centers x_cent and y_cent. This function always produces a positive
	radius.'''
    return np.sqrt((x_pos - x_cent) ** 2 + (y_pos - y_cent) ** 2)


def polar_coor_data_frame(data_frame, x_cent, y_cent):
    '''This function takes a data frame and adds columns for the radius and theta of
	each position based on x_cent and y_cent of the image'''
    data_frame['theta'] = calc_angle(data_frame['x pos'], data_frame['y pos'], x_cent, y_cent)
    data_frame['r'] = calc_radius(data_frame['x pos'], data_frame['y pos'], x_cent, y_cent)

def calc_x_from_polar(r, theta, x_cent):
    '''This function will calculate the x coordinate from r and theta. This 
    will reverse my custom polar transform'''
    return x_cent + r * np.sin((theta-180.0)*np.pi/180.0)

def calc_y_from_polar(r, theta, y_cent):
    '''This function will calculate the y coordinate from r and theta. This 
    will reverse my custom polar transform'''
    return y_cent - r * np.cos((theta-180.0)*np.pi/180.0)

def calc_cent_from_polar(x, y, r, theta):
    '''This function will calculate the center of the fitted circle of my custom
    polar transform'''
    x_cent = x + r * np.sin(theta * np.pi / 180.0)
    y_cent = y - r * np.cos(theta * np.pi / 180.0)
    return x_cent, y_cent

def nn_distance_angle_seperation(data_frame, number_of_bins, x_cent, y_cent):
    '''Seperates the nearest neighbor distances into angular bins based
	on the peak location. This function only works if you are picking peaks
	from the images and not using tracked data.
	
	:pram peaks: List of peaks from func find_peak_locations
	:pram nn_dist: List of nearest neighbors from func find_nn_peaks
	:pram number_of_bins: int.
	'''
    if 360 % number_of_bins != 0:
        print "Error: Number of bins must be divide 360 without remainders"
        return None
    if number_of_bins / 2 % 2 != 0:
        print "Error: Half the number of bins must be even"
        return None
    bin_size = 360 / number_of_bins
    bin_limits = []
    for i in range(number_of_bins):
        if i == 0:
            bin_limits.append([0 + bin_size / 2.0, 360 - bin_size / 2.0])
        else:
            bin_limits.append([(i * bin_size) + bin_size / 2.0, (i * bin_size) - bin_size / 2.0])
    radial_bins = []
    for step in range(len(bin_limits)):
        # Define the less than and greater than criterea
        less_than = calc_angle(data_frame['x pos'], data_frame['y pos'], x_cent, y_cent) <= bin_limits[step][0]
        greater_than = calc_angle(data_frame['x pos'], data_frame['y pos'], x_cent, y_cent) >= bin_limits[step][1]
        if step == 0:
            radial_bins.append(np.concatenate([[v[1]] for v in data_frame[less_than | greater_than]['nn_dist'].values]))
        # radial_bins.append(np.concatenate([nn_dist[v][np.logical_or(peaks[v]<=bin_limits[step][0],peaks[v]>=bin_limits[step][1])] for v in range(len(peaks))]))
        else:
            radial_bins.append(np.concatenate([[v[1]] for v in data_frame[less_than & greater_than]['nn_dist'].values]))
            # radial_bins.append(np.concatenate([nn_dist[v][np.logical_and(peaks[v]<=bin_limits[step][0],peaks[v]>=bin_limits[step][1])] for v in range(len(peaks))]))
    return radial_bins, bin_limits


def plot_polar_transform(polar_transform_output):
    '''This function plots the polar transform from the output from
	polar_transform_points'''
    plt.scatter(polar_transform_output[:, 0], polar_transform_output[:, 1])
    plt.xlabel('Radius (px)')
    plt.ylabel('Degrees')
    # plt.axis('equal')
    plt.show()


def import_diatrack_txt(file_path=None):
    '''This function imports the data from a Diatrack text
	file into a pandas dataframe'''
    if file_path == None:
        root = Tkinter.Tk()
        root.withdraw()
        file_path = tkFileDialog.askopenfilename(filetypes=[('diatrack files', '.txt')])
    track_id = 0
    imported_data = []
    for line in open(file_path, 'r').readlines():
        if line.startswith('f'):
            continue
        else:
            track_id += 1
            templine = string.split(line)
            templine = [float(i) for i in templine]
            frame = int(templine[0])
            for pos in range(1, len(templine) - 1, 3):
                imported_data.append([frame, track_id, templine[pos], templine[pos + 1]])
                frame += 1
    imported_data = pd.DataFrame(imported_data, columns=['frame', 'track id', 'x pos', 'y pos'])
    imported_data = imported_data.sort(['frame', 'track id'])
    imported_data.index = range(len(imported_data))
    print file_path
    return imported_data, file_path


def find_longest_traj(data_frame):
    '''In a data frame created from import_matlab_gui() or
	import_diatrack_txt and will return a data frame that contains only
	the longest trajectory 'track id' values'''
    grouped = data_frame.groupby('track id')
    max_len, max_name = [0, 0]
    for name, grp in grouped:
        if len(grp) > max_len:
            max_name = name
            max_len = len(grp)
        else:
            continue
    return data_frame[data_frame['track id'] == max_name]


def import_mosaic_trajectories(file_path=None):
    '''This function imports the data from the Mosaic trajectories table file
	into a pandas dataframe'''
    if file_path == None:
        root = Tkinter.Tk()
        root.withdraw()
        file_path = tkFileDialog.askopenfilename(filetypes=[('mosaic trajectories', '.xls')])
    imported_data = pd.read_csv(file_path, delim_whitespace=True, usecols=[2, 1, 3, 4])
    imported_data.columns = ['track id', 'frame', 'y pos', 'x pos']
    imported_data = imported_data[['frame', 'track id', 'x pos', 'y pos']]
    imported_data = imported_data.sort(['frame', 'track id']).reset_index(drop=True)
    imported_data['frame'] = imported_data['frame'] + 1
    print file_path
    return imported_data, file_path


def find_nn_theta(grp):
    '''This function will find all the nearest neighbors for theta and add 
    columns to the data frame that include the nn number (1st, 2nd,...), the
    particle id of the nn, and the distance of that nn. The function respects
    boundary conditions of the data such that it is periodic at 360 degrees.
    The way to use this function is as such:
    
    df=df.groupby('frame', group_keys=False).apply(find_nn_theta).reset_index()
    df=df[['frame','track id','x pos','y pos','nn_part','nn_dist','theta','r','passing_event','theta_nn_num','theta_nn_id','theta_nn_dist']]
    
    The group_keys kwrg prevents a redundant frames column. Reseting the index
    will give the data frame a regular integer index like before the function
    is applied. The second line rearranges the columns to the correct order.
    '''
    from periodic_kdtree import PeriodicCKDTree
    bounds = np.array([360])
    data = grp['theta'].values
    data = np.reshape(data, [len(data), 1])
    tree = PeriodicCKDTree(bounds, data)
    d, i = tree.query(data, k=len(data))
    if len(d) == 1:  # If only one particle return the group
        # grp['theta_nn_num'],grp['theta_nn_id'],grp['theta_nn_dist']=[np.nan,np.nan,np.nan]
        return grp
    # Create particle id column
    particle_ids = grp['track id'].values
    track_ids = np.tile(particle_ids, (len(particle_ids), 1))
    track_ids = track_ids.T[:, 1:].flatten()
    nn_ids = particle_ids[i]
    # Create nn number column (1st, 2nd, etc)
    nn_num = np.arange(len(particle_ids))
    nn_num = np.tile(nn_num, (len(particle_ids), 1))[:, 1:]
    nn_num = nn_num.flatten()
    # Create corresponding nn track id
    nn_ids = nn_ids[:, 1:].flatten()
    nn_dist = d[:, 1:].flatten()
    # Merge with current group
    nn_df = pd.DataFrame(np.vstack((track_ids, nn_num, nn_ids, nn_dist)).T,
                         columns=['track id', 'theta_nn_num', 'theta_nn_id', 'theta_nn_dist'])
    new_df = pd.merge(grp, nn_df, left_on='track id', right_on='track id')
    return new_df


def nn_distance_angle_seperation_ver_2(data_frame, number_of_bins, nn_num_limit):
    '''Seperates the nearest neighbor distances into angular bins based
    on their theta values and their nn_dist. This function only works 
    if you use the find_nn_ver_2 where each nn has its own line in the data_frame.
    
    :pram data_frame: data_frame of the particle tracks you want to seperate
    :pram number_of_bins: number of angular bins you want
    :pram (int) nn_num_limit: Up to which nn_num you want to find the distances for 
    '''
    # Make sure the number of bins is even and divides 360 w/0 remainders
    if 360 % number_of_bins != 0:
        print "Error: Number of bins must be divide 360 without remainders"
        return None
    if number_of_bins / 2 % 2 != 0:
        print "Error: Half the number of bins must be even"
        return None

    # Create the bin limits with bin centers at 12, 3, 6, and 9 o'clock
    bin_size = 360 / number_of_bins
    bin_limits = []
    for i in range(number_of_bins):
        if i == 0:
            bin_limits.append([0 + bin_size / 2.0, 360 - bin_size / 2.0])
        else:
            bin_limits.append([(i * bin_size) + bin_size / 2.0, (i * bin_size) - bin_size / 2.0])
    radial_bins = []

    # Find the theta values of each of the nearest neighbors. This is done by reindexing w.r.t. frame, nn_id
    # track id and pick from that index (frame, track id, nn_id). This is used to check that the nearest
    # neighbor is in the same theta bin as the other particle
    nn_theta_data_frame = data_frame.copy()
    theta_nn_ids = nn_theta_data_frame[['frame', 'nn_id', 'track id']]
    frame_particle_id_reindex = nn_theta_data_frame.set_index(['frame', 'track id', 'nn_id'])
    theta_values_nn = frame_particle_id_reindex.ix[[tuple(v) for v in theta_nn_ids.values]]['theta'].values
    nn_theta_data_frame['nn_theta_val'] = theta_values_nn

    for step in range(len(bin_limits)):
        if step == 0:
            # Need a special step for the first bin in order to satisfy the 360 to 0 transition.
            # Find the particles that have a theta in the first bin, also drop NaNs for nn_id
            valid_positions = nn_theta_data_frame[((nn_theta_data_frame.theta <= bin_limits[step][0]) |
                                                   (nn_theta_data_frame.theta >= bin_limits[step][1])) &
                                                  -(nn_theta_data_frame.nn_id == np.nan)]
            # Drop duplicate nearest neighbor pairs (using i<j indexing)
            valid_positions = valid_positions[valid_positions['track id'] < valid_positions.nn_id]
            # Choose only nearest neighbors in less than or equal to nn_num_limit (e.g. nn_num_limit=1 then
            # only use first nearest neighbor)
            valid_positions = valid_positions[valid_positions.nn_num <= nn_num_limit]
            # Only use nearest neighbors that have both particles in the same bin.
            valid_positions = valid_positions[(valid_positions.nn_theta_val <= bin_limits[step][0]) |
                                              (valid_positions.theta >= bin_limits[step][1])]
            radial_bins.append(valid_positions.nn_dist.values)
        else:
            # Every other angular bin uses this block of code. The processes is the same as above
            valid_positions = nn_theta_data_frame[((nn_theta_data_frame.theta <= bin_limits[step][0]) &
                                                   (nn_theta_data_frame.theta >= bin_limits[step][1])) &
                                                  -(nn_theta_data_frame.nn_id == np.nan)]
            valid_positions = valid_positions[valid_positions['track id'] < valid_positions.nn_id]
            valid_positions = valid_positions[valid_positions.nn_num <= nn_num_limit]
            valid_positions = valid_positions[(valid_positions.nn_theta_val <= bin_limits[step][0]) &
                                              (valid_positions.theta >= bin_limits[step][1])]
            radial_bins.append(valid_positions.nn_dist.values)
    return radial_bins, bin_limits


def restructure_nn_data_frame(data_frame):
    '''Restructures the data frame after using the function find_nn such that 
    each nn has its on line in the dataframe. If there are no nearest neighbors
    the first one is shown as NaN particle id and inf distance. This function
    returns the new dataframe with a reset index.
    '''
    index_list = data_frame.index
    nn_part = np.concatenate(data_frame['nn_part'].values)
    nn_dist = np.concatenate(data_frame['nn_dist'].values)
    nn_num = np.tile(np.arange(12), len(nn_part) / len(data_frame['nn_part'].iloc[0]))
    index_new_nn_part = np.repeat(data_frame.index.values, len(data_frame['nn_part'].iloc[0]))
    new_nn_part = pd.DataFrame(np.vstack([nn_num, nn_part, nn_dist]).T,
                               columns=['nn_num2', 'nn_part2', 'nn_dist2'],
                               index=index_new_nn_part)
    # Filter out nn_num=0
    new_nn_part = new_nn_part.query('nn_num2>0')
    # Filter all NaNs and infs above nn_num=2
    new_nn_part = new_nn_part[(new_nn_part.nn_num2 >= 1) &
                              -((new_nn_part.nn_num2 >= 2) &
                                (new_nn_part.nn_part2.isnull()) &
                                (new_nn_part.nn_dist2 == np.inf))]
    return data_frame.join(new_nn_part).reset_index(drop=True)


def matlab_gui_to_mosaic_text_files(file_path=None, save_directory=None):
    '''
    Turns .mat particle postions from the Raghu code into individual files that can be imported
    into the Mosaic plugin for linking.
    
    The funciton takes a file path for the .mat file you want to open and takes a file path
    (which will become a directory) to save the individual files to
    
    :params file_path: The file path of the .mat file. If none, can pick with GUI
    :params save_directory: File path to save the files. A directory will be created 
    with this name in the parent directory. If none, can pick with GUI
    '''
    if file_path == None:
        root = Tkinter.Tk()
        root.focus()
        file_path = tkFileDialog.askopenfilename(parent=root,
                                                 initialdir="C:\Users\Scherer Lab E\Downloads\TrackingGUI_and_associated_files_20July2014 My Version")
        root.destroy()
    matlab = import_matlab_gui(file_path)
    data = matlab_gui_to_data_frame(matlab)
    if save_directory == None:
        root = Tkinter.Tk()
        root.focus()
        save_directory = tkFileDialog.asksaveasfilename(parent=root,
                                                        initialfile=os.path.splitext(os.path.split(file_path)[1])[0],
                                                        initialdir=os.path.split(file_path)[0])
        root.destroy()
        try:
            os.makedirs(save_directory)
        except:
            print 'error'
    for name, group in data.groupby('frame'):
        data_file = open(os.path.join(save_directory, 'frame_' + str(int(name) - 1)), 'w')
        data_file.write('frame ' + str(int(name) - 1) + '\n')
        for particle in range(len(group)):
            data_file.write("{0:.6f}".format(group.iloc[particle]['x pos']) + ' '
                            + "{0:.6f}".format(group.iloc[particle]['y pos']) + ' 0.000000\n')
        data_file.close()


def pandas_data_frame_to_mosaic_text_files(file_path=None, save_directory=None):
    '''
    Turns .pandas particle postions in a pickled data frame into individual files that can be imported
    into the Mosaic plugin for linking.
    
    The funciton takes a file path for the .pandas file you want to open and takes a file path
    (which will become a directory) to save the individual files to
    
    :params file_path: The file path of the .pandas file. If none, can pick with GUI
    :params save_directory: File path to save the files. A directory will be created 
    with this name in the parent directory. If none, can pick with GUI
    '''
    if file_path == None:
        root = Tkinter.Tk()
        root.focus()
        file_path = tkFileDialog.askopenfilename(parent=root,
                                                 initialdir="C:\Users\Scherer Lab E\Downloads\TrackingGUI_and_associated_files_20July2014 My Version")
        root.destroy()
    data = pd.read_pickle(file_path)
    if save_directory == None:
        root = Tkinter.Tk()
        root.focus()
        save_directory = tkFileDialog.asksaveasfilename(parent=root,
                                                        initialfile=os.path.splitext(os.path.split(file_path)[1])[0],
                                                        initialdir=os.path.split(file_path)[0])
        root.destroy()
    try:
        os.makedirs(save_directory)
    except:
        print 'error'
    cur_frame = 0
    for name, group in data.groupby('frame'):
        if name - 1 != cur_frame:
            while name - 1 != cur_frame:
                data_file = open(os.path.join(save_directory, 'frame_' + str(int(cur_frame))), 'w')
                data_file.write('frame ' + str(int(cur_frame)) + '\n')
                data_file.close()
                cur_frame += 1
        data_file = open(os.path.join(save_directory, 'frame_' + str(int(name) - 1)), 'w')
        data_file.write('frame ' + str(int(name) - 1) + '\n')
        for particle in range(len(group)):
            data_file.write("{0:.6f}".format(group.iloc[particle]['x pos']) + ' '
                            + "{0:.6f}".format(group.iloc[particle]['y pos']) + ' 0.000000\n')
        cur_frame = name
        data_file.close()

def view_trajectories(data_frame, particle_size=6.0, tail_length=10, image_size=(380,380)):
    """Visualize all the trajectories in the data frame
    
    :param data_frame: the data_frame that contains the trajectory information
    :param particle_size: the size of each particle to draw (in pixels)
    :param tail_length: the number of previous points to draw the trajectory tail to
    :param image_size: Tuple containing the size of the original image
    """
    
    image_frames = []
    traj_palette_count = 0
    # Colors to use for Trajectories
    trajectory_palette = ['yellow', 'firebrick', 'lime', 'cyan', 'peachpuff', 'mediumaquamarine', 'lavenderblush', 'plum', 'turquoise', 'wheat', 'palevioletred']
    track_color = {}
    data_frame[['frame','track id']] = data_frame[['frame','track id']].astype(int)
    re_index_df = data_frame.set_index(['track id', 'frame'])
    for frame, grp in data_frame.groupby('frame'):
        img_frame = Image.new('RGB', image_size, 'black')
        img_draw = ImageDraw.Draw(img_frame)
        for idx, particle in grp.iterrows():
            
            # Draw a circle for a particle in each frame
            ellipse_bb = [tuple(particle['x pos':'y pos'].values-particle_size),
                          tuple(particle['x pos':'y pos'].values+particle_size)]
            img_draw.ellipse(ellipse_bb, outline='red')
            track_num = particle['track id']
            
            # Determine to color to use for the trajectory
            try:
                track_color[str(track_num)]
            except KeyError:
                track_color[str(track_num)] = trajectory_palette[traj_palette_count % len(trajectory_palette)]
                traj_palette_count += 1
            
            # Determine the line segments to draw based on tail length
            line_segments = data_frame[(data_frame['track id']==track_num) & 
                                       ((int(frame)-tail_length <= data_frame['frame']) &
                                        (data_frame['frame'] <= frame))]
            line_segments = line_segments[['x pos','y pos']].values/1 
            line_segments = line_segments.flatten()
            if len(line_segments)/2 == 1:
                img_draw.point(list(line_segments), fill=track_color[str(track_num)])
            elif len(line_segments)/2 > 1:
                img_draw.line(list(line_segments), fill=track_color[str(track_num)])
        
        # Convert image to array and display it
        image_frame_array = np.array(img_frame)
        img_frame.close()
        image_frame_array = image_frame_array/255.0
        image_frames.append(image_frame_array)
        
    skimage.viewer.CollectionViewer(image_frames).show()
    gc.collect()
    
def view_trajectories_new_particles(data_frame, particle_size=6.0, frame_window=5, tail_length=10, image_size=(380,380)):
    """Visualize all the trajectories in the data frame where another particle appears
    or disappears.
    
    :param data_frame: the data_frame that contains the trajectory information
    :param particle_size: the size of each particle to draw (in pixels)
    :param frame_window: The number of frames to display before and after a particle change event
    :param tail_length: the number of previous points to draw the trajectory tail to
    :param image_size: Tuple containing the size of the original image
    """
    image_frames = []
    traj_palette_count = 0
    trajectory_palette = ['yellow', 'firebrick', 'lime', 'cyan', 'peachpuff', 'mediumaquamarine', 'lavenderblush', 'plum', 'turquoise', 'wheat', 'palevioletred']
    track_color = {}
    prev_tracks = set()
    frame_check = -1
    for frame, grp in data_frame.groupby('frame'):
        # Check if a particle change occurred
        current_tracks = set(list(grp['track id']))
        if prev_tracks != current_tracks and frame != frame_check+1:
            changing_particle =  current_tracks - prev_tracks
            # if len(prev_tracks) < len(current_tracks):
            #     changing_particle =  current_tracks - prev_tracks
            if len(prev_tracks) > len(current_tracks):
                changing_particle =  prev_tracks - current_tracks
            prev_tracks = set(grp['track id'])
            frame_check = frame
            # Determine which frames to draw
            frames_to_draw = data_frame[((int(frame)-frame_window <= data_frame['frame']) &
                                        (data_frame['frame'] <= frame_window + frame))]
            new_particle_group = []
            # Generate an image for each frame in the group
            for frame_draw, frame_grp in frames_to_draw.groupby('frame'):
                img_frame = Image.new('RGB', image_size, 'black')
                img_draw = ImageDraw.Draw(img_frame)
                img_draw.text((3,3), 'Frame '+str(frame_draw),  fill='gray')
                # For each particle in the frame, draw the particle and trajectory tail
                for idx, particle in frame_grp.iterrows():
                    if particle['track id'] in changing_particle:
                        ellipse_bb = [tuple(particle['x pos':'y pos'].values-particle_size),
                                      tuple(particle['x pos':'y pos'].values+particle_size)]
                        img_draw.ellipse(ellipse_bb, outline='green')
                    else:
                        ellipse_bb = [tuple(particle['x pos':'y pos'].values-particle_size),
                                      tuple(particle['x pos':'y pos'].values+particle_size)]
                        img_draw.ellipse(ellipse_bb, outline='red')
                    track_num = particle['track id']
                    try:
                        track_color[str(track_num)]
                    except KeyError:
                        track_color[str(track_num)] = trajectory_palette[traj_palette_count % len(trajectory_palette)]
                        traj_palette_count += 1
                    line_segments = data_frame.loc[(data_frame['track id']==track_num) & 
                                               ((int(frame_draw)-tail_length <= data_frame['frame']) &
                                                (data_frame['frame'] <= frame_draw))]
                    line_segments = line_segments[['x pos', 'y pos']].values/1 
                    line_segments = line_segments.flatten()
                    if len(line_segments)/2 == 1:
                        img_draw.point(list(line_segments), fill=track_color[str(track_num)])
                    elif len(line_segments)/2 > 1:
                        img_draw.line(list(line_segments), fill=track_color[str(track_num)])
                
                # Convert the image to a numpy array 
                image_frame_array = np.array(img_frame)
                img_frame.close()
                image_frame_array = image_frame_array/255.0
                new_particle_group.append(image_frame_array)
                
            image_frames += new_particle_group
        #frame_check = frame
    skimage.viewer.CollectionViewer(image_frames).show()
    gc.collect()

def trackpy_rot_motion_linker(data_frame, search_range, rot_velocity=0.0, memory=0, theta_lim_bias=[0,360], **kwargs):
    '''A wrapper for trackpy linking that includes a predictor for rotating particles

    :params data_frame: DataFrame containing all the particle position information
    :params search_range: Max distance a particle can move between frames
    :params rot_velocity: The bias (in degrees) that a candidate particle should be
    found at. This is positive for positive L's
    :params memory: The number of frames a particle can disappear for and still be 
    considered the same particle.
    :params theta_lim_bias: The limits in degrees theta to apply the rotational bias.
    If a particle is outside this range then no bias is applied.
    :params kwrgs: Additional keyword arguments passed to trackpy.link_df
    '''

    # Find the particle locations in polar coords
    xf, yf, rf = least_sq_fit_circle(data_frame)
    polar_coor_data_frame(data_frame, xf, yf)

    # Generate the predictor function
    @trackpy.predict.predictor
    def predict(t1, particle):
        theta = calc_angle(particle.pos[0], particle.pos[1], xf, yf)
        r = calc_radius(particle.pos[0], particle.pos[1], xf, yf)
        if theta_lim_bias[0] < theta < theta_lim_bias[1]:
            new_theta = theta + rot_velocity * (t1 - particle.t)
            new_theta %= 360.0
            new_x = calc_x_from_polar(r, new_theta, xf)
            new_y = calc_y_from_polar(r, new_theta, yf)
            return np.array((new_x,new_y))
        else:
            return np.array((particle.pos[0], particle.pos[1]))
    
    # Track the data and restructure the resulting DataFrame
    trackpy.link_df(data_frame, search_range, memory=memory, pos_columns=['x pos', 'y pos'],
                    retain_index=True, link_strategy='numba', predictor=predict, **kwargs)
    data_frame['track id'] = data_frame['particle']
    del data_frame['particle']

def notebook_title_info(notebook_full_path):
    import re
    notebook_name = os.path.split(notebook_full_path)[-1]
    serial_number = re.search('(Ana_.* - )', notebook_name)
    serial_number = serial_number.groups()[0][:-3]
    subtitle = re.search('( - .*)', notebook_name)
    subtitle = subtitle.groups()[0][3:]
    return [notebook_name, serial_number, subtitle]

def add_list_of_dfs_to_hdf(hdf_obj, dfs_list, experiment_name, movie_names_list, selector_metadata, individual_metadata=None):
    '''A function for adding DataFrames to an HDF5 file and making entries
    to the index table used to select based on metadata
    
    :param hdf_obj: an HDFStore pandas object to store the data in
    :param dfs_list: A list of the data frames you want to add to the HDF file
    :param (str) experiment_name: A string to tag the experiment that all 
    DataFrames in dfs_list falls under
    :param (list) movie_names_list: A list of strings to name each of the 
    DataFrames by in dfs_list. Must be the same length as dfs_list
    :param (dict) selector_metadata: A dictionary of keys and values to describe
    the data that will be put in the indexer table. They can contain single
    values or a list of length equal to dfs_list.
    '''
    selector_df = pd.DataFrame(selector_metadata)
    #selector_df['mov_index'] = selector_df.index.values+ 1 
    for num,i in enumerate(dfs_list):
        movie_name = re.search('(Mov_[0123456789]{8})',movie_names_list[num])
        key = experiment_name+'/'+movie_name.groups()[0]
        selector_df.loc[num, 'mov_index'] = int(movie_name.groups()[0][-2:])
        selector_df.loc[num, 'key'] = key
        hdf_obj.put(key, i)
    # Pandas will make the mov_index column Floats because appending will make it
    # have NaNs.
    selector_df.loc[:,'mov_index'] = selector_df.loc[:,'mov_index'].astype('float64')
    try:
        hdf_obj.get('index')
    except KeyError:
        hdf_obj.put('index', selector_df)
    if selector_df.isin(hdf_obj.index.reset_index()).all().all() == True:
        print "Values already in index!"
        return
    elif selector_df.isin(hdf_obj.index.reset_index()).any().any() == False:
        print "Some values of the index match while others don't"
        print selector_df[selector_df.isin(hdf_obj.index)]
        return 
    else:
        hdf_obj.append('index', selector_df)
        hdf_obj.put('index', hdf_obj.index.reset_index(drop=True))
        return

def find_k_nn(grp, num_nn_k=None):
    '''This function will find the k nearest neighbors in x-y and add 
    columns to the data frame that include the nn number (1st, 2nd,...), the
    particle id of the nn, and the distance of that nn.
    The way to use this function is as such:
    
    df.groupby('frame',group_keys=False).apply(find_k_nn, num_nn_k=X)
    
    The group_keys kwarg prevents a redundant frames column. Reseting the index
    will give the data frame a regular integer index like before the function
    is applied. The second line rearranges the columns to the correct order.
    '''
    xy_data = grp[['x pos','y pos']]
    tree = scipy.spatial.KDTree(xy_data)
    if num_nn_k == None:
        num_nn_k = len(xy_data)
    elif num_nn_k > len(xy_data)-1:
        num_nn_k = len(xy_data)
    elif num_nn_k != None:
        num_nn_k += 1
    nn_d, nn_i = tree.query(xy_data, k=num_nn_k)
    if len(nn_d)==1: # If only one particle return the group
        grp['nn_num'] = np.nan
        grp['nn_id'] = np.nan
        grp['nn_dist'] = np.nan
        return grp
    particle_ids = grp['track id'].values
    track_ids = np.tile(particle_ids, (num_nn_k,1))
    track_ids = track_ids.T[:,1:].flatten()
    nn_ids = particle_ids[nn_i][:,1:].flatten()
    nn_dist = nn_d[:,1:].flatten()
    # Create nn number column (1st, 2nd, etc)
    nn_num = np.arange(len(nn_d[0]))
    nn_num = np.tile(nn_num,(len(particle_ids),1))[:,1:]
    nn_num = nn_num.flatten()
    # Merge with current group
    nn_df = pd.DataFrame(np.vstack((track_ids, nn_num, nn_ids, nn_dist)).T, columns=['track id','nn_num','nn_id','nn_dist'])
    new_df = pd.merge(grp, nn_df, left_on='track id', right_on='track id')
    return new_df

def find_k_nn_theta(grp, num_nn_k=None):
    '''This function will find the k nearest neighbors for theta and add 
    columns to the data frame that include the nn number (1st, 2nd,...), the
    particle id of the nn, and the distance of that nn. The function respects
    boundary conditions of the data such that it is periodic at 360 degrees.
    The way to use this function is as such:
    
    df=df.groupby('frame', group_keys=False).apply(find_k_nn_theta, num_nn_k=x)
    
    The group_keys kwrg prevents a redundant frames column. Reseting the index
    will give the data frame a regular integer index like before the function
    is applied. The second line rearranges the columns to the correct order.
    '''
    from periodic_kdtree import PeriodicCKDTree
    bounds = np.array([360])
    data = grp['theta'].values
    data = np.reshape(data, [len(data), 1])
    tree = PeriodicCKDTree(bounds, data)
    if num_nn_k == None:
        num_nn_k = len(data)
    elif num_nn_k > len(data)-1:
        num_nn_k = len(data)
    elif num_nn_k != None:
        num_nn_k += 1
    d, i = tree.query(data, k=num_nn_k)
    if len(d) == 1:  # If only one particle return the group
        grp['theta_nn_num'] = np.nan
        grp['theta_nn_id'] = np.nan
        grp['theta_nn_dist'] = np.nan        
        return grp
    # Create particle id column
    particle_ids = grp['track id'].values
    track_ids = np.tile(particle_ids, (num_nn_k, 1))
    track_ids = track_ids.T[:, 1:].flatten()
    nn_ids = particle_ids[i]
    # Create nn number column (1st, 2nd, etc)
    nn_num = np.arange(num_nn_k)
    nn_num = np.tile(nn_num, (len(particle_ids), 1))[:, 1:]
    nn_num = nn_num.flatten()
    # Create corresponding nn track id
    nn_ids = nn_ids[:, 1:].flatten()
    nn_dist = d[:, 1:].flatten()
    # Merge with current group
    nn_df = pd.DataFrame(np.vstack((track_ids, nn_num, nn_ids, nn_dist)).T,
                         columns=['track id', 'theta_nn_num', 'theta_nn_id', 'theta_nn_dist'])
    new_df = pd.merge(grp, nn_df, left_on='track id', right_on='track id')
    return new_df

def find_k_nn_xy_and_theta(grp, num_nn_k=None):
    '''This function will find the k nearest neighbors in xy and also in theta
    and add them to the DataFrame. The way to use this function is:
    
    df=df.groupby('frame', group_keys=False).apply(find_k_nn_xy_and_theta, num_nn_k=x)
    '''
    grp_copy = grp.copy()
    xy_grp = find_k_nn(grp_copy, num_nn_k)
    theta_grp = find_k_nn_theta(grp_copy, num_nn_k)
    new_df = pd.concat([xy_grp, theta_grp.loc[:,['theta_nn_num', 'theta_nn_id', 'theta_nn_dist']]], axis=1)
    return new_df

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