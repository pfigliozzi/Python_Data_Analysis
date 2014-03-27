import string
import scipy
import Tkinter, tkFileDialog
import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist
import matplotlib.pyplot as plt
import os
import sys

def import_matlab_gui(file_path=None):
	'''This function will import data files that are generated from 
	Rahgu's tracking GUI. Pick a matlab file and the output will be
	an array'''
	import scipy.io
	if file_path==None:
		import Tkinter, tkFileDialog
		root = Tkinter.Tk()
		root.withdraw()
		file_path = tkFileDialog.askopenfilename(filetypes=[('matlab files','.mat')])
	print 'Loading '+file_path
	m=scipy.io.loadmat(str(file_path))['objs_link']
	if len(m)==0:
		m=scipy.io.loadmat(str(file_path))['objs']
	return m

def matlab_gui_to_data_frame(m):
	import pandas
	return pandas.DataFrame({'x pos' : m[0][:], 'y pos' : m[1][:], 'frame' : m[4][:], 'particle id' : m[3][:], 'track id' : m[5][:]})

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
	tree=scipy.spatial.KDTree(grp[['x pos','y pos']])
	nn=tree.query(grp[['x pos','y pos']], k=12)
	reindex_value=grp.index[0]
	grp['nn_part']=[grp['track id'][v].values for v in nn[1]+reindex_value]
	grp['nn_dist']=[v for v in nn[0]]
	return grp

def plot_particle_positions(data_frame):
	'''If the particle data is in a data frame (from matlab_gui_to_data_frame)
	then this function will plot the position of all the particles.
	Keep in mind that the y axis is flipped after images are run through
	tracking software'''
	plt.scatter(data_frame['x pos'], data_frame['y pos'])
	plt.axis('equal')
	plt.xlabel('Pixels')
	plt.ylabel('Pixels')
	plt.show()
	
def polar_transform_points(data_frame, origin=None):
	'''Polar transforms the particle positions in all frames around
	a defined center (origin).'''
	xy=data_frame[['x pos','y pos']].values
	if origin is None:
		origin = np.mean(xy, axis=0)
	output= np.zeros(xy.shape)
	#This finds the radius of each point fromt the center
	output[:,0]= np.sqrt((origin[0]-xy[:,0])**2 +(origin[1]-xy[:,1])**2)
	#This finds the angle of each point starting with 0 at 12 o'clock
	#and going counter clockwise
	output[:,1]= 180+180*np.arctan2(xy[:,0]-origin[0],xy[:,1]-origin[1])/np.pi
	return output
	
def y_axis_flip(data_frame, img_height):
	'''The Raghu code spits out particle positions with the a shift
	of the origin from the top left to the bottom left resulting in
	a flip in the y-axis for all the data. This function flips it 
	back to values that match the image. It requires knowledge of
	the height of the original data.'''
	data_frame['y pos']=abs(data_frame['y pos']-img_height)

def least_sq_fit_circle(data_frame):
	'''Calculate the best fit circle by least squares for the particle
	locations for all frames in a data frame. This funtion finds
	the center via center of mass'''
	from scipy.spatial.distance import cdist
	import scipy.optimize
	xy=data_frame[['x pos','y pos']].values
	#Find center via center of mass
	center = np.mean(xy, axis=0).reshape(1,2)
	radius = cdist(center,xy).mean()
	#Define the error function to minimize
	def err_function((x_center, y_center,radius)):
		err = [np.linalg.norm([x-x_center,y-y_center])-radius for x,y in xy]
		return (np.array(err)**2).sum()
	xf,yf,rf=scipy.optimize.fmin(err_function,[center[0,0],center[0,1],radius])
	return xf,yf,rf

def filter_data_frame_by_radius(data_frame, x_cent, y_cent, r_in, r_out):
	'''A simple function that will create a data frame that will
	exclude a set of points not within two set radii'''
	#outer_box = data_frame[['x pos','y pos']].applymap(lambda x: cdist([x_cent,y_cent],[x,x])>r_out)
	outer_box = np.sqrt((data_frame['x pos']-x_cent)**2 + (data_frame['y pos']-y_cent)**2)<r_out
	inner_box = np.sqrt((data_frame['x pos']-x_cent)**2 + (data_frame['y pos']-y_cent)**2)>r_in
	return data_frame[outer_box & inner_box]

def calc_angle(x_pos,y_pos,x_cent,y_cent):
	'''This function calculates the angle of the positions with respect 
	to the centers x_cent and y_cent. This function is designed so that
	the origin is at the 12 o'clock position and increases going counter-
	clockwise
	'''
    return 180+180*np.arctan2(x_pos-x_cent,y_cent-y_pos)/np.pi
	
def nn_distance_angle_seperation(data_frame, number_of_bins, x_cent, y_cent):
	'''Seperates the nearest neighbor distances into angular bins based
	on the peak location. This function only works if you are picking peaks
	from the images and not using tracked data.
	
	:pram peaks: List of peaks from func find_peak_locations
	:pram nn_dist: List of nearest neighbors from func find_nn_peaks
	:pram number_of_bins: int.
	'''
	if 360%number_of_bins!=0:
		print "Error: Number of bins must be divide 360 without remainders"
		return None
	if number_of_bins/2%2!=0:
		print "Error: Half the number of bins must be even"
		return None
	bin_size=360/number_of_bins
	bin_limits=[]
	for i in range(number_of_bins):
		if i==0:
			bin_limits.append([0+bin_size/2,360-bin_size/2])
		else:
			bin_limits.append([(i*bin_size)+bin_size/2,(i*bin_size)-bin_size/2])
	radial_bins=[]
	for step in range(len(bin_limits)):
		#Define the less than and greater than criterea
		less_than=calc_angle(data_frame['x pos'],data_frame['y pos'],x_cent,y_cent)<=bin_limits[step][0]
		greater_than=calc_angle(data_frame['x pos'],data_frame['y pos'],x_cent,y_cent)>=bin_limits[step][1]
		if step==0:
			radial_bins.append(np.concatenate([[v[1]] for v in data_frame[less_than | greater_than]['nn_dist'].values]))
			#radial_bins.append(np.concatenate([nn_dist[v][np.logical_or(peaks[v]<=bin_limits[step][0],peaks[v]>=bin_limits[step][1])] for v in range(len(peaks))]))
		else:
			radial_bins.append(np.concatenate([[v[1]] for v in data_frame[less_than & greater_than]['nn_dist'].values]))
			#radial_bins.append(np.concatenate([nn_dist[v][np.logical_and(peaks[v]<=bin_limits[step][0],peaks[v]>=bin_limits[step][1])] for v in range(len(peaks))]))
	return radial_bins,bin_limits

def plot_polar_transform(polar_transform_output):
	'''This function plots the polar transform from the output from
	polar_transform_points'''
	plt.scatter(polar_transform_output[:,0],polar_transform_output[:,1])
	plt.xlabel('Radius (px)')
	plt.ylabel('Degrees')
	#plt.axis('equal')
	plt.show()
	
def import_diatrack_txt(file_path=None):
	'''This function imports the data from a Diatrack text
	file into a pandas dataframe'''
	if file_path==None:
		root = Tkinter.Tk()
		root.withdraw()
		file_path=tkFileDialog.askopenfilename(filetypes=[('diatrack files','.txt')])
	track_id=0
	imported_data=[]
	for line in open(file_path,'r').readlines():
		if line.startswith('f'):
			continue
		else:
			track_id+=1
			templine=string.split(line)
			templine=[float(i) for i in templine]
			frame=int(templine[0])
			for pos in range(1,len(templine)-1,3):
				imported_data.append([frame,track_id,templine[pos],templine[pos+1]])
				frame+=1
	imported_data=pd.DataFrame(imported_data,columns=['frame','track id','x pos','y pos'])
	imported_data= imported_data.sort(['frame','track id'])
	imported_data.index = range(len(imported_data))
	print file_path
	return  imported_data,file_path
	
def find_longest_traj(data_frame):
	'''In a data frame created from import_matlab_gui() or 
	import_diatrack_txt and will return a data frame that contains only
	the longest trajectory 'track id' values'''
	grouped=data_frame.groupby('track id')
	max_len,max_name=[0,0]
	for name,grp in grouped:
		if len(grp)>max_len:
			max_name=name
			max_len=len(grp)
		else:
			continue
	return data_frame[data_frame['track id']==max_name]