import string
import Tkinter, tkFileDialog
import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist
import matplotlib.pyplot as plt
import os

def import_diatrack_txt():
	'''This function imports the data from a Diatrack text
	file into a pandas dataframe'''
	root = Tkinter.Tk()
	root.withdraw()
	track_id=0
	imported_data=[]
	file_path=tkFileDialog.askopenfilename(filetypes=[('diatrack files','.txt')])
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

def calc_disp(data_frame,waiting_time=1):
	'''Insert the data frame and the waiting time for displacement 
	calculation and it will return the data frame with an extra column
	for the displacement'''
	track=data_frame.sort(['track id','frame'])
	#Take the diff of the x's  and y's by their waiting_time
	delta_x=track['x pos'].diff(periods=waiting_time)
	delta_y=track['y pos'].diff(periods=waiting_time)
	disp=np.sqrt(delta_x**2 + delta_y**2)
	disp.name = 'disp '+str(waiting_time)
	#Join the displacement series with the original df and shift it
	#such that the NaN fall on the end of a particle track
	try:
		track=track.join(disp.shift(-waiting_time))
	except:
		pass
	#Fix boundaries so that particles near the end of their track 
	#have a displacement of NaN
	track['disp '+str(waiting_time)][track['track id'] != track['track id'].shift(-waiting_time)] = np.nan
	return track.sort(['frame','track id'])

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

def plot_polar_transform(polar_transform_output):
	'''This function plots the polar transform from the output from
	polar_transform_points'''
	plt.scatter(polar_transform_output[:,0],polar_transform_output[:,1])
	plt.xlabel('Radius (px)')
	plt.ylabel('Degrees')
	#plt.axis('equal')
	plt.show()

diatrack_list,file_path=import_diatrack_txt()
plot_particle_positions(diatrack_list)
#Ask if you want to cut out a radial part of the data
looks_good=raw_input('Does this data need cutting? (y or n) ')
#This while loop will allow you to iterate and refine the radial
#cutting. The Least squares fitting could use some refinement
#in order to better fit the circle
while looks_good=='y':
	xf,yf,rf=least_sq_fit_circle(diatrack_list)
	transform=polar_transform_points(diatrack_list)
	plot_polar_transform(transform)
	radial_bins=input('Define radial bounds (r1,r2) ')
	tmp=filter_data_frame_by_radius(diatrack_list,xf,yf,radial_bins[0],radial_bins[1])
	plot_particle_positions(tmp)
	looks_good=raw_input('Does this data need more cutting? (y or n)')
	diatrack_list=tmp
	if looks_good=='n':
		diatrack_list=tmp
		break
#disp_lag_time=input('What is displacement lag time? ')
#Hard coded lag time values
disp_lag_time=[1,3,5,10,20]
#Calculate the displacements for each lag time
for lag_time in disp_lag_time:
	diatrack_list=calc_disp(diatrack_list,waiting_time=lag_time)
frame_rate=input('What is the frame rate? (in fps): ')
l_value=input('What is L for this movie? ')
over_plate=raw_input('Is this over a nanoplate? (y or n) ')
if over_plate=='y':
	over_plate=' Over Plate '
else:
	over_plate=' Over Glass '

um_conv=6.5/60/1.6/2 #Conversion factor between px and um
disp_call_list=['disp '+str(v) for v in disp_lag_time]
#disp_in_nm_and_sec=diatrack_list['disp '+str(lag_time)]*um_conv*frame_rate/disp_lag_time
disp_mean=[np.sqrt(diatrack_list['disp '+str(v)].mean()*um_conv*frame_rate/v) for v in disp_lag_time]
disp_call_list=[disp_call_list[v]+' mean='+str('{0:.2f}'.format(disp_mean[v])) for v in range(len(disp_mean))]

#Plotting function for plotting several displacements for each movie
plt.hist([np.sqrt(diatrack_list['disp '+str(v)].dropna()*um_conv*frame_rate/v) for v in disp_lag_time],normed=True,range=[0,100],bins=50,histtype='step',label=disp_call_list)
plt.legend()
plt.xlabel('Speed (um/s)')
plt.ylabel('PDF')
plt.title(str(os.path.split(file_path)[1])+' L='+str(l_value)+over_plate+'Lag Time '+str(disp_lag_time))
plt.show()

def eliminate_b_box_from_data_frame(data_frame):
	'''This function allows you to eliminate points within a bounding box
	from upper left to lower right between two points.'''
	plot_particle_positions(data_frame)
	#Ask if you want to cut out a radial part of the data
	looks_good=raw_input('Does this data need cutting? (y or n) ')
	#This while loop will allow you to iterate and refine the radial
	#cutting. The Least squares fitting could use some refinement
	#in order to better fit the circle
	while looks_good=='y':
		print "Add a bounding box to eliminate (x1,y1,x2,y2) from upper left to lower right: "
		plot_particle_positions(data_frame)
		b_box=input('Bounding Box to eliminate (x1,y1,x2,y2): ')
		cut_out_x1=data_frame['x pos']<b_box[0]
		cut_out_x2=data_frame['x pos']>b_box[2]
		cut_out_y1=data_frame['y pos']>b_box[1]
		cut_out_y2=data_frame['y pos']<b_box[3]
		data_frame=data_frame[cut_out_x1 | cut_out_x2 | cut_out_y1 | cut_out_y2]
		plot_particle_positions(data_frame)
		looks_good=raw_input('Does this data need more cutting? (y or n) ')
		if looks_good=='n':
			break
	return data_frame
	
def plot_ang_vs_time(data_frame,x_cent,y_cent):
	track_id_grp=data_frame.groupby('track id')
	for name,grp in track_id_grp:
		theta=[]
		for value in range(len(grp)):
			theta.append(180+180*np.arctan2(data_frame['x pos'][value]-x_cent,data_frame['y pos'][value]-y_cent)/np.pi)
		plt.plot(theta)
		
		