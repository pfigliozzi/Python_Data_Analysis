def import_matlab_gui():
	'''This function will import data files that are generated from 
	Rahgu's tracking GUI. Pick a matlab file and the output will be
	an array'''
	import scipy.io
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
	grp['nn_part']=[grp['particle id'][v].values for v in nn[1]+reindex_value]
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

#def angular_range_of_data_frame(data_frame, 
	
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
		less_than=(180+180*np.arctan2(data_frame['x pos']-x_cent,data_frame['y pos']-y_cent)/np.pi)<=bin_limits[step][0]
		greater_than=(180+180*np.arctan2(data_frame['x pos']-x_cent,data_frame['y pos']-y_cent)/np.pi)>=bin_limits[step][1]
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

'''This is the begining of the script portion of the analysis.
This essential runs all the functions above in the right order to 
perform the analysis properly.'''
track_data=import_matlab_gui()
track_data=matlab_gui_to_data_frame(track_data)
img_height=input('Height of raw data (px)')
y_axis_flip(track_data,img_height)
xf,yf,rf=least_sq_fit_circle(track_data)
polar_transform=polar_transform_points(track_data, origin=[xf,yf])
plot_polar_transform(polar_transform)
radial_bins=input('Decide 4 radial boundaries (1,2,3,4)')

'''This section seperates the data into radial bins'''
figure_data=[]
for bound in range(len(radial_bins)-1):
	rad_sep_df=filter_data_frame_by_radius(track_data, xf, yf, radial_bins[bound], radial_bins[bound+1])
	ang_sep,ang_limits=nn_distance_angle_seperation(rad_sep_df,4,xf,yf)
	figure_data.append([[radial_bins[bound],radial_bins[bound+1]],ang_sep,ang_limits])

nm_conv=6.5/60/1.6/2 #Conversion factor between px and nm
fig_data_lab_units=[]
for rad_bin in range(len(figure_data)): #this loop seperates the data into angular bins
	fig_data_lab_units.append(np.copy(figure_data[rad_bin]))
	fig_data_lab_units[rad_bin][1]=[v*nm_conv for v in fig_data_lab_units[rad_bin][1]]

'''This section allows you to plot the three radial bins nn histograms
next to eachother'''
fg1,axarr=plt.subplots(1,3,sharex=True, sharey=True)
for figure in range(len(fig_data_lab_units)):
	labels=[str(v[1])+'-'+str(v[0]) for v in fig_data_lab_units[figure][2]]
	axarr[figure].hist(fig_data_lab_units[figure][1],normed=True,histtype='step',label=labels,range=(0,2),bins=60)
	axarr[figure].set_title('Radius bounds'+str(fig_data_lab_units[figure][0][0])+'-'+str(fig_data_lab_units[figure][0][1]))
plt.xlabel('1st NN Sep (um)')
plt.ylabel('Prob Density')
#fg1.subplots_adjust(vspace=0)
axarr[2].legend()

'''This small section of the code allows you to combine the output of 
nn_distance_angle_seperation() from 4 sections to 2 where the two 
parts are the parallel and perpendicular parts of the ring. The
last two lines convert the units from pixels to um'''
ang_sep=np.array((np.concatenate((ang_sep[0],ang_sep[2])),np.concatenate((ang_sep[1],ang_sep[3]))))
um_conv=6.5/60/1.6/2
ang_sep_conv=ang_sep*um_conv

def norm(x,mean,std,a):
	return a*(1/(std*np.sqrt(2*np.pi)))*np.exp(-((x-mean)**2)/(2*std**2))

def plot_par_and_perp_nn(ang_sep,plot_title):
	labels=['Parellel (top/bottom)','Perpendicular (left/right)']
	n,b,p=plt.hist([v for v in ang_sep],bins=100,range=(0,3),normed=True,histtype='step',label=labels)
	plt.legend()
	plt.title(plot_title)
	plt.xlabel('Seperation (um)')
	plt.ylabel('Probability Density')
	plt.show()
