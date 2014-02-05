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

def refined_least_sq_fit_circle(data_frame, x_center, y_center, r_inner, r_outer):
	

	
track_data=import_matlab_gui()
track_data=matlab_gui_to_data_frame(m)
img_height=input('Height of raw data (px)')
y_axis_flip(track_data,img_height)

