
def import_matlab_gui():
	'''This function will import data files that are generated from 
	Rahgu's tracking GUI. Pick a matlab file and the output will be
	an array'''
	import scipy.io
	import Tkinter, tkFileDialog
	root = Tkinter.Tk()
	root.withdraw()
	file_path = tkFileDialog.askopenfilename(filetypes=[('matlab files','.mat')])
	m=scipy.io.loadmat(str(file_path))['objs_link']
	if len(m)==0:
		m=scipy.io.loadmat(str(file_path))['objs']
	return m

def matlab_gui_to_data_frame(m):
	import pandas
	return pandas.DataFrame({'x pos' : m[0][:], 'y pos' : m[1][:], 'frame' : m[4][:], 'particle id' : m[3][:], 'track id' : m[5][:]})

def find_nn(grp):
	'''Used to find the nearest neighbors in a data frame generated from matlab_gui_to_data_frame.
	This adds both the nearest neighbor track ids and the distances to each particle in each frame.
	Use this by using df.groupby('frame').apply(find_nn) where df is the data frame generated from matlab_gui_to_data_frame'''
	tree=scipy.spatial.KDTree(grp[['x pos','y pos']])
	nn=tree.query(grp[['x pos','y pos']], k=12)
	reindex_value=grp.index[0]
	grp['nn_part']=[grp['track id'][v].values for v in nn[1]+reindex_value]
	grp['nn_dist']=[v for v in nn[0]]
	return grp
	
def add_nn_to_df(grouped):
	for name, group in grouped:
		tree=scipy.spatial.KDTree(group[['x pos','y pos']])
		nn=tree.query(group[['x pos','y pos']], k=12)
		grp['nn_part']=[grp['track id'][v].values for v in nn[1]]
		grp['nn_dist']=[v for v in nn[0]]
	return grouped
		
		
def plot_traj(df,track_id):
	trackdf=df[df['track id']==track_id]
	trackdf.plot(x='x pos', y='y pos')
	plt.axis('equal')
	plt.show()
	
def open_tif_series():
	'''Opens a tif series as a giant numpy array. Only works for files
	whose files names are numbers. Running the function will allow you
	to select the first image and it will open the rest from there'''
	import Tkinter, tkFileDialog
	import os
	from PIL import Image
	root = Tkinter.Tk()
	root.withdraw()
	file_path = tkFileDialog.askopenfilename(filetypes=[('tif files','.tif')])
	f_path,f_name=os.path.split(file_path)
	f_base,f_ext=os.path.splitext(f_name)
	im=Image.open(file_path)
	listnames=[]
	for file in os.listdir(f_path):
		try:
			listnames.append(int(os.path.splitext(file)[0]))
		except:
			continue
	maxfile=len(listnames)
	array=np.zeros((im.size[1],im.size[0],maxfile))
	array[:,:,0]=np.array(im.getdata()).reshape(im.size[1],im.size[0])
	start_index=os.listdir(f_path).index(f_name)+1
	index=1
	for file in os.listdir(f_path)[start_index:]:
		if f_ext==os.path.splitext(file)[1]:
			im=Image.open(os.path.join(f_path,file))
			array[:,:,index]=np.array(im.getdata()).reshape(im.size[1],im.size[0])
			index+=1
			print 'finished frame'+str(index)
		else:
			continue
	return array

def determine_threshold(avg_image, block_size=100, offset=-2):
	'''Used to determine the right threshold value for the segmentation. The input
	is the average of the image stack of particles on the ring. Play with the offset
	and block size to make a clear ring with minimal background noise. Negative values
	of offset should reduce background noise. This functions returns the thresholded array
	in addition to showing what it looks like.'''
	from skimage.filter import threshold_adaptive
	threshold=threshold_adaptive(avg_image, block_size, offset=offset)
	import matplotlib.pyplot as plt
	plt.imshow(threshold)
	plt.show()
	return threshold

def obj_segmentation(threshold_image, largest_obj=1):
	'''Segments a thresholded image and filters out small objects except for the
	largest 'n' objects.
	
	:pram threshold_image: thresholded image
	:type threshold_image: bool. np array.
	:pram largest_obj: the number of largest objects to keep (default=1)
	:type largest_obj: int.
	
	:returns: bool. array of just the largest object(s)
	'''
	label_objs, num_labels = scipy.ndimage.label(threshold_image)
	sizes=np.bincount(np.int64(label_objs.ravel()))
	sizes[0]=0
	mask_sizes=sizes>=sizes[sizes.argsort()[-largest_obj]]
	segmented=mask_sizes[label_objs]
	plt.imshow(segmented)
	plt.show()
	return segmented
	
def polar_transform(image, origin=None):
	'''Does a polar transform of an image. Note, this is hard coded to not work
	for images larger than 512x512.
	
	:pram image: image
	:type image: float. np array.
	:pram origin: origin of the transform
	:type origin: tuple.
	'''
	if origin is None:
		origin = (image.shape[0]//2,image.shape[1]//2)
	maxr=np.sqrt((origin[0])**2+origin[1]**2)
	thetabins=360
	def cart2polar(outputcoords):
		#r=np.sqrt((origin[0]-outputcoords[0])**2 +(origin[1]-outputcoords[1])**2)
		#theta=180+180*np.arctan2(outputcoords[1]-origin[1],outputcoords[0]-origin[0])/np.pi
		x=outputcoords[0]*np.cos((outputcoords[1]-180)*np.pi/180)+origin[0]
		y=outputcoords[0]*np.sin((outputcoords[1]-180)*np.pi/180)+origin[1]
		#theta_index=np.round((theta*360/np.pi)/(np.pi*2))
		return (x,y)
	return scipy.ndimage.geometric_transform(image, cart2polar, output_shape=(400,360))
	
def transform_images(array, origin=None):
	transform=np.zeros((400,360,len(array[1,1,:])))
	for img in range(len(array[1,1,:])):
		transform[:,:,img]=polar_transform(array[:,:,img],origin=origin)
		print "Frame "+str(img)+" complete!"
	return transform
	
def profile_and_image_plot(frames):
	import Tkinter, tkFileDialog
	import os
	root = Tkinter.Tk()
	root.withdraw()
	file_path = tkFileDialog.askdirectory()
	for frame in range(len(frames[1,1,:])):
		fig=plt.figure()
		ax2=plt.subplot(212)
		plt.xlabel('Angle (degrees)')
		plt.ylabel('Sum Intensity')
		plt.plot(np.average(frames[:,:,frame], axis=0))
		plt.xlim((0,360))
		plt.ylim((100,130))
		ax1=plt.subplot(211)
		plt.ylabel('Radius')
		plt.xticks
		plt.setp(ax1.get_xticklabels(), visible=False)
		plt.ScalarFormatter(useOffset=True)
		fmt=plt.ScalarFormatter(useOffset=True)
		fmt.set_useOffset(-100)
		plt.gca().yaxis.set_major_formatter(fmt)
		plt.imshow(frames[:,:,frame],cmap=plt.cm.Greys_r)
		plt.tight_layout()
		plt.savefig(os.path.join(file_path,str(frame+1)+'.png'), format='png')
		del fig
		plt.clf()
		plt.close()
		print 'Frame '+str(frame)+' is complete!'

def subtract_baseline(one_trace, subtract):
	new_trace=np.copy(one_trace)-subtract
	for x in range(len(new_trace)):
		if new_trace[x]<0:
			new_trace[x]=0
	return new_trace

def subtract_baseline_2D(array, subtract):
	new_array=np.copy(array)
	for step in range(len(new_array[1,:])):
		new_array[:,step]=subtract_baseline(new_array[:,step],subtract)
	return new_array

def find_peak_locations(array, baseline, order=3):
	temp_array=subtract_baseline_2D(array,baseline)
	peaks=[]
	for step in range(len(temp_array[1,:])):
		peaks.append(scipy.signal.argrelmax(temp_array[:,step],order=order,mode='wrap')[0])
	return peaks

def find_nn_peaks(peaks):
	nn_dist=[]
	for step in peaks:
		tree=scipy.spatial.KDTree(np.transpose(np.array([step,np.zeros(len(step))])))
		nn=tree.query(np.transpose(np.array([step,np.zeros(len(step))])),k=2)
		if np.abs(step[0]+np.abs(step[-1]-360))<nn[0][0,1]:
			nn[0][0,1]=np.abs(step[0]+np.abs(step[-1]-360))
		if np.abs(step[-1]-360)+step[0]<nn[0][-1,1]:
			nn[0][-1,1]=np.abs(step[-1]-360)+step[0]
		nn_dist.append(nn[0][:,1])
	return nn_dist
		
def nn_distance_angle_seperation(peaks,nn_dist, number_of_bins):
	'''Seperates the nearest neighbor distances into radial bins based
	on the peak location
	
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
		if step==0:
			radial_bins.append(np.concatenate([nn_dist[v][np.logical_or(peaks[v]<=bin_limits[step][0],peaks[v]>=bin_limits[step][1])] for v in range(len(peaks))]))
		else:
			radial_bins.append(np.concatenate([nn_dist[v][np.logical_and(peaks[v]<=bin_limits[step][0],peaks[v]>=bin_limits[step][1])] for v in range(len(peaks))]))
	return radial_bins,bin_limits

def create_hist_seperated_angles(deg_bins,bin_limits):
	#hist_data=[np.concatenate(deg_1),np.concatenate(deg_2),np.concatenate(deg_3),np.concatenate(deg_4),np.concatenate(deg_5),np.concatenate(deg_6)]
	labels=[str(v[1])+'-'+str(v[0]) for v in bin_limits]
	plt.hist(deg_bins,bins=60,range=(0,60),normed=True,histtype='step',label=labels)
	plt.legend()
	plt.xlabel('Degrees')
	plt.ylabel('Probability Density')
	plt.show()


