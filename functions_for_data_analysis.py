
def import_matlab_gui():
	import scipy.io
	import Tkinter, tkFileDialog
	root = Tkinter.Tk()
	root.withdraw()
	file_path = tkFileDialog.askopenfilename(filetypes=[('matlab files','.mat')])
	m=scipy.io.loadmat(str(file_path))['objs_link']
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
	
	'''
	:returns: bool. array of just the largest object(s)
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
	if origin==None:
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
