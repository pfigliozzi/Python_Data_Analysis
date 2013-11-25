
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
		
		

[v for v in nn[1]]