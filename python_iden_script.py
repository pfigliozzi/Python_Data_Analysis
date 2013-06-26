from __future__ import division
import matplotlib
import os
from os import listdir
import trackpy.tracking as pt
import trackpy.identification as tid
import matplotlib.pyplot as plt
import random
import numpy as np
import string, Tkinter, tkFileDialog
from decimal import *
import Image
import numpy as np

def open_16bit_img(FileName):
	'''Simple function to open 16bit .tiff images. Matplotlib apparently handles this very well'''
	return plt.imread(FileName) 

def iden_img_sequence(FileName, p_rad, hwhm, d_rad, threshold, mask_rad):
	'''Identifies all the particles in a image series where each image is it's one file in a specific folder'''
	f_path,f_name=os.path.split(FileName)
	f_base,f_ext=os.path.splitext(f_name)
	start_index=listdir(f_path).index(f_name)
	prams=[p_rad, hwhm, d_rad, threshold, mask_rad]
	locs=[]
	for file in listdir(f_path)[start_index:]:
		if f_ext==os.path.splitext(file)[1]:
			im=open_16bit_img(os.path.join(f_path,file))
			locs.append(iden_image(im, *prams)[0])
			print file
		else:
			print 'continue'
			continue
	return locs

def track_data(points_list):
	data=[]
	for frame in range(len(points_list)):
		frame_data=[]
		data.append(frame_data)
		for pts in range(len(points_list[frame][0,:])):
			data[frame].append(pt.PointND(float(frame+1), np.asarray((points_list[frame][0,pts],points_list[frame][1,pts]))))
	hash_generator = lambda: pt.Hash_table((20, 20), .5)
	t=pt.link(data,2)
	return t
	'''fig = plt.figure()
	ax = fig.gca()
	for trk in t:
		x, y = zip(*[p.pos for p in trk.points])
		ax.plot(x, y)
	plt.axis('equal')
	plt.show()'''
	
def iden_image(img_array, p_rad, hwhm, d_rad, threshold, mask_rad):
	'''Takes an image array and all the parameters necessary to do the identification'''
	bp_img = tid.band_pass(img_array, p_rad, hwhm)
	res_lm = tid.find_local_max(bp_img, d_rad, threshold)
	locs, mass, r2 = tid.subpixel_centroid(bp_img, res_lm, mask_rad)
	return locs, mass, r2
			

root = Tkinter.Tk()
root.withdraw()
FileName = tkFileDialog.askopenfilename()