from __future__ import division
import matplotlib

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
	return plt.imread(Filename) 

def iden_img_sequence(FileName, p_rad, hwhm, d_rad, threshold, mask_rad):
	f_path,f_name=os.path.split(FileName)
	f_base,f_ext=os.path.splitext(f_name)
	start_index=listdir(f_path).index(f_name)
	prams=[p_rad, hwhm, d_rad, threshold, mask_rad]
	for file in listdir(f_path)[start_index:]
		if f_ext==os.path.splitext(file)[1]:
			im=open_16bit_img(os.path.join(f_path,file))
			iden_image(im, *prams)
	
def iden_image(img_array, p_rad, hwhm, d_rad, threshold, mask_rad):
	'''Takes an image array and all the parameters necessary to do the identification'''
	bp_img = tid.band_pass(img_array, p_rad, hwhm)
	res_lm = tid.find_local_max(bp_img, d_rad, threshold)
	locs, mass, r2 = tid.subpixel_centroid(bp_img, res_lm, mask_rad)
	return locs, mass, r2
			

root = Tkinter.Tk()
root.withdraw()
FileName = os.path.split(tkFileDialog.askopenfilename())