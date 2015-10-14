import string
import scipy
import Tkinter, tkFileDialog
import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist
import matplotlib.pyplot as plt
import os
import sys
import skimage.viewer
from PIL import Image, ImageDraw, ImageColor
import trackpy
import trackpy.predict

def plot_nn_dist_par_perp(data_frame,fit_params,plot_name):
    ang_sep,ang_bins=nn_distance_angle_seperation_ver_2(data_frame,8,1)
    print len(ang_sep)
    ang_sep=np.array([np.concatenate((ang_sep[v],ang_sep[v+4])) for v in range(len(ang_sep)/2)])
    um_conv=6.5/60/1.6/2
    print len(ang_sep)
    ang_sep_conv=ang_sep*um_conv
    #labels=['Parellel','Perpendicular']
    print ang_bins
    #print ang_sep
    labels=[str(v[0])+'-'+str(v[1]) for v in ang_bins]
    plt.figure(figsize=[4,3])
    n,b,p=plt.hist([v for v in ang_sep_conv],bins=100,
                   range=(0,3),normed=True,histtype='step',label=labels)
    plt.legend()
    plt.xlabel('Seperation (um)')
    plt.ylabel('Probability Density')
    plt.title(plot_name)
    #fig_dir="J:\Pat's Projects\Dynamical Phase Transition\Figures for Paper\\nn_histogram_comparison\\"
    #plt.savefig(fig_dir+plot_name[:-4]+".png", dpi=200)
    plt.show()