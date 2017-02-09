import string
import copy
import scipy
import Tkinter, tkFileDialog
import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist
import matplotlib.pyplot as plt
import os
import re
import sys
import cPickle
import glob
from PIL import Image, ImageDraw, ImageColor
sys.path.append(os.path.abspath("C:\Users\Scherer Lab E\Documents\GitHub\Python_Data_Analysis"))
import half_nanoplate_functions as hnf
import common_functions


'''Restoring from checkpoint'''
data_dir="C:\Users\Scherer Lab E\Downloads\TrackingGUI_and_associated_files_20July2014 My Version"
os.chdir(data_dir)
file_list = glob.glob('Mov_012014*processed_linked.pandas')
final_positions = [pd.read_pickle(i) for i in file_list]
# data_list = [common_functions.matlab_gui_to_data_frame(i) for i in data_list]

common_functions.view_trajectories_new_particles(final_positions[5].query('100 < theta < 310').query('frame < 1000'), frame_window=1, image_size=[390,390])