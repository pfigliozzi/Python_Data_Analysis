import string, scipy
from decimal import *
import numpy as np
from matplotlib import pyplot as plt
import os
from os import listdir

def read_diatrack_file(file_path):
	'''Reads the Diatrack text File into a dictionary
	
	This should eventually be rewritten to import as numpy array, it would be much faster'''
    particle=0
    trajdb={}
    for line in open(file_path,'r').readlines():
        if line.startswith('f'):
            continue
        else:
            particle+=1
            templine=string.split(line)
            intline=[float(i) for i in templine]
            spt=[(particle,intline[0],(len(intline)-1)/3)]
            for k in range(1,len(intline),3):
                spt.append((intline[k],intline[k+1]))
            trajdb[str(particle)]=spt
            #print 'Particle '+str(particle)+' is done!'
    trajdb['range']=[1,particle]
    return trajdb
	
def displacement(x1,x2,y1=0,y2=0,z1=0,z2=0):
    """Calculates the displacement in up to 2 dimensions. If y1, y2, z1, and z2 are undefined it will calculate the 1 dimension displacement"""
    return ((x1-x2)**2+(y1-y2)**2+(z1-z2)**2)**0.5

def calc_single_traj_fft(traj_dict,traj_number,len_fft):
    """Calculated the fft of the tau=1 displacements of a single trajectory"""
    single_traj=traj_dict[str(traj_number)]
    single_traj_displacements=[]
    NFFT=2**nextpow2(len_fft)
    for pos in range(1,len(single_traj[1:])):
        single_traj_displacements.append(displacement(single_traj[pos][0],single_traj[pos+1][0],single_traj[pos][1],single_traj[pos+1][1]))
    print single_traj_displacements[:4],single_traj[:4],len_fft
    return scipy.fft(single_traj_displacements,n=int(NFFT))/len_fft
	
def add_two_unequal_lists(lista,listb):
    """A function that adds the first n values of a list to the second n values of a list where n is the lenght of the shorter list."""
    c=[]
    #lista=lista.tolist()
    #listb=listb.tolist()
    if len(lista) < len(listb):
        for x in range(len(lista)):
            c.append(lista[x]+listb[x])
        c.extend(listb[len(lista):])
    else:
        for x in range(len(listb)):
            c.append(listb[x]+lista[x])
        c.extend(lista[len(listb):])
    return c
	
def avg_all_trajs_fft(traj_dict):
    '''Sums up the FFTs of all the trajectories in a trajectory dictionary'''
    fft_sum=[]
    max_frame=0
    for particle in range(1,traj_dict['range'][1]):
        if len(traj_dict[str(particle)][1:])>max_frame:
	    max_frame=len(traj_dict[str(particle)][1:])
    len_fft=max_frame-1 #Because displacements are FFTed the number of points is number of frames -1
    for particle in range(1,traj_dict['range'][1]):
        if len(traj_dict[str(particle)][1:])<max_frame:
            continue
        else:    
            current_traj_fft=calc_single_traj_fft(traj_dict,particle,len_fft)
            plt.plot(current_traj_fft)
            fft_sum=add_two_unequal_lists(fft_sum,current_traj_fft.tolist())
    NFFT=2**nextpow2(len_fft)
    fft_sum=2*np.array(fft_sum)/float(traj_dict['range'][1])
    return abs(fft_sum[:(NFFT/2)+1])**2,(500/2)*np.linspace(0,1,(NFFT/2)+1)

"""def calc_single_disp(traj_dict):
    for traj in range(1,traj_dict['range'][1]):
        singletraj=trajdb[str(traj)]
    for pos in range(1,len(singletraj[1:0])):
                displacement(singletraj[pos][0],singletraj[pos][1],singletraj[pos+1][0],singletraj[pos+1][1])"""
def nextpow2(n):
    """Finds the power of 2 (the exponent) that is greater than or equal to n"""
    m_f = np.log2(n)
    m_i = np.ceil(m_f)
    return m_i

def openimgseq16(FileName):
	'''Opens a series of individual tiff files as 16bit array. It starts at the Filename as the first entry'''
    f_path,f_name=os.path.split(FileName)
	f_base,f_ext=os.path.splitext(f_name)
	start_index=listdir(f_path).index(f_name)
	for file in listdir(f_path)[start_index:]
		if f_ext==os.path.splitext(file)[1]:
			im=Image.open(os.path.join(f_path,file))
			
			
    start_num_string=''    
    for index in iter_loc:
        start_num_string+=FileName[index]
    starting_number=int(start_num_string)
    image=Image.open(FileName)
    return np.array(image.getdata(),np.uint16).reshape(image.size[1],image.size[0])

def discover_iterator(FileName):
	iter_loc=[]
    for x in range(len(FileName)):
        digits=0
        if FileName[x] in string.digits:
            iter_loc.append(x)
            digits=1
        elif digits==1 and FileName[x] not in string.digits:
            digits=0
        else:
            continue
    return iter_loc
    
