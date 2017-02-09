'''Takes a Diatrack File and converts it into a database (shelf). The shelf is then immediately used to construct the TAMSD'''

import numpy, cPickle, fileinput, string, sys, heapq, os, shelve, string, Tkinter, tkFileDialog

#Formation
#Step 1: Read in each line of data file and parse into integers.



def Number_of_Frames(filename):
    lines=[]
    f=open(filename,'r')
    hiframe=0
    for line in f:
        if line.startswith('f'):
            continue
        temp=string.split(line)
        lines.append(temp)
        if (float(temp[0])+(len(temp[1:])/3)-1)>hiframe:
            hiframe=float(temp[0])+(len(temp[1:])/3)-1
    heap=heapq.nlargest(1,lines,len)
    if (len(heap[0])-1)/3>hiframe:
        return (len(heap[0])-1)/3
    else:
        return int(hiframe)

def Parse_Line(f):
    for line in f:
        templine= string.split(line)
        intLine=[map(int,x) for x in templine]
        intLine=[item for sublist in intLine for item in sublist]
    return intLine

root = Tkinter.Tk()
root.withdraw()
file_path = tkFileDialog.askopenfilename()
print "Opened "+file_path
output=tkFileDialog.asksaveasfilename(initialfile=file_path)
output=output+'.txt'
print "Will save as "+output
#print Parse_Line(f)
#frameCount=2 #frameCount is going to scan over all the frames in a single movie file = (length of longest line-1)/3


f=open(file_path,'r')
trajdb=shelve.open(file_path+'_scratch')
particle=0
#print Number_of_Frames(filename)
for line in f.readlines():
    if line.startswith('f'):
        continue
    else:
        particle+=1
        templine= string.split(line)
        #print templine
        #intline=[map(float,x) for x in templine]
        #intline=[item for sublist in intline for item in sublist]
        intline=[float(i) for i in templine]
        #print intline
        #print "This is line1", line, "This is line2"
        #if frameCount>=intline[0] and frameCount<(intline[0]+(len(intline)-1)/3):
        spt=[(particle,intline[0],(len(intline)-1)/3)]
            #first coord of the first entry refer current position in list w/ 0 index
            # second coord is what frame the particle first appeared, index starting with 1
            # third coord is number of frames the particle formed a trajectory. 
        for k in range(1,len(intline),3):
            spt.append((intline[k],intline[k+1]))
        trajdb[str(particle)]=spt
        print 'Particle '+str(particle)+' is done!'
            #frame.append(spt)
trajdb['range']=[1,particle]

'''This is where the MSD calculation begins'''

import numpy, pylab, string, fileinput, PIL, Image, ImageDraw, cPickle, math, shelve
import sets, sys

def displacement(x1,y1,x2,y2):
    return ((x1-x2)**2+(y1-y2)**2)

def findmsd(numberoftraj,trajs):
    masterlist=[]
    for dt in range(1, numberoftraj):
        #print dt,trajs
        total=0
        numtrajs=0
        for currentpoint in range(len(trajs)-dt):
            #print trajs[currentpoint][0],trajs[currentpoint][1],trajs[currentpoint+dt][0],trajs[currentpoint+dt][1]
            total+=displacement(trajs[currentpoint][0],trajs[currentpoint][1],trajs[currentpoint+dt][0],trajs[currentpoint+dt][1])
            numtrajs+=1
        avg=total/numtrajs
        masterlist.append([dt,avg])
    return masterlist

print str(trajdb['range'][1])+' particles!'

dat=open(output,'w')
dat.close()
dat=open(output,'a')
import numpy as np
first=1
for particles in range(trajdb['range'][0],trajdb['range'][1]):
	templist=findmsd(trajdb[str(particles)][0][2],trajdb[str(particles)][1:])
    #print templist
	if first==1:
		master=np.array(templist)
		first=2
	else:
		master=np.append(master,templist,axis=0)
	for points in templist:
		dat.write(str(points[0])+'     '+str(points[1])+'\n')
	print 'Particle '+str(particles)+' is done!'
dat.close()

import numpy as np

"""dat=open(output,'r')
first=1
print "reading in TAMSD file"
for line in dat.readlines():
	if first==1:
		templine= string.split(line)
		master=np.array([[float(i) for i in templine]])
		first=2
	else:
		templine= string.split(line)
		master=np.append(master,[[float(i) for i in templine]],axis=0)

dat.close()"""
print "beginning to process EAMSD"
maxdt=max(master[:,0])
mindt=min(master[:,0])
templist=[]
dat=open(output[:-4]+'_EA.txt','w')
dat.close()
dat=open(output[:-4]+'_EA.txt','a')
for i in range(int(mindt),int(maxdt)):
	templist.append([i,np.average(master[np.where(master[:,0]==i),1])])
	print "Tau = "+str(i)+" of "+str(maxdt)+" is done!"
for points in templist:
	dat.write(str(points[0])+'     '+str(points[1])+'\n')
dat.close()
	