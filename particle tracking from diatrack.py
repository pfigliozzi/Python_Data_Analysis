#Copyright 2012 Thomas A Caswell
#tcaswell@uchicago.edu
#http://jfi.uchicago.edu/~tcaswell
#
#This program is free software; you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation; either version 3 of the License, or (at
#your option) any later version.
#
#This program is distributed in the hope that it will be useful, but
#WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
#General Public License for more details.
#
#You should have received a copy of the GNU General Public License
#along with this program; if not, see <http://www.gnu.org/licenses>.

from __future__ import division
#import matplotlib
#matplotlib.use('qt4agg')
import trackpy.tracking as pt
import matplotlib.pyplot as plt
import numpy as np
import sys, string

def numframes(xcoor,ycoor):
	for x in range(2):
		xline=string.split(xcoor.readline())
	for y in range(2):
		yline=string.split(ycoor.readline())
	if len(xline)==len(yline):
		return len(xline)
	else:
		print "X and Y coordinates not the same number of frames!"
		
def makecoorlist(numframes):
	mainlist=[]
	for fr in range(numframes):
		mainlist.append([])
	return mainlist
	

# generate fake data
levels = []
# 15 planes
for i in range(15):
    level = []
	# add the current level to the list of levels
    levels.append(level)
    #print level
	# a 15 by 15 grid
    for j in range(15):
        for k in range(15):
			# displace the location from the grid by a guassian with width 1/10
			level.append(pt.PointND(0,np.asarray((j+2,k+2))+np.random.randn(2)/10))
			#print np.random.randn(2)/10
			


xfil=open(str(sys.argv[1]),'r')
yfil=open(str(sys.argv[2]),'r')
stop=0
data=makecoorlist(numframes(open(str(sys.argv[1]),'r'),open(str(sys.argv[2]),'r')))
linenumber=1
for (linex,liney) in zip(xfil.readlines(),yfil.readlines()):
	if linex.startswith('f') or liney.startswith('f'):
		continue
	else:
		templinex=string.split(linex)
		templiney=string.split(liney)
		intlinex=[float(i) for i in templinex]
		intliney=[float(i) for i in templiney]
		for index in range(len(intlinex)):
			if intlinex[index]==0.00 and intliney[index]==0.00:
				continue
			else:
				data[index].append(pt.PointND(index+1,np.asarray((intlinex[index],intliney[index]))))
		print "Line "+str(linenumber)+" is done!"
		linenumber+=1
		#print [intlinex[0],intliney[0]]
		
#print data[0]
print "Tracking started!"
#print levels
# do the tracking
t = pt.link_full(data[:],(1930,1930),3,pt.Hash_table)

print "Linking is done!"
# plot tracks
fig = plt.figure()
ax = fig.gca()
out=open('trackeddata.txt','w')
out=open('trackeddata.txt','a')
for track in t:
     pos=track.points
     out.write(str(pos[0].t)+'    ')
     for p in pos:
         out.write(str(p.pos[1])+'    '+str(p.pos[0])+'    '+str(1.00)+'    ')
     out.write('\n')
out.close()
    
#plt.show()
