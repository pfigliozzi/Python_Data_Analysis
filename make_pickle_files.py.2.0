import numpy, cPickle, fileinput, string, sys, heapq, os, shelve

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

filename = sys.argv[1]
#print Parse_Line(f)
#frameCount=2 #frameCount is going to scan over all the frames in a single movie file = (length of longest line-1)/3
mkdir=str.split(filename,'/')
try:
    os.mkdir(mkdir[0]+'/Pickles')
except OSError:
    pass
f=open(str(sys.argv[1]),'r')
output=str(mkdir[0]+"/Pickles/"+mkdir[1]+".trajectories")
trajdb=shelve.open(str(output))
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
trajdb.close()
#print frame
#output=str(mkdir[0]+"/Pickles/"+mkdir[1]+".pickle."+str(frameCount))
#output=open(mkdir[0]+"/Pickles/"+mkdir[1]+".pickle."+str(frameCount),'w')
#cPickle.dump(frame,output)
#print "Frame "+str(frameCount)+" is done!"
    


    


#Step 2: scan each line/list indexed by frame number. If frame number is between (first number of list) + (length of list) write nested list containing [[frame#-initial frame,[all trajectories]].
#Step 3: add control characters [( for example to separate particle trajectories in each line
#Step 4: cPickle the underlying data structure for each line
#Step 5: I think you would just add each line as a separate cPickle write call, but THIS SHOULD BE CHECKED


#Extraction from file
#Step 1: if memory is not an issue, unpickle the whole file. Alternatively, cpickle.load() line by line by making the first line the number of lines (frames) in the file.
# Step2: Once you have unpickled the line, do whatever test to it, run PIL, form the frame, etc...

