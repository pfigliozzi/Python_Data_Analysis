import string, Tkinter, tkFileDialog
from decimal import *

root = Tkinter.Tk()
root.withdraw()
file_path = tkFileDialog.askopenfilename()
print "Opened "+file_path
particle=0
trajdb={}
file_save=tkFileDialog.asksaveasfilename(initialfile=file_path)
print "Will save as "+file_save
pxcutoff=input('Pixel cutoff around edges=')
imgsize=input('Pixel image size (x,y)')
for line in open(file_path,'r').readlines():
    if line.startswith('f'):
        continue
    else:
        particle+=1
        templine=string.split(line)
        intline=[Decimal(i) for i in templine]
        spt=[(particle,intline[0],(len(intline)-1)/3)]
        for k in range(1,len(intline),3):
            spt.append((intline[k],intline[k+1]))
        trajdb[str(particle)]=spt
        #print 'Particle '+str(particle)+' is done!'
trajdb['range']=[1,particle]
print str(particle)+" trajectories analyzed"
innertrajs=[]
for traj in range(1,trajdb['range'][1]):
    singletraj=trajdb[str(traj)]
    trajswitch=0
    for pos in singletraj[1:]:
        if pos[0]<pxcutoff or pos[1]<pxcutoff or pos[0]>(imgsize[0]-pxcutoff) or pos[1]>(imgsize[1]-pxcutoff):
            print 'Particle '+str(singletraj[0][0])+'/'+str(particle)+' is too close to edge!'
            trajswitch=1            
            break
        else:
            continue
    if trajswitch==1:
        continue
    else:
        innertrajs.append(singletraj)
        #print 'Particle '+str(singletraj[0][0])+' is okay'
out=open(file_save,'w')
out=open(file_save,'a')
for traj in innertrajs:
    out.write(str.rjust(str(traj[0][1]),6))    
    for pos in traj[1:]:
        #line=repr(pos[0]).rjust(13),repr(pos[1]).rjust(13),repr(1.00).rjust(13)
        out.write(str.rjust(str(pos[0]),13)+str.rjust(str(pos[1]),13)+str.rjust(str(1.00),13))
    out.write('\n')
    print 'Traj '+str(traj[0][0])+'/'+str(particle)+' is done'
out.close()