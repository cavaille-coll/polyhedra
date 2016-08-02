from __future__ import print_function
import numpy as np
from MMTK import *
from MMTK.Trajectory import Trajectory, SnapshotGenerator, TrajectoryOutput
from Library import *

# "de divina proportione"
phi=(1+np.sqrt(5))/2 
#particle number
n=13 
# step number
s=1000
 # step's length
dt=0.05
#pparticle parameters
x=np.zeros(n,float)
y=np.zeros(n,float)
z=np.zeros(n,float)
vx=np.zeros(n,float)  
vy=np.zeros(n,float)  
vz=np.zeros(n,float)  
m=np.zeros(n,float)
#lennard jones parameters
e=10.
d=1.
A=4.*e*12.*d**12
B=4.*e*6.*d**6
Alj=A/12.
Blj=B/6.
Eout = open("E.dat", 'w')
Eminout = open("Emin.dat", 'w')
coef=(np.sqrt(2)-1)/(phi-1)
edge=2.**(1./6.)*d
#transition steps
ts=0.01 
#sampling frequence
sa=2
K=0

def E(x,y,z,d,e,n,Alj, Blj,c):
    vlj=0.
    
    for a in  range(n): 
        for b in range(a+1,n):
            r2= ((x[a]-x[b])**2)+((y[a]-y[b])**2+(z[a]-z[b])**2)  
            r4=r2*r2
            r6=r4*r2
            r12=r6*r6
            vlj+=Alj/r12-Blj/r6
    return vlj

universe = InfiniteUniverse()

for a in range(n):
    m[a]=1
    universe.addObject(Atom('Ar', position=Vector(x[a],y[a],z[a])))

trajectory = Trajectory(universe, "poly.nc", "w", "many particles")

snapshot = SnapshotGenerator(universe,
                             actions = [TrajectoryOutput(trajectory,
                                                         ["all"], 0, None, 1)])

fout = open("Pl.xyz", 'w')
fout.write('13''\n')
fout.write('Polyhedra'+'\n')
cindex=-1
for c in np.arange(0.0,(phi-1+ts),ts):
    cindex+=1
    # if (c%((phi-1)/sa)) == 0:
    w=[0,0,0]
    #use of the list's being cyclic backward
    index=0
    #radius
    r=edge*(np.sqrt(2)-(c*coef))/2.
    for i in range(3):
        w[i]=0
        for j in range(2):
            w[i-1]=((-1)**j)*r
            for k in range(2):
                w[i-2]=(1+c)*((-1)**k)*r
                
                x[index]=w[0]
                y[index]=w[1]
                z[index]=w[2]                   
                index+=1
    x[12]=0
    y[12]=0
    z[12]=0


    Eout.write("{} {}\n".format(c,E(x,y,z,d,e,n,Alj,Blj,c)))
    Eout.flush()

    fx,fy,fz = Flj(x,y,z,A,B)

    # loop over time steps
    for i in range(s):
    #   calculate forces
    #velocity verlet
        for a in range(n):
            universe.atomList()[a].setPosition(Vector(x[a],y[a],z[a]))
        snapshot()
        Eminout.write("{} {} {} {}\n".format(c,i*dt,cindex*s+i,E(x,y,z,d,e,n,Alj,Blj,c)))
        Eminout.flush()    
        for a in range(n):
            # 1/2 velocity
            vx[a]=vx[a]+.5*(fx[a]/m[a])*dt
            vy[a]=vy[a]+.5*(fy[a]/m[a])*dt
            vz[a]=vy[a]+.5*(fz[a]/m[a])*dt
            x[a]=vx[a]*dt+x[a]
            y[a]=vy[a]*dt+y[a]
            z[a]=vx[a]*dt+z[a]
            universe.atomList()[a].setPosition(Vector(x[a],y[a],z[a]))

        fx,fy,fz = Flj(x,y,z,A,B)
        for a in range(n):    
            vx[a]=0
            vy[a]=0
            vz[a]=0.

            #coordinates's writing
            fout.write("{} {} {} {} {}".format(i*dt, x[a], y[a], vx[a], vy[a]))
        

        fout.write('\n') 

fout.close()
Eout.close()
Eminout.close()