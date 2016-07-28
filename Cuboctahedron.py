from __future__ import print_function
import numpy as np

def Cuboctahedron(e):
    # e= edge length
    #according to Wikipedia, the coordinates of a cuboctahedron of edge length  are
    m=e*np.sqrt(2)/2.

    w=[0,0,0]
    #use of the list's being cyclic backward
    index=0
    for i in range(3):
        w[i]=0
        for j in range(2):
            w[i-1]=((-1)**j)*m      
            for k in range(2):
                w[i-2]=((-1)**k)*m
                # fout.write('Ar'+str(w)+'\n')
                # print(str(w)+' \n')
                x[index]=w[0]
                y[index]=w[1]
                z[index]=w[2]
                index+=1
    
    for a in range(n):
            x[a]=x[a]*np.srt(2)/2*e
            y[a]=y[a]*np.srt(2)/2*e
            z[a]=z[a]*np.srt(2)/2*e

    for a in range(n):
        fout.write("Ar {} {} {} \n".format(x[a],y[a],z[a]))
    fout.close()
    return

def E(x,y,z,d,e,n,Alj, Blj):
    ## # kinetic energy
    ## ke=.5*m*(np.dot(vx,vx)+np.dot(vy,vy))
    #Lennard-Jones constants
    A=4.*e*12.*d**12
    B=4.*e*6.*d**6
    Alj=A/12.
    Blj=B/6.
    #Lennard-Jones potential
    vlj=0.
    for a in  range(n): 
        for b in range(a+1,n):
            #r= distance between 2 particles
            #rn=r**n
            r2= ((x[a]-x[b])**2)+((y[a]-y[b])**2+(z[a]-z[b])**2)  
            r4=r2*r2
            r6=r4*r2
            r12=r6*r6
            vlj+=Alj/r12-Blj/r6
 
    print(vlj)
    #
# coordinates file definition
fout = open("Co.xyz", 'w')
    fout.write('12''\n')
    fout.write('Cuboctahedron'+'\n')
#Total energy file definition
fout = open("harout.dat", 'w')
    # Eout = open("E.dat", 'w')

#number of vertices/particles
n=12
#coordinates arrays
x=np.zeros(n,float)
y=np.zeros(n,float)
z=np.zeros(n,float)

#lennard jones parameters
e=10.
d=1.


edge=2.**(1./6.)*d
Cuboctahedron(edge)
E(x,y,z,d,e,n,Alj,Blj)