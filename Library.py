from __future__ import print_function
import numpy as np
import random
from sys import argv

#harmonic oscillator's force
def f_ext(x,K):
    f = -K*x
    return f
#Lennard-Jones' potential
def Flj(x,y,z,A,B):
    #lennard jones constant
    A = 4.*e*12.*d**12
    B = 4.*e*6.*d**6
    fx = np.zeros(len(x))
    fy = np.zeros(len(x))
    # fz = np.zeros(len(x))
    for a in range(len(x)):
        for b in range(len(x)):
            if a < b:      
                r2 = ((x[a]-x[b])**2)+((y[a]-y[b])**2)
                # +((z[a]-z[b])**2)    
                r4 = r2*r2
                r8 = r4*r4
                r14 = r2*r4*r8
                fx[a] += (A/r14-B/r8)*(x[a]-x[b])
                fx[b] += -fx[a]
                fy[a] += +(A/r14-B/r8)*(y[a]-y[b])
                fy[b] += -fy[a]
                # fz[a] += +(A/r14-B/r8)*(z[a]-z[b])
                # fz[b] += -fz[a]
    return fx,fy,
    #fz
# Total energy
def E(x,y,z,vx,vy,vz,K,m,n,Alj, Blj):
    #lennard jones constant
    A = 4.*e*12.*d**12
    B = 4.*e*6.*d**6
    # Energy calculus constants
    Alj = A/12.
    Blj = B/6.
    # kinetic energy
    ke = (.5*(np.dot(vx*m,vx)+np.dot(vy*m,vy)))
    # Harmonic Oscillator Potential Energy
    hope = 0.5*K*(np.dot(x,x)+np.dot(y,y))
    #Lennard Jones
    vlj = 0.
    for a in  range(n): 
        for b in range(a+1,n):
            r2 = ((x[a]-x[b])**2)+((y[a]-y[b])**2)
            #+((z[a]-z[b])**2)  
            r4 = r2*r2
            r6 = r4*r2
            r12 = r6*r6
            vlj += (Alj/r12-Blj/r6)
    #total energy
    Et = ke+hope+vlj
    #ptential energy
    pe = (hope+vlj)
    return ke,hope,vlj,Et

def D(ra,rb):
    #distances computation with [x,y,z] arrays
    di = np.sqrt((ra[0]-rb[0])**2+(ra[1]-rb[1])**2+(ra[2]-rb[2])**2)
    return di

def M(a,b):
    #distances computation from cartesian center
    A = np.sqrt(a**2+b**2)
    return A

def checkL(l,x):
        # if a given value x is not in list l
        #then true
        return_val = True
        for i in l:
                if i == x:
                        return_val = False
        return return_val

def checkH(l,x):
    #if 2 list have less than 3 common values
    #then true
    a = 0    
    return_val = True
    for i in l:
        for j in x :   
            if i == j:
                a += 1
    if a > 2:
        return_val = False
    return return_val

def tpD(n,x,y,z):
# test pair distances
    rab_list = []
    for a in range(n):
      for b in range(a+1,n):
          rab = np.sqrt((x[a]-x[b])**2+(y[a]-y[b])**2+(z[a]-z[b])**2)
          rab_list.append(rab)
    rab_list.sort()
    rab_list = [ '%.2f' % elem for elem in rab_list ]
    print (rab_list)
    return
