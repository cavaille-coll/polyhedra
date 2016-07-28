from __future__ import print_function
import numpy as np
import random
from sys import argv

#harmonic oscillator's force
def f_ext(x,K):
    f=-K*x
    return f
#Lennard-Jones' potential
def Flj(x,y,z,A,B):
    fx = np.zeros(len(x))
    fy = np.zeros(len(x))
    fz = np.zeros(len(x))
    for a in range(len(x)):
        for b in range(len(x)):
            if a < b:      
                r2= ((x[a]-x[b])**2)+((y[a]-y[b])**2)+((z[a]-z[b])**2)    
                r4=r2*r2
                r8=r4*r4
                r14=r2*r4*r8
                fx[a] += (A/r14-B/r8)*(x[a]-x[b])
                fx[b] += -fx[a]
                fy[a] += +(A/r14-B/r8)*(y[a]-y[b])
                fy[b] += -fy[a]
                fz[a] += +(A/r14-B/r8)*(z[a]-z[b])
                fz[b] += -fz[a]
    return fx,fy,fz
# Total energy
def E(x,y,z,vx,vy,vz,K,m,n,Alj, Blj):
    # kinetic energy
    ke=(.5*(np.dot(vx*m,vx)+np.dot(vy*m,vy)))
    # Harmonic Oscillator Potential Energy
    hope=0.5*K*(np.dot(x,x)+np.dot(y,y))
    #Lennard Jones
    vlj=0.
    for a in  range(n): 
        for b in range(a+1,n):
            r2= ((x[a]-x[b])**2)+((y[a]-y[b])**2)+((z[a]-z[b])**2)  
            r4=r2*r2
            r6=r4*r2
            r12=r6*r6
            vlj+=(Alj/r12-Blj/r6)
    #total energy
    Et=ke+hope+vlj
    #ptential energy
    pe=(hope+vlj)
    return ke,hope,vlj,Et

#distances computation
def D(ra,rb):
    di=np.sqrt((ra[0]-rb[0])**2+(ra[1]-rb[1])**2+(ra[2]-rb[2])**2)
    return di

def mirror_reflection(r,a1,b1,index):
    for a in range(a1,n):  
        if(all(b1)!=a):
            for b in range(len(x60)):
                if((a!=b) and D(r[0][a],r[0][b]))==edge:                          
                    for c in range(len(x60)):            
                        if((a!=c)and(b!=c)and D(r[0][a],r[0][c]))==edge:  
                            for d in range(len(x60)):            
                                if((a!=d)and(b!=d)and(c!=d)and D(r[0][a],r[0][b]))==edge:  
                                    ab=r[0][a]-r[0][b]
                                    ac=r[0][a]-r[0][c]
                                    ad=r[0][a]-r[0][d]
                                    cosine_angle1 = np.dot(ab, ac) / (np.linalg.norm(ab) * np.linalg.norm(ac))
                                    BAC = np.arccos(cosine_angle1)
                                    cosine_angle2 = np.dot(ab, ad) / (np.linalg.norm(ab) * np.linalg.norm(ad))
                                    BAD = np.arccos(cosine_angle2)                   
                                
                                    if BAC==BAD:
                                        a1=a
                                        b1[index]=b
                                        index+=1

                                        x=(x60[a]-x60[b])/2
                                        y=(y60[a]-y60[b])/2
                                        z=(z60[a]-z60[b])/2
                                        for e in (range(n)):
                                            phi=np.arctan2(y,x)
                                            nu=np.arctan2(y,z)
                                            alpha=np.arctan2(y60[e],x60[e])
                                            beta=np.arctan2(y60[e],z60[e])
                                            
                                            xt[e]=M(x60[e],y60[e])*np.cos(alpha+2*(phi-alpha))
                                            yt[e]=M(x60[e],y60[e])*np.sin(alpha+2*(phi-alpha))
                                            zt[e]=M(z60[e],y60[e])*np.cos(beta+2*(nu-beta))
                                           
                                            xf[e]=M(xt[e],yt[e])*np.cos(alpha)
                                            yf[e]=M(xt[e],yt[e])*np.sin(alpha)
                                            zf[e]=M(zt[e],yt[e])*np.cos(beta)
                                        for e in (range(n)):
                                            r[index+1][e]=[xf[e],yf[e],zf[e]]
                                        print(index)
                                        return r,a1,b1,index,xf,yf,zf
def M(a,b):
    A=np.sqrt(a**2+b**2)
    return A