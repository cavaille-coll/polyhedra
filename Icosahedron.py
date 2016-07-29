import numpy as np
def icosahedron(e):
    #triplet list: [0]=x,[1]=y,[2]=z
    w=[0,0,0]
    #use of the list's being cyclic backward
    index=0
    for i in range(3):
        w[i]=0
        for j in range(2):
            w[i-1]=((-1)**j)
            for k in range(2):
                w[i-2]=phi*((-1)**k)
                x[index]=w[0]
                y[index]=w[1]
                z[index]=w[2]
                index+=1
    for a in range(n):
        x[a]=x[a]/2*e
        y[a]=y[a]/2*e
        z[a]=z[a]/2*e
    return x,y,z
    
def E(x,y,z,d,e,n):
    #Lennard-Jones formula constants
    A=4.*e*12.*d**12
    B=4.*e*6.*d**6
    Alj=A/12.
    Blj=B/6.
    #Lennard-Jones potential
    vlj=0.
    for a in  range(n): 
        for b in range(a+1,n):
            r2= ((x[a]-x[b])**2)+((y[a]-y[b])**2+(z[a]-z[b])**2)  
            r4=r2*r2
            r6=r4*r2
            r12=r6*r6
            vlj+=Alj/r12-Blj/r6
    return vlj
#file definition
fout = open("Ih.xyz", 'w')
fout.write('12''\n')
fout.write('Icosahedron'+'\n')
#golden mean
phi=(1+np.sqrt(5))/2

#vertices/particles number
n=12

#vertices/particles coordinates array 
x=np.zeros(n,float)
y=np.zeros(n,float)
z=np.zeros(n,float)

#lennard jones parameters
e=10.
d=1.

#supposed edge length of minimal energy
edge=2.**(1./6.)*d
x,y,z=icosahedron(edge)
vlj=E(x,y,z,d,e,n)

for a in range(n):
       fout.write("Ar {} {} {} \n".format(x[a],y[a],z[a]))
fout.close()
