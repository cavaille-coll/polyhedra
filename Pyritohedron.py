import numpy as np

# pyritohedron = half irregular pentagonal dodecahedron 
# 2 groups of equal length edges :24 and 6 (30 total) 

def pyritohedron(e,h):
#Accordibg to Wikipedia
#The vertices of an pyritohedron with edge-length 2, centered at the origin,
# are described by all the cyclic permutations of:
#(0, +/-(1 + h), +/-(1 − h**2)
#and (+/-1, +/-1, +/-1), the cube's coordinates

# if h= 1 then rhombic dodecahedron
# if h= 0 then cube
# if h= 1/phi then regular dodecahedron
# e = edge length
    
    w = [0,0,0]
    #determination of coordinates
    #use of the list's being periodic backward
    index = 0
    for i in range(3):
        w[i] = 0
        for j in range(2):
            w[i-1] = (1+h)*((-1)**j)      
            for k in range(2):
                w[i-2] = (1-h**2)*((-1)**k)
                x[index] = w[0]
                y[index] = w[1]
                z[index] = w[2]
                index+=1
    #cubic coordinates
    for i in range(2):
        w[0] = ((-1)**i)  
        for j in range(2):
            w[1] = ((-1)**j)      
            for k in range(2):
                w[2] = ((-1)**k)
                x[index] = w[0]
                y[index] = w[1]
                z[index] = w[2]
                index += 1
    #output
    return x,y,z

#system's energy
def E(x,y,z,d,e,n,Alj, Blj):
    #lennard jones constants
    A = 4.*e*12.*d**12
    B = 4.*e*6.*d**6
    
    Alj = A/12.
    Blj = B/6.
    vlj = 0.
    for a in  range(n): 
        for b in range(a+1,n):
            r2 = ((x[a]-x[b])**2)+((y[a]-y[b])**2+(z[a]-z[b])**2)  
            r4 = r2*r2
            r6 = r4*r2
            r12 = r6*r6
            vlj += Alj/r12-Blj/r6
 
    return vlj
#golden mean
phi = (1+np.sqrt(5))/2
#vertices/particles number
n = 20
#vertices/particles coordinates array 
x = np.zeros(n,float)
y = np.zeros(n,float)
z = np.zeros(n,float)

#lennard jones parameters
e = 10.
d = 1.

edge = 2.**(1./6.)*d
#file definition
fout = open("Dh.xyz", 'w')
    fout.write('20''\n')
    fout.write('Pyritohedron'+'\n')

x,y,z = pyritohedron(1,1/phi)
vlj = E(x,y,z,d,e,n,Alj, Blj)
#file outpout
for a in range(n):
    fout.write("Ar {} {} {} \n".format(x[a],y[a],z[a]))
fout.close()