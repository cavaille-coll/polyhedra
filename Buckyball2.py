# this code's aim is the study of the C60 Buckminsterfullerene,
# a.k.a. Buckyball 

# the Buckyball function 
# generate the coordinates of the truncated icosahedron's 60 vertices
# each vertices is equvalent to a carbon

# the E function giv the whole system's Lennard-Jones potential
# the following functions compute icosahedral symmetry group's permutations

import numpy as np
from Library import *
#e=edge length
def Buckyball(e):
#according to Wikipedia,
#the coordinates of a edge length 2 truncated icosahedron are
#(0, +/-1, +/-3phi)
# (+/-2, +/-(1 + 2phi), +/-phi)
# (+/-1, +/-(2 + phi), +/-2phi)
#permuted 

    #triplet list: [0]=x,[1]=y,[2]=z
    w = [0,0,0]
    #use of the list's being periodic backward
    #e.g: w[-1]= w[2]
    
    #index= particle identifier
    index = 0
    # computing (0, +/-1, +/-3phi)
    #positionning the first number
    for i in range(3):
        w[i] = 0
        #2nd number position and value 
        #j=0/1, hence -1**j = 1/-1
        for j in range(2):
            w[i-1] = 3*phi*((-1)**j)  
            #ibidem
            for k in range(2):
                w[i-2] = ((-1)**k)
                x60[index]=w[0]
                y60[index]=w[1]
                z60[index]=w[2]
                index += 1
    #computing (+/-2, +/-(1 + 2phi), +/-phi)           
    for i in range(3):
        for h in range(2):
            w[i] = phi*((-1)**h)
            for j in range(2):
                w[i-1] = (1+2*phi)*((-1)**j)
                for k in range(2):
                    w[i-2] = 2*((-1)**k)
                    x60[index] = w[0]
                    y60[index] = w[1]
                    z60[index] = w[2]
                    index += 1
    #computing(+/-1, +/-(2 + phi), +/-2phi)
    for i in range(3):
        for h in range(2):
            w[i] = ((-1)**h)
            for j in range(2):
                w[i-1] = 2*phi*((-1)**j)
                for k in range(2):
                    w[i-2] = (2+phi)*((-1)**k)
                    x60[index] = w[0]
                    y60[index] = w[1]
                    z60[index] = w[2]
                    index+=1
    # adjusting to an edge of length 1
    # e become a coefficient
    for a in range(n):
        x60[a] = x60[a]/2*e
        y60[a] = y60[a]/2*e
        z60[a] = z60[a]/2*e
    return x60,y60,z60
    
#computing Lennard-Jones potential
def E(x,y,z,d,e,n):
    # groupment of the Lennard-Jones formula's constants
    A = 4.*e*12.*d**12
    B = 4.*e*6.*d**6
    Alj = A/12.
    Blj = B/6.
    # defining Lennard-Jones potential
    vlj = 0.
    for a in  range(n): 
        for b in range(a+1,n):
            #r= distance a-b
            #rn=r**n
            r2 = ((x60[a]-x60[b])**2)+((y60[a]-y60[b])**2+(z60[a]-z60[b])**2)
            r4 = r2*r2
            r6 = r4*r2
            r12 = r6*r6
            vlj += Alj/r12-Blj/r6
    ##if you want
    # print(vlj)
def mirror_reflection(r,b1,index):
#according to the icosahedral symmetry
#each Buckyball's mirror pass trough the middle of 2 opposite edges
# wich belong to 2 hexagon
#the plane is also perpendicular to the edges
#exploiting this property, we search an egdge a-b belonging to 2 hexagons
    d1 = 0
    #searching 4 points forming a central vertex (a) and its surrounfing
    for a in range(n):  
        #retroaction: avoid to compute for a then b as central vertex
        if all(b1 != a):
            # print('a    ',a)
            for b in range(n):
                # D function = distance a-b
                if((a != b) and D(r[0][a],r[0][b])) == edge: 
                    # print('b',b)
                    for c in range(n):            
                        if((a != c) and (b != c) and (d1 != c) 
                        and D(r[0][a],r[0][c])) == edge:  
                            # print('c',c)
                            for d in range(n):            
                                if((a != d) and (b != d) and (c != d) 
                                and D(r[0][a],r[0][d])) == edge:  
                                    # print('d',d)
                                    #normalize the following vectors
                                    ab=r[0][a]-r[0][b]
                                    ac=r[0][a]-r[0][c]
                                    ad=r[0][a]-r[0][d]
                                    #compute angle BAC
                                    cosine_angle1 = (np.dot(ab, ac) /
                                                     (np.linalg.norm(ab) *
                                                     np.linalg.norm(ac)))
                                    BAC = np.arccos(cosine_angle1)
                                    #compute angle BAd
                                    cosine_angle2 = (np.dot(ab, ad) /
                                                     (np.linalg.norm(ab) *
                                                     np.linalg.norm(ad)))
                                    BAD = np.arccos(cosine_angle2)                   
                                
                                    if BAC == BAD:
                                    #then ac and ad belong to a pentagon 
                                    #and consequently ab to 2 hexagons
                                        
                                    #assigns current a,b values to retroaction
                                        d1 = d
                                        b1[index] = b
                                        index += 1
                                        # print('yes')

                                        #compute edge's middle point
                                        x = (x60[a]+x60[b])/2
                                        y = (y60[a]+y60[b])/2
                                        z = (z60[a]+z60[b])/2
                                        
                                        #proceeding to the reflexion
                                        for e in (range(n)):
                                            #plane' polar coordinate 
                                            phi = np.arctan2(y,x)
                                            #plane' azimutal coordinate 
                                            nu = np.arctan2(y,z)
                                            #point's coordinates
                                            alpha = np.arctan2(y60[e],x60[e])
                                            beta = np.arctan2(y60[e],z60[e])
                                            #rotating forward to get an horizontal plane
                                            #reflecting the point
                                            xt[e] = M(x60[e],y60[e])*np.cos(alpha+2*(phi-alpha))
                                            yt[e] = M(x60[e],y60[e])*np.sin(alpha+2*(phi-alpha))
                                            zt[e] = M(z60[e],y60[e])*np.cos(beta+2*(nu-beta))
                                            ydz[e] = M(z60[e],y60[e])*np.sin(beta+2*(nu-beta))
                                            #rotating backward to get the actual reflection
                                            xf[e] = M(xt[e],yt[e])*np.cos(alpha)
                                            yf[e] = M(xt[e],yt[e])*np.sin(alpha)
                                            zf[e] = M(zt[e],ydz[e])*np.cos(beta)

                                        #assigning coordinate set to the 3d array
                                        for e in (range(n)):
                                            r[index+1][e] = [xf[e],yf[e],zf[e]]                                                       
    #returning coordinates values
    return r,b1,index,r,xf,yf,zf
 
def vertex_rotation(r,b1,index):
#the Buckyball rotation aroud its pentagons centers (icosahedron's vertices)
#is symmetric for each 2pi/5 giration (4 possible)

#searching a pentagon as a suite of 5 points 
#in which each one is adjacent to the next
#and the last to the first


    #like b1 previously, avoid to reuse the same pentagon
    points = []


    for a in range(n):  
        if checkL(points,a):
            for b in range(n):
                if((a != b) and checkL(points,b)
                and D(r[0][a],r[0][b])) == edge:
                    for c in range(n):
                        if((a != c) and (b != c) and checkL(points,c)
                        and D(r[0][b],r[0][c])) == edge:  
                            for d in range(n):
                                if((a != d) and (b != d) and checkL(points,d)
                                and (c != d) and D(r[0][c],r[0][d])) == edge:  
                                    for e in range(n):
                                        if((a != e) and (b != e) and (c != e)
                                        and checkL(points,e) and ( d!=e )
                                        and (D(r[0][d],r[0][e])) ==
                                             D(r[0][e],r[0][a]) == edge):
                                            
                                            points.append(a)
                                            points.append(b)
                                            points.append(c)
                                            points.append(d)
                                            points.append(e)
                                            
                                            for g in range(1,5):  
                                                index+=1
                                                #defining pentagon's center
                                                x = (x60[a]+x60[b]+x60[c]+x60[d]+x60[e])/5
                                                y = (y60[a]+y60[b]+y60[c]+y60[d]+y60[e])/5
                                                z = (z60[a]+z60[b]+y60[c]+y60[d]+y60[e])/5
                                                    #record the pentagon for retroaction
                                                    
                                                #proceeding to the giration
                                                for f in (range(n)):
                                                    #like mirror_reflection
                                                    phi = np.arctan2(y,x)
                                                    nu = np.arctan2(y,z)
                                                    #giration angle
                                                    delta = 2*np.pi/5
                                                    alpha = np.arctan2(y60[f],x60[f])
                                                    beta = np.arctan2(y60[f],z60[f])

                                                    xt[f] = M(x60[f],y60[f])*np.cos(alpha-phi+delta*g)
                                                    yt[f] = M(x60[f],y60[f])*np.sin(alpha-phi+delta*g)
                                                    zt[f] = -M(z60[f],y60[f])*np.cos(beta-nu+delta*g)
                                                    ydz[f]= M(z60[f],y60[f])*np.sin(beta-nu+delta*g)

                                                    xf[f] = M(xt[f],yt[f])*np.cos(alpha)
                                                    yf[f] = M(xt[f],yt[f])*np.sin(alpha)
                                                    zf[f] = M(zt[f],ydz[f])*np.cos(beta)
                                                for e in (range(n)):
                                                    r[index+1][e] = [xf[e],yf[e],zf[e]]
    return r,b1,index,g,xf,yf,zf

def face_rotation(r,b1,index):
#the Buckyball rotation aroud its hexagons centers (icosahedron'sface center)
#is symmetric for each 2pi/3 giration (+/-2 possible)

#searching a hexagon as a suite of 6 points 
#in which each one is adjacent to the next
#and the last to the first
    
    points = []

    for a in range(n):
        print 'a  ', a
        for b in range(n):
            if((a != b) and D(r[0][a],r[0][b])) == edge:
                print'b', b
                for c in range(n):            
                    if((a != c) and (b != c) and D(r[0][b],r[0][c])) == edge:
                        print'c',c
                        for d in range(n):            
                            if((a != d) and (b != d) and (c != d)
                            and D(r[0][c],r[0][d])) == edge:  
                                print 'd', d
                                for e in range(n):            
                                    if((a != e) and (b != e) and (c != e)
                                    and (d != e) and (D(r[0][d],r[0][e])) == 
                                                      edge):
                                        print 'e',e
                                        for f in range(n):            
                                            if((a != f) and (b != f) 
                                            and (c != f) and (d != f)
                                            and (e != f)
                                            and (D(r[0][e],r[0][f])) == 
                                                 D(r[0][f],r[0][a]) == edge):
                                                print'f', f
                                                #defining hexagon's center
                                                x = (x60[a]+x60[b]+x60[c]+x60[d]+x60[e]+x60[f])/6
                                                y = (y60[a]+y60[b]+y60[c]+y60[d]+y60[e]+y60[f])/6
                                                z = (z60[a]+z60[b]+y60[c]+y60[d]+y60[e]+z60[f])/6
                                                
                                                w = [a,b,c,d,e,f]
                                                if checkH(points,w):
                                                    print w
                                                    points.append(a)
                                                    points.append(b)
                                                    points.append(c)
                                                    points.append(d)
                                                    points.append(e)
                                                    points.append(f)
                                                    
                                                    for g in range(1,3):
                                                        index += 1
                                                        #proceeding to the giration
                                                        for h in (range(n)):
                                                            #like mirror_reflection
                                                            phi = np.arctan2(y,x)
                                                            nu = np.arctan2(y,z)
                                                            #giration angle
                                                            delta = 2*np.pi/3
                                                            alpha = np.arctan2(y60[h],x60[h])
                                                            beta = np.arctan2(y60[h],z60[h])
                                                            

                                                            xt[h] = M(x60[h],y60[h])*np.cos(alpha-phi+delta*g)
                                                            yt[h] = M(x60[h],y60[h])*np.sin(alpha-phi+delta*g)
                                                            zt[h] = M(z60[h],y60[h])*np.cos(beta-nu+delta*g)
                                                            ydz[h] = M(z60[h],y60[h])*np.sin(beta-nu+delta*g)

                                                            xf[h] = M(xt[h],yt[h])*np.cos(alpha)
                                                            yf[h] = M(xt[h],yt[h])*np.sin(alpha)
                                                            zf[h] = M(zt[h],ydz[h])*np.cos(beta)
                                                        for e in (range(n)):
                                                            r[index+1][e] = [xf[e],yf[e],zf[e]]
                                                        print 'yes'

    return r,b1,index,g,xf,yf,zf

# "de divina proportione", phi, or golden ratio
phi = (1+np.sqrt(5))/2 
#particle number
n = 60 
#3d array containing each set of 60 coordinate triplets(xyz)
r = np.zeros((130,n,3,),float)

x60 = np.zeros(n,float)
y60 = np.zeros(n,float)
z60 = np.zeros(n,float)

# x axis reflection
xf = np.zeros(n,float)
yf = np.zeros(n,float)
zf = np.zeros(n,float)

b1 = np.zeros(60,float)

xt = np.zeros(n,float)
yt = np.zeros(n,float)
zt = np.zeros(n,float)
ydz = np.zeros(n,float)
 
#defining the xyz file
fout = open("Bl.xyz", 'w')
fout.write('60''\n')
fout.write('Buckyball'+'\n')

#lennard jones parameters
e = 10. #epsilon
d = 1. #sigma

#supposed edge length of minimal energy
edge = 2.**(1./6.)*d 
#radius length as a function of edge
radius = edge*np.sqrt(1+9*phi**2)/2 
#function call
x60,y60,z60 = Buckyball(edge)

#E(x60,y60,z60,d,e,n)

# for a in range(n):
#    fout.write("Ar {} {} {} \n".format(x60[a],y60[a],z60[a]))
# fout.close()

#coordinates inversion
xrxyz =- x60
yrxyz =- y60
zrxyz =- z60

#assigning the coordinate set to the 3d array
for a in range(n):
    r[0][a] = [x60[a],y60[a],z60[a]]
for a in range(n):
    r[1][a] = [xrxyz[a],yrxyz[a],zrxyz[a]]

index = 0
#previous "a" value

r,b1,index,r,xf,yf,zf = mirror_reflection(r,b1,index)    

print(index)
index = 0
r,b1,index,g,xf,yf,zf = vertex_rotation(r,b1,index)

print(index)
index = 0
r,b1,index,g,xf,yf,zf = face_rotation(r,b1,index)
print(index)
for a in range(n):
    fout.write("Ar {} {} {} \n".format(xf[a],yf[a],zf[a]))
fout.close()