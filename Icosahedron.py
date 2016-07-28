import numpy as np
phi=(1+np.sqrt(5))/2
n=12
x=np.zeros(n,float)
y=np.zeros(n,float)
z=np.zeros(n,float)
def icosahedron(e):
	m=e/2
	fout = open("Ih.xyz", 'w')
	fout.write('12''\n')
	fout.write('Icosahedron'+'\n')

	w=[0,0,0]
	#use of the list's being cyclic backward
	index=0
	for i in range(3):
		w[i]=0
		for j in range(2):
			w[i-1]=((-1)**j)
			for k in range(2):
				w[i-2]=phi*((-1)**k)
				# fout.write('Ar'+str(w)+'\n')
				# print(str(w)+' \n')
				x[index]=w[0]
				y[index]=w[1]
				z[index]=w[2]
				index+=1
	for a in range(n):
        x[a]=x[a]/2*e
        y[a]=y[a]/2*e
        z[a]=z[a]/2*e

	for a in range(n):
	 	fout.write("Ar {} {} {} \n".format(x[a],y[a],z[a]))
	fout.close()

	# # test pair distances
	# rab_list=[]
	# for a in range(12):
	# 	for b in range(a+1,12):
	# 		rab=np.sqrt((x[a]-x[b])**2+(y[a]-y[b])**2+(z[a]-z[b])**2)
	# 		rab_list.append(rab)
	# rab_list.sort()
	# rab_list= [ '%.2f' % elem for elem in rab_list ]
	# print rab_list
	# return
def E(x,y,z,d,e,n,Alj, Blj):
    vlj=0.
    for a in  range(n): 
        for b in range(a+1,n):
            r2= ((x[a]-x[b])**2)+((y[a]-y[b])**2+(z[a]-z[b])**2)  
            r4=r2*r2
            r6=r4*r2
            r12=r6*r6
            vlj+=Alj/r12-Blj/r6
 
    print(vlj)

#lennard jones parameters
e=10.
d=1.
A=4.*e*12.*d**12
B=4.*e*6.*d**6
Alj=A/12.
Blj=B/6.

edge=2.**(1./6.)*d
icosahedron(edge)
E(x,y,z,d,e,n,Alj, Blj)