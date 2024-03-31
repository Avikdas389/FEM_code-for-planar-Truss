# -*- coding: utf-8 -*-
"""
Created on Sun Dec  3 15:07:24 2023

@author: user
"""

import numpy as np
import matplotlib.pyplot as plt

NNodes=np.array([[0,0],[80,0],[70,20],[60,40],[50,60],[30,60],[20,40],[10,20]]) 
NElems=np.array([[1,2],[2,3],[3,4],[4,5],[5,6],[6,7],[7,0],[0,2],[2,6],[4,6],[1,7],[3,7],[3,5]]) 
# near and far node numbers for elements
#Node= 8
#Elem= 13
Node= int(input("Enter the number of Nodes :")) 
Elem= int(input("Enter the number of elements :")) 
K=np.zeros([2*Node,2*Node]) 
A = float(input("Enter the area of the cross section :"))
E = float(input("Enter the Elastic Modulus :"))
F=np.array([0,0,0,0,0,0,0,0,0,0,7000,0,8000,0,9000,0])#*****Loading*****#
t=len(F)
print(t)
#******PLotting of Structure *********
x=np.linspace(0,0,Node)
y=np.linspace(0,0,Node)

for i in range(0,Node):
  x[i] = NNodes[i,0]
  y[i] = NNodes[i,1]


node0 = [x[0],y[0]]
node1 = [x[1],y[1]]
node2=[x[2],y[2]]

for i in range(1,Node-1):
    x01 = [x[i], x[i+1]]
    y01 = [y[i], y[i+1]]
    plt.plot(x01, y01, 'bo', linestyle="-")

x01 = [x[2], x[0]]
y01 = [y[2], y[0]]
plt.plot(x01, y01, 'bo', linestyle="-")

x01 = [x[2], x[0]]
y01 = [y[2], y[0]]
plt.plot(x01, y01, 'bo', linestyle="-")

x01 = [x[2], x[6]]
y01 = [y[2], y[6]]
plt.plot(x01, y01, 'bo', linestyle="-")

x01 = [x[4], x[6]]
y01 = [y[4], y[6]]
plt.plot(x01, y01, 'bo', linestyle="-")

x01 = [x[1], x[7]]
y01 = [y[1], y[7]]
plt.plot(x01, y01, 'bo', linestyle="-")

x01 = [x[3], x[7]]
y01 = [y[3], y[7]]
plt.plot(x01, y01, 'bo', linestyle="-")

x01 = [x[3], x[5]]
y01 = [y[3], y[5]]
plt.plot(x01, y01, 'bo', linestyle="-")

x01 = [x[7], x[0]]
y01 = [y[7], y[0]]
plt.plot(x01, y01, 'bo', linestyle="-")


#******Matrices********
c=np.linspace(0,0,Elem)
s=np.linspace(0,0,Elem)
l=np.linspace(0,0,Elem)
for i in range(Elem):
  n=NElems[i,0]
  f=NElems[i,1]
  l[i]=np.sqrt((NNodes[f,0]-NNodes[n,0])**2+(NNodes[f,1]-NNodes[n,1])**2)
  c[i]=(NNodes[f,0]-NNodes[n,0])/l[i]
  s[i]=(NNodes[f,1]-NNodes[n,1])/l[i]
  # the dofs for any node with node number 'n' is 2n and 2n+1, n=0,1,2...
  K[2*n,2*n]=K[2*n,2*n]+c[i]**2/l[i]
  K[2*n,2*n+1]=K[2*n,2*n+1]+c[i]*s[i]/l[i]
  K[2*n,2*f]=K[2*n,2*f]-c[i]**2/l[i]
  K[2*n,2*f+1]=K[2*n,2*f+1]-c[i]*s[i]/l[i]

  K[2*n+1,2*n+1]=K[2*n+1,2*n+1]+s[i]**2/l[i]
  K[2*n+1,2*f]=K[2*n+1,2*f]-c[i]*s[i]/l[i]
  K[2*n+1,2*f+1]=K[2*n+1,2*f+1]-s[i]**2/l[i]
  K[2*n+1,2*n]=K[2*n+1,2*n]+c[i]*s[i]/l[i]


  K[2*f,2*f]=K[2*f,2*f]+c[i]**2/l[i]
  K[2*f,2*f+1]=K[2*f,2*f+1]+c[i]*s[i]/l[i]
  K[2*f,2*n]=K[2*f,2*n]-c[i]**2/l[i]
  K[2*f,2*n+1]=K[2*f,2*n+1]-c[i]*s[i]/l[i]


  K[2*f+1,2*f+1]=K[2*f+1,2*f+1]+s[i]**2/l[i]
  K[2*f+1,2*n]=K[2*f+1,2*n]-c[i]*s[i]/l[i]
  K[2*f+1,2*n+1]=K[2*f+1,2*n+1]-s[i]**2/l[i]
  K[2*f+1,2*f]=K[2*f+1,2*f]+c[i]*s[i]/l[i]

print(K*A*E)

#******Boundary Conditions ********
user_input = input("Enter Degree of freedom where Displacent is zero : ") #0 1 2 3
d = user_input.split()
d = [int(element) for element in d]
N_d= len(d)

for i in range(N_d):
    s= d[i]-i
    
    K=np.delete(K, s, 0)  
    K=np.delete(K, s, 1)
    F=np.delete(F,s,0)


K_modified= K*A*E
print (K_modified)

#*********Solution************
All_DOF= np.zeros(2*Node)
Q= np.linalg.solve(K_modified, F)
Qg=np.linspace(0,2*Node-1,2*Node) # all the dof indices
All_DOF=np.linspace(0,0,2*Node)# displacements at all the dofs
NZQ=np.setdiff1d(Qg,d)  # the array containing the indices of the free dofs
ZQ=np.setdiff1d(Qg,NZQ)  # the array containing the indices of the constrained dofs

nNZDOF= 2*Node - N_d

for i in range(nNZDOF):
  j=int(NZQ[i])
  All_DOF[j]=Q[i]
print('displacement at all dofs')
print(All_DOF)

#********Plotting deformed Stucture**
p=np.linspace(0,0,Node)
q=np.linspace(0,0,Node)

for i in range(0,Node):
  for j in range(i*2,i*2+2):
   if j % 2 == 0:
      p[i] = NNodes[i,0]+ All_DOF[j]
   else:
     q[i] = NNodes[i,1]+ All_DOF[j]

#
for i in range(1,Node-1):
    p01 = [p[i], p[i+1]]
    q01 = [q[i], q[i+1]]
    plt.plot(p01, q01, 'red', linestyle="-.")


p01 = [p[2], p[0]]
q01 = [q[2], q[0]]
plt.plot(p01, q01, 'red', linestyle="-.")

p01 = [p[2], p[6]]
q01 = [q[2], q[6]]
plt.plot(p01, q01, 'red', linestyle="-.")


p01 = [p[4], p[6]]
q01 = [q[4], q[6]]
plt.plot(p01, q01, 'red', linestyle="-.")

p01 = [p[1], p[7]]
q01 = [q[1], q[7]]
plt.plot(p01, q01, 'red', linestyle="-.")

p01 = [p[3], p[7]]
q01 = [q[3], q[7]]
plt.plot(p01, q01, 'red', linestyle="-.")

p01 = [p[3], p[5]]
q01 = [q[3], q[5]]
plt.plot(p01, q01, 'red', linestyle="-.")

p01 = [p[7], p[0]]
q01 = [q[7], q[0]]
plt.plot(p01, q01, 'red', linestyle="-.")

# **Add text annotations at specific data points*** #


for i in range(0,Node):
    plt.text(p[i], q[i], f' {i}', ha='center', fontsize=18, color='black')
    
plt.show()    
#******End of code*******#
