# -*- coding: utf-8 -*-
"""
Created on Sun Dec  3 08:07:15 2023

@author: user
"""
import numpy as np
import matplotlib.pyplot as plt

NNodes=np.array([[0,0],[40,0],[20,40]]) 
NElems=np.array([[0,1],[1,2],[0,2]]) 
# near and far node numbers for elements
#Node=3
#Elem=3
Node= int(input("Enter the number of Nodes :"))
Elem= int(input("Enter the number of elements :"))
K=np.zeros([2*Node,2*Node]) 
A = float(input("Enter the area of the cross section m^2 :"))
E = float(input("Enter the Elastic Modulus in Pascal :"))
F=np.array([0,0,0,0,40000,-50000])#*****Loading*****

Reac = [0]*Elem
#******PLotting of Structure*********
x=np.linspace(0,0,Node)
y=np.linspace(0,0,Node)

for i in range(0,Node):
  x[i] = NNodes[i,0]
  y[i] = NNodes[i,1]


node0 = [x[0],y[0]]
node1 = [x[1],y[1]]
node2=[x[2],y[2]]

for i in range(0,Node-1):
    x01 = [x[i], x[i+1]]
    y01 = [y[i], y[i+1]]
    plt.plot(x01, y01, 'bo', linestyle="-")




x01 = [x[2], x[0]]
y01 = [y[2], y[0]]
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

K0=K*A*E
#******Boundary Conditions ********
user_input = input("Enter Degree of freedom where Displacent is zero : ") #0 1 2 3
d = user_input.split()
d = [int(element) for element in d]
N_d= len(d)

for i in range(N_d):
    p= d[i]-i
    
    K=np.delete(K, p, 0)  
    K=np.delete(K, p, 1)
    F=np.delete(F,p,0)


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
prod=np.matmul(K0,All_DOF.T)
print(All_DOF)

#******Finding out the Reaction****#
for i in range(nNZDOF):
  j=int(NZQ[i])
  Reac[i] = prod[j]

print(Reac)    
#******Solving for axial force*****#

Ax=np.zeros([Elem,1]) # axial force vector

print(np.shape(K0))

print(type(c))  
print(type(s)) 

for i in range(Elem):
 B= np.array([-c[i],-s[i],c[i],s[i]])
 n=NElems[i,0]
 f=NElems[i,1]

 q=np.array([All_DOF[2*n],All_DOF[2*n+1],All_DOF[2*f],All_DOF[2*f+1]])
 Ax[i]=np.dot(B,q.T)/l[i]*A*E
print('Axial force',Ax)
print('Axial force')



#********Plotting deformed Stucture**
p=np.linspace(0,0,Node)
q=np.linspace(0,0,Node)

for i in range(0,Node):
  for j in range(i*2,i*2+2):
   if j % 2 == 0:
      p[i] = NNodes[i,0]+ All_DOF[j]
   else:
     q[i] = NNodes[i,1]+ All_DOF[j]


for i in range(0,Node-1):
    p01 = [p[i], p[i+1]]
    q01 = [q[i], q[i+1]]
    plt.plot(p01, q01, 'red', linestyle="-.")


p01 = [p[2], p[0]]
q01 = [q[2], q[0]]
plt.plot(p01, q01, 'red', linestyle="-.")

# **Add text annotations at specific data points*** #


for i in range(0,Node):
    plt.text(x[i], y[i], f' {i}', ha='center', fontsize=18, color='black')

#******End of code*******#


