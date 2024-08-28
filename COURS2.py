# -*- coding: utf-8 -*-
"""
Created on Tue Mar 14 13:55:03 2023

@author: lhoes
"""


import numpy as np
import matplotlib.pyplot as plt

h=1 #day
G=(2.95824*10**-4)
A=2
mo=1

t=np.arange(0,10000+h,h)
u=np.zeros((t.shape[0],6))
u[0]=np.array([A, 0, 0, 0, (G/A)**(1/2),0])

a=np.zeros((t.shape[0]))
e=np.zeros((t.shape[0]))
j=np.zeros((t.shape[0]))
a[0]=A
e[0]=0
j[0]=0

def f(u,t):
    return np.array([u[3], u[4], u[5], -G*mo/((u[0]**2+u[1]**2+u[2]**2)**(3/2))*u[0], -G*mo/((u[0]**2+u[1]**2+u[2]**2)**(3/2))*u[1], -G*mo/((u[0]**2+u[1]**2+u[2]**2)**(3/2))*u[2]])

ju=np.array([0, 0, 0]) #Jupiter traj
mj=1/1047.348625
T=4331.797849

def ff(u,t):
    ju=np.array([5.2*np.cos(2*np.pi*t/T), 5.2*np.sin(2*np.pi*t/T), 0])
    return np.array([u[3], u[4], u[5], 
                     -G*mo/((u[0]**2+u[1]**2+u[2]**2)**(3/2))*u[0]-G*mj*((u[0]-ju[0])/(((u[0]-ju[0])**2+(u[1]-ju[1])**2+(u[2]-ju[2])**2)**(3/2))+ju[0]/(ju[0]**2+ju[1]**2+ju[2]**2)**(3/2)), 
                     -G*mo/((u[0]**2+u[1]**2+u[2]**2)**(3/2))*u[1]-G*mj*((u[1]-ju[1])/(((u[0]-ju[0])**2+(u[1]-ju[1])**2+(u[2]-ju[2])**2)**(3/2))+ju[1]/(ju[0]**2+ju[1]**2+ju[2]**2)**(3/2)),  
                     -G*mo/((u[0]**2+u[1]**2+u[2]**2)**(3/2))*u[2]-G*mj*((u[2]-ju[2])/(((u[0]-ju[0])**2+(u[1]-ju[1])**2+(u[2]-ju[2])**2)**(3/2))+ju[2]/(ju[0]**2+ju[1]**2+ju[2]**2)**(3/2))])

f=ff
#exact=np.exp(-t)
#error=np.zeros(t.shape)

for i in range(t.shape[0]-1):
   k1=f(u[i],t[i])
   k2=f(u[i]+h/2*k1,t[i])
   k3=f(u[i]+h/2*k2,t[i])
   k4=f(u[i]+h*k3,t[i])
   u[i+1]=u[i]+h/6*(k1+2*k2+2*k3+k4)
   
   r=u[i,:3]
   rp=u[i,3:]
   nr=np.linalg.norm(r)
   nrp=np.linalg.norm(rp)
   a[i+1]=(2/nr-(nrp**2)/(G*mo))**(-1)
   et=np.cross(rp,np.cross(r,rp))/(G*mo)-r/nr
   e[i+1]=np.linalg.norm(et)
   k=np.cross(r,rp)/(nr*nrp)
   j[i+1]=np.arccos(k[2])*(180/np.pi)
   
   #error[i]=(np.abs(exact[i]-y[i]))/y[i]

#error[-1]=(np.abs(exact[-1]-y[-1]))/y[-1]
#plt.plot(t,error)

plt.figure('plot')
plt.plot(u[:,0],u[:,1])
lim=[-5,5]
plt.xlim(lim)
plt.ylim(lim)
plt.grid()
#%%

dist=(u[:,0]**2+u[:,1]**2)**(1/2)
plt.plot(t,dist)
plt.grid()

#%%
'''
r=u[:,:3]
rp=u[:,3:]
nr=np.linalg.norm(r)
nrp=np.linalg.norm(rp)
a=(2/nr-(nrp**2)/(G*mo))**(-1)
e=np.cross(rp,np.cross(r,rp))/(G*mo)-r/nr
e=np.linalg.norm(e)
k=np.cross(r,rp)/(nr*nrp)
i=np.arccos(k[:,2])*180/np.pi
'''
plt.figure('a')
plt.plot(t[1:],a[1:])
plt.figure('e')
plt.plot(t[1:],e[1:])
plt.figure('i')
plt.plot(t,j)











