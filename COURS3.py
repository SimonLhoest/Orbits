# -*- coding: utf-8 -*-
"""
Created on Tue Mar 28 14:04:37 2023
ALLER RETOUR JUPITER CIRULAIRE
@author: lhoes
"""

import numpy as np
import matplotlib.pyplot as plt

h=1 #day 0.1 c'est ultra good, mais long #euuuuuuuuuuuhhhhhhhhhh
G=(2.95824*10**-4)
A=3.27
mo=1

t=np.arange(0,36500+h,h)
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
   k2=f(u[i]+h/2*k1,t[i]+h/2)
   k3=f(u[i]+h/2*k2,t[i]+h/2)
   k4=f(u[i]+h*k3,t[i]+h)
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
u2=u.copy()
a2=a.copy()
e2=e.copy()

h=-h
print(h)
for i in range(t.shape[0]-1,0,-1):
    
    k1=f(u[i],t[i])
    k2=f(u[i]+h/2*k1,t[i]+h/2)
    k3=f(u[i]+h/2*k2,t[i]+h/2)
    k4=f(u[i]+h*k3,t[i]+h)
    u[i-1]=u[i]+h/6*(k1+2*k2+2*k3+k4)

    r=u[i,:3]
    rp=u[i,3:]
    nr=np.linalg.norm(r)
    nrp=np.linalg.norm(rp)
    a[i-1]=(2/nr-(nrp**2)/(G*mo))**(-1)
    et=np.cross(rp,np.cross(r,rp))/(G*mo)-r/nr
    e[i-1]=np.linalg.norm(et)
    #k=np.cross(r,rp)/(nr*nrp)
    #j[i-1]=np.arccos(k[2])*(180/np.pi)
        
#error[-1]=(np.abs(exact[-1]-y[-1]))/y[-1]
#plt.plot(t,error)
'''
plt.figure('plot')
plt.plot(u[:,0],u[:,1])
lim=[-5,5]
plt.xlim(lim)
plt.ylim(lim)
plt.grid()
'''
#%% ellipse
fig,axs=plt.subplots(2,1)
plt.grid(1)
axs[0].plot(u2[:,0],u2[:,1])
plt.grid(1)
axs[1].plot(u[:,0],u[:,1])
plt.grid(1)


#%% semi major et excentricity
fig, axs = plt.subplots(2,2)

axs[0,0].plot(t,a2)
axs[1,0].plot(t,e2)
axs[0,1].plot(t,a)
plt.plot(t,e)
#plt.plot(t,j)

#%% LOUPED
plt.grid(1)
plt.plot(t,a-a2)

'''
for 2 au (h,error)
0.1: <1.10**(-6)
1  : <8.10**(-6)  <--
10 : <10.10**(-5)

for 3.27: 5.10**(-5)
for 3.5 : 8.10**(-5)
'''

#%% NORME ERREUR position 0
print(np.linalg.norm(u2[0,:3]-u[0,:3]))

''' euhhhhhhhhhhhhhh le -h était inversé et c'est mieux
1au=150.10e6 km
for 2au (h,error)
10 : 0.0033198564512326523 3s~
1  : 3.495190251120644e-05 15s~ <--- 5.10e3 km
0.1: 2.430907514781819e-05 3min~

for 3.27: 0.00021103155320902287
for 3.5 : 0.00023643770062880856

'''


'''
1au=150.10e6 km
for 2au (h,error)
10 : 3.8899241058412355e-05 3s~   (0.0034069005562069685 ????????????)
1  : 0.00034601083704638983 15s~ 
0.5: 0.00017302224038868397 30s <--- 26k km
0.1: 3.460471077430032e-05 2-3min~

for 3.27: 0.0004710784192234184
for 3.5 : 0.0019754600696854038 (0.012691948492063827)

'''







