# -*- coding: utf-8 -*-
"""
Created on Tue Mar 28 14:04:37 2023
COMPARE JUPITER CIRCULAIRE ET ELLIPTIQUE
2007VW266
@author: lhoes
"""

import numpy as np
import matplotlib.pyplot as plt

h=1 #day 0.1 c'est ultra good, mais long #euuuuuuuuuuuhhhhhhhhhh
G=(2.95824*10**-4)
k=G**(1/2)
mo=1
epoch=2456600.5
t=np.arange(0+epoch,epoch+36500+h,h)
u=np.zeros((t.shape[0],6))

def R1(O):
    return np.array([[1, 0, 0],
                     [0, np.cos(O), np.sin(O)],
                     [0, -np.sin(O), np.cos(O)]])
def R3(O):
    return np.array([[np.cos(O), np.sin(O), 0],
                     [-np.sin(O), np.cos(O), 0],
                     [0, 0, 1]])

#Orbital elements of Jupiter
aj=5.202575
ecc=0.048908
inc=1.3038*np.pi/180
Omega=100.5145*np.pi/180
omega=274.8752*np.pi/180
Mo=80.0392*np.pi/180
mj=1/1047.348625
T=4331.797849

ju=np.array([0, 0, 0]) #Jupiter traj
n=k/((aj**3)**(1/2))
r3Omega=R3(-Omega)
r1inc=R1(-inc)
r3omega=R3(-omega)

#Orbital elements of asteroid
A=5.454
ecca=0.3896
inca=108.358
Omegaa=276.509*np.pi/180
omegaa=226.107*np.pi/180
Moa=146.88*np.pi/180
TT=(4*np.pi**2*A**3/(G*mo))**(1/2)

#Initial position of asteroids
na=k/((A**3)**(1/2))
Ma=Moa+na*(epoch-epoch)
Ea=Moa
while abs(Ea-ecca*np.sin(Ea)-Ma)>10**(-6):
    Ea=Ea-(Ea-ecca*np.sin(Ea)-Ma)/(1-ecca*np.cos(Ea))

ra=A*(1-ecca*np.cos(Ea))
Xa=A*(np.cos(Ea)-ecca)
Ya=A*(1-ecca**2)**(1/2)*np.sin(Ea)
Xpa=-na*A**2/ra*np.sin(Ea)
Ypa=na*A**2/ra*(1-ecca**2)**(1/2)*np.cos(Ea)

ar3Omega=R3(-Omegaa)
ar1inc=R1(-inca)
ar3omega=R3(-omegaa)
u[0,:3]=(ar3Omega@ar1inc@ar3omega@np.array([[Xa],[Ya],[0]])).reshape((3,))
u[0,3:]=(ar3Omega@ar1inc@ar3omega@np.array([[Xpa],[Ypa],[0]])).reshape((3,))



a=np.zeros((t.shape[0]))
e=np.zeros((t.shape[0]))
j=np.zeros((t.shape[0]))
a[0]=A
e[0]=0.3896
j[0]=108.358




def ff(u,t):
    M=Mo+n*(t-epoch)
    E=Mo
    while abs(E-ecc*np.sin(E)-M)>10**(-6):
        E=E-(E-ecc*np.sin(E)-M)/(1-ecc*np.cos(E))
    
    X=aj*(np.cos(E)-ecc)
    Y=aj*(1-ecc**2)**(1/2)*np.sin(E)
    
    jusec=r3Omega@r1inc@r3omega@np.array([[X],[Y],[0]])
    ju[0]=jusec[0]
    ju[1]=jusec[1]
    ju[2]=jusec[2]
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

#%% ellipse
fig,axs=plt.subplots(2,1)
plt.grid(1)
plt.grid(1)
axs[1].plot(u[:,0],u[:,1])
plt.grid(1)


#%% semi major et excentricity

plt.plot(t,a)

#%%
plt.plot(t,e)
#plt.plot(t,j)

#%% Calculus Jupiter traj
def Jup(t):
    M=Mo+n*(t-epoch)
    E=Mo
    while abs(E-ecc*np.sin(E)-M)>10**(-6):
        E=E-(E-ecc*np.sin(E)-M)/(1-ecc*np.cos(E))
    
    X=aj*(np.cos(E)-ecc)
    Y=aj*(1-ecc**2)**(1/2)*np.sin(E)
    
    return (r3Omega@r1inc@r3omega@np.array([[X],[Y],[0]])).reshape((3,))

ju=np.zeros((t.shape[0],3))
for i in range(t.shape[0]):
    ju[i]=Jup(t[i])

#%% 3D
ax = plt.figure().add_subplot(projection='3d')

ax.plot(u[:,0],u[:,1],u[:,2])
ax.plot(ju[:,0],ju[:,1],ju[:,2])
plt.show()

#%% DIST AST-JUP

dist=((u[:,0]-ju[:,0])**2+(u[:,1]-ju[:,1])**2+(u[:,2]-ju[:,2])**2)**(1/2)
plt.plot(t,dist)
plt.grid()





