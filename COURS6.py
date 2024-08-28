# -*- coding: utf-8 -*-
"""
Created on Tue Mar 28 14:04:37 2023
COMPARE JUPITER CIRCULAIRE ET ELLIPTIQUE
2007VW266
@author: lhoes
"""

import numpy as np
import matplotlib.pyplot as plt
import random as rnd

h=1 #day 0.1 c'est ultra good, mais long #euuuuuuuuuuuhhhhhhhhhh
G=(2.95824*10**-4)
k=G**(1/2)
mo=1
epoch=2456600.5
t=np.arange(0+epoch,epoch+365000+h,h)
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

dA=0.0156
decca=0.00170
dinca=0.0261
dOmegaa=0.001144*np.pi/180
domegaa=0.0501*np.pi/180
dMoa=0.604*np.pi/180


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


#Initiate arrays for a e inc
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
   K=np.cross(r,rp)/(nr*nrp)
   j[i+1]=np.arccos(K[2])*(180/np.pi)

Nclones=4
uclones=np.zeros((t.shape[0],6,Nclones))
#rnd.seed(time.time())
rnd.seed(42)
for clone in range(Nclones):
    A=rnd.uniform(A-dA,A+dA)
    ecca=rnd.uniform(ecca-decca,ecca+decca)
    inca=rnd.uniform(inca-dinca,inca+dinca)
    Omegaa=rnd.uniform(Omegaa-dOmegaa,Omegaa+dOmegaa)
    omegaa=rnd.uniform(omegaa-domegaa,omegaa+domegaa)
    Moa=rnd.uniform(Moa-dMoa,Moa+dMoa)
    
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

    #Initiate arrays for a e inc
    a=np.zeros((t.shape[0]))
    e=np.zeros((t.shape[0]))
    j=np.zeros((t.shape[0]))
    a[0]=A
    e[0]=ecca
    j[0]=inca
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
       K=np.cross(r,rp)/(nr*nrp)
       j[i+1]=np.arccos(K[2])*(180/np.pi)
    uclones[:,:3,clone]=u[:,:3]
    uclones[:,3,clone]=a
    uclones[:,4,clone]=e
    uclones[:,5,clone]=i
    


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
ax = plt.figure("3d clones").add_subplot(projection='3d')
for i in range (Nclones):
    ax.plot(uclones[:,0,i],uclones[:,1,i],uclones[:,2,i],label="clone"+str(i))

ax.plot(u[:,0],u[:,1],u[:,2],label="og")
ax.plot(ju[:,0],ju[:,1],ju[:,2],label="Jupiter")
plt.xlabel("x")
plt.ylabel("y")
ax.legend()
plt.show()


#%% 2D clones

ax = plt.figure("2d clones").add_subplot()
for i in range (Nclones):
    ax.plot(uclones[:,1,i],uclones[:,2,i],label="clone"+str(i))

ax.plot(u[:,1],u[:,2],label="og")
plt.ylabel("y")
plt.xlabel("z")
plt.grid()
ax.legend()
plt.show()

#%% Semi maj axis
#t=(t-epoch)/365.25
option=3
ax = plt.figure("Semi-major axis").add_subplot()
for i in range (Nclones):
    ax.plot(t,uclones[:,option,i],label="clone"+str(i))

ax.plot(t,a,label="og")
plt.ylabel("Semi-major axis (au)")
plt.xlabel("t (year)")
plt.grid()
ax.legend()
plt.show()
#%%ecc
option=4
ax = plt.figure("Eccentricity").add_subplot()
for i in range (Nclones):
    ax.plot(t,uclones[:,option,i],label="clone"+str(i))

ax.plot(t,e,label="og")
plt.ylabel("Eccentricity (au)")
plt.xlabel("t (year)")
plt.grid()
ax.legend()
plt.show()







