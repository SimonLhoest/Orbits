# -*- coding: utf-8 -*-
"""
Created on Tue Mar 28 14:04:37 2023
COMPARE JUPITER CIRCULAIRE ET ELLIPTIQUE
@author: lhoes
"""

import numpy as np
import matplotlib.pyplot as plt

h=1 #day 0.1 c'est ultra good, mais long #euuuuuuuuuuuhhhhhhhhhh
G=(2.95824*10**-4)
k=G**(1/2)
A=3.27
mo=1

#Orbital elements of Jupiter
aj=5.202575
ecc=0.048908
inc=1.3038*np.pi/180
Omega=100.5145*np.pi/180
omega=274.8752*np.pi/180
Mo=80.0392*np.pi/180
epoch=2456600.5
mj=1/1047.348625
T=4331.797849
ju=np.array([0, 0, 0]) #Jupiter traj
n=k/((aj**3)**(1/2))

''' MUST START AT T0, CHANGE IT'''
t=np.arange(0+epoch,epoch+36500+h,h)
u=np.zeros((t.shape[0],6))
u[0]=np.array([A, 0, 0, 0, (G/A)**(1/2),0])

a=np.zeros((t.shape[0]))
e=np.zeros((t.shape[0]))
j=np.zeros((t.shape[0]))
a[0]=A
e[0]=0
j[0]=0

def R1(O):
    return np.array([[1, 0, 0],
                     [0, np.cos(O), np.sin(O)],
                     [0, -np.sin(O), np.cos(O)]])
def R3(O):
    return np.array([[np.cos(O), np.sin(O), 0],
                     [-np.sin(O), np.cos(O), 0],
                     [0, 0, 1]])
r3Omega=R3(-Omega)
r1inc=R1(-inc)
r3omega=R3(-omega)



def ff(u,t):
    ''' RECALCULER N ??????????????????????'''
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
   #k=np.cross(r,rp)/(nr*nrp)
   #j[i+1]=np.arccos(k[2])*(180/np.pi)
   
   #error[i]=(np.abs(exact[i]-y[i]))/y[i]
u2=u.copy()
a2=a.copy()
e2=e.copy()

def ff(u,t):
    ju=np.array([5.2*np.cos(2*np.pi*t/T), 5.2*np.sin(2*np.pi*t/T), 0])
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

plt.plot(t,a2)
plt.plot(t,a)
'''
AHHHHHH MAIS JE VEUX CIRCULAR VS ELLIPTIQUE PAS DIFFERENTE AU VS AU
'''
#%%
plt.plot(t,e2)
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







