# -*- coding: utf-8 -*-
"""
Created on Thu Apr 20 07:38:38 2023

@author: eshan_user
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress
from scipy.signal import find_peaks


k = [5,10,10,3,1,5,25,1]
Uf = 4
Af = 12
Afm = 4
E = 0.09

K = [150,50,50]

h = [3,3,3]


delta_t = 0.001
maxtime = 100
x_0 = [1]
y_0 = [1]
z_0 = [1]


def parkin(k,K,h,E,Uf,Af,Afm,x,y,z):
    k1,k2,k3,k4,k5,k6,k7,k8 = k
    K1,K2,K3 = K
    h1,h2,h3 = h    
    dxdt = k1*E*Uf - k3*E*Af*x - (k2*x*Uf**h1)/(K1+Uf**h1) + k4*E*Uf*y - k5*E*Afm*x + k6*E*Uf*z+ (k7*y*x**h2)/(K2+x**h2) + (k8*y*x**h3)/(K3+x**h3)
    dydt = k3*E*Af*x - (k7*y*x**h2)/(K2+x**h2)
    dzdt = k5*E*Afm*x - (k8*z*x**h3)/(K3+x**h3)
    return dxdt, dydt, dzdt

E = np.linspace(0.85,1,1000)

E_rec = []

x_rec = []
y_rec = []
z_rec = []
n = 300

for i in range(len(E)):
    x0 =1
    x = x0
    y0 =1 
    y=y0
    z0 = 1
    z=z0
    
    for j in range(n):
        x = parkin(k,K,h,E[i],Uf,Af,Afm,x,y,z)[0]*10**-2
        y = parkin(k,K,h,E[i],Uf,Af,Afm,x,y,z)[1]*10**-2
        z = parkin(k,K,h,E[i],Uf,Af,Afm,x,y,z)[2]*10**-2
    
    for j in range(n):
        x = parkin(k,K,h,E[i],Uf,Af,Afm,x,y,z)[0]*10**-2
        y = parkin(k,K,h,E[i],Uf,Af,Afm,x,y,z)[1]*10**-2
        z = parkin(k,K,h,E[i],Uf,Af,Afm,x,y,z)[2]*10**-2
        
        E_rec.append(E[i])
        x_rec.append(x)
        y_rec.append(y)
        z_rec.append(z)
        
plt.scatter(E_rec,x_rec,s = 0.1)
plt.xlabel('Pink1 Concentration (E)')
plt.ylabel('Mitochondrial Ubiquitin level')
plt.ylim(-10,10)
#plt.xlim(0.05,2.5)
#plt.yscale('log')
plt.title('Orbit Diagram of Mitochondrial Ubiquitin')
plt.show()
plt.scatter(E_rec,y_rec,s=0.1)
plt.xlabel('Pink1 Concentration (E)')
plt.ylabel('Active Parkin level')
#plt.xlim(0.05,2.5)
plt.ylim(-10,10)
plt.title('Orbit Diagram of Active Parkin')
plt.show()
plt.scatter(E_rec,z_rec,s=0.1)
plt.xlabel('Pink1 Concentration (E)')
plt.ylabel('Active Mutant Parkin level')
plt.title('Orbit Diagram of Active Mutant Parkin')
#plt.xlim(0.05,2.5)
plt.ylim(-10,10)
plt.show()

plt.figure()
ax = plt.axes(projection = '3d')

ax.plot3D(x_rec,y_rec,z_rec)
ax.view_init(20,120)
plt.show()
plt.figure()
ax = plt.axes(projection='3d')
ax.set_title('Ub - Am plane')
ax.plot3D(x_rec,y_rec,E_rec)
ax.view_init(-120,-40)
plt.show()

plt.figure()
ax = plt.axes(projection='3d')
ax.set_title('Ub - Amut plane')
ax.plot3D(x_rec,E_rec,z_rec)
ax.view_init(-120,120)
plt.show()
plt.figure()
ax = plt.axes(projection='3d')
ax.set_title('Am - Amut plane')
ax.plot3D(E_rec,y_rec,z_rec)
ax.view_init(-120,120)
plt.show()
        

        
        
