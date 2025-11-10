# -*- coding: utf-8 -*-
"""
Created on Wed Apr 19 20:50:13 2023

@author: eshan_user
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress



k = [5,10,10,3,1,5,25,1]
Uf = 4
Af = 12
Afm = 4
E = 0.05

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

fig = plt.figure()
ax = plt.axes(projection = '3d')
for i in range(len(x_0)):
    x0 = x_0[i]
            
    y0 = y_0[i]
                
    z0 = z_0[i]
    x=x0
    y=y0
    z=z0
    t= 0
    count = 1
    x_rec = np.array([x])
    y_rec = np.array([y])
    z_rec = np.array([z])
    t_rec = np.array([t])
    while (t < maxtime):
         x = x + parkin(k,K,h,E,Uf,Af,Afm,x,y,z)[0]*delta_t
         y = y + parkin(k,K,h,E,Uf,Af,Afm,x,y,z)[1]*delta_t
         z = z + parkin(k,K,h,E,Uf,Af,Afm,x,y,z)[2]*delta_t
         t = t + delta_t
         if np.isnan(x):
             break
         elif np.isnan(y):
             break
         elif np.isnan(z):
             break
         elif np.isnan(t):
             break         
                        
         count += 1 
         x_rec = np.append(x_rec,x)
         y_rec = np.append(y_rec,y)
         z_rec = np.append(z_rec,z)
         t_rec = np.append(t_rec,t)
    ax.plot3D(x_rec,y_rec,z_rec)
    #ax.view_init(-5,20)
ax.set_xlabel('Activated Ubiquitin')
ax.set_ylabel('Activated Parkin')
ax.set_zlabel('Activated Mutant Parkin')
ax.set_title('E = ' + str(E))
plt.show()

plt.figure()
plt.xlabel('Mitochondrial Ubiquitin')
plt.ylabel('Activated Mutant Parkin')
plt.plot(x_rec,z_rec)

plt.figure()
plt.xlabel('Mitochondrial Ubiquitin')
plt.ylabel('Activated Parkin')
plt.plot(x_rec,y_rec)

plt.figure()
plt.xlabel('Activated Parkin')
plt.ylabel('Activated Mutant Parkin')
plt.plot(y_rec,z_rec)
