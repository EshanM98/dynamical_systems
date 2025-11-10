# -*- coding: utf-8 -*-
"""
Created on Thu Mar  9 22:49:51 2023

@author: eshan_user
"""

import numpy as np
import matplotlib.pyplot as plt
from intersect import intersection
v=1

alpha1 = 1
beta1 = 200
h1 = 4
k1 = 30
alpha2 = 1
beta2 = 10
h2 = 4
k2 = 1
omega = np.linspace(0,5,100)
x = np.linspace(0,5,5)
y_quiver= np.linspace(0,5,5)
y = np.linspace(0,5,100)
c= 0.8
b= 500/140
k = 405/14


    
def prot_act(x):
    x1, y1 = x
    dx = alpha1*(1-x1) - (beta1*x1*(v*y1)**h1)/(k1+(v*y1)**h1)
    dy = alpha2*(1-y1) - (beta2*y1*x1**h2)/(k2 + x1**h2)
    return np.array([dx, dy])
 


def phase_plane_plotter(f, range_x = (0,1), range_y = None,
                     num_grid_points = 25, show = False):
    
    if range_y is None:
        range_y = range_x
    x = np.linspace(range_x[0], range_x[1], num_grid_points)                                                             
    y = np.linspace(range_y[0], range_y[1], num_grid_points)                                                             

    grid = np.meshgrid(x, y)
    
    df_fill = np.zeros((num_grid_points, num_grid_points, 2))
    for m in range(num_grid_points):
        for n in range(num_grid_points):
            df = f([grid[0][m,n], grid[1][m,n]])
            df_fill[m, n, 0] = df[0]
            df_fill[m, n, 1] = df[1]
            
    nullcline_x = alpha1/(alpha1 + (beta1*(v*y)**h1)/(k1 + (v*y)**h1))

    nullcline_y =  alpha2/(alpha2 + (beta2*x**h2)/(k2+x**h2))
    
    plt.plot(x,nullcline_y, 'r--', label = 'dy/dt = 0')
    plt.plot(nullcline_x,y, 'g--', label = 'dx/dt = 0')
    #null_1,null_2 = intersection(nullcline_x,y,x,nullcline_y)
    
    #plt.plot(null_1,null_2, 'mo', markersize = 12)
    plt.streamplot(grid[0],grid[1],df_fill[:,:,0],df_fill[:,:,1],linewidth = 1)
    plt.legend()
    plt.xlabel('Active Cdc2-cyclin B')
    plt.ylabel('Active Wee1')
    plt.title('Phase plane plot of Cdc2-cyclin B/Wee1 system with v = ' +str(v))
    if show:
        plt.show()

if __name__ == "__main__":
    phase_plane_plotter(prot_act, range_x = (0, 1), show = True)
    
def monotonicity(x1,y1):
    return [x1*(-x1+y1), 3*y1*((-x1+c+(b*y1**4)/(k+y1**4)))]

x_lim = np.linspace(0,10,100)
y_lim= np.linspace(0,20,100)
X,Y = np.meshgrid(x_lim,y_lim)
u,v = np.zeros(X.shape), np.zeros(Y.shape)
    
ix, iy = X.shape
    
for i in range(ix):
    for j in range(iy):
        x_phase = X[i,j]
        y_phase = Y[i,j]
        system = monotonicity(x_phase,y_phase)
        u[i,j] = system[0]
        v[i,j] = system[1]
plt.quiver(X,Y, u,v, color = '0.8')
plt.streamplot(X,Y,u,v,linewidth = 1,density = 2)
plt.xlabel('x')
plt.ylabel('y')
plt.ylim(0,20)
plt.xlim(0,10)
plt.title('Phase Plane for Steady State System')
plt.plot()
plt.show()

def steady_state(omega):
    return c + (b*omega**4)/(k+omega**4)

line1 = intersection(omega,steady_state(omega),omega,y)
plt.plot(line1[0][0],line1[1][0], 'ro')
plt.plot(line1[0][1],line1[1][1], 'ro', mfc = 'none')
plt.plot(line1[0][2],line1[1][2],'ro')
plt.quiver(x,y_quiver, np.sign(steady_state(x)-y_quiver), np.sign(steady_state(x)-y_quiver),angles = 'xy',scale_units ='xy')
plt.plot(omega,steady_state(omega),label = 'Steady State')
plt.plot(omega,y)
plt.xlabel('Stimulus(\u03C9)')
plt.ylabel('Response (\u03B7)')
plt.title('Static Characteristic Curve of Input/Output Steady State')
plt.legend()
plt.show()