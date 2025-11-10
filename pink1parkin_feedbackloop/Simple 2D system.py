import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress
from intersect import intersection


k = [0.5 ,1.5 ,0.4 ,1 ,0.6 ,10]
Uf = 5
Af = 7
K = [3,1]
h = [4,8]
E = [0.05,5]

range_x = np. linspace(-5,5,1000)
range_y = np. linspace(-5,5,1000)


def x_null_deg(k,h,K,E,Af,Uf,y):
    k1 ,k2 ,k3 ,k4 ,k5 ,k6 = k
    h1 ,h2 = h
    K1,K2 = K
    return (k1*E*Uf + k5*Uf*y)/((k2*Uf**h1)/(K1 + Uf**h1) + k3*E*Af)


def y_null_deg(k,h,K,E,Af,Uf,x):
    k1 ,k2 ,k3 ,k4 ,k5 ,k6 = k
    h1 ,h2 = h
    K1,K2 = K
    return (k3*E*Af*x)/(k4+(k6*x**h2)/(K2 + x**h2))
 
def parkin_deg(k,K,h,Uf,Af,E,x,y):
    dxdt = k[0]*E*Uf - (k[1]*x*Uf**h[0])/(K[0]+Uf**h[0]) - k[2]*E*Af*x +k[3]*y +k[4]*Uf*y
    dydt = k[2]*E*Af*x - k[3]*y - (k[5]*y*x**h[1])/(K[1]+x**h[1])
    return dxdt , dydt


def phase_plane_plotter_deg(f, range_x ,range_y ,show = False ):
    for m in range(len(E)):
        x_lim = range_x
        y_lim = range_y
    
        X,Y = np.meshgrid(x_lim,y_lim)
        x1, y1 = np.zeros(X.shape) , np.zeros(Y.shape)
        
        ix, iy = X.shape
 
 
        for i in range(ix):
            for j in range(iy):
                x = X[i, j]
                y = Y[i, j]

                system = f(k,K,h,Uf,Af,E[m] ,x,y)
                x1[i, j] = system[0]
                y1[i, j] = system[1]
        
        null_y = y_null_deg(k,h,K,E[m],Af,Uf,x_lim)
        null_x = x_null_deg(k,h,K,E[m],Af,Uf,y_lim)
        null1= intersection(null_x ,y_lim,x_lim, null_y)
        #print ( null1 , null2 )
 
        plt.xlabel ( 'Mitochondrial Ubiquitin')
        plt.ylabel ( 'Mitochondria Parkin')
        plt.streamplot(X,Y,x1 ,y1 , density =2)
        plt.plot (x_lim , null_y , label = 'dy/dt = 0')
        plt.plot (null_x ,y_lim , label = 'dx/dt=0')
        plt.plot (null1[0][0] , null1[1][0] , 'mo' , markersize= 12)
        plt.legend( title = 'E = '+str(E[m]))
        plt.xlim(-5,5)
        plt.ylim(-5,5)
        #plt . title ( ’ Ubiquitin/Parkin levels with E = ’ + str(E[m]))
        if show:
            plt.show()
 
 
#phase_plane_plotter(parkin ,range_x ,range_y ,show=True)
phase_plane_plotter_deg(parkin_deg,range_x ,range_y ,show=True)
 
 
x_0 = [1]
y_0 = [1]


delta_t = 0.0001
maxtime = 10
E_test = 3.6

for i in range(len(x_0)):
    x0 = x_0[i]
    y0 = y_0[i]
 
    x=x0
    y=y0
 
    t= 0
    count = 1
    
    x_rec = np.array ([x])
    y_rec = np.array ([y])
 
    t_rec = np.array ([t])
    while (t < maxtime):
        x = x + parkin_deg(k,K,h,Uf,Af,E_test ,x,y)[0]*delta_t
        y = y + parkin_deg(k,K,h,Uf,Af,E_test ,x,y)[1]*delta_t
        t = t + delta_t
        if np.isnan(x):
            break
        elif np.isnan(y):
            break
        elif np.isnan(t):
            break
 
        count += 1
        x_rec = np.append(x_rec ,x)
        y_rec = np.append(y_rec ,y)
        
        t_rec = np.append(t_rec ,t)

plt.plot (t_rec ,x_rec , label = 'Mitochondrial Ubiquitin')
plt.plot (t_rec ,y_rec , label = 'Mitochondrial Parkin')
plt.xlabel ( 'Time')
plt.ylabel ( 'Protein levels')
#plt.title ( ' Protein ')
#plt.xlim (2.5 ,3)
plt.legend( title = 'E = '+str(E_test))
plt.show()