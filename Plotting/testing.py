import numpy as np
import matplotlib.pyplot as plt
import math


def F(x,y,epsilon):
    xout = x + epsilon*(y + epsilon*x*(1.0-x));
    yout = y + epsilon*x*(1.0-x);
    return xout , yout

def computeOrbit(x0 , y0 , eps , n):
    x = np.zeros(n)
    y = np.zeros(n)
    x[0] = x0
    y[0] = y0
    for i in range(1 , n):
        x[i] , y[i] = F(x[i-1] , y[i-1] , eps)
    return x , y


def DF(x,y,epsilon):
    A = np.zeros((2,2))
    A[0,0] = 1 + epsilon*epsilon*(1.0 - 2.0*x)
    A[1,0] = epsilon*(1.0 - 2.0*x)
    A[0,1] = epsilon
    A[1,1] = 1
    return A

n = 50
x0 = np.linspace(-2 , 2 , n)

eps = 1.9
fig = plt.figure()
ax = fig.add_subplot(1,1,1)
for i in range(0 , n):
    for j in range(0,n):
        orbx , orby = computeOrbit(x0[i] , x0[j] , eps , 100)
        ax.plot(orbx[1:] , orby[1:] , '.r' , markersize=0.5)
ax.set_xlim(-0.5 , 3)
ax.set_ylim(-2 , 2)

n = 500
x0 = np.linspace(-2 , 2 , n)

fig2 = plt.figure()
ax2 = fig2.add_subplot(1,1,1)
for i in range(0 , n):
    orbx , orby = computeOrbit(x0[i] , 0 , eps , n)
    ax2.plot(orbx[1:] , orby[1:] , '.r' , markersize=0.5)
for i in range(0 , n):
    orbx , orby = computeOrbit(0,x0[i], eps , n)
    ax2.plot(orbx[1:] , orby[1:] , '.r' , markersize=0.5)
ax2.set_xlim(-0.5 , 3)
ax2.set_ylim(-2 , 2)


plt.show()












