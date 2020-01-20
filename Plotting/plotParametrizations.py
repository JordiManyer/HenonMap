import numpy as np
import matplotlib.pyplot as plt
from math import *
import csv


# orbit class
class orbit:
    def __init__(self , m , x , y):
        self.m = m
        self.x = x
        self.y = y


# we read the orbits
filename = "../cmake-build-debug/Stable.out"
with open(filename, 'r') as csvfile:
    myreader = csv.reader(csvfile, delimiter=' ', quotechar='|', skipinitialspace=True)

    row = next(myreader)
    m = int(row[0]) # Number of points in the orbit

    x = np.zeros(m)
    y = np.zeros(m)

    for j in range(m):
        row = next(myreader)
        x[j] = float(row[0])
        y[j] = float(row[1])
    parS = orbit(m,x,y)

filename = "../cmake-build-debug/Unstable.out"
with open(filename, 'r') as csvfile:
    myreader = csv.reader(csvfile, delimiter=' ', quotechar='|', skipinitialspace=True)

    row = next(myreader)
    m = int(row[0]) # Number of points in the orbit

    x = np.zeros(m)
    y = np.zeros(m)

    for j in range(m):
        row = next(myreader)
        x[j] = float(row[0])
        y[j] = float(row[1])
    parU = orbit(m,x,y)


orbitsS = []
filename = "../cmake-build-debug/Stable_orbits.out"
with open(filename, 'r') as csvfile:
    myreader = csv.reader(csvfile, delimiter=' ', quotechar='|', skipinitialspace=True)

    row = next(myreader)
    ns = int(row[0]) # Number of orbits

    for i in range(0,ns):
        row = next(myreader)
        m = int(row[0]) # Number of points in the orbit

        x = np.zeros(m)
        y = np.zeros(m)

        for j in range(m):
            row = next(myreader)
            x[j] = float(row[0])
            y[j] = float(row[1])
        orbitsS.append(orbit(m,x,y))


orbitsU = []
filename = "../cmake-build-debug/Unstable_orbits.out"
with open(filename, 'r') as csvfile:
    myreader = csv.reader(csvfile, delimiter=' ', quotechar='|', skipinitialspace=True)

    row = next(myreader)
    nu = int(row[0]) # Number of orbits

    for i in range(0,nu):
        row = next(myreader)
        m = int(row[0]) # Number of points in the orbit

        x = np.zeros(m)
        y = np.zeros(m)

        for j in range(m):
            row = next(myreader)
            x[j] = float(row[0])
            y[j] = float(row[1])
        orbitsU.append(orbit(m,x,y))




fig1 = plt.figure()
ax1 = fig1.add_subplot(1,1,1)
ax1.plot(parS.x , parS.y)
ax1.plot(parU.x , parU.y)
ax1.set_xlim(-1,1)
ax1.set_ylim(-1,1)

fig2 = plt.figure()
ax2 = fig2.add_subplot(1,1,1)
ax2.plot(parS.x , parS.y , 'b')
ax2.plot(parU.x , parU.y , 'r')
for i in range(0,ns):
    ax2.plot(orbitsS[i].x , orbitsS[i].y , '.b')
for i in range(0,nu):
    ax2.plot(orbitsU[i].x , orbitsU[i].y , '.r')
ax2.set_xlim(-1,1)
ax2.set_ylim(-1,1)


plt.show()



