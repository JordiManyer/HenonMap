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
def readParams(i1 , i2):
    filename = outdir+"Stable_" + str(i1) + "_" + str(i2) + ".out"
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

    filename = outdir+"Unstable_" + str(i1) + "_" + str(i2) + ".out"
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

    return parS , parU


def readOrbits(i1 , i2):
    orbitsS = []
    filename = outdir+"StableOrbits_" + str(i1) + "_" + str(i2) + ".out"
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
    filename = outdir+"UnstableOrbits_" + str(i1) + "_" + str(i2) + ".out"
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
    return orbitsS , orbitsU


def readEps():
    filename = outdir+"epsilonVec.out"
    with open(filename, 'r') as csvfile:
        myreader = csv.reader(csvfile, delimiter=' ', quotechar='|', skipinitialspace=True)

        row = next(myreader)
        n = int(row[0])

        epsvec = []
        for i in range(0,n):
            row = next(myreader)
            epsvec.append(float(row[0]))
        return np.array(epsvec)



##################################################################3
outdir = "../Outputs6/"
epsVec = readEps()

neps = 1
epsSelect = np.array([1.5])

nmax = 100
fig1 = plt.figure()
fig2 = plt.figure()
axes1 = []
axes2 = []
for i in range(0 , neps):
    eps = epsSelect[i]
    index = np.where( np.abs(eps - epsVec) == np.min(np.abs(eps - epsVec)))[0][0]
    print(epsVec[index])

    parS0 , parU0 = readParams(index , 0)
    orbitsS0 , orbitsU0 = readOrbits(index , 0)
    ns0 = len(orbitsS0)
    nu0 = len(orbitsU0)

    if np.abs(eps) > 2.0:
        parS1 , parU1 = readParams(index , 1)
        orbitsS1 , orbitsU1 = readOrbits(index , 1)
        ns1 = len(orbitsS1)
        nu1 = len(orbitsU1)

    axes1.append(fig1.add_subplot(neps,1,i+1))
    axes1[i].plot(parS0.x , parS0.y , 'b')
    axes1[i].plot(parU0.x , parU0.y , 'r')
    if np.abs(eps) > 2.0:
        axes1[i].plot(parS1.x , parS1.y , 'g')
        axes1[i].plot(parU1.x , parU1.y , 'k')
    axes1[i].set_xlim(-2,3)
    axes1[i].set_ylim(-2,2)

    axes2.append(fig2.add_subplot(neps,1,i+1))
    axes2[i].plot(parS0.x , parS0.y , 'b')
    axes2[i].plot(parU0.x , parU0.y , 'r')
    if np.abs(eps) > 2.0:
        axes2[i].plot(parS1.x , parS1.y , 'g')
        axes2[i].plot(parU1.x , parU1.y , 'k')
    for j in range(0,ns0):
        axes2[i].plot(orbitsS0[j].x[:nmax] , orbitsS0[j].y[:nmax] , '.b' , markersize=0.7)
    if np.abs(eps) > 2.0:
        for j in range(0,ns1):
            axes2[i].plot(orbitsS1[j].x[:nmax] , orbitsS1[j].y[:nmax] , '.g' ,markersize=0.7)
    for j in range(0,nu0):
        axes2[i].plot(orbitsU0[j].x[:nmax] , orbitsU0[j].y[:nmax] , '.r' , markersize=0.7)
    if np.abs(eps) > 2.0:
        for j in range(0,nu1):
            axes2[i].plot(orbitsU1[j].x[:nmax] , orbitsU1[j].y[:nmax] , '.k' , markersize=0.7)
    axes2[i].set_xlim(-2,3)
    axes2[i].set_ylim(-2,2)


plt.show()



