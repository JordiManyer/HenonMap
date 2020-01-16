import numpy as np
import matplotlib.pyplot as plt
import math




def DF(x,y,epsilon):
    A = np.zeros((2,2))
    A[0,0] = 1 + epsilon*epsilon*(1.0 - 2.0*x)
    A[1,0] = epsilon*(1.0 - 2.0*x)
    A[0,1] = epsilon
    A[1,1] = 1
    return A



eps = 0.5
D = DF(0,0,eps)

w , v = np.linalg.eig(D)

print("Eigenvalues: ")
print(w)
print("Eigenvectors: ")
print(v)











