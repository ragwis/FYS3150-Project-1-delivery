# -*- coding: utf-8 -*-
"""
Created on Mon Sep  7 18:21:38 2020

@author: rawis
"""
import numpy as np
import matplotlib.pyplot as plt

"""
make the toeplitz matrix with lower diagonal a, diagonal b and upper diagonal c, 
assumes that |a|+1 = |b| = |c|+1
def toeplitz_matrix(a,b,c):
    n = len(b)
    A = np.zeros([n,n])
    indexes = np.arange(0,n) #same as range, but an array(can use +-*/ etc) #array from 0 to n-1
    
    A[indexes, indexes] = b #diagonal
    A[indexes[:-1]+1, indexes[:-1]]=a #subdiagonal
    A[indexes[:-1], indexes[:-1]+1]=c #superdiagonal
    
    return A
"""

"""
solve Au=g for u, where A is a toeplitz matrix with off-diagonals e and diagonal d.
e and d are vectors with not neccesarilly equal entries.
"""
def solve_toeplitz(e,d,g):
    n=len(d)
    d2=np.zeros([n])
    g2=np.zeros([n])
    u =np.zeros([n])
    
    #boundary conditions(?)
    u[0]  = 0
    u[-1] = 0
    
    #forward subs
    d2[0] = d[0]
    g2[0] = g[0]
    for i in range(1,n-1):
        d2[i] = d[i] - e[i]**2/d2[i-1]
        g2[i] = g[i] - e[i-1]*g2[i-1]/d2[i-1]
    print(d2)
    print(g2)
    #backward subs
    for i in range(n-2,0,-1):
        u[i] = (g2[i]-e[i]*u[i+1])/d2[i]
    
    return u


#run by selecting all and F9:
n=1000
x = np.linspace(0,1,n)
e = np.zeros(n-1)+1
d = np.zeros(n)-2
g = -(1/n)**2 * 100*np.exp(-10*x)
u_numerical = solve_toeplitz(e,d,g)
u_exact     = 1-(1-np.exp(-10))*x-np.exp(-10*x)
plt.plot(x, u_numerical, label='$u_{computed}$')
plt.plot(x, u_exact, label='$u_{exact}$')
plt.title(f'n = {n}')
plt.xlabel('x')
plt.ylabel('u(x)')
plt.legend()


"""
Then you should code the above algorithm and solve the problem for matrices
of the size 10 × 10, 100 × 100 and 1000 × 1000. That means that you select
n = 10, n = 100 and n = 1000 grid points.

Compare your results (make plots) with the closed-form solution for the
different number of grid points in the interval x ∈ (0, 1). The different number
of grid points corresponds to different step lengths h.
"""

def closed_form_solution_u(x):
    u = 1-(1-np.exp(-10))*x -np.exp(-10*x)
    return u

def f(x):
    return 100*np.exp(-10*x)

n_1 = 10
n_2 = 100
n_3 = 1000
x_0 = 0
x_end = 1
x_1 = np.linspace(x_0, x_end, n_1)
x_2 = np.linspace(x_0, x_end, n_2)
x_3 = np.linspace(x_0, x_end, n_3)
d_1 = -2*np.ones(n_1)
d_2 = -2*np.ones(n_2)
d_3 = -2*np.ones(n_3)



"""
c)
Use thereafter the fact that the matrix has identical matrix
elements along the diagonal and identical (but different) values for the nondiagonal elements. Specialize your algorithm to the special case and find the
number of floating point operations for this specific tri-diagonal matrix. Compare
the CPU time with the general algorithm from the previous point for matrices
up to n = 10**6 grid points.
"""
























