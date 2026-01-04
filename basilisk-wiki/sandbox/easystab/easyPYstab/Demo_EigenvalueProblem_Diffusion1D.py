#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Demonstration of how to solve eigenvalue problems with easyPYstab
This program considers the 1D diffusion equation for T(x,t)

d T / d t = mu d^2 T / d x^2  on the interval x in [0, L]
with BC T(x=0,t) = T(x=L,t) = 0

We look for eigenvalue/eigenmode solutions 

T(x,t) = \hat{T}(x)_n e^{lambda_n t} 

(\hat{T}_n, lambda_n) are the solutions of 

lambda_n \hat{T}_n = A \hat{T}_n

where A is the discretized version of the second derivative operator


Created on Thu Dec 11 08:17:37 2025
(starting from a previous program in Matlab/Octave)

@author: fabred
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import numpy as np
import math
import scipy.linalg as la
import matplotlib.pyplot as plt

import EasyStab as ES



#### PART 0 : parameters and theoretical solution ####
plt.close('all')
N = 20         # number of gridpoints
L = np.pi      # domain length
mu = 1.0       # diffusion coefficient

# theoretical eigenvalues
n = np.arange(1, N+1)
s_theory = - mu * (np.pi**2) / (L**2) * n**2



#### PART I : Finite differences ####
print("\nPart I : resolution with Finite Differences \n ")

# discretization
x,dx,dxx,w = ES.dif1D(Dtype='fd',N=N,parameters = {"a":0,"b":L})

# build matrices A and B
I = np.eye(N+1)
A = mu * dxx.copy()
B = I.copy()
# Modify fist/last rows to impose BC
for ib in [0, N]:
    B[ib, :] = 0.0
    A[ib, :] = 0.0
    A[ib, ib] = -1.0

# solve generalized eigenvalue problem
s,U = la.eig(A, B)
# sort by decreasing real part (i.e. least negative eigenvalues first)
idx = np.argsort(-np.real(s))
s = s[idx]
U = U[:, idx]


# Prepare a figure with two subplots 
fig, axes = plt.subplots(2, 1)

# Plot theoretical vs computed eigenvalues
axes[0].plot(np.real(s_theory), np.zeros_like(s_theory), 'ro', label='theory')
axes[0].plot(np.real(s), np.imag(s), 'g+',markeredgewidth=2, label='computed (FD)')
axes[0].set_xlim([-25.001, 1])
axes[0].set_ylim([-1, 5])
axes[0].set_xlabel('Re(λ)')
axes[0].set_ylabel('Im(λ)')
axes[0].set_title('Eigenvalues spectrum')
axes[0].legend()
axes[0].grid(True)

# Plot first few eigenmodes
for ind in range(min(5, U.shape[1])):
    axes[1].plot(x, np.real(U[:, ind]), '+-', label=f'mode {ind} (real part)')
    if np.any(np.imag(U[:, ind]) != 0):
        axes[1].plot(x, np.imag(U[:, ind]), '--', label=f'mode {ind} (imag part)')
axes[1].set_xlabel('x')
axes[1].set_ylabel('T(x)')
axes[1].set_title('First eigenmodes')
axes[1].grid(True)
axes[1].legend()
fig.tight_layout()
fig.show()

#### PART II : using Chebyshev discretization ####
plt.pause(0.001)
print("\nPart II : use Chebychev to compare \n ")
input("Click 'enter' to continue...")

# build Cnhebyshev discretization, then redo everything
x,dx,dxx,w = ES.dif1D(Dtype="Cheb",N=N,parameters = {"a":0,"b":L})

# build matrices A and B
I = np.eye(N+1)
A = mu * dxx.copy()
B = I.copy()
for ib in [0, N]:
    B[ib, :] = 0.0
    A[ib, :] = 0.0
    A[ib, ib] = -1.0
# solve generalized eigenvalue problem, sort and filter
s,U = la.eig(A, B)
idx = np.argsort(-np.real(s))
s = s[idx]
U = U[:, idx]


# Add results on the previous figure
axes[0].plot(np.real(s), np.imag(s), 'bx',markeredgewidth=2, label='computed (Cheb)')
axes[0].legend()
for ind in range(min(5, U.shape[1])):
    axes[1].plot(x, np.real(U[:, ind]), 'x--')
fig.show()



#### PART III : time loop with initial condition ####
plt.pause(0.001)
print("\nPart III : construction of the solution of initial value problem with projection on eigenmodes \n ")
input("Click 'enter' to continue...")

# Time evolution from initial condition (random)
n_modes = min(10, U.shape[1])
a = np.random.randn(n_modes)*.4
a[2] = 1.0

t_list = np.linspace(0, 1.0, 51)
for t in t_list:
    q = np.zeros(N+1)
    for j in range (n_modes): 
      if math.isfinite(s[j]):
        q += np.real(U[:,j] * (a[j] * np.exp(s[j] * t)))
    plt.clf()
    plt.plot(x, q,'b-', label='T(x,t)')
    # optionally plot each component
    for im in range(n_modes):
        plt.plot(x, (U[:, im] * a[im] * np.exp(s[im] * t)), 'r--', alpha=0.3)
    if t==0: 
        qrange = [min(q),max(q)]
        plt.ylim(qrange)
        plt.xlabel('x')
        plt.ylabel('T(x,0)')
        plt.title("Initial condition and Fourier components \n (type enter to start time loop)")
        plt.pause(0.01)
        input("Click 'enter' to start simulation...")
    plt.title(f"t = {t:.3f}")
    plt.ylabel('T(x,t)')
    plt.grid(True)
    plt.ylim(qrange)
    plt.pause(0.1)
plt.show()


