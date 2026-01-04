#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

This is the library containing the basic bricks of the easyPYstab project:
    - function "dif1D" to construct mesh, derivative matrices and integration weights

Created on Fri Jan 31 18:58:21 2025

@author: fabred

This program belongs to the easyPYstab project protected under CC licence
"""
import matplotlib.pyplot as plt           # Pour tracer
import numpy as np                        # Pour calculer
from numpy import*                        # Pour construire les tableaux
from matplotlib.pyplot import*            # Pour construire le graphe
"""
# Function Dif1D (central tool of the easystab tookbox)
"""

def dif1D(Dtype='fd',N=50,parameters = {"a":0,"b":1}):
# This function builds the basic bricks : mesh and derivative matrices
# at order one and two. Currently implemented are :
#     "fd"    -> Finite differences on interval [a b]
#     "Cheb"  -> Chebychev method on interval [a b]
 
  d  = np.zeros((N+1,N+1));
  dd = np.zeros((N+1,N+1));
  x = np.zeros(N+1); w = np.zeros(N+1);
  
  """
  ## Finite differences
  """
  if Dtype=="fd":  
        dx = (parameters["b"]-parameters["a"])/N; 
        x = np.linspace(parameters["a"],parameters["b"],N+1)
        for i in range(1,N): 
            # centered formulas leading to a tridiagonal matrix
            d[i,i] = 0;
            d[i,i+1] =  1/(2*dx); 
            d[i,i-1] = -1/(2*dx);
            dd[i,i] = -2/dx**2;
            dd[i,i+1] =  1/(dx**2); 
            dd[i,i-1] =  1/(dx**2);        
        # uncentred formulas for first and last grid points 
        d[0,0] = -3/(2*dx); d[0,1] = 4/(2*dx) ; d[0,2] = -1/(2*dx);
        d[N,N-2] =  1/(2*dx); d[N,N-1] = -4/(2*dx); d[N,N]   =  3/(2*dx);
        dd[0,0] = 2/dx**2; dd[0,1] = -5/dx**2; dd[0,2] = 4/dx**2; dd[0,3] = -1/dx**2; 
        dd[N,N] = 2/dx**2; dd[N,N-1] = -5/dx**2; dd[N,N-2] = 4/dx**2; dd[N,N-3] = -1/dx**2; 
        # "weight" to compute integrals using the trapeze rule
        for i in range(1,N): 
            w[i] = dx;
        w[1] = dx/2; w[N] = dx/2;
  """
  ## Chebyshev
  """      
  elif Dtype=="Cheb":
        # maillage sur l'interalle [-1,1] ordonné de 1 à -1
        s = np.cos(np.pi * np.arange(N + 1) / N)
        # Calcul des coefficients pour la matrice de dérivée première
        for i in range(N + 1):
          for j in range(N + 1):
            if i != j:
                c_i = 2 if i == 0 or i == N else 1
                c_j = 2 if j == 0 or j == N else 1
                d[i, j] = (c_i / c_j) * ((-1) ** (i + j)) / (s[i] - s[j])
        # Remplissage de la diagonale de la matrice de dérivée première
        for i in range(N + 1):
            d[i, i] = -np.sum(d[i, :])

        # Calcul de la matrice de dérivée seconde
        dd = np.dot(d, d)
        x = s
        
        # rescaling si intervalle different de [1,-1]
        if parameters["a"]!=1 and parameters["a"]!=-1:
            scale = (parameters["b"]-parameters["a"])/2
            x = (parameters["b"]+parameters["a"])/2+scale*x
            d = 1/scale*d
            dd = 1/scale**2*dd
  else:
     x = N      
  return x,d,dd,w       


"""
# A few auxiliary functions
"""


"""
## Fonction myinput (input with default value)
""" 


def myinput(question,default):
    answer = input(f"{question} [default : {default}] : " )
    if not answer:
        return default 
    else:
        return answer

        

