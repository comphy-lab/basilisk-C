#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Demonstration of how to solve a generalized eigenvalue problem with scipy

This program corresponds to the example treated in lecture 3 of 
the "Introdution to instabilities course" of M2 DET (Toulouse)"

Created on Fri Jan 31 18:58:21 2025

@author: fabred

This program belongs to the easyPYstab project 
"""


import numpy as np
from scipy import linalg as LA
from scipy.sparse import linalg as LAS 

# construction of matrices
a = np.array([[1., 0.,0.], [0., 2.,0.],[1.,1.,1.]])
b = np.array([[1., 0.,0.], [0., 1.,0.],[0.,0.,0.]])

# resolution with direct method
print("\n Resolution with eig (direct method)\n")
eigenvalues,eigenvectors = LA.eig(a,b=b)
print ("Eigenvalues  : \n"); print(eigenvalues)
print ("Eigenvectors : \n"); print(eigenvectors)

# resolution with iterative method
print("\n Resolution with eigs (iterative method)\n")
eigenvalues,eigenvectors = LAS.eigs(a,M=b)
print ("Eigenvalues  : \n"); print(eigenvalues)
print ("Eigenvectors : \n"); print(eigenvectors)