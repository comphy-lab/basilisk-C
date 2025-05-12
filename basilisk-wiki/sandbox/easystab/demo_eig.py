# Demonstration of how to solve a generalized eigenvalue problem with scipy
# This program corresponds to the example treated in lecture 3.

import numpy as np
from scipy import linalg as LA
from scipy.sparse import linalg as LAS 

# construction of matrices
a = np.array([[1., 0.,0.], [0., 2.,0.],[1.,1.,1.]])
b = np.array([[1., 0.,0.], [0., 1.,0.],[0.,0.,0.]])

# resolution with direct method
eigenvalues,eigenvectors = LA.eig(a,b=b)
print(eigenvalues)
print(eigenvectors)

# resolution with iterative method
eigenvalues,eigenvectors = LAS.eigs(a,M=b)
print(eigenvalues)
print(eigenvectors)