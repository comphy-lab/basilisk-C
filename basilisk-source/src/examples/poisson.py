"""
# Poisson problem with Python

In this example we access lower-level Basilisk functions from Python.

As in the [2D turbulence example](example.py), we import
matplotlib, numpy and the *stream* Basilisk module. We also use math
functions."""

import matplotlib.pyplot as plt
import numpy as np
import stream as bas
from math import *

"""
We initialise the $256^2$ regular grid."""

N = 256
bas.init_grid(N)

"""
We allocate two new scalar fields *a* and *b*."""

a = bas.scalar()
b = bas.scalar()

"""
And initialize them with simple arithmetic functions."""

a.f = lambda x,y: 0.
b.f = lambda x,y: sin(2.*pi*x)*cos(2.*pi*y)

"""
We then call the [multigrid Poisson solver](/src/poisson.h) to solve
$$
\nabla^2 a = b
$$
"""

bas.poisson(a,b)

"""
And use matplotlib to display the solution."""

x = np.linspace(0, 1, N)
y = np.linspace(0, 1, N)
X,Y = np.meshgrid(x,y)
plt.imshow(a.f(X,Y))
plt.show()

"""
Note that we could reuse the same *stream.py* module because the
streamfunction--vorticity solver (which we don't use in this example)
also exports the interface for the Poisson solver (which we do use in
this example)."""
