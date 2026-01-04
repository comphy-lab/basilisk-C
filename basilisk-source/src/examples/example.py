"""
# Python interface

Basilisk also includes a high-level interface for the [Python
programming
language](http://en.wikipedia.org/wiki/Python_%28programming_language%29).

This example documents how to use this interface to setup the
[decaying two-dimensional turbulence](turbulence.c) simulation.

First check that Basilisk is [configured
properly](/src/INSTALL#using-basilisk-with-python) for use with
Python.

## Basilisk python module

The first step is to create a python module importing the necessary
functions from Basilisk. For this example we just need to create a
file called *stream.c* with the following lines

~~~literatec
#include "grid/multigrid.h"
#include "navier-stokes/stream.h"
~~~

i.e. we will use the streamfunction--vorticity Navier--Stokes solver
with a multigrid discretisation.

To create the corresponding python module, assuming that a [Makefile has
been setup properly](/Tutorial#using-makefiles), we can just do

~~~bash
make stream.py
~~~

## Python program

The next step is to create the python program itself, let's call it
*example.py*.

We first import [matplotib](http://matplotlib.org/),
[numpy](http://www.numpy.org/) and our Basilisk module.

"""

from __future__ import print_function
import matplotlib.pyplot as plt
import numpy as np
import stream as bas

"""
We setup a *numpy* regular grid which will hold the values
interpolated from the Basilisk field."""

N = 256
x = np.linspace(0, 1, N)
y = np.linspace(0, 1, N)
X,Y = np.meshgrid(x,y)

"""
We define a function which sets initial conditions for the vorticity
$\omega$: a white noise."""

def init(i,t):
    bas.omega.f = bas.noise

"""
We define a function which uses *matplotlib* to generate a graph
of the $\omega$ field."""

def graph(i,t):
    print ("t=",t)
    Z = bas.omega.f(X,Y)
    plt.cla()
    plt.imshow(Z)
    plt.pause(0.0001)
    
"""
We setup Basilisk's resolution to match that of our *numpy* grid."""

bas.init_grid(N)

"""
We register the initial condition (at *t=0*) and graph function
(at *t = 0,10,20,...,1000*)."""

bas.event(init, t = 0.)
plt.ion()
plt.figure()
bas.event(graph, t = range(0,1000,10))

"""
We run the simulation."""

bas.run()
