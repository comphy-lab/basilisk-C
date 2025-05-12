#!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np
import pypg as bas

plt.ion()

N = 4
nl = 4
x = np.linspace(0, 1, N)
y = np.linspace(0, 1, N)
X,Y = np.meshgrid(x,y)
ne = N*N*nl

# arrays are 
b0 = np.zeros((nl,N,N))
b0_lt = np.zeros((nl,N,N))
db0 = np.zeros((nl,N,N))


for l in range(0,nl):
  b0[l,:,:] = (l+1)*np.sin(2*np.pi*X)*np.cos(2*np.pi*X)

b0_lt[:,:,:] = 1e-5*(np.random.rand(nl,N,N) - 0.5)

def graph(i,t):
    print ("t=",t)
    bas.pyset_field(1,b0.reshape(ne))
#    print(bas.fonh.x.f(X,Y))

bas.init_grid(N)
bas.pyinit_vertgrid(nl)
bas.pyset_vars()

#bas.event(graph, i = 2)

b0_sav = 1.*b0

db1 = 0*db0
b1 = b0 + b0_lt
bas.pystep(b0.reshape(ne),db0.reshape(ne))
bas.pystep(b1.reshape(ne),db1.reshape(ne))

db6 = db1 - db0
db2 = db6/b0_lt
db3 = 0*db2

bas.pystep_lt(b0.reshape(ne),db3.reshape(ne),b0_lt.reshape(ne))

#db4 = np.where(db3 != 0, (db1 - db0)/db3, 0.)
db4 = db6/db3
db5 = db6 -db3
print (db5) # should be only zeros
print (db4) # should be only ones

bas.pytrash_vars()

#bas.pyget_field(1,b0.reshape(ne))

# bas.pyset_field(1,b0.reshape(ne))
# bas.pystep(b0.reshape(ne),db0.reshape(ne))
# bas.pystep(b0.reshape(ne),db0.reshape(ne))
# bas.pystep(b0.reshape(ne),db0.reshape(ne))


#bas.trash_vars()
#bas.run()

#plt.contourf(X,Y,db0[-1,:,:].T); plt.colorbar()
