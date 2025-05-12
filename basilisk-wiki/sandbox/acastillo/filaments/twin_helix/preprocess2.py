#!/usr/bin/env python
# coding: utf-8

# In[193]:


import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib import cm
from scipy.io import loadmat

mat = loadmat('Parameters_02_24mm_14deg_FOV_B_option3_pos2.mat')

Gamma_star = 0.95
n_div = 48
alpha = 4./3.

Rs1 = 7.5
Hs1 = 5.25
Rs2 = 6.5
Hs2 = alpha * Hs1
a   = 0.1

R = Rs2 + 1.0/(1 + Gamma_star);
tau2 = (1 + (2*np.pi*R/Hs1)**2)**(-1)
bmax = Hs2/tau2


z_max = alpha / np.abs(1 - alpha) * Hs1
t_max = z_max * (2*np.pi)
N = np.max([1., alpha]) / np.abs(1 - alpha) * n_div;

t0 = mat['new_b'][:]/bmax

r1 = mat['r1'][:]
t1 = mat['t1'][:]
z1 = mat['z1'][:]

r2 = mat['r2'][:]
t2 = mat['t2'][:]
z2 = mat['z2'][:]

x1 = r1 * np.cos(t1)
y1 = r1 * np.cos(t1)

x2 = r2 * np.cos(t2)
y2 = r2 * np.cos(t2)


dr1 = mat['fd_r1'][:]
dt1 = mat['fd_t1'][:]
dz1 = mat['fd_z1'][:]
dr2 = mat['fd_r2'][:]
dt2 = mat['fd_t2'][:]
dz2 = mat['fd_z2'][:]

d2r1 = mat['fd2_r1'][:]
d2t1 = mat['fd2_t1'][:]
d2z1 = mat['fd2_z1'][:]
d2r2 = mat['fd2_r2'][:]
d2t2 = mat['fd2_t2'][:]
d2z2 = mat['fd2_z2'][:]


# In[194]:


beta = 3/4
ntimes = 18
shift_t = 2*np.pi/beta
shift_z = Hs2

import numpy.matlib
r1_long = np.matlib.repmat(r1, ntimes, 1)
t1_long = np.matlib.repmat(t1, ntimes, 1) + np.kron(np.arange(ntimes).T*shift_t, np.ones_like(t1).T).T
z1_long = np.matlib.repmat(z1, ntimes, 1) + np.kron(np.arange(ntimes).T*shift_z, np.ones_like(z1).T).T - 42

dr1_long = np.matlib.repmat(dr1, ntimes, 1)
dt1_long = np.matlib.repmat(dt1, ntimes, 1)
dz1_long = np.matlib.repmat(dz1, ntimes, 1)
d2r1_long = np.matlib.repmat(d2r1, ntimes, 1)
d2t1_long = np.matlib.repmat(d2t1, ntimes, 1)
d2z1_long = np.matlib.repmat(d2z1, ntimes, 1)

r2_long = np.matlib.repmat(r2, ntimes, 1)
t2_long = np.matlib.repmat(t2, ntimes, 1) + np.kron(np.arange(ntimes).T*shift_t, np.ones_like(t2).T).T
z2_long = np.matlib.repmat(z2, ntimes, 1) + np.kron(np.arange(ntimes).T*shift_z, np.ones_like(z2).T).T - 42

dr2_long = np.matlib.repmat(dr2, ntimes, 1)
dt2_long = np.matlib.repmat(dt2, ntimes, 1)
dz2_long = np.matlib.repmat(dz2, ntimes, 1)
d2r2_long = np.matlib.repmat(d2r2, ntimes, 1)
d2t2_long = np.matlib.repmat(d2t2, ntimes, 1)
d2z2_long = np.matlib.repmat(d2z2, ntimes, 1)

t0_long = np.matlib.repmat(t0, ntimes, 1) + np.kron(np.arange(ntimes).T*1, np.ones_like(t0).T).T - 6
t0_long = t0_long*2*np.pi

x1_long = r1_long*np.cos(t1_long)
x2_long = r2_long*np.cos(t2_long)
y1_long = r1_long*np.sin(t1_long)
y2_long = r2_long*np.sin(t2_long)

dx1_long = dr1_long * np.cos(t1_long) - dt1_long * r1_long * np.sin(t1_long)
dx2_long = dr2_long * np.cos(t2_long) - dt2_long * r2_long * np.sin(t2_long)

dy1_long = dr1_long * np.sin(t1_long) + dt1_long * r1_long * np.cos(t1_long)
dy2_long = dr2_long * np.sin(t2_long) + dt2_long * r2_long * np.cos(t2_long)

d2x1_long = np.cos(t1_long) * (d2r1_long - r1_long * dt1_long * dt1_long) - np.sin(t1_long) * (2 * dr1_long * dt1_long + r1_long * d2t1_long)
d2x2_long = np.cos(t2_long) * (d2r2_long - r2_long * dt2_long * dt2_long) - np.sin(t2_long) * (2 * dr2_long * dt2_long + r2_long * d2t2_long)

d2y1_long = np.sin(t1_long) * (d2r1_long - r1_long * dt1_long * dt1_long) + np.cos(t1_long) * (2 * dr1_long * dt1_long + r1_long * d2t1_long)
d2y2_long = np.sin(t2_long) * (d2r2_long - r2_long * dt2_long * dt2_long) + np.cos(t2_long) * (2 * dr2_long * dt2_long + r2_long * d2t2_long)

n_seg1 = np.shape(r1_long)[0]
n_seg2 = np.shape(r2_long)[0]

plt.plot(dx1_long)
plt.plot(d2x1_long)


# In[195]:


import csv
with open('hext.asc', 'w', newline='') as csvfile:
    spamwriter = csv.writer(csvfile, delimiter=' ',
                            quotechar='|', quoting=csv.QUOTE_MINIMAL)
    spamwriter.writerow([int(n_seg1)])
    spamwriter.writerow([Rs1])
    spamwriter.writerow([a])
    spamwriter.writerow([Hs1])
    for ti, xi, yi, zi, dxi, dyi, dzi, d2xi, d2yi, d2zi in zip(t0_long.flatten(), x1_long.flatten(), y1_long.flatten(), z1_long.flatten(), dx1_long.flatten(), dy1_long.flatten(), dz1_long.flatten(), d2x1_long.flatten(), d2y1_long.flatten(), d2z1_long.flatten()):
        spamwriter.writerow([ti, xi, yi, zi, dxi, dyi, dzi, d2xi, d2yi, d2zi])

import csv
with open('hint.asc', 'w', newline='') as csvfile:
    spamwriter = csv.writer(csvfile, delimiter=' ',
                            quotechar='|', quoting=csv.QUOTE_MINIMAL)
    spamwriter.writerow([int(n_seg2)])
    spamwriter.writerow([Rs2])
    spamwriter.writerow([a])
    spamwriter.writerow([Hs2])
    for ti, xi, yi, zi, dxi, dyi, dzi, d2xi, d2yi, d2zi in zip(t0_long.flatten(), x2_long.flatten(), y2_long.flatten(), z2_long.flatten(), dx2_long.flatten(), dy2_long.flatten(), dz2_long.flatten(), d2x2_long.flatten(), d2y2_long.flatten(), d2z2_long.flatten()):
        spamwriter.writerow([ti, xi, yi, zi, dxi, dyi, dzi, d2xi, d2yi, d2zi])


# In[196]:


plt.plot(t0_long, r1_long*np.cos(t1_long))
plt.plot(t0_long, r2_long*np.cos(t2_long), '--')
plt.show()
plt.plot(t0_long, r1_long*np.sin(t1_long))
plt.plot(t0_long, r2_long*np.sin(t2_long), '--')
plt.show()
plt.plot(t0_long, z1_long)
plt.plot(t0_long, z2_long, '--')
plt.show()
plt.plot(t0_long, r1_long)
plt.plot(t0_long, r2_long, '--')
plt.show()
plt.plot(t0_long, t1_long/(2*np.pi))
plt.plot(t0_long, t2_long/(2*np.pi), '--')

