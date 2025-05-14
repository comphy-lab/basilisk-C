import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib import cm
from scipy.io import loadmat
mat = loadmat('save_R0.5_h1_alpha1.25_epsilon0.03c.mat')

n_div = 32
alpha = 1.25
Rs1 = 1
Hs1 = 1
Rs2 = 0.5
Hs2 = alpha * Hs1
a   = 0.03
z_max = alpha / np.abs(1 - alpha) * Hs1
t_max = z_max * (2*np.pi)
N = np.max([1., alpha]) / np.abs(1 - alpha) * n_div;

t0 = np.linspace(0, t_max, N)

r1 = mat['hext']['r'][0][0].flatten()
t1 = mat['hext']['t'][0][0].flatten()
z1 = mat['hext']['z'][0][0].flatten()

r2 = mat['hint']['r'][0][0].flatten()
t2 = mat['hint']['t'][0][0].flatten()
z2 = mat['hint']['z'][0][0].flatten()

t0_long = np.hstack([t0[0:-1] - t_max, t0, t0[1:] + t_max])
r1_long = np.hstack([r1[0:-1],         r1, r1[1:]        ])
t1_long = np.hstack([t1[0:-1] - t_max, t1, t1[1:] + t_max])
z1_long = np.hstack([z1[0:-1] - z_max, z1, z1[1:] + z_max])

r2_long = np.hstack([r2[0:-1],         r2, r2[1:]        ])
t2_long = np.hstack([t2[0:-1] - t_max, t2, t2[1:] + t_max])
z2_long = np.hstack([z2[0:-1] - z_max, z2, z2[1:] + z_max])

x1_long = r1_long * np.cos(t1_long)
x2_long = r2_long * np.cos(t2_long)
y1_long = r1_long * np.sin(t1_long)
y2_long = r2_long * np.sin(t2_long)

n_seg1 = np.shape(r1_long)[0]
n_seg2 = np.shape(r2_long)[0]


import csv
with open('hext.asc', 'w', newline='') as csvfile:
    spamwriter = csv.writer(csvfile, delimiter=' ',
                            quotechar='|', quoting=csv.QUOTE_MINIMAL)
    spamwriter.writerow([int(n_seg1)])
    spamwriter.writerow([Rs1])
    spamwriter.writerow([a])
    spamwriter.writerow([Hs1])
    for ti, xi, yi, zi in zip(t0_long, x1_long, y1_long, z1_long):
        spamwriter.writerow([ti, xi, yi, zi])

import csv
with open('hint.asc', 'w', newline='') as csvfile:
    spamwriter = csv.writer(csvfile, delimiter=' ',
                            quotechar='|', quoting=csv.QUOTE_MINIMAL)
    spamwriter.writerow([int(n_seg2)])
    spamwriter.writerow([Rs2])
    spamwriter.writerow([a])
    spamwriter.writerow([Hs2])
    for ti, xi, yi, zi in zip(t0_long, x2_long, y2_long, z2_long):
        spamwriter.writerow([ti, xi, yi, zi])
