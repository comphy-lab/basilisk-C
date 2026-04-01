import matplotlib.pyplot as plt
import numpy as np

# loading data
raw = np.loadtxt("T_profile.dat")

# get number of timesteps and layers
dt=0.
Nl=0
while dt==0:
    Nl=Nl+1
    dt=raw[Nl][0]
Nt = int(raw[-1][0]/dt)+1

data = np.zeros((Nt,Nl))
time = np.zeros(Nt)
layers = np.arange(Nl)
for i in range(Nt):
    time[i] = raw[i*Nl][0]
    for k in range(Nl):
        index = i*Nl + k
        data[i,k]=raw[index][-1]

# plot
fig, ax = plt.subplots(1,1,figsize = (7,3),constrained_layout=True,dpi=100)
s=ax.pcolormesh(time,layers,data.T, cmap='magma')
plt.colorbar(s,ax=ax, label='T (°C)')
ax.set_xlabel('time (s)')
ax.set_ylabel('layer')
fig.savefig('T_profile.png')
plt.show()
