import matplotlib.pyplot as plt
import numpy as np
import traceback

# loading data
rawT = np.loadtxt("T_profile.dat")
rawU = np.loadtxt("u_profile.dat")


def get_Nl_Nt(rawfile):
    dt=0.
    Nl=0
    while dt==0:
        Nl=Nl+1
        dt=rawfile[Nl][0]
    Nt = int(rawfile[-1][0]/dt)+1
    return Nl, Nt

def read_raw(raw, Ny, Nx):
    data = np.zeros((Ny,Nx))
    X = np.zeros(Nx)
    Y = np.zeros(Ny)
    for i in range(Nx):
        X[i] = raw[i*Ny][0]
        for k in range(Ny):
            index = i*Ny + k
            Y[k] = raw[k][1]
            data[k,i]=raw[index][-1]
    return X, Y, data

# > T_profile
Nl, Nt = get_Nl_Nt(rawT)
time, layers, dataT = read_raw(rawT, Nl, Nt)
T0 = dataT[-1,0]
# plot
fig, ax = plt.subplots(1,1,figsize = (7,3),constrained_layout=True,dpi=100)
# s = ax.pcolormesh(time, layers, dataT/T0, shading='nearest', cmap='magma')
# plt.colorbar(s,ax=ax, label='$T/Ts_{0}$')
s = ax.pcolormesh(time, layers, dataT, shading='nearest', cmap='magma')
plt.colorbar(s,ax=ax, label='$T (°C)$')
ax.set_xlabel('time (s)')
ax.set_ylabel('layer')
fig.savefig('T_profile.png')

# > wT_profile
Nl, Nt = get_Nl_Nt(rawU)
time, layers, dataU = read_raw(rawU, Nl, Nt)
#flx0 = 500/(1025*4.2e3)
# plot
fig, ax = plt.subplots(1,1,figsize = (7,3),constrained_layout=True,dpi=100)
s=ax.pcolormesh(time,layers,dataU, shading="nearest", cmap='magma')
plt.colorbar(s,ax=ax, label=r"U ($m.s^{-1}$)")
ax.set_xlabel('time (s)')
ax.set_ylabel('layer')
fig.savefig('u_profile.png')

plt.show()

