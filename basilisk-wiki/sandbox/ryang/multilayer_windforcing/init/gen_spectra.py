#
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
from specgen import spectrum_PM, spectrum_gen_linear, eta_random
#
L = 200;    # domain size
P = 0.0125;        # energy level (estimated so that kpHs is reasonable)
kp = 2.*np.pi/40.
#P = 0.3          # energy level (estimated so that kpHs is reasonable)
#N = 4;
#L = 2.0*np.pi    # domain size
#kp = 2*np.pi/(L/N);    # peak wavenumber
N_mode = 64       # number of nodes
N_power = 2       # directional spreading coeff
#
def shape (kmod):
    ''' Choose values here '''
    global P, kp
    F_kmod = spectrum_PM (P=P, kp=kp, kmod=kmod)
    return F_kmod

kmod, F_kmod, kx, ky, F_kxky_tile = spectrum_gen_linear(shape, N_mode=N_mode, L=L, N_power=N_power)
#
''' Generate a grid in x-y to visualize random eta '''
N_grid = 256; #L = 200
x = np.linspace(-L/2,L/2,N_grid); y = np.linspace(-L/2,L/2,N_grid)
x_tile, y_tile = np.meshgrid(x, y)
kx_tile, ky_tile = np.meshgrid(kx,ky)
t = 0
eta_tile, phase_tile = eta_random(t, kx_tile, ky_tile, F_kxky_tile, x_tile, y_tile)
#print('kpHs = %g' %(kp*np.std(eta_tile)*4));
#''' Visualization '''
#plt.imshow(eta_tile, cmap='RdBu_r', vmax=eta_tile.max(), vmin=-eta_tile.max())
#plt.colorbar();
#plt.savefig('init.png',dpi=200)
#plt.close()
#
''' When happy, output the spectra '''
print('Spectrum array shape:', F_kxky_tile.shape)
#path = "/projects/DEIKE/jiarongw/multilayer/" # path of choice
#path = "/home/nn8802/Documents/simulation/basilisk_dotc_ml/spectra_nicolo/"
path = "./";
#
file = "F_kxky_P%g_L%gm_N%g" %(P,L,N_power)
fF = open(path + file, 'bw')
print('Written to', path+file)
F_output = F_kxky_tile.astype('float32'); F_output.tofile(fF)
#
file = "kx_P%g_L%gm_N%g" %(P,L,N_power)
fF = open(path + file, 'bw')
print('Written to', path+file)
F_output = kx.astype('float32'); F_output.tofile(fF)
#
file = "ky_P%g_L%gm_N%g" %(P,L,N_power)
fF = open(path + file, 'bw')
print('Written to', path+file)
F_output = ky.astype('float32'); F_output.tofile(fF)

print(np.std(eta_tile)*4*kp,np.std(eta_tile))
#
''' This is a 0.5 meter box case '''
#kp = 2*np.pi/0.1 
#P = 0.0012
#L = 0.5
#N_mode = 32       
#N_power = 5    

''' This is a 500 meter box case '''
#kp = 2*np.pi/100 
#P = 0.08
#L = 0.5
#N_mode = 32      
#N_power = 5     