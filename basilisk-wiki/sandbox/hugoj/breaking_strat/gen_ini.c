
#define g_ 9.81       // before including spectrum.h
#include "hugoj/lib/spectrum.h" // Initial conditions generation
#include "bderembl/libs/extra.h"      // parameters from namlist


// Default parameters
char namlist[80] = "namelist.toml";    // file name of namlist
double P = 0.2;         // energy level 
int coeff_kpL0 = 10 []; // kpL0 = coeff_kpL0 * pi
int N_mode = 32 [];     // number of modes in wavenumber space
int N_power = 5 [];     // directional spreading coeff
double L = 200.0;       // size of the domain
double kp=1.0;          // peak wave length



int main(int argc, char *argv[]){
  
  // read parameters from namelist
  params = array_new();
  add_param("P", &P, "double");
  add_param("coeff_kpL0", &coeff_kpL0, "int");
  add_param("N_mode", &N_mode, "int");
  add_param("N_power", &N_power, "int");

  // Search for the configuration file with a given path or read params.in
  if (argc == 2)
    strcpy(file_param, argv[1]);
  else
    strcpy(file_param, namlist);
  read_params(file_param);
  
  // setting values from namelist
  kp = PI * coeff_kpL0 / L; 

  // generate spectrum F(kx,ky)
  T_Spectrum spectrum;
  spectrum = spectrum_gen_linear(N_mode, N_power, L, P, kp);

  // save files
  FILE *fptr1 = fopen("F_kxky", "wb");
  fwrite(spectrum.F_kxky, sizeof(double), N_mode*(N_mode+1), fptr1);
  fclose(fptr1);

  FILE *fptr2 = fopen("kx", "wb");
  fwrite(spectrum.kx, sizeof(double), N_mode, fptr2);
  fclose(fptr2);

  FILE *fptr3 = fopen("ky", "wb");
  fwrite(spectrum.ky, sizeof(double), N_mode+1, fptr3);
  fclose(fptr3);

  FILE *fptr4 = fopen("omega", "wb");
  fwrite(spectrum.omega, sizeof(double), N_mode*(N_mode+1), fptr4);
  fclose(fptr4);

  FILE *fptr5 = fopen("phase", "wb");
  fwrite(spectrum.phase, sizeof(double), N_mode*(N_mode+1), fptr5);
  fclose(fptr5);
  
}

/**
## Results

~~~pythonplot Initial 2D spectrum
import numpy as np
import matplotlib.pyplot as plt 

kx = np.fromfile("kx")
ky = np.fromfile("ky")
raw = np.fromfile("F_kxky")


datapy = np.zeros((len(kx),len(ky)))
for row in range(len(kx)):
    for col in range(len(ky)):
        index = row*len(kx) + col
        #index = col*len(kx_py) + row
        datapy[row, col] = raw[index]

fig, ax = plt.subplots(1,1,figsize = (7,5),constrained_layout=True,dpi=100)
s=ax.pcolormesh(kx,ky,datapy.T, shading='nearest', vmin=0,vmax=3)
ax.set_xlabel('kx')
ax.set_ylabel('ky')
plt.colorbar(s,ax=ax)
fig.savefig('F_kxky.png')
~~~

**/
