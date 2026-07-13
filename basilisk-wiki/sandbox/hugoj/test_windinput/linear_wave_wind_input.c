/**
# Wave growth: from linear to breaking

## Analytical solution

A linear surface gravity wave will decay in a viscous fluid. The rate of this
decay is $E(t)=E_0 e^{-4\nu k^2 t}$ (Lamb 1932)

## Observations

Re=40000
- N=1024 nl=15 n'est pas convergé.
- 1 mode, N=2048 nl=15 fait apparaitre des vitesses verticale très fortes ...
- 1 mode, N=4096 pareil
Instabilité à lambda / 4 ?? Pour lambda=1, L=1, instabilité à -0.25, 0, 0.25

- 2 modes, N=2048 ok !
- Diminuer Re fait augmenter les instabilités ?? + message "warning: could not conserve barotropic flux"


Cas qui fonctionne mais trop dissipatif:
k=2pi L=1 ak=0.05 g=9.81 theta_H=0.6 Re=40000 N=1024 nl=15

## References

Horace Lamb, Hydrodynamics (6th ed., 1932), Chapter XI, Article 348, "Effect of
Viscosity on Water-Waves," pp. 623–625

et Article 349 pour une dérivation plus sérieuse


*/

#include "grid/multigrid1D.h"
#include "layered/hydro.h"
#include "layered/nh.h"
#include "layered/remap.h"
#include "layered/perfs.h"
#include "bderembl/libs/netcdf_bas.h"

#if STOKES
#include "stokes.h" // link from src/test/stokes.h
double k_ = 2.*pi, h_ = 0.5, g_ = 1.0, ak = 0.35; 
double RE = 40000.;
#else
double k_ = 2.*pi, h_ = 0.5, g_ = 9.81, ak = 0.05; 
double RE = 40000.;
#endif

#define T0  (2.*pi/sqrt(g_*k_))
double lam, p0;
double base_gpe=0.;
double NT0 = 100.;
scalar ps[], Fp[]; 



int main()
{
  L0=1.;
  origin (-L0/2.);
  periodic (right);
  N = 512;
  nl = 15; // 60
  G = g_;
  lam = 2*pi/k_*L0;
  nu = sqrt(g_*lam*lam*lam)/RE; // sqrt(g_*lam*lam*lam)
  CFL_H = 1.0;
  //max_slope = 1.; // a bit less dissipative
  theta_H = 0.55; // >0.5 or else energy increases even without forcing ? 0.51
  p0=1.; //0.001; //4*1*nu*k_*sqrt(9.81*k_)/k_;
  base_gpe=0.5*g_*h_*h_;
  //printf("%f", nu);
  run();
}



event init (i = 0)
{
  foreach() {
    zb[] = -h_*lam;
#if STOKES
#warning "STOKES branch active"
    double H =  wave(x, 0) - zb[];  
#else
    double H =  ak/k_*cos(k_*x) - zb[];
#endif
    double z = zb[];
    foreach_layer() {
      h[] = H/nl;
      z += h[]/2.;
#if STOKES
      u.x[] = u_x(x, z);  
      w[] = u_y(x, z);  
#else
      u.x[] = ak/k_*sqrt(g_*k_)*exp(k_*z)*cos(k_*x); // 
      w[] = ak/k_*sqrt(g_*k_)*exp(k_*z)*sin(k_*x);  // 
#endif
      z += h[]/2.;
    }
  }
  create_nc({zb, eta, h, u.x, w, Fp}, "out.nc");
}

// event profiles (t += T0/4.) { 
//   foreach (serial) {
//     double H = zb[];
//     foreach_layer()
//       H += h[];
//     fprintf (stderr, "%g %g\n", x, H);
//   }
//   fprintf (stderr, "\n\n");
// }

#if !STOKES
event viscous_term (i++) {
  // vertical diffusion is done in diffusion.h (u.x) and nh.h (w)
  horizontal_diffusion ({u.x, w}, nu, dt);
}
#endif


#if 1
#if FORCING
#warning "FORCING accel branch active"
event acceleration(i++, last){
  double detadx = 0.;
  foreach() {
    detadx = (eta[1] - eta[-1])/(2*Delta);
    ha.x[0,0,nl-1] -= hf.x[]*(p0 * detadx * sin(atan(detadx))); //   
  }
}
#endif
#endif

#if 0
#if FORCING
#warning "FORCING accel branch active"
event acceleration(i++, last){
  double detadx = 0.;
  foreach_face(x) {
    detadx = (eta[] - eta[-1])/Delta;
    ha.x[0,0,nl-1] += hf.x[]*(p0 * detadx * sin(atan(detadx))); //   
  }
}
#endif
#endif

#if 0
#if FORCING
#warning "FORCING diffusion branch active"
event forcing(i++, last){
  foreach(){
    double detadx = (eta[1] - eta[-1])/(2*Delta);
    ps[] = p0 * detadx;
    Fp[] = ps[]* sin(atan(detadx));
    vertical_diffusion (point, h, u.x, dt, nu,
			                  Fp[],      // my forcing
                        u_b.x[], lambda_b.x[]); // keep the default values
  }
}
#endif
#endif 

event logfile (i++; t <= NT0*T0) // target: at least 300*T0
{
  double ke = 0., gpe = 0.;
  foreach (reduction(+:ke) reduction(+:gpe)) {
    foreach_layer() {
      double norm2 = sq(w[]);
      foreach_dimension()
	norm2 += sq(u.x[]);
      ke += norm2*h[]*dv();
    }
    gpe += sq(eta[])*dv();
  }
  printf ("%g %g %g\n", t/T0, ke/2., g_*gpe/2.);
}


event movie (i+=3)
{
  static FILE * fp = popen ("gnuplot", "w");
  if (i == 0)
    fprintf (fp, "set term pngcairo font ',9' size 800,250;"
      "set size ratio -1\n");  
  fprintf (fp,
    "set output 'plot%04d.png'\n"
    "set title 't = %.2f'\n"
    "p [%g:%g][-0.1:0.15]'-' u 1:(-1):2 w filledcu lc 3 t ''",
    i/3, t/(k_/sqrt(g_*k_)), X0, X0 + L0);
  fprintf (fp, "\n");
  foreach (serial) {
    double H = 0.;
    foreach_layer()
      H += h[];
    fprintf (fp, "%g %g %g", x, zb[] + H, zb[]);
    fprintf (fp, "\n");
  }
  fprintf (fp, "e\n\n");
  fflush (fp);
  write_nc();
}

event moviemaker (t = end)
{
  system ("rm -f movie.mp4 && "
	  "ffmpeg -r 25 -f image2 -i plot%04d.png "
	  "-vcodec libx264 -vf format=yuv420p -movflags +faststart "
	  "movie.mp4 2> /dev/null && "
	  "rm -f plot*.png");
}

/**
~~~pythonplot Wave energy
import numpy as np
import matplotlib.pyplot as plt

g = 9.81
L0 = 1
ak = 0.05
k = 2*np.pi
lam = 2*np.pi/k*L0
Re = 40000

nu = np.sqrt(g*lam**3)/Re
T0 = (2*np.pi/np.sqrt(g*k))

def E_linwave(E0,nu,ak,k,t):
  print('\nWave decay for linear wave (theory)')
  print(f'nu={nu},ak={ak},k={k/np.pi}pi\n')
  return E0*np.exp(-4*nu*k**2*t)

data = np.loadtxt("out",skiprows=1)
data_ref = np.loadtxt('../no_forcing/out',skiprows=1)
time = data[:,0]
ke  = data[:,1]
gpe = data[:,2]
E = ke + gpe
timeref = data_ref[:,0]
Eref = data_ref[:,1]+data_ref[:,2]
E0 = E[0]
Eth = E_linwave(E0,nu,ak,k,time*T0)
fig, ax = plt.subplots(figsize=(5, 5))
ax.plot(timeref,Eref/Eref[0],c='k',ls='--',label='no forcing')
ax.plot(time,2*ke/E[0],color='b', label='2*ke')
ax.plot(time,2*gpe/E[0],color='g', label='2*gpe')
ax.plot(time,E/E[0], color='k', label='E')
ax.plot(time, Eth/E[0], color='r', label=r'$E(t)=E_0 e^{-4 \nu k^2 t}$')
ax.set_xlabel("t/T0")
ax.set_ylabel("E/E0")
ax.legend(loc="lower left")
plt.tight_layout()
plt.savefig("energy.png", dpi=150)
plt.show()
~~~
*/
