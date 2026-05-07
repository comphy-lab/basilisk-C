#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "../soluble-surfactant.h"
#include "view.h"

double R0 = 0.2;
double K = 1.;
double Pe = 16.;
int LEVEL = 6;

int main() {
  origin (-L0/2., -L0/2.);
  N = 256;
  reinit_skip_steps = 1; 
  advect_diff_phase_field = 1; 
  Di0 = sq(2.*R0)*K/Pe;
  Db0 = 0.;
  ra = 0.;
  rd = 0.;

  for (LEVEL = 5; LEVEL <= 7; LEVEL++) {
    N = 1 << LEVEL;
    counter = 0;
    run();
  }
}

event init(i = 0) {
  refine((sq(0.3) - pow(x, 2) - pow(y, 2)) > 0 && level < LEVEL);
  fraction(f, -(sq(0.2) - pow(x, 2) - pow(y, 2)));
  event("properties2");
  
  foreach() { 
    double deltas = (pfield[]*(1. - pfield[]))/EPSILON;
    double gamma0 = 1.;
    ci[] = gamma0*deltas;
    cb[] = 0.;
  }
  if (grid->maxdepth == 7){
    draw_vof (c = "f", lw = 2);
    squares (color = "ci", spread = -1, linear = true);
    vectors (u = "u", scale = 1);
    
    save ("fields-initial.png");
  }
}

event stability(i++) {
  trash({u, uf});
  foreach_face(y) uf.x[] = K*(sqrt(x*x + y*y))*cos(atan2(y,x));
  foreach_face(x) uf.y[] = K*(sqrt(x*x + y*y))*sin(atan2(y,x));
  
  foreach() {
    u.x[] = K*(sqrt(x*x + y*y))*cos(atan2(y,x));
    u.y[] = K*(sqrt(x*x + y*y))*sin(atan2(y,x));
  }
}

int count = 0;
event logfile(t = 0; t += 0.01; t <= 0.7) {
  double cnet = 0.; //net surfactant
  scalar deltas[];

  foreach() {
    deltas[] = pfield[]*(1. - pfield[])/EPSILON;
    cnet += ci[]*dv();
    gamma2[] = ci[]*4.*EPSILON;
  }
    
  char filename[200];
  sprintf(filename,"expanding-circle-test-%d",(1 << grid->maxdepth));
  static FILE * fp2 = fopen(filename,"w");
  fprintf (fp2, "%.17g %.17g %.17g\n", t, statsf(ci).max/statsf(deltas).max, cnet);
  fflush (fp2);
  count++;
}

event end( t = 0.7){
  if (grid->maxdepth == 7){
    draw_vof (c = "f", lw = 2);
    squares (color = "gamma2", spread = -1., min = 0, max = 1, linear = true);
    save ("fields-end.png");
  }
}

/**
# Surfactant evolution over the circle
![Surfactant concentration and Interface at t = 0](expanding_circle/fields-initial.png){ width="20%" }
![Surfactant concentration and Interface at t = 0.8](expanding_circle/fields-end.png){ width="20%" }
![Time evolution of Surfactant concentration](expanding_circle/expanding_circle_test.png){ width="50%" }
![Time evolution of Surfactant concentration](expanding_circle/expanding_circle_test_mc.png){ width="50%" }
![Time evolution of Surfactant concentration](expanding_circle/expanding_circle_test_order.png){ width="50%" }
~~~pythonplot

import numpy as np
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 15})	
import matplotlib.ticker as mticker


plt.figure()


ts32, gammat32,cnet32 = np.loadtxt('expanding-circle-test-32',delimiter=' ',unpack=True)
ts64, gammat64,cnet64 = np.loadtxt('expanding-circle-test-64',delimiter=' ',unpack=True)
ts128, gammat128,cnet128 = np.loadtxt('expanding-circle-test-128',delimiter=' ',unpack=True)

anexp = np.exp(-ts128);

plt.plot(ts128,anexp,'r',label='$exp(-U_r t)$');

plt.plot(ts32[1:-1:2],gammat32[1:-1:2],'bx',label='$1/6.4 R_0$');

plt.plot(ts64[1:-1:2],gammat64[1:-1:2],'gx',label='$1/12.8 R_0$');

plt.plot(ts128[1:-1:2],gammat128[1:-1:2],'kx',label='$1/25.6 R_0$');


plt.ylim(0,1.2)
plt.xlim(0,0.8)

plt.legend();
plt.xlabel(r'$tU_r/R_0$')
plt.ylabel(r'$\Gamma(t)$')
plt.tight_layout()

plt.savefig('expanding_circle_test.png')
plt.figure()


plt.semilogy(ts32[1:-1:2]*5,cnet32[1:-1:2]/cnet32[1],'bx',label='$1/6.4 R_0$');

plt.semilogy(ts64[2:-1:2]*5,cnet64[2:-1:2]/cnet64[1],'gx',label='$1/12.8 R_0$');

plt.semilogy(ts128[3:-1:2]*5,cnet128[3:-1:2]/cnet128[1],'kx',label='$1/25.6 R_0$');
plt.legend();

#plt.ylim(0.99,1.01)
plt.xlim(0,3.5)
plt.xlabel(r'$tU_r/R_0$')
plt.ylabel(r'$\int c\; dv$')

plt.savefig('expanding_circle_test_mc.png')
plt.figure()

error32 = np.max(np.abs((anexp - gammat32) / anexp))
error64 = np.max(np.abs((anexp - gammat64) / anexp))
error128 = np.max(np.abs((anexp - gammat128) / anexp))
dx = [6.4, 12.8, 25.6]
erv = [error32, error64, error128]
ax = plt.loglog(dx, erv, 'r-x')
order = -2;
yl =  [x ** order for x in dx]
ylp = [x * 1  for x in yl]

plt.loglog(dx, ylp, 'k')
plt.axis([5, 30, 1e-3, 1])
plt.legend(['error', '2nd order'])
plt.xlabel('$R_0/\Delta x$',fontsize = 15)

ax = plt.gca()

ax.xaxis.set_major_formatter(mticker.ScalarFormatter())
ax.xaxis.get_major_formatter().set_scientific(False)
ax.xaxis.get_major_formatter().set_useOffset(False)

ax.xaxis.set_minor_formatter(mticker.ScalarFormatter())
ax.xaxis.get_minor_formatter().set_scientific(False)
ax.xaxis.get_minor_formatter().set_useOffset(False)



plt.ylabel('$||(\Gamma_c - \Gamma_{exact})/\Gamma_{exact}||_{\infty}$',fontsize=15)
plt.tight_layout()

plt.savefig('expanding_circle_test_order.png')
plt.figure()

~~~

## References

~~~bib

@article{stone1990simple,
  title={A simple derivation of the time-dependent convective-diffusion equation for surfactant transport along a deforming interface},
  author={Stone, HA},
  journal={Physics of Fluids A: Fluid Dynamics},
  volume={2},
  number={1},
  pages={111--112},
  year={1990},
  publisher={American Institute of Physics}
}

@article{atasi2018influence,
  title={Influence of soluble surfactants and deformation on the dynamics of centered bubbles in cylindrical microchannels},
  author={Atasi, Omer and Haut, Benoit and Pedrono, Annaig and Scheid, Benoit and Legendre, Dominique},
  journal={Langmuir},
  volume={34},
  number={34},
  pages={10048--10062},
  year={2018},
  publisher={ACS Publications}
}
~~~

*/
