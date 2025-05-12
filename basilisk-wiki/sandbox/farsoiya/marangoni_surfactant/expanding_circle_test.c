/**
# Expansion of circle

The uniformly distributed surfactant initially redistributes as the circle expands for a given radial velocity. 
[Stone, 1990](#stone1990simple) and [Atasi et. al., 2018](#atasi2018influence)
*/

#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "surfactant-transport.h"
#include "view.h"


double R0 = 0.2;
double K = 1.;
double Pe = 16.;
int LEVEL = 6;
int main (int argc, char * argv[])
{
  
  /**
  The domain is a square box centered on the origin and of length $L0=1$ */
   
  origin (-L0/2., -L0/2.);

   //Redistance at each comp step
  reinit_skip_steps = 1; 
  //If redistancing not done every step then use transport equation for phase-field
  advect_diff_phase_field = 1; 
  D_s = sq(2.*R0)*K/Pe;
  // Iterate over 3 resolutions
  for (LEVEL = 5; LEVEL <= 7; LEVEL++) {
    N = 1 << LEVEL;
    counter = 0;
    run();
  }

}



/**
We initialize a circle */
event init (i = 0)
{
     refine( (sq(0.3)- pow(x,2) - pow(y,2)) > 0 && level < LEVEL);

    fraction(f, -(sq(0.2)- pow(x,2) - pow(y,2)) );
    event ("properties2");
  
 // Initialize surfactant 
  foreach(){ 
  //  d2[] = (0.2) - sqrt(pow(x,2) + pow(y,2)) ;
  //  pfield[] = 0.5*(1. - tanh((d2[])/2./EPSILON));
  //  pfield[] = clamp(pfield[], 0., 1.);
    double deltas = (pfield[]*(1. - pfield[]))/EPSILON;
    double gamma0 = 1.;
    c1[] = gamma0*deltas;
  }

  if (grid->maxdepth == 7){
    //  view (fov = 30, near = 0.01, far = 1000,
	  // tx = 0.009, ty = -0.076, tz = -0.291,
	  // width = 1239, height = 575);
    draw_vof (c = "f", lw = 2);
    squares (color = "c1", spread = -1, linear = true);
    vectors (u = "u", scale = 1);
    
    save ("fields-initial.png");
  }
 
}

// Radial Velocity 
event stability (i++) {

	trash ({u,uf});

	foreach_face(y) uf.x[] = K*(sqrt(x*x + y*y))*cos(atan2(y,x));
	foreach_face(x) uf.y[] = K*(sqrt(x*x + y*y))*sin(atan2(y,x));

	foreach(){
		u.x[] = K*(sqrt(x*x + y*y))*cos(atan2(y,x));
		u.y[] = K*(sqrt(x*x + y*y))*sin(atan2(y,x));
	}

}

/*
event adapt (i++) {

	scalar pfield2[];
	foreach(){
     pfield2[] = 0.5*(1. - tanh((d2[])/4./EPSILON));
     //c2[] = clamp(pfield[], 0., 1.);	
 }

 double uxemax = 0.1*statsf(u.x).max;
 double uyemax = 0.1*statsf(u.y).max;
// double uzemax = 0.1*statsf(u.z).max;
 double c1emax = 0.001*statsf(c1).max;
 double pfemax = 0.001*statsf(pfield).max;
 double pf2emax = 0.001*statsf(pfield2).max;
 double d2emax = 0.1*statsf(d2).max;

  adapt_wavelet ({c1,pfield2}, (double[]){c1emax,pf2emax}, LEVEL); 
//  printf("\n%d %g",i,t); fflush(stdout);
}
*/

int count = 0;
// event logfile(i += 1){
event logfile(t = 0; t += 0.01; t <= 0.7){
double cnet = 0.; //net surfactant
scalar deltas[];
  foreach(){
   // pfield[] = 0.5*(1. - tanh((d2[])/2./EPSILON));
   // pfield[] = clamp(pfield[], 0., 1.);
    deltas[] = pfield[]*(1. - pfield[])/EPSILON;
    cnet += c1[]*dv();
    gamma2[] = c1[]*4.*EPSILON;
  }
  char filename[200];
  sprintf(filename,"expanding-circle-test-%d",(1 << grid->maxdepth));
  static FILE * fp2 = fopen(filename,"w");
  //double theta = M_PI/2.;
  fprintf (fp2, "%.17g %.17g %.17g\n", t, statsf(c1).max/statsf(deltas).max, cnet);
  fflush (fp2);

	
  count++;

}

event end( t = 0.7){

  if (grid->maxdepth == 7){
    //  view (fov = 30, near = 0.01, far = 1000,
	  // tx = 0.009, ty = -0.076, tz = -0.291,
	  // width = 1239, height = 575);
    draw_vof (c = "f", lw = 2);
    squares (color = "gamma2", spread = -1., min = 0, max = 1, linear = true);
    // vectors (u = "u", scale = 1);
    
    save ("fields-end.png");
  }
  // save();
}

/**
# Surfactant evolution over the circle
![Surfactant concentration and Interface at t = 0](expanding_circle_test/fields-initial.png){ width="20%" }
![Surfactant concentration and Interface at t = 0.8](expanding_circle_test/fields-end.png){ width="20%" }
![Time evolution of Surfactant concentration](expanding_circle_test/expanding_circle_test.png){ width="50%" }
![Time evolution of Surfactant concentration](expanding_circle_test/expanding_circle_test_mc.png){ width="50%" }
![Time evolution of Surfactant concentration](expanding_circle_test/expanding_circle_test_order.png){ width="50%" }
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


plt.plot(ts32[1:-1:2]*5,cnet32[1:-1:2]/cnet32[1],'bx',label='$1/6.4 R_0$');

plt.plot(ts64[2:-1:2]*5,cnet64[2:-1:2]/cnet64[1],'gx',label='$1/12.8 R_0$');

plt.plot(ts128[3:-1:2]*5,cnet128[3:-1:2]/cnet128[1],'kx',label='$1/25.6 R_0$');
plt.legend();

plt.ylim(0.99,1.01)
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
