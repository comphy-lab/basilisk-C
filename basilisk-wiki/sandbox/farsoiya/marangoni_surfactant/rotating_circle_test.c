/**
# Rotation of circle

The uniformly distributed surfactant initially distributed as sinusoid. The solid body rotation of circle rotates the field and diffuses along
 [Atasi et. al., 2018](#atasi2018influence)
*/

#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "surfactant-transport.h"
#include "view.h"


   //Boundary Conditions
c1[top] = dirichlet(0.);
c1[bottom] = dirichlet(0.);
c1[left] = dirichlet(0.);
c1[right] = dirichlet(0.);

   //Boundary Conditions
pfield[top] = dirichlet(0.);
pfield[bottom] = dirichlet(0.);
pfield[left] = dirichlet(0.);
pfield[right] = dirichlet(0.);

int counter2 = 0;
double db = 0.4 [*];
double Pe = 160. [*];
double omega = 10. [*];

int main (int argc, char * argv[])
{
  
  /**
  The domain is a square box centered on the origin and of length $L0=1$ */
   
   //Redistance at each comp step
  reinit_skip_steps = 10; 
  //If redistancing not done every step then use transport equation for phase-field
  advect_diff_phase_field = 1; 

  //Peclet number
  
  D_s = sq(db)*omega/Pe;

  //Do not advect velocity (since it is reinitialised at each timestep)
  stokes = true;
  
  //Use strict minmod limiting when advecting c1
  theta = 1.;
  c1.gradient = minmod2;
 
  // Iterate over 3 resolutions
  for (int LEVEL = 5; LEVEL <= 7; LEVEL++) {

    N = 1 << LEVEL;
    origin (-L0/2. + L0/N/2., -L0/2 + L0/N/2.);
    counter2 = 0;
    counter = 0;
    run();
  }

}



/**
We initialize a circle */

event init (i = 0)
{
    fraction(f, -(sq(db/2.)- pow(x,2) - pow(y,2)) );
    event ("properties2");

  
 // Thickness of the phase-field over which surfactant is distributed
  // Initialize surfactant 
  foreach(){ 
    d2[] = (db/2.) - sqrt(pow(x,2) + pow(y,2)) ;
    pfield[] = 0.5*(1. - tanh((d2[])/2./EPSILON));
    pfield[] = clamp(pfield[], 0., 1.);
    double deltas = (pfield[]*(1. - pfield[]))/EPSILON;
    double gamma0 = 2. + sin(atan2(y,x));
    c1[] = gamma0*deltas;
    gamma2[] = c1[]*4.*EPSILON;
  }
  event("stability");
  if (grid->maxdepth == 7){
    //  view (fov = 30, near = 0.01, far = 1000,
	  // tx = 0.009, ty = -0.076, tz = -0.291,
	  // width = 1239, height = 575);
    draw_vof (c = "f", lw = 2);
     squares (color = "gamma2", min = 0, max = 3, spread = -1., linear = true);
  
    
    save ("rotating_fields_initial.png");
  }

}

// Radial Velocity 
event stability (i++) {
  foreach_face (x)
    uf.x[] = omega*y;
  foreach_face (y)
    uf.y[] = - omega*x;
  foreach(){
    u.x[] =   omega*y;
    u.y[] = - omega*x;
  }
}

// event logfile(i += 1){
event logfile(t = 0; t += 0.01; t <= 2.){
  double cnet = 0;
  scalar deltas[];
  foreach() {
    deltas[] = pfield[]*(1. - pfield[])/EPSILON;
    cnet += c1[]*dv();
  }
  
  char filename[200];
  sprintf(filename,"rotating-circle-test-%d",(1 << grid->maxdepth));
  static FILE * fp2 = fopen(filename,"w");
  double theta = M_PI/2.;
  double c_piby2 = interpolate(c1, db/2.*cos(theta), db/2.*sin(theta), 0.);
  double deltas_piby2 = interpolate(deltas, db/2.*cos(theta), db/2.*sin(theta), 0.);
  fprintf (fp2, "%.17g %.17g %.17g\n", t, c_piby2/deltas_piby2,cnet);
  //  fprintf (fp2, "%.17g %.17g %.17g\n", t, gamma2_theta,cnet);
  fflush (fp2);



  counter2++;
//  exit(0);

}

event end( t = 2.){

  if (grid->maxdepth == 7){
    //  view (fov = 30, near = 0.01, far = 1000,
	  // tx = 0.009, ty = -0.076, tz = -0.291,
	  // width = 1239, height = 575);
    draw_vof (c = "f", lw = 2);
    squares (color = "gamma2", min = 0, max = 3, spread = -1., linear = true);
    // vectors (u = "u", scale = 1);
    save ("rotating_fields_end.png");
  }
  // save();
}

/**
# Surfactant evolution over the circle


![ t = 0](rotating_circle_test/rotating_fields_initial.png){ width="40%" }
![ t = 3.2](rotating_circle_test/rotating_fields_end.png){ width="40%" }
![ Surfactant concentration](rotating_circle_test/rotating_circle_test.svg)
![ Relative error mass](rotating_circle_test/rotating_circle_test_mc.svg)

~~~pythonplot

import numpy as np
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 10})	
import matplotlib.ticker as mticker
	
plt.figure()

omega = 10.;
tn = 2.*np.pi/omega;

ts32, gammat32, cnet32 = np.loadtxt('rotating-circle-test-32',delimiter=' ',unpack=True)
ts64, gammat64, cnet64 = np.loadtxt('rotating-circle-test-64',delimiter=' ',unpack=True)
ts128, gammat128, cnet128 = np.loadtxt('rotating-circle-test-128',delimiter=' ',unpack=True)


Pe = 160;
anexp = 2. + np.cos(omega*ts64)*np.exp(-4.*omega*ts64/Pe);
plt.plot(ts64[0:-1:1]/tn,anexp[0:-1:1],'r',label='Analytical');

plt.plot(ts32[0:-1:1]/tn,gammat32[0:-1:1],'b',label='$1/6.4 d_b$');

plt.plot(ts64[0:-1:1]/tn,gammat64[0:-1:1],'g',label='$1/12.8 d_b$');

plt.plot(ts128[0:-1:1]/tn,gammat128[0:-1:1],'k',label='$1/25/6 d_b$');




plt.ylim(2.9,3.1)
plt.xlim(0,0.1)

plt.legend(fontsize=10);
plt.xlabel(r'$t \omega /(2\pi)$')
plt.ylabel(r'$\Gamma(t)$')
plt.tight_layout()

plt.savefig('rotating_circle_test.svg')

plt.figure()


plt.figure()

plt.semilogy(ts32[0:-1:2]*5,1 - cnet32[0:-1:2]/cnet32[0],'bx',label='$1/12.8 d_b$');

plt.semilogy(ts64[0:-1:2]*5,1 - cnet64[0:-1:2]/cnet64[0],'gx',label='$1/25.6 d_b$');

plt.semilogy(ts128[0:-1:2]*5,1 - cnet128[0:-1:2]/cnet128[0],'kx',label='$1/51.2 d_b$');
plt.legend();

#plt.ylim(-1e-3,1e-3)
#plt.xlim(0,3.5)
plt.xlabel(r'$t\omega/(2\pi)$')
plt.ylabel(r'$\int c\; dv$')

plt.savefig('rotating_circle_test_mc.svg')

plt.figure()

error32 = np.max(np.abs((anexp - gammat32) / anexp))
error64 = np.max(np.abs((anexp - gammat64) / anexp))
error128 = np.max(np.abs((anexp - gammat128) / anexp))

d0bydx = [6.4*2, 12.8*2, 25.6*2] #number of grid points in bubble diameter
erv = [error32, error64, error128]
ax = plt.loglog(d0bydx, erv, 'r-x')
order = -2;
d0bydx_for_line = [9*2,10*2,20*2]
yl =  [x ** order for x in d0bydx_for_line]
ylp = [x * 2  for x in yl]

plt.loglog(d0bydx_for_line, ylp, 'k')
#plt.axis([1, 80, 0.01, 0.2])
plt.legend(['error', '2nd order'])
plt.xlabel('$d_0/\Delta x$',fontsize = 15)

ax = plt.gca()

ax.xaxis.set_major_formatter(mticker.ScalarFormatter())
ax.xaxis.get_major_formatter().set_scientific(False)
ax.xaxis.get_major_formatter().set_useOffset(False)

ax.xaxis.set_minor_formatter(mticker.ScalarFormatter())
ax.xaxis.get_minor_formatter().set_scientific(False)
ax.xaxis.get_minor_formatter().set_useOffset(False)



plt.ylabel('$||(\Gamma_c - \Gamma_{exact})/\Gamma_{exact}||_{\infty}$',fontsize=15)

plt.savefig('rotating_circle_test_order.svg')

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
