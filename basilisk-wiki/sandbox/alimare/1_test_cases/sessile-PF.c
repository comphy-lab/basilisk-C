/**
# Sessile drop

A sessile drop is a drop of liquid at rest on a solid surface. In the
absence of gravity, the shape of the drop is controlled by surface
tension only. An important parameter is the "contact angle" $\theta$ between
the solid surface and the interface. In the absence of gravity, the
drop is hemispherical and it is easy to show that the relation between
the radius of the drop $R$ and its volume $V$ is (for two-dimensional
drops)
$$
V = R^2 (\theta - \sin\theta\cos\theta)
$$

To test this relation, a drop is initialised as a half-disk (i.e. the
initial contact angle is 90$^\circ$) and the contact angle is varied
between 15$^\circ$ and 165$^\circ$. The drop oscillates and eventually relaxes
to its equilibrium position. This equilibrium is exact to within
machine accuracy. The curvature along the interface is constant.

![Animation of the density](sessile-PF/rho.mp4)


Here we plot the internal energy of during our calculation:

$$
e(\rho,\theta) = 4\theta\rho - 3 \rho^{2}
$$

The kinetic energy 
$$
ev(\rho,v) = \frac12 \rho v^2
$$

and the interfacial energy:

$$
ei(\rho,lambda) = 0.5 \lambda |\nabla \rho|^2
$$


~~~gnuplot Internal energy
set xlabel "time"
set ylabel "energy"
plot 'log' u 1:2 w l t 'e (sum)'
~~~

~~~gnuplot Kinetic energy
set xlabel "time"
set ylabel "energy"
plot 'log' u 1:3 w l t 'ev'
~~~

~~~gnuplot Interfacial energy
set xlabel "time"
set ylabel "energy"
plot 'log' u 1:4 w l t 'ev'
~~~

First we reach a minimum and then converge toward a plateau.

*/
double THETAE =  70.*M_PI/180.;
#define MU 1e-1
#include "../vdw.h"
#include "view.h"
#include "../diverging_cool_warm.h"
#define LEVEL 8

#define LAMBDA 3.e-2
#define rho_c 1 
#define theta 0.95
#define p_c 1

#include "../weno5.h"

double P0(double x)
{
  double rhop;
  rhop=x/rho_c;
  return p_c*rhop*theta*((8/(3-rhop) - 3*rhop/theta));
  // return 0.5*sq(x);
}

double myContactAngle(Point point, scalar rho, scalar normGradRho, 
  double THETAE){
  return max(0.,exp(-20*sq(rho[]- rho_c))*cos(THETAE)*normGradRho[]-1.e-3);
}

int main()
{
  // periodic(top);
  init_grid (1 << LEVEL);
  // L0 = 2.e-4;
  L0 = 2;
  origin (-L0/2.);
  DT=0.5e-3;
  mgu.nrelax = 20;
  lambda=LAMBDA;
  TOLERANCE = 1.e-6;
  periodic(right);
  run();
}

/**
The initial drop is a quarter of a circle. */

event init (t = 0)
{
    foreach()
      {
        double dist = (sq(x) + sq(y) - sq(0.10*L0)); /*sgined squared distance*/
        foreach_dimension()
          mu.x[] = MU;
        rho[] = rho_c *(1- 0.4*tanh(360/sq(L0)*dist));
        q.x[] = q.y[] = 0.;
      }
    boundary(all);

  rho[bottom]  = neumann(myContactAngle(point,rho,normGradRho,THETAE));
  rhop[bottom] = neumann(myContactAngle(point,rhop,normGradRho,THETAE));
  u.n[bottom] = dirichlet(0.);
  u.t[bottom] = dirichlet(0.);

}


event adapt(i++,last){
  adapt_wavelet ({rho,q.x,q.y},
    (double[]){1.e-3,1.e-2,1.e-2},LEVEL, 3);
}


event movie(t+=6.e-2,last){
  view(ty = -0.5);
  isoline("rho", rho_c);

  scalar logGrad[];

  foreach()
    logGrad[] = log(normGradRho[]);

  squares("logGrad", min = -3, max = 2,  map = mycoolwarm);
  // squares("u.x*u.x+u.y*u.y");
  stats s = statsf(rho);
  save("logGrad.mp4");
  fprintf(stderr, "##%g %g %g\n", t,s.min, s.max);
  isoline("rho", rho_c);
  s = statsf(u.x);
  squares("u.x", min = 0.5*s.min, max = 0.5*s.max,  map = mycoolwarm);
  save("ux.mp4");
  s = statsf(rho);
  isoline("rho", rho_c);
  squares("rho", min = 0.6, max = 1.4,  map = mycoolwarm);
  save("rho.mp4");
}

event logfile(i+=20){


/**
We output the internal energy of our system.
*/
  scalar e[],ev[],ePF[];
  vector gradrho[];
  foreach()
    foreach_dimension()
      gradrho.x[] = (WENO3_x (point, rho, -1) - WENO3_x (point, rho, 1));
  boundary ((scalar* ){gradrho});

  foreach(){
    e[] = 4*theta*rho[] - 3*sq(rho[]);
    ev[] = 0.;
    ePF[] = 0.;
    foreach_dimension(){
      ev[] += 0.5*rho[]*sq(u.x[]);
      ePF[] += 0.5*LAMBDA*sq(gradrho.x[]);
    }
  }

  stats s = statsf(e), s2 = statsf(ev), s3 = statsf(ePF);
  fprintf(stderr, "%g %g %g %g\n", t, s.sum,s2.sum,s3.sum);
}

event end (t = 100){}

/**

## See also

* [Similar test with
   Gerris](http://gerris.dalembert.upmc.fr/gerris/tests/tests/sessile.html)
*/
