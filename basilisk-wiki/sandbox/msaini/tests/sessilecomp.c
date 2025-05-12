 /**
# Sessile bubble

A sessile bubble is a bubble of liquid at rest on a solid surface. In the
absence of gravity, the shape of the bubble is controlled by surface
tension only. An important parameter is the "contact angle" $\theta$ between
the solid surface and the interface. In the absence of gravity, the
bubble is hemispherical and it is easy to show that the relation between
the radius of the bubble $R$ and its volume $V$ is (for two-dimensional
bubbles)
$$
V = R^2 (\theta - \sin\theta\cos\theta)
$$

To test this relation, a bubble is initialised as a half-disk (i.e. the
initial contact angle is 90$^\circ$) and the contact angle is varied
between 15$^\circ$ and 165$^\circ$. The bubble oscillates and eventually relaxes
to its equilibrium position. This equilibrium is exact to within
machine accuracy. The curvature along the interface is constant.

Note that shallower angles are [not accessible yet](/src/contact.h).

~~~gnuplot Equilibrium shapes for $15^\circ \leq \theta \leq 165^\circ$
set term push
set term @SVG size 640,180
set size ratio -1
unset key
unset xtics
unset ytics
unset border
set xrange [-1.6:1.6]
set yrange [0:]
plot 'out' u 1:2 w l, '' u (-1.*$1):2 w l lt 1, 0 lt -1
set term pop
~~~
*/

//#include "axi.h"
#include "../src/tpccontact.h"
#include "../src/compressible-tension.h"
//#include "view.h"

/**
To set the contact angle, we allocate a [height-function
field](/src/heights.h) and set the contact angle boundary condition on
its tangential component. */

vector h[];
double theta0 = 120.;
h.t[bottom] = contact_angle (theta0*pi/180.);
double p0,pg0;
double CFLac = 2.;

int main()
{
  size (2);

  /**
     We use a constant viscosity. */
  
  mu1 = 5.;
  mu2 = 5.;
  gamma1 = 1.4;
  gamma2 = 7.14;
  PI2 = 300;

  /**
     We must associate the height function field with the VOF tracer, so
     that it is used by the relevant functions (curvature calculation in
     particular). */
  
  f.height = h;

  /**
     We set the surface tension coefficient and run for the range of
     contact angles. */
  
  f.sigma = 1.;
  N = 128;
  for (theta0 = 15; theta0 <= 165; theta0 += 30)
    run();
}

event init (t = 0)
{
  fraction (f, (sq(x) + sq(y) - sq(0.5)));
  p0 = 1.;
  pg0 = 1. + 2.*f.sigma/0.5;
  foreach(){
    frho1[] = 1.*f[];
    frho2[] = 1.*(1. - f[]);
    p[] = pg0*f[] + p0*(1. - f[]);
    q.x[] = 0.;
    q.y[] = 0.;
    fE1[] = p[]/(gamma1 - 1.)*f[];
    fE2[] = (1. - f[])*(p[] + gamma2*PI2)/(gamma2 - 1.);
  }
  boundary({frho1,frho2,p,q,fE1,fE2});
}

event stability(i++)
{
  
  foreach()
    {
      double cspeed;
       if(f[] <= 0.0001)
  	     cspeed = sqrt(gamma2*(p[]+PI2)/frho2[]*(1.-f[]));
  	   else if(f[] >= 0.9999)
  	     cspeed = sqrt(gamma1*(p[]+PI1)/frho1[]*f[]);
  	   else
  	     cspeed = max(sqrt(gamma2*(p0+PI2)),sqrt(gamma1*(pg0)));
  	     double dtmaxac = CFLac*Delta/cspeed;
  	   dtmax = min(dtmax,dtmaxac);
    }
}

/* int count = 0.; */
/* event movie (t += 2./40){ */
/* view (fov = 22.5796, quat = {0,0,0,1}, tx = -0.475266, ty = -0.463264, bg = {1,1,1}, width = 600, height = 600, samples = 1); */
/*   clear(); */
/*   draw_vof ("f", lw = 2); */
/*   squares ("p", linear = true); */
/*   box (notics = true); */
/*   //cells(); */
/*   /\* mirror ({0,1}) { *\/ */
/*   /\*   draw_vof ("f", lw = 2); *\/ */
/*   /\*   squares ("unorm", linear = true); *\/ */
/*   /\*   box (notics = true); *\/ */
/*   /\* } *\/ */
    
/*   char name_exp[50]; */
/*   sprintf(name_exp,"movie-%5.5d.mp4",(int)theta0); */
/*   save (name_exp); */
/*   count ++ */
/* } */
 
/**
At equilibrium (t = 10 seems sufficient), we output the interface
shape and compute the (constant) curvature. */

event end (t = 10.)
{
  output_facets (f, stdout);
}

/**

## See also

* [Similar test with Gerris](http://gerris.dalembert.upmc.fr/gerris/tests/tests/sessile.html)
* [Similar test with incompressible solver](http://basilisk.fr/src/test/sessile.c)
*/
