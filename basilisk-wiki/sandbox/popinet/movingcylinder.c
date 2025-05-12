/**
# An example of moving solid

This is a fairly rough way of approximating boundary conditions on a
moving (or stationary) solid.

A cylinder moves through a liquid with a velocity corresponding to a
Reynolds number of 160 i.e. this is equivalent in principle to the
[von Karman vortex street example](/src/examples/karman.c) but in a
different reference frame. */

#include "navier-stokes/centered.h"
#include "fractions.h"

int main(){
  L0 = 8.;
  origin (-0.5, -L0/2.);
  N = 512;
  const face vector muc[] = {0.00078125,0.00078125};
  mu = muc;
  run();
}

/**
We just make an (empty) long channel. */

event init (t = 0) {
  mask(y > 0.5 ? top: y < -0.5 ? bottom : none);
}

/**
The moving cylinder is created by first initialising a volume fraction
field at each timestep corresponding to the position and shape of the
moving cylinder. Note that *cylinder* should really be a temporary
field inside *moving_cylinder()*. We just keep it as a global field to
be able to display it in the movie below. */

scalar cylinder[];

event moving_cylinder (i++) {
  coord vc = {1.,0.}; // the velocity of the cylinder
  fraction (cylinder, - (sq(x - vc.x*t) + sq(y - vc.y*t) - sq(0.0625)));

  /**
  We then use this (solid) volume fraction field to impose the
  corresponding velocity field inside the solid (as a volume-weighted 
  average). */
  
  foreach()
    foreach_dimension()
      u.x[] = cylinder[]*vc.x + (1. - cylinder[])*u.x[];
}

/**
We output some convergence statistics... */

event logfile (i++)
  fprintf (stderr, "%d %g %d %d %d %d\n", i, t,
	   mgp.i, mgp.nrelax, mgu.i, mgu.nrelax);

/**
... and a movie of vorticity. 

![Animation of the vorticity field.](movingcylinder/vort.gif)
*/

event images (t += 0.1; t <= 10.0) {
  static FILE * fp = popen ("ppm2gif > vort.gif", "w");
  scalar omega[];
  vorticity (u, omega);

  /**
  Cells for which *m* is negative will be black in the movie. */

  scalar m[];
  foreach()
    m[] = 0.5 - cylinder[];

  output_ppm (omega, fp, box = {{-0.5,-0.5},{7.5,0.5}}, mask = m,
	      min=-10, max=10, linear=true);
}

/**
We adapt according to the error on the velocity field as usual. */

event adapt (i++) {
  adapt_wavelet ((scalar *){u}, (double[]){3e-2,3e-2}, 9, 4);
}

/**
## See also

This method was used in Gerris in particular for these papers:

~~~bib
@InProceedings{wu2007,
  author =  {C. J. Wu and L. Wang},
  title =  {Direct numerical simulation of self-propelled 
swimming of 3D bionic fish school},
  booktitle =  {Computational Mechanics, Proceedings of ISCM 2007},
  year =  {2007},
  pages =  {1-17},
  url =  {http://gerris.dalembert.upmc.fr/Wu-Wang-2007-Comp_Mech.pdf},
}

@Article{Lin-Lin2016,
author={Lin-Lin, Zhu and Hui, Guan and Chui-Jie, Wu},
title={Three-dimensional numerical simulation of a bird 
model in unsteady flight},
journal={Computational Mechanics},
year={2016},
pages={1--11},
issn={1432-0924},
doi={10.1007/s00466-015-1233-3},
url={http://dx.doi.org/10.1007/s00466-015-1233-3}
}
~~~
*/
