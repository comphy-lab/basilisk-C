/**
# Granular shear cell with cylinder

We consider a fixed cylinder in a shear cell geometry filled with granular media with a confinement pressure at the top and gravity in the bulk.

![Domain and mesh](particle_granular2/schematic.png)

Note that, I have used a modified form of granular.h as granular_novof.h

## My notes:

24.5.24: 

- VVI Note: muv.x[] = fm.x[]*1.; is needed for embedded geometry.
- Granular media: using granular_novof.h: pass 
- Check Granular with G=0: Pass

25.5.24: 

- Check Granular with G=10, sym BC: Pass. drag-lift asymptote to a steady value (both negative !)
- Granular with G=10, asym. BC [Up=0.1]: Pass - cyl at sub-yield
- Granular with G=10, asym. BC [Up=2, maxlevel=12]: Code blows up @t=1.2
- Granular with G=10, asym. BC [Up=3.5, maxlevel=14]: Code blows up @t=18, noisy drag-lift

30.5.24: 

- changed interpolation scheme for muv.x[]: Code converges ! noisy drag-lift
- restriction on muI[]= cs[} ? (mus + dmu*In[]/(I0 + In[])) : 0. ==> no change. 
- double-projection.h needed to remove noisy drags force (see test/starting.c)  
 
 */

#include "embed.h"
#include "granular_novof.h"
#include "utils.h"

//Media properties in SI units.
#define RHO1 1000
#define G 10. 

//Reference scales in SI units.
#define Lref 7.5e-2
#define Pref 1600
#define Uref sqrt(Pref/RHO1)

// Dimensional input parameters in SI units
#define Upt 1.8
#define Upb 0.2
#define Pc Pref 


// Dimensionless input parameters
#define ust (Upt/Uref)			// top plate
#define usb (-Upb/Uref)			// bottom plate
#define cp (Pc/Pref)
#define FrSqInv (RHO1*G*Lref/Pc)	// 1/Froude^2

#define tend 10.

scalar shear[];
int maxlevel = 10;


int main() {
  periodic(right);			// Periodic BC left-right
  L0 = 1.;
  origin (-0.5, -L0/2.);
  N = 64;
  mu = muv;
  // Gravity as body force    
  const face vector g[] = {0.,-1.*FrSqInv};
  a = g;
  //DT=1e-3;
  run(); 
}


/**
Shear cell BCs

The velocities at the top and botom are chosen such that $u=0$ at the cylinder center, so that the drag force, $F_D\sim 0$
*/
//Velocity at the top
  u.n[top] = dirichlet(0);
  u.t[top] = dirichlet(ust);
  
//Confining pressure at the top
  p[top] = dirichlet(cp);
  pf[top] = dirichlet(cp);

  u.n[bottom] = dirichlet(0);
  u.t[bottom] = dirichlet(usb);


/**
The cylinder is no-slip. 
*/

u.n[embed] = dirichlet(0.);
u.t[embed] = dirichlet(0.);


event init (t = 0) {
  
  solid (cs, fs, sq(x) + sq(y) - sq(0.125/2.));
 
  /**
  We set the initial velocity field. */
  foreach()
    u.x[] = cs[] ? 0. : 0.;
}

/** ASCII output file
*/
event profiles (t += 1.)
{
    FILE * fp = fopen("zprof", "w");
    foreach()
    shear[] = (u.x[0,1] - u.x[0,-1])/(2.*Delta);
    boundary ({shear});
    for (double y = -0.5; y < 0.5; y += 1./128)
        fprintf (fp, "%g %g %g %g %g\n", y, interpolate (u.x, 0, y), interpolate (u.x, 0.3, y), interpolate (shear, 0, y), interpolate (shear, 0.3, y));
    fclose (fp);
}


/**
We adapt according to the error on the embedded geometry
*/
event adapt (i++) {
  adapt_wavelet ({cs}, (double[]){1e-4}, maxlevel, 5);
}

/**
Logging the drag and lift forces
*/

event logfile (i += 10)
{
  coord Fp, Fmu;
  embed_force (p, u, mu, &Fp, &Fmu);
  fprintf (stderr, "%d %g %g %.6f %.6f %.6f %.6f\n",
	   i, t, dt, Fp.x, Fp.y, Fmu.x, Fmu.y);
}


/**
Output files for contour maps
*/
event contours (t+=1., t<=tend){
  char s[80];
  sprintf (s, "map_pressure_shear_t%g",t);
  FILE *  fp = fopen(s, "w");
  output_field ({p,shear,In,muI}, fp, n=512);
  fclose(fp);
}


/** 
Details of the cylinder
*/
void cpout (FILE * fp)
{
  foreach (serial)
    if (cs[] > 0. && cs[] < 1.) {
      coord b, n;
      double area = embed_geometry (point, &b, &n);
      x += b.x*Delta, y += b.y*Delta;

      fprintf (fp, "%g %g %g %g %g %g\n",
	       x, // 1
	       y, // 2
	       atan2(y, x), // 3
	       embed_interpolate (point, p, b), // 4
	       area*Delta, // 5
	       embed_vorticity (point, u, b, n)); // 6
    }
}

event snapshot (i += 10) {
  FILE * fp = fopen ("cp", "w");
  cpout (fp);
  fclose (fp);
}

/**

## Results
~~~gnuplot Velocity profiles
set xlabel "u/U" font "Times-Roman,18"
set ylabel "z/H" font "Times-Roman,18"
set key right bottom
p[][-.5:.5] 'zprof' u 2:1 w l lw 2 t 'particle','' u 3:1 w l lw 2 t 'bulk'

~~~

~~~gnuplot Map for Inertial number
set pm3d map; 
set palette rgbformulae 22,13,-31;
set size ratio -1
set xlabel 'x' font "Times-Roman,18"
set ylabel 'y' font "Times-Roman,18"
set cblabel 'I' font "Times-Roman,18"
splot[][] 'map_pressure_shear_t10' u 1:2:5,'cp' u 1:2:(0.0) w p pt 7 ps .2 lc rgb "black" notitle
reset

~~~

~~~gnuplot Drag, lift forces and its components
set xlabel 't' font "Times-Roman,18"
set ylabel 'F' font "Times-Roman,18"
set key right bottom
set grid
p 'log' u 2:($4+$6) w l lw 4 lc 1 lt 1 t 'drag(Total)','' u 2:($5+$7) w l lw 4 lc 2 lt 2 t 'lift(Total)','' u 2:($4) w l lw 2 lc 1 lt 1 t 'drag(p)','' u 2:5 w l lw 2 lc 2 lt 2 t 'lift(p)','' u 2:6 w p lw 2 lc 1 pt 2 t 'drag(vis.)','' u 2:7 w p lw 2 lc 2 pt 2 t 'lift(vis.)'
reset

~~~

*/
