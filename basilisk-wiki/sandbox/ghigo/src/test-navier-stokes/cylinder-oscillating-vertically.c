/**
# Flow past a vertically oscillating cylinder for $Re=185$ and $f=0.8f_0$

This test case is inspired from the experimental work of [Gu et al.,
1994](#gu1994) and the numerical simulations of [Guilmineau et al.,
2002](#guilmineau2002). The authors studied the vortex switching
induced by a harmonic forcing of a cylinder oscillating tranverselly
to the freestream direction. Depending on the ratio $f/f_0$, the
vortex shedding frequency should synchronize with the oscillation
frequency after a few periods.

This test case was also reproduced numerically by [Uhlman,
2005](#uhlman2005), [Kim et al., 2006](#kim2006), [Yang et al.,
2009](#yang2009), [Young et al., 2009](#young2009), [Schneiders et
al., 2013](#schneiders2013), [Meinke et al., 2013](#meinke2013) and
[Xin et al., 2018](#xin2018), using the following numerical
parameters:

                        $f_0$ (for $U=1,\, D=1$) $f/f_0$
----------------------- ------------------------ --------------------------- -----------------------------------------------------------------------
Guilmineau et al., 2002 0.195                    0.8, 0.9, 1, 1.1, 1.12, 1.2 $\frac{U_{\infty}\Delta t}{D} = 2e^{-3}$, $\frac{D}{\Delta} = ?$
Uhlman, 2005            ?                        0.8                         $\frac{U_{\infty}\Delta t}{D} = 1e^{-2}$, $\frac{D}{\Delta} = 38$
Kim et al., 2006        ?                        0.8, 0.9, 1, 1.1, 1.12, 1.2 $\frac{U_{\infty}\Delta t}{D} = ?$, $\frac{D}{\Delta} = 30$
Yang et al, 2009        0.156                    0.8                         $\frac{U_{\infty}\Delta t}{D} = 2e^{-3}$, $\frac{D}{\Delta} = 50$
Young et al, 2009       0.197/0.199              0.8, 0.9, 1, 1.1, 1.12, 1.2 $\frac{U_{\infty}\Delta t}{D} = ?$, $\frac{D}{\Delta} = 15-20$
Schneiders et al., 2013 0.195                    0.8                         $\frac{U_{\infty}\Delta t}{D} = ?\, (CFL=0.5)$, $\frac{D}{\Delta} = 50$
Meinke et al,. 2013     0.195                    0.8, 0.9, 1, 1.1, 1.12, 1.2 $\frac{U_{\infty}\Delta t}{D} = ?\, (CFL=0.5)$, $\frac{D}{\Delta} = 33$
Xin et al., 2018        0.195                    0.8, 0.9, 1, 1.1, 1.12, 1.2 $\frac{U_{\infty}\Delta t}{D} = 2e^{-2}$, $\frac{D}{\Delta} = 16$

We solve here the Navier-Stokes equations and add the cylinder using
an [embedded boundary](/src/embed.h) for $Re=185$ and a forcing
frequency $f_e=0.8 f_0$, where $f_0=0.195$ is the natural shedding
frequency of the cylinder. Note that this is a relatively challenging
test case as the motion of the cylinder is perpendicular to the
direction of the flow, which results in difficulties in imposing the
values in emerged cells. */

#include "grid/quadtree.h"
#include "../myembed.h"
#include "../mycentered.h"
#include "../myembed-moving.h"
#include "../myperfs.h"
#include "view.h"

/**
## Reference solution */

#define d    (1.)
#define Re   (185.) // Particle Reynolds number
#define p_f0 (0.195) // Natural vortex shedding frequency
#define p_f  (0.8*(p_f0)) // Oscillating frequency
#define uref (max (1., 0.4*M_PI*(p_f))) // Reference velocity, uref
#define tref (min ((d)/(uref), 1./(p_f))) // Reference time, tref

/**
We define the cylinder's imposed motion. */

# define p_acceleration(t,f) (0.8*sq(M_PI*(f))*cos(2.*M_PI*(f)*(t))) // Particle acceleration
# define p_velocity(t,f)     (0.4*M_PI*(f)*sin(2.*M_PI*(f)*(t))) // Particle velocity
# define p_displacement(t,f) (-0.2*cos(2.*M_PI*(f)*(t))) // Particle displacement

/**
We also define the shape of the domain. */

#define cylinder(x,y) (sq (x) + sq (y) - sq ((d)/2.))

void p_shape (scalar c, face vector f, coord p)
{
  vertex scalar phi[];
  foreach_vertex()
    phi[] = (cylinder ((x - p.x), (y - p.y)));
  boundary ({phi});
  fractions (phi, c, f);
  fractions_cleanup (c, f,
		     smin = 1.e-14, cmin = 1.e-14);
}

/**
## Setup

We need a field for viscosity so that the embedded boundary metric can
be taken into account. */

face vector muv[];

/**
We define the mesh adaptation parameters. */

#define lmin (7)  // Min mesh refinement level (l=7 is 2pt/d)
#define lmax (12) // Max mesh refinement level (l=12 is 64pt/d)
#define cmax (1.e-2*(uref)) // Absolute refinement criteria for the velocity field

/**
We finally define a useful function that allows us to define a
rectangular block mesh refinement (BMR) around the cylinder. */

#define rectangle(x,hx,y,hy) (union (					\
				     union ((x) - (hx)/2., -(x) - (hx)/2.), \
				     union ((y) - (hy)/2., -(y) - (hy)/2.)) \
			      )

int main ()
{  
  /**
  The domain is $64\times 64$. */

  L0 = 64.;
  size (L0);
  origin (-L0/2., -L0/2.);

  /**
  We set the maximum timestep. */

  DT = 1.e-2*(tref);
 
  /**
  We set the tolerance of the Poisson solver. */

  TOLERANCE    = 1.e-6;
  TOLERANCE_MU = 1.e-6*(uref);  

  /**
  We initialize the grid. */
  
  N = 1 << (lmin);
  init_grid (N);
  
  run();
}

/**
## Boundary conditions

We use inlet boundary conditions. */

u.n[left] = dirichlet ((uref));
u.t[left] = dirichlet (0);
p[left]   = neumann (0);
pf[left]  = neumann (0);

u.n[right] = neumann (0);
u.t[right] = neumann (0);
p[right]   = dirichlet (0);
pf[right]  = dirichlet (0);

/**
We give boundary conditions for the face velocity to "potentially"
improve the convergence of the multigrid Poisson solver. */

uf.n[left]   = (uref);
uf.n[bottom] = 0;
uf.n[top]    = 0;

/**
## Properties */

event properties (i++)
{
  foreach_face()
    muv.x[] = (uref)*(d)/(Re)*fm.x[];
  boundary ((scalar *) {muv});
}

/**
## Initial conditions */

event init (i = 0)
{
  /**
  We set the viscosity field in the event *properties*. */

  mu = muv;
  
  /**
  We use "third-order" [face flux interpolation](/src/embed.h). */
    
#if ORDER2
  for (scalar s in {u, p, pf})
    s.third = false;
#else
  for (scalar s in {u, p, pf})
    s.third = true;
#endif // ORDER2

  /**
  We use a slope-limiter to reduce the errors made in small-cells. */

#if SLOPELIMITER
  for (scalar s in {u, p, pf}) {
    s.gradient = minmod2;
  }
#endif // SLOPELIMITER
  
#if TREE
  /**
  When using *TREE* and in the presence of embedded boundaries, we
  should also define the gradient of *u* at the cell center of
  cut-cells. */
#endif // TREE
  
  /**
  We initialize the embedded boundary. */

  /**
  We first define the particle's initial position. */

  p_p.y = (p_displacement (0., (p_f)));
  
#if TREE
#if BMR
  int lvl = (lmin);
  double wxmin = (L0)/2.*(d), wxmax = 3.*(d);
  double xmin = wxmin/2. - 4.*(d), xmax = wxmax/2. - 1.*(d);
  double wymin = 8.*(d), wymax = 3.*(d);
  double ymin = 0., ymax = 0.;
  while (lvl < (lmax)) {
    refine (
	    (rectangle ((x - (p_p.x + (xmin + (xmax - xmin)/((lmax) - (lmin))*((lvl + 1) - (lmin))))),
			wxmin + (wxmax - wxmin)/((lmax) - (lmin))*((lvl + 1) - (lmin)),
			(y - (p_p.y + (ymin + (ymax - ymin)/((lmax) - (lmin))*((lvl + 1) - (lmin))))),
			wymin + (wymax - wymin)/((lmax) - (lmin))*((lvl + 1) - (lmin)))) <= 0. &&
	    level < (lvl + 1)
	    );
    p_shape (cs, fs, p_p);
    lvl ++;
  }

  /**
  In the case of the moving cylinder, we don't unrefine inside the
  cylinder. */
  
  unrefine ((rectangle ((x - (p_p.x + xmin)), wxmin,
			(y - (p_p.y + ymin)), wymin)) > 0. &&
	    level > 1);
#else // AMR
  /**
  When using *TREE*, we refine the mesh around the embedded
  boundary. */
  
  astats ss;
  int ic = 0;
  do {
    ic++;
    p_shape (cs, fs, p_p);
    ss = adapt_wavelet ({cs}, (double[]) {1.e-30},
			maxlevel = (lmax), minlevel = (1));
  } while ((ss.nf || ss.nc) && ic < 100);
#endif // BMR
#endif // TREE
  
  p_shape (cs, fs, p_p);

  /**
  We initialize the particle's speed and accelerating. */
  
  p_au.y = (p_acceleration (0., (p_f)));
  p_u.y  = (p_velocity     (0., (p_f)));

  /**
  We finally initialize the velocity. */

  foreach() 
    u.x[] = cs[]*(uref);
  boundary ((scalar *) {u});
}

/**
## Embedded boundaries 

The cylinder's position is advanced to time $t + \Delta t$. */

event advection_term (i++)
{
  p_au.y = (p_acceleration (t + dt, (p_f)));
  p_u.y  = (p_velocity     (t + dt, (p_f)));
  p_p.y  = (p_displacement (t + dt, (p_f)));
}

/**
## Adaptive mesh refinement */

#if TREE && !BMR
event adapt (i++)
{
  adapt_wavelet ({cs,u}, (double[]) {1.e-2,(cmax),(cmax)},
  		 maxlevel = (lmax), minlevel = (1));

  /**
  We do not need here to reset the embedded fractions to avoid
  interpolation errors on the geometry as the is already done when
  moving the embedded boundaries. It might be necessary to do this
  however if surface forces are computed around the embedded
  boundaries. */
}
#endif // TREE && !BMR

/**
## Outputs */

event logfile (i++; t <= 20./(p_f))
{
  coord Fp, Fmu;
  embed_force (p, u, mu, &Fp, &Fmu);

  double CD = (Fp.x + Fmu.x)/(0.5*sq ((uref))*(d));
  double CL = (Fp.y + Fmu.y)/(0.5*sq ((uref))*(d));

  fprintf (stderr, "%d %g %g %d %d %d %d %d %d %g %g %g %g %g %g %g %g %g %g %g %g\n",
	   i, t/(1./(p_f)), dt/(1./(p_f)),
	   mgp.i, mgp.nrelax, mgp.minlevel,
	   mgu.i, mgu.nrelax, mgu.minlevel,
	   mgp.resb, mgp.resa,
	   mgu.resb, mgu.resa,
	   p_p.x/(d), p_p.y/(d),
	   CD, CL,
	   Fp.x, Fp.y, Fmu.x, Fmu.z);
  fflush (stderr);
}

/**
## Snapshots 

We also compute the distribution of the pressure coefficient $C_p$ and
vorticity $\omega$ at the surface of the cylinder as it reaches: (i)
the extreme upper position and (ii) the center position while moving
downwards. */

void cpout (FILE * fp)
{
  foreach(serial)
    if (cs[] > 0. && cs[] < 1.) {

      coord b, n;
      embed_geometry (point, &b, &n);
      double xe = x + b.x*Delta, ye = y + b.y*Delta;

      fprintf (fp, "%g %g %g\n",
	       M_PI + atan2(ye - p_p.y, xe - p_p.x), // 1
	       embed_interpolate (point, p, b)/(0.5*sq (uref)*(d)), // 2
	       embed_vorticity (point, u, b, n) //3
	       );
      fflush (fp);
    }
}

event snapshots (t = {0., 19.5/(p_f), 19.75/(p_f)})
{
  int phase = max(0, floorf ((t/(1./(p_f)) - 19.)*360));

  char name2[80];
  sprintf (name2, "cp-phi-%d-pid-%d", phase, pid());
  FILE * fp = fopen (name2, "w");
  cpout (fp);
  fclose (fp);
  
  scalar omega[];
  vorticity (u, omega);
  
  /**
  We first plot the entire domain. */
 
  view (fov = 20, camera = "front",
	tx = 0., ty = 0.,
	bg = {1,1,1},
	width = 800, height = 800);

  draw_vof ("cs", "fs", lw = 5);
  cells ();
  sprintf (name2, "mesh-phi-%d.png", phase);
  save (name2);
    
  draw_vof ("cs", "fs", filled = -1, lw = 5, fc = {1,1,1});
  squares ("u.x", map = cool_warm);
  sprintf (name2, "ux-phi-%d.png", phase);
  save (name2);

  draw_vof ("cs", "fs", filled = -1, lw = 5, fc = {1,1,1});
  squares ("u.y", map = cool_warm);
  sprintf (name2, "uy-phi-%d.png", phase);
  save (name2);

  draw_vof ("cs", "fs", filled = -1, lw = 5, fc = {1,1,1});
  squares ("p", map = cool_warm);
  sprintf (name2, "p-phi-%d.png", phase);
  save (name2);

  draw_vof ("cs", "fs", filled = -1, lw = 5, fc = {1,1,1});
  squares ("omega", map = cool_warm);
  sprintf (name2, "omega-phi-%d.png", phase);
  save (name2);

  /**
  We then zoom on the cylinder. */

  view (fov = 2.5, camera = "front",
	tx = -(p_p.x + 3*(d))/(L0), ty = 0., bg = {1,1,1},
	width = 800, height = 400);

  draw_vof ("cs", "fs", lw = 5);
  cells ();
  sprintf (name2, "mesh-zoom-phi-%d.png", phase);
  save (name2);

  draw_vof ("cs", "fs", filled = -1, lw = 5, fc = {1,1,1});
  squares ("u.x", map = cool_warm);
  sprintf (name2, "ux-zoom-phi-%d.png", phase);
  save (name2);

  draw_vof ("cs", "fs", filled = -1, lw = 5, fc = {1,1,1});
  squares ("u.y", map = cool_warm);
  sprintf (name2, "uy-zoom-phi-%d.png", phase);
  save (name2);

  draw_vof ("cs", "fs", filled = -1, lw = 5, fc = {1,1,1});
  squares ("p", map = cool_warm);
  sprintf (name2, "p-zoom-phi-%d.png", phase);
  save (name2);

  draw_vof ("cs", "fs", filled = -1, lw = 5, fc = {1,1,1});
  squares ("omega", map = cool_warm);
  sprintf (name2, "omega-zoom-phi-%d.png", phase);
  save (name2);
}

/**
## Results

#### Pressure and vorticity for level = 12

<p>

<img src="cylinder-oscillating/p-zoom-level-11-phi-180.png" alt="drawing" width="25%"/>
<img src="cylinder-oscillating/omega-zoom-level-11-phi-180.png" alt="drawing" width="25%"/> 
<br />
<img src="cylinder-oscillating/p-zoom-level-11-phi-270.png" alt="drawing" width="25%"/>
<img src="cylinder-oscillating/omega-zoom-level-11-phi-270.png" alt="drawing" width="25%"/>

</p>

#### Periodicity

We first plot the time evolution of the drag and lift coefficients and
verify that the flow has reached a periodic state after 9 cycles.

~~~gnuplot Time evolution of the drag coefficient $C_D$
set terminal svg font ",16"
set key top right spacing 1.1
set xlabel 't/T'
set ylabel 'C_D'
set xrange [0:20]
set yrange [0.5:2.5]
plot 'log' u 2:16 w l lw 1.5 lc rgb "black" t "Basilisk, l=12"
~~~

~~~gnuplot Time evolution of the lift coefficient $C_L$
set ylabel 'C_L'
set yrange[-1.5:2]
plot 'log' u 2:17 w l lw 1.5 lc rgb "black" t "Basilisk, l=12"
~~~

~~~gnuplot Time evolution of the timestep $\Delta t$
set ylabel "dt/T"
set yrange [1.e-4:1e0]
set logscale y
plo 'log' u 2:3 w l lw 1.5 lc rgb "black" t "Basilisk, l=12"
~~~

#### Comparison with the results of fig. 9 from [Schneiders et al., 2013](#schneiders2013).

We next plot the drag coefficent $C_D$ as a function of the vertical
displacement $y/D$ over the last period. We compare the results with
those of fig. 9 from [Schneiders et al., 2013](#schneiders2013).

~~~gnuplot Drag coefficient $C_D$ as a function of the displacement $y/D$
set xlabel 'y/D'
set ylabel 'C_D'
set xrange[-0.25:0.25]
set yrange[1.1:1.5]
unset logscale
plot '../data/Schneiders2013/Schneiders2013-fig9.csv' u 1:2 w p ps 0.7 pt 6 lc rgb "black" \
     t "fig. 9, Schneiders et al., 2013", \
     '< sort -k1,1 log | awk -f ../data/cylinder-oscillating-vertically/lastperiod.awk' w l lw 2 lc rgb "black" \
     t "Basilisk, l=12"
~~~

#### Surface pressure coefficient $C_p$ around the cylinder

We now plot the distribution of the pressure coefficent $C_p$ at the
surface of the cylinder at (i) the extreme upper position and (ii) the
center position while goind downwards. We compare the results with
those of fig. 17 from [Guilmineau et al., 2002](#guilmineau2002) and
fig. 12b from [Schneiders et al., 2013](#schneiders2013).

~~~gnuplot (i) Pressure coefficient $C_p$ at the extreme upper position
set key top left
set xlabel 'theta'
set ylabel 'C_p'
set xrange[0:360]
set yrange[-2:2]
plot '../data/Guilmineau2002/Guilmineau2002-fig17-fsf0-0p8.csv' u 1:2 w p ps 0.7 pt 6 lc rgb "black" \
     t "fig. 17, Guilmineau et al., 2002", \
     '< cat cp-phi-180-pid-* | sort -k1,1' u (-(180./3.14*$1 - 360.)):($2) w l lw 2 lc rgb "black" \
     t 'Basilisk, l=12'
~~~

~~~gnuplot (ii) Pressure coefficient $C_p$ at the center position while going downwards
plot '../data/Schneiders2013/Schneiders2013-fig12b.csv' u 1:2 w p ps 1. pt 6 lc rgb "black" \
     t "fig. 12b, Schneiders et al., 2013", \
     '< cat cp-phi-270-pid-* | sort -k1,1' u (-(180./3.14*$1 - 360.)):($2) w l lw 2 lc rgb "black" \
     t 'Basilisk, l=12'
~~~

#### Surface vorticity $\omega$ around the cylinder

We finally plot the distribution of the vorticity $\omega$ at the
surface of the cylinder at the extreme upper position. We compare the
results with those of fig. 13 from [Guilmineau et al.,
2002](#guilmineau2002).

~~~gnuplot Surface vorticity $\omega$ at the extreme upper position
set key bottom right
set ylabel 'omega'
set yrange[-40:30]
plot '../data/Guilmineau2002/Guilmineau2002-fig13-fsf0-0p8.csv' u 1:2 w p ps 1 pt 6 lc rgb "black" \
     t "fig. 13, Guilmineau et al., 2002", \
     '< cat cp-phi-180-pid-* | sort -k1,1' u (-(180./3.14*$1 - 360.)):($3) w l lw 2 lc rgb "black" \
     t 'Basilisk, l=12'
~~~

## References

~~~bib
@article{gu1994,
  title={Timing of vortex formation from an oscillating cylinder},
  author={Gu, W. and Chyu, C. and Rockwell, D.},
  journal={Physics of Fluid},
  volume={6},
  pages={3677--3682},
  year={1994},
  publisher={AIP Publishing}
}

@article{guilmineau2002,
  title={A NUMERICAL SIMULATION OF VORTEX SHEDDING FROM AN OSCILLATING CIRCULAR CYLINDER},
  author={Guilmineau, E. and Queutey, P.},
  journal={Journal of Fluids and Structures},
  volume={16},
  pages={773--794},
  year={2002},
  publisher={Elsevier}
}

@article{uhlman2005,
  title={An immersed boundary method with direct forcing for the simulation of particulate flows},
  author={Uhlman, M.},
  journal={Journal of Computational Physics},
  volume={209},
  pages={448--476},
  year={2005},
  publisher={Elsevier}
}

@article{kim2006,
  title={Immersed boundary method for flow around an arbitrarily moving body},
  author={Kim, D. and Choi, H.},
  journal={Journal of Computational Physics},
  volume={212},
  pages={662--680},
  year={2006},
  publisher={Elsevier}
}

@article{yang2009,
  title={A smoothing technique for discrete delta functions with application to immersed boundary method in moving boundary simulations},
  author={Yang, X. and Zhang, X. and Li, Z. and He, G.-W.},
  journal={Journal of Computational Physics},
  volume={228},
  pages={7821--7836},
  year={2009},
  publisher={Elsevier}
}

@article{young2009,
  title={A novel immersed boundary procedure for flow and heat simulations with moving boundary},
  author={Young, D.L. and Jan, Y.J. and Chiu, C.L.},
  journal={Computers \& Fluids},
  volume={38},
  pages={1145--1159},
  year={2009},
  publisher={Elsevier}
}

@article{schneiders2013,
  title={An accurate moving boundary formulation in cut-cell methods},
  author={Schneiders, L. and Hartmann, D. and Meinke, M. and Schroder, W.},
  journal={Journal of Computational Physics},
  volume={235},
  pages={786--809},
  year={2013},
  publisher={Elsevier}
}

@article{meinke2013,
  title={A cut-cell method for sharp moving boundaries in Cartesian grids},
  author={Meinke, M. and Schneiders, L. and Gunther, C. and Schroder, W.},
  journal={Computers \& Fluids},
  volume={85},
  pages={135--142},
  year={2013},
  publisher={Elsevier}
}

@article{xin2018,
  title={A radial basis function based ghost cell method with improved mass conservation for complex moving boundary flows},
  author={Xin, J. and Shi, F. and Jin, Q. and Lin, C.},
  journal={Computers \& Fluids},
  volume={176},
  pages={210--225},
  year={2018},
  publisher={Elsevier}
}
~~~
*/
