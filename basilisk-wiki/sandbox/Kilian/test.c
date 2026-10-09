/**
# Two-dimensional dry granular collapse with selectable basal friction

## Introduction

<p align="justify">
This code is used in my first-year Master's internship work on the numerical
modelling of granular landslides. It reproduces the collapse of a finite-width
dry granular pile released on an inclined plane followed by a horizontal flat
bed. It extends the basic 1D dry configuration studied by Poulain *et al.*
(2023) to a 2D horizontal domain in order to reveal the lateral spreading of
the pile during its collapse.
</p>

<p align="justify">
The flow is computed numerically with the hydrostatic `layered/hydro.h` solver
using a single layer. The unknowns are the vertical thickness $h(x,y,t)$ and
the depth-averaged horizontal velocity $\mathbf{u}=(u_x,u_y)$. The friction
source term is selected between a constant Coulomb coefficient (`muc`) and the
variable Pouliquen--Forterre law (`muv`), which is the default choice.
</p>

## Geometry sketch

<p align="justify">
Figure 1 shows a longitudinal sketch of the granular pile on the slope,
together with the notation used below. Since the code is 2D, the triangular
pile extends over a finite transverse width $W$. Although the sketch shows a
polynomial transition, the code uses a sharp break at $x=x_{\rm break}$.
</p>

<p align="center">
  <img src="https://raw.githubusercontent.com/kilian-gthub/basilisk-ressources/main/ressources/pente_droite.png"
       width="50%">
</p>

<p align="center">
*Figure 1 - Longitudinal sketch of the initial granular pile and the inclined
bed followed by the horizontal flat bed.*
</p>

## Solved equations

<p align="justify">
With $z_b(x,y)$ the bed elevation and $G$ the gravity parameter, the
single-layer hydrostatic model can be written as
</p>

$$
\frac{\partial h}{\partial t}
+\boldsymbol{\nabla}\boldsymbol{\cdot}(h\mathbf{u})=0,
$$

$$
\frac{\partial(h\mathbf{u})}{\partial t}
+\boldsymbol{\nabla}\boldsymbol{\cdot}
\left(h\mathbf{u}\otimes\mathbf{u}+\frac{Gh^2}{2}\mathbf{I}\right)
=-Gh\boldsymbol{\nabla}z_b
-Gh\mu\frac{\mathbf{u}}{|\mathbf{u}|}.
$$

<p align="justify">
The different terms in the equation are handled automatically by Basilisk.
Only the last term is applied explicitly by the `friction` event, without
changing the direction of the velocity vector.
</p>

# Code

The fixed grid used here contains $N\times N$ cells. The hydrostatic multilayer
solver is used with a single layer (`nl = 1`).
*/

#include "grid/multigrid.h"
#include "layered/hydro.h"

#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/**
## Parameters and friction choice

The code allows us to choose between variable and constant friction. This
choice can be made without editing the source by passing `muv` or `muc` as the
second command-line argument, after the grid resolution.

The constant value is $\mu_c=\tan(34^\circ)$. The variable law uses the three
friction angles and the same values of $\beta$ and $\xi$ as the
one-dimensional dry-collapse version studied by Poulain *et al.*
*/

bool use_constant_mu = false;

double tmax;
double pile_length, pile_height, pile_width;
double slope, xbreak;

double mu1, mu2, mu3;
double mu_constant;
double beta_pf, xi_pf;
double grain_length;

// The Gnuplot pipe remains open between outputs and is closed at tmax.
FILE * animation_fp = NULL;
const char * animation_filename = "dry_landslide_2D_3D.gif";

/**
## Bed and initial pile

The bed is uniform in the transverse direction and is defined by

$$
z_b(x)=
\begin{cases}
-x\tan\theta, & x<x_{\rm break},\\
-x_{\rm break}\tan\theta, & x\geq x_{\rm break}.
\end{cases}
$$

The initial thickness is triangular in $x$ and uniform over the finite width
$W$:

$$
h(x,y,0)=
\begin{cases}
H\dfrac{x+L_{\rm pile}}{L_{\rm pile}},
&-L_{\rm pile}\leq x\leq0,
\quad |y|<W/2,\\
0,&\text{elsewhere}.
\end{cases}
$$
*/

double bed_elevation (double xp)
{
  if (xp < xbreak)
    return -slope*xp;

  return -slope*xbreak;
}

double initial_pile (double xp, double yp)
{
  if (xp < -pile_length || xp > 0. || fabs(yp) >= pile_width/2.)
    return 0.;

  return pile_height*(xp + pile_length)/pile_length;
}

/**
## Pouliquen--Forterre friction law

The granular Froude number and the magnitude of the driving surface slope are

$$
Fr=\frac{|\mathbf{u}|}{\sqrt{Gh\cos\theta_{\rm loc}}},
\qquad
S=\left|\boldsymbol{\nabla}(z_b+h)\right|,
\qquad
\cos\theta_{\rm loc}=
\frac{1}{\sqrt{1+|\boldsymbol{\nabla}z_b|^2}}.
$$

The complete three-regime law implemented below is

$$
\mu(h,Fr,S)=
\begin{cases}
\displaystyle
\min\!\left(\mu_3+\frac{\mu_2-\mu_1}{1+h/L},S\right),
&Fr=0,\\[1.1em]
\displaystyle
\mu_3+\frac{\mu_2-\mu_1}{1+h/L}
+\left(\frac{Fr}{\beta}\right)^\xi(\mu_1-\mu_3),
&0<Fr<\beta,\\[1.1em]
\displaystyle
\mu_1+\frac{\mu_2-\mu_1}{1+(h/L)\beta/Fr},
&Fr\geq\beta,
\end{cases}
$$

where $\mu_i=\tan\delta_i$.
*/

double pouliquen_forterre (double Fr, double h_over_L,
                           double driving_slope)
{
  if (Fr <= 0.)
    return min(mu3 + (mu2 - mu1)/(1. + h_over_L), driving_slope);

  if (Fr < beta_pf)
    return mu3 + (mu2 - mu1)/(1. + h_over_L)
      + pow(Fr/beta_pf, xi_pf)*(mu1 - mu3);

  return mu1 + (mu2 - mu1)/(1. + h_over_L*beta_pf/Fr);
}

/**
## Main program

A $128^2$-cell grid is used by default, but the resolution can be changed by
passing $N$ as the first command-line argument.
*/

int main (int argc, char * argv[])
{
  slope = tan(35.*pi/180.);

  pile_length = 0.131;
  pile_height = slope*pile_length;
  pile_width = 0.30;
  xbreak = 0.25;

  mu1 = tan(22.1*pi/180.);
  mu2 = tan(31.8*pi/180.);
  mu3 = tan(23.3*pi/180.);
  mu_constant = tan(34.*pi/180.);

  beta_pf = 0.136;
  xi_pf = 1e-3;
  grain_length = 1e-3;

  L0 = 1.10;
  origin (-0.16, -L0/2.);
  N = argc > 1 ? atoi(argv[1]) : 128;

  if (argc > 2) {
    if (!strcmp(argv[2], "muc"))
      use_constant_mu = true;
    else if (!strcmp(argv[2], "muv"))
      use_constant_mu = false;
    else {
      fprintf (stderr, "ERROR: friction must be 'muv' or 'muc'.\n");
      return 1;
    }
  }

  nl = 1;
  G = 1.;
  DT = 1e-4;
  tmax = 8.;

  run();
}

/**
## Initial condition

The bed and the pile are assigned at the centre of each $(x,y)$ cell. Both
velocity components are initially zero.
*/

event init (i = 0)
{
  foreach() {
    zb[] = bed_elevation(x);
    h[] = initial_pile(x, y);
    u.x[] = 0.;
    u.y[] = 0.;
  }
}

/**
## Basal friction in both horizontal directions

Centered differences provide both components of the bed and thickness
gradients. For `muv`, these quantities define $Fr$, the local bed angle and the
static driving slope. For `muc`, the update is applied with $\mu=\mu_c$. The
current velocity magnitude is used at each iteration so that the source term
cannot reverse the velocity.
*/

event friction (i++)
{
  foreach() {
    if (h[] > dry) {
      double dhdx = (h[1,0] - h[-1,0])/(2.*Delta);
      double dhdy = (h[0,1] - h[0,-1])/(2.*Delta);
      double dzbdx = (zb[1,0] - zb[-1,0])/(2.*Delta);
      double dzbdy = (zb[0,1] - zb[0,-1])/(2.*Delta);

      double driving_x = -(dzbdx + dhdx);
      double driving_y = -(dzbdy + dhdy);
      double driving_slope = sqrt(sq(driving_x) + sq(driving_y));

      double velocity = sqrt(sq(u.x[]) + sq(u.y[]));
      double bottom_slope = sqrt(sq(dzbdx) + sq(dzbdy));
      double local_cosine = 1./sqrt(1. + sq(bottom_slope));
      double Fr = velocity > 1e-12 ?
        velocity/sqrt(G*h[]*local_cosine) : 0.;

      double mu_local = use_constant_mu ? mu_constant :
        pouliquen_forterre(Fr, h[]/grain_length, driving_slope);
      double du = G*mu_local*dt;

      if (velocity <= du) {
        u.x[] = 0.;
        u.y[] = 0.;
      }
      else {
        double factor = 1. - du/velocity;
        u.x[] *= factor;
        u.y[] *= factor;
      }
    }
  }
}

/**
## Three-dimensional animation

*/

event animatedplot (t = 0.; t <= tmax; t += 0.05)
{
  if (!animation_fp) {
    animation_fp = popen ("gnuplot 2> dry_landslide_2D_gnuplot.log", "w");
    if (!animation_fp) {
      fprintf (stderr, "WARNING: unable to start Gnuplot.\n");
      return 0;
    }

    fprintf (animation_fp,
             "set encoding utf8\n"
             "set term gif animate delay 4 size 1600,1200 enhanced "
             "font 'Liberation Sans,18'\n"
             "set output '%s'\n"
             "set view 62,150,1.15,1.05\n"
             "unset hidden3d\n"
             "set pm3d scansautomatic ftriangles border retrace\n"
             "unset colorbox\n"
             "set xyplane at -0.40\n"
             "set xrange [%g:%g]\n"
             "set yrange [%g:%g]\n"
             "set zrange [-0.40:0.14]\n"
             "set xlabel 'x' font 'Liberation Sans,22'\n"
             "set ylabel 'y' font 'Liberation Sans,22'\n"
             "set zlabel 'z' font 'Liberation Sans,22'\n"
             "unset key\n"
             "set object 101 rect from screen 0.12,0.012 "
             "to screen 0.145,0.028 front fc rgb '#C7AD9A' "
             "fs solid 1.0 border lc rgb '#6A6058'\n"
             "set label 101 'topography' at screen 0.152,0.019 left "
             "front font 'Liberation Sans,18'\n"
             "set object 102 rect from screen 0.41,0.012 "
             "to screen 0.435,0.028 front fc rgb '#FFF1B8' "
             "fs solid 1.0 border lc rgb '#A89540'\n"
             "set label 102 'granular material' at screen 0.442,0.019 left "
             "front font 'Liberation Sans,18'\n"
             "set object 103 rect from screen 0.72,0.018 "
             "to screen 0.75,0.022 front fc rgb '#6A6058' "
             "fs solid 1.0 noborder\n"
             "set label 103 'surface mesh' at screen 0.757,0.019 left "
             "front font 'Liberation Sans,18'\n",
             animation_filename, X0, X0 + L0, Y0, Y0 + L0);
  }

  fprintf (animation_fp,
           "set title 'Dry 2D collapse - %s - t = %.2f' "
           "font 'Liberation Sans,26'\n"
           "splot '-' using 1:2:3 with pm3d fillcolor rgb '#C7AD9A' "
           "notitle, "
           "'-' using 1:2:3 with pm3d fillcolor rgb '#C7AD9A' notitle, "
           "'-' using 1:2:3 with pm3d fillcolor rgb '#C7AD9A' notitle, "
           "'-' using 1:2:3 with pm3d fillcolor rgb '#C7AD9A' notitle, "
           "'-' using 1:2:3 with pm3d fillcolor rgb '#C7AD9A' notitle, "
           "'-' using 1:2:3 with lines lc rgb '#6A6058' lw 1.0 notitle, "
           "'-' using 1:2:3 with lines lc rgb '#6A6058' lw 1.0 notitle, "
           "'-' using 1:2:3 with pm3d fillcolor rgb '#FFF1B8' notitle, "
           "'-' using 1:2:3 with lines lc rgb '#6A6058' lw 1.0 notitle, "
           "'-' using 1:2:3 with lines lc rgb '#6A6058' lw 1.0 notitle\n",
           use_constant_mu ? "constant mu" : "variable mu (PF)", t);

  int np = 128;
  int mesh_stride = 2;
  double eps = 1e-9*L0;
  double xmin = X0 + eps, xmax = X0 + L0 - eps;
  double ymin = Y0 + eps, ymax = Y0 + L0 - eps;
  double ds = (xmax - xmin)/np;
  double zbase = -0.40;
  double mesh_offset = 2e-4;

  // Filled upper surface of the ground block.
  for (int j = 0; j <= np; j++) {
    double yp = ymin + j*ds;
    for (int i = 0; i <= np; i++) {
      double xp = xmin + i*ds;
      fprintf (animation_fp, "%g %g %g\n",
               xp, yp, interpolate (zb, xp, yp));
    }
    fputc ('\n', animation_fp);
  }
  fprintf (animation_fp, "e\n");

  // Four vertical sides close the ground block down to zbase.
  for (int i = 0; i <= np; i++)
    fprintf (animation_fp, "%g %g %g\n", xmin + i*ds, ymin, zbase);
  fputc ('\n', animation_fp);
  for (int i = 0; i <= np; i++) {
    double xp = xmin + i*ds;
    fprintf (animation_fp, "%g %g %g\n",
             xp, ymin, interpolate (zb, xp, ymin));
  }
  fprintf (animation_fp, "\ne\n");

  for (int i = 0; i <= np; i++)
    fprintf (animation_fp, "%g %g %g\n", xmin + i*ds, ymax, zbase);
  fputc ('\n', animation_fp);
  for (int i = 0; i <= np; i++) {
    double xp = xmin + i*ds;
    fprintf (animation_fp, "%g %g %g\n",
             xp, ymax, interpolate (zb, xp, ymax));
  }
  fprintf (animation_fp, "\ne\n");

  for (int j = 0; j <= np; j++)
    fprintf (animation_fp, "%g %g %g\n", xmin, ymin + j*ds, zbase);
  fputc ('\n', animation_fp);
  for (int j = 0; j <= np; j++) {
    double yp = ymin + j*ds;
    fprintf (animation_fp, "%g %g %g\n",
             xmin, yp, interpolate (zb, xmin, yp));
  }
  fprintf (animation_fp, "\ne\n");

  for (int j = 0; j <= np; j++)
    fprintf (animation_fp, "%g %g %g\n", xmax, ymin + j*ds, zbase);
  fputc ('\n', animation_fp);
  for (int j = 0; j <= np; j++) {
    double yp = ymin + j*ds;
    fprintf (animation_fp, "%g %g %g\n",
             xmax, yp, interpolate (zb, xmax, yp));
  }
  fprintf (animation_fp, "\ne\n");

  // Ground mesh in both directions, slightly raised to remain visible with
  // the Gnuplot version used by the Basilisk server.
  for (int j = 0; j <= np; j += mesh_stride) {
    double yp = ymin + j*ds;
    for (int i = 0; i <= np; i += mesh_stride) {
      double xp = xmin + i*ds;
      fprintf (animation_fp, "%g %g %g\n", xp, yp,
               interpolate (zb, xp, yp) + mesh_offset);
    }
    fputc ('\n', animation_fp);
  }
  fprintf (animation_fp, "e\n");

  for (int i = 0; i <= np; i += mesh_stride) {
    double xp = xmin + i*ds;
    for (int j = 0; j <= np; j += mesh_stride) {
      double yp = ymin + j*ds;
      fprintf (animation_fp, "%g %g %g\n", xp, yp,
               interpolate (zb, xp, yp) + mesh_offset);
    }
    fputc ('\n', animation_fp);
  }
  fprintf (animation_fp, "e\n");

  // Filled granular free surface; dry points are hidden with NaN.
  for (int j = 0; j <= np; j++) {
    double yp = ymin + j*ds;
    for (int i = 0; i <= np; i++) {
      double xp = xmin + i*ds;
      double hp = interpolate (h, xp, yp);
      if (hp > 1e-4)
        fprintf (animation_fp, "%g %g %g\n",
                 xp, yp, interpolate (zb, xp, yp) + hp);
      else
        fprintf (animation_fp, "%g %g NaN\n", xp, yp);
    }
    fputc ('\n', animation_fp);
  }
  fprintf (animation_fp, "e\n");

  // Granular mesh in both directions, also raised by a negligible display
  // offset so that the filled surface cannot hide it.
  for (int j = 0; j <= np; j += mesh_stride) {
    double yp = ymin + j*ds;
    for (int i = 0; i <= np; i += mesh_stride) {
      double xp = xmin + i*ds;
      double hp = interpolate (h, xp, yp);
      if (hp > 1e-4)
        fprintf (animation_fp, "%g %g %g\n",
                 xp, yp, interpolate (zb, xp, yp) + hp + mesh_offset);
      else
        fprintf (animation_fp, "%g %g NaN\n", xp, yp);
    }
    fputc ('\n', animation_fp);
  }
  fprintf (animation_fp, "e\n");

  for (int i = 0; i <= np; i += mesh_stride) {
    double xp = xmin + i*ds;
    for (int j = 0; j <= np; j += mesh_stride) {
      double yp = ymin + j*ds;
      double hp = interpolate (h, xp, yp);
      if (hp > 1e-4)
        fprintf (animation_fp, "%g %g %g\n",
                 xp, yp, interpolate (zb, xp, yp) + hp + mesh_offset);
      else
        fprintf (animation_fp, "%g %g NaN\n", xp, yp);
    }
    fputc ('\n', animation_fp);
  }
  fprintf (animation_fp, "e\n");
  fflush (animation_fp);
}

/**
## Clean termination

Closing the Gnuplot pipe finalises all the frames of the GIF.
*/

event stop (t = tmax)
{
  if (animation_fp) {
    fprintf (animation_fp, "unset output\n");
    pclose (animation_fp);
    animation_fp = NULL;
  }

  return 1;
}

/**
# Main result

## Three-dimensional collapse animation

<p align="justify">
The animation shows the release of the initial pile, its acceleration down the
inclined plane, and its longitudinal and transverse spreading. The 3D view
highlights the collapse dynamics and the final deposit footprint on the flat
bed.
</p>

<p align="center">
  <img src="test/dry_landslide_2D_3D.gif"
       alt="Three-dimensional animation of the two-dimensional dry granular collapse"
       width="70%">
</p>

<p align="center">
*Figure 2 - Three-dimensional animation of the two-dimensional dry granular collapse.*
</p>

# References

* POLGE POULICHET, K. (2026). *Numerical modelling of landslides under different mechanical conditions*. First-year Master's internship report, Sorbonne Université, Institut Jean le Rond d'Alembert. [Report PDF](https://raw.githubusercontent.com/kilian-gthub/basilisk-ressources/main/ressources/Rapport_POLGE_POULICHET_Kilian.pdf)
* [P. Poulain *et al.*, *Performance and limits of a shallow-water model for landslide-generated tsunamis: from laboratory experiments to simulations of flank collapses at Montagne Pelée (Martinique)*, Geophysical Journal International, 233, 796--825, 2023.](https://www.ipgp.fr/~mangeney/poulain-etal_gji-2023.pdf)
* [O. Pouliquen & Y. Forterre, *Friction law for dense granular flows: application to the motion of a mass down a rough inclined plane*, Journal of Fluid Mechanics, 453, 133--151, 2002.](https://yoelforterre.wordpress.com/wp-content/uploads/2016/09/jfmmoire02.pdf)
* [S. Popinet, *A vertically-Lagrangian, non-hydrostatic, multilayered model for multiscale free-surface flows*, Journal of Computational Physics, 418, 109609, 2020.](https://doi.org/10.1016/j.jcp.2020.109609)
*/
