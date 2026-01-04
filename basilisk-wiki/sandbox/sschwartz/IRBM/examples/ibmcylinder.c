#include "navier-stokes/centered.h"
#include "../ibm.h"
#include "view.h"

#define L0 15.
#define D 0.5

#define LEVEL 10
#define MIN_LEVEL 5

int Re;
double U0 =  1.;               // inlet velocity
double t_end = 50;
coord ci = {L0/4., L0/2.};     // initial coordinates of cylinder

face vector muv[];

/**
Boundary conditions: left=inlet, right=outlet.*/

u.n[left] = dirichlet (U0);
u.t[left] = dirichlet (0);
p[left]   = neumann (0);
pf[left]  = neumann (0);

u.n[right] = neumann (0);
u.t[right] = neumann (0);
p[right]   = dirichlet (0);
pf[right]  = dirichlet (0);

u.n[top] = neumann (0);
p[top] = neumann (0);
pf[top] = neumann (0);

u.n[bottom] = neumann (0);
p[bottom] = neumann (0);
pf[bottom] = neumann (0);

int main() 
{
  size(L0);
  init_grid (1 << (LEVEL - 2));
  mu = muv;
  TOLERANCE = 1.e-6; 
  CFL = 0.5;
  spreading_func = 4; // spreading function size (4 cells)
  Ni = 10;             // # of MDF iterations

  Re = 100;
  run();
}

event init (t = 0) {
  vc.x = 0, vc.y = 0;
  solid (ibm, ibmf, sq(x - ci.x - vc.x) + sq(y - ci.y - vc.y) - sq(D/2));
  refine (ibm[] < 1 && ibm[] > 0 && level < LEVEL);
  solid (ibm, ibmf, sq(x - ci.x - vc.x) + sq(y - ci.y - vc.y) - sq(D/2));

  foreach()
    u.x[] = U0 * ibm[];
}

event properties (i++) 
{
  foreach_face()
    muv.x[] = fm.x[]*(U0)*(D)/(Re);
   boundary ((scalar *) {muv});
}

SurfacePoints surf;
scalar cf[];

event logfile (i++; t <= t_end) 
{
  coord F = ibm_force();
  double CD = (F.x)/(0.5*sq(U0)*(D));
  double CL = (F.y)/(0.5*sq(U0)*(D));

  coord Fp, Fmu;
  ibm_force_int (p, u, mu, &Fp, &Fmu);

  double cftotal = ibm_skin_friction (u, mu, cf);

  surf = get_surfacepoints (u);
  if (i == 0)
    fprintf (stderr, "0:Ni | 1:i | 2:t | 3:Re | 4:mgp.i | 5:mgp.nr | 6:mgu.i | 7:mgu.nr | 8:CD | 9:CL | 10:Fpx | 11:Fpy | 12:Fmux | 13:Fmuy | 14:ul2 | 15:umin | 16:umax | 17:vl2 | 18:vmin | 19:vmax | 20: nump | 21:cf\n");
  fprintf (stderr, "%d %d %g %d %d %d %d %d %g %g %g %g %g %g %g %g %g %g %g %g %d %g\n",
          Ni, i, t, Re, mgp.i, mgp.nrelax, mgu.i, mgu.nrelax, CD, CL, Fp.x, Fp.y, Fmu.x, Fmu.y, 
          surf.ul2, surf.umin, surf.umax, surf.vl2, surf.vmin, surf.vmax, surf.total, cftotal);
}


event flux_profile (t = end)
{
    scalar omega[];
    vorticity(u, omega);

    char name[80];
    sprintf (name, "%d-uerror-%d.png", Re, Ni);
    view (fov = 1, tx = -0.25, ty = -0.5,
          width = 3000, height = 1500); 
    clear();
    draw_vof ("ibm", "ibmf", lw = 5, lc = {0,0,0});
    double umax = fabs(surf.umax - vc.x);
    squares ("fabs(verror)", min = -umax, max = umax, map = cool_warm, cbar = true);
    save (name);

    sprintf (name, "%d-verror-%d.png", Re, Ni);
    clear();
    draw_vof ("ibm", "ibmf", lw = 5, lc = {0,0,0});
    double vmax = fabs(surf.vmax - vc.y);
    squares ("fabs(uerror)", min = -vmax, max = vmax, map = cool_warm, cbar = true);
    save (name);


    sprintf (name, "%d-fields-%d", Re, Ni);
    FILE * fp1 = fopen (name, "w");
    output_field ({all}, fp = fp1, n = 2*pow(2,LEVEL), linear = true);
    fclose(fp1);
}

/**
We output the velocity and pressure profiles at certain "cuts" across the domain. */

event profile (t = t_end) {
  double delta = L0/(pow(2,LEVEL));
  char name[80];

  sprintf (name, "%d-vprofx1-%d", Re, Ni); // x = center of cylinder
  FILE * fv1 = fopen(name, "w");
  for(double i = 0; i <= L0; i += delta) {
    foreach_point (ci.x, i) {
      fprintf (fv1, "%g %g %g %g %g\n", x, y, u.x[], u.y[], p[]);
    }
  }
  fflush (fv1);
  fclose (fv1);

  sprintf (name, "%d-vprofx2-%d", Re, Ni); // x = outlet
  FILE * fv2 = fopen(name, "w");
  for(double i = 0; i <= L0; i += delta) {
    foreach_point (L0-delta, i) {
      fprintf (fv2, "%g %g %g %g %g\n", x, y, u.x[], u.y[], p[]); 
    }
  }
  fflush (fv2);
  fclose (fv2);

  sprintf (name, "%d-vprofy1-%d", Re, Ni); // y = L0/2
  FILE * fv3 = fopen(name, "w");
  for(double i = 0; i <= L0; i += delta) {
    foreach_point (i, ci.y) {
      fprintf (fv3, "%g %g %g %g %g\n", x, y, u.x[], u.y[], p[]);
    }
  }
  fflush (fv3);
  fclose (fv3);
}

/**
Output flow field images (if told to). */
#if MOVIE
event snapshot (t = t_end) {
  scalar omega[];
  vorticity (u, omega);

  char name[80];
  view (fov = 2, tx = -0.26, ty = -0.5,
        width = 3000, height = 1500); 
  sprintf (name, "%d-pressure-%d.png", Re, Ni);
  clear();
  draw_vof ("ibm", "ibmf", lw = 5, lc = {0,0,0});
  squares ("p", min = -0.5, max = 0.5, map = cool_warm);
  save (name);

  sprintf (name, "%d-vort-%d.png", Re, Ni);
  clear();
  squares ("omega", min = -3, max = 3, map = cool_warm);
  draw_vof ("ibm", "ibmf", lw = 5, lc = {0,0,0});
  save (name);

  sprintf (name, "%d-pressureiso-%d.png", Re, Ni);
  clear();
  draw_vof ("ibm", "ibmf", lw = 5, lc = {0,0,0});
  isoline ("p", n = 20, min = -0.5, max = 0.5, lc = {0,0,0});
  squares ("p", min = -0.5, max = 0.5, map = cool_warm);
  save (name);

  sprintf (name, "%d-vortiso-%d.png", Re, Ni);
  clear();
  isoline ("omega", n = 15, min = -3, max = 3, lc = {0,0,0});
  squares ("omega", min = -3, max = 3, map = cool_warm);
  draw_vof ("ibm", "ibmf", lw = 5,  lc = {0,0,0});
  save (name);
}
#endif

event surface_values (t = end)
{
    char name[80];
    sprintf (name, "%d-surfval-%d", Re, Ni);
    FILE * fp = fopen (name, "w");
    fprintf(fp, "1:x | 2:y | 3:theta | 4:area | 5:cp | 6:cf | 7:unf | 8:utf | 9:uflux\n");

    foreach() {
        if (ibm[] > 0 && ibm[] < 1) {
            coord n = facet_normal (point, ibm, ibmf), midPoint, b;
            double alpha = plane_alpha (ibm[], n);
            double area = plane_area_center (n, alpha, &b) * pow(Delta, dimension-1);

            coord cellCenter = {x,y,z};
            foreach_dimension() {
                midPoint.x = cellCenter.x + b.x*Delta;
                n.x *= -1;
            }
            normalize(&n);

            // calculate pressure coefficient
            double boundaryPressure = extrapolate_scalar (point, ibm, midPoint, n, p);
            double cp = boundaryPressure / (0.5 * sq(U0));

            // calculate skin friction coefficient
            double cfi = cf[] / (0.5 * sq(U0));
            double theta = atan2 ((y - ci.y), (x - ci.x)) * (180/M_PI);

            fprintf(fp, "%g %g %g %g %g %g %g %g %g\n", 
                         x, y, theta, area, cp, cfi, unf[], utf[], uflux[]);
        }
    }

    fclose(fp);
}

/**
We adapt based on the smeared boundary force and velocity. The solid surface
should be at the maximum level of refinement (for now).*/

event adapt (i++) {
  adapt_wavelet ({forceTotal, u}, (double[]){3e-3, 3e-3, 3e-3, 3e-3},
		 maxlevel = LEVEL, minlevel = MIN_LEVEL);
}

