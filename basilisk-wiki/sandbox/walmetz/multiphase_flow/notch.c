/**
# Flow around a 2D circular cylinder with a square notch, containing a gas bubble. 
The gas bubble as several states (to reproduce an hypothetical external pressure). It has been first written by Francesco Picella.
It lacks an experimental reference... Enzo Santoromito (from PRISME) is doing it in its PhD :)

This case has no direct application. It has been created to study the effect of external pressure on the fluid-gas interface.
*/

#include "embed.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "navier-stokes/conserving.h"
#include "tension.h"
#include "popinet/contact/contact-embed.h"
#include "view.h"

#define T 70.
#define tMovie 0.1

double R0 = 0.5;
double XC = 0.;
double YC = 0.;
double NS = 0.1; // notch size
double THETA = pi/2.; // notch tangential location


/**
Solid Geometry. A cylinder with a square notch.*/

// Geometry ==========================================================================
static inline double cylinder_notch (double x, double y,
                                      double x_cylinder, double y_cylinder,
                                      double cylinder_radius,
                                      double notch_size,
                                      double notch_location)
{
  double ct = cos(notch_location), st = sin(notch_location);

  // Square center on the original circle boundary
  double xs = x_cylinder + (cylinder_radius-notch_size/2.)*ct;
  double ys = y_cylinder + (cylinder_radius-notch_size/2.)*st;

  // Local square coordinates:
  // X is radial, Y is tangential.
  double X =  (x - xs)*ct + (y - ys)*st;
  double Y = -(x - xs)*st + (y - ys)*ct;

  // Negative inside the circular solid
  double phi_circle =
    sq(x - x_cylinder) + sq(y - y_cylinder) - sq(cylinder_radius);

  // Negative inside the square notch
  double phi_square =
    max(fabs(X), fabs(Y)) - notch_size/2.;

  // Cylinder minus square
  return max(phi_circle, -phi_square);
}
/**
Gas geometry. A square notch within a cylinder.*/
static inline double notch_bubble (double x, double y,
                                   double x_cylinder, double y_cylinder,
                                   double cylinder_radius,
                                   double notch_size,
                                   double notch_location)
{
  double ct = cos(notch_location), st = sin(notch_location);

  // Le centre du carré reste inchangé pour s'aligner parfaitement avec le solide
  double xs = x_cylinder + (cylinder_radius - notch_size/2.)*ct;
  double ys = y_cylinder + (cylinder_radius - notch_size/2.)*st;

  // Local square coordinates: X est radial, Y est tangentiel
  double X =  (x - xs)*ct + (y - ys)*st;
  double Y = -(x - xs)*st + (y - ys)*ct;

  // 1. Frontière du cercle (inchangée)
  double phi_inside_circle =
    sq(cylinder_radius) - sq(x - x_cylinder) - sq(y - y_cylinder);

  // --- MODIFICATION POUR LA PRESSION et coller à Enzo---
  // delta est la distance dont la bulle est repoussée vers le fond
  double delta = notch_size * 0.25;

  // Le front extérieur est repoussé : X < (notch_size/2 - delta)
  double limite_X = max(X + delta, -X);

  double phi_inside_square =
    notch_size/2. - max(limite_X, fabs(Y));

  // Bubble = intersection
  return min(phi_inside_circle, phi_inside_square);
}
//=================================================================================



/**
Also, to avoid unphysical transient at startups, I implement
a simple ramp for external velocity.*/
static inline double ramp (double t, double t_ramp, double value)
{
  return value*min(t/t_ramp, 1.);
}


int maxlevel = 12;
int minlevel = 8;
double Reynolds = 1900.;


int main() {
  size (32.);

  /**
  We set the origin */

  origin (-L0/2., -L0/2.);

  init_grid (1 << (maxlevel-4));


  /**
  We use a constant viscosity. */

        mu2 = 1.*R0*2./Reynolds; // outer fluid viscosity 
        /**
        We set plastron viscosity equal to 1/100 of the main fluid one.*/
  mu1 = 0.018*mu2; // plastron viscosity 
  /**
  We set the surface tension coefficient. */

  f.sigma = 0.073;

        /**
        Set a constant contact angle.*/
        const scalar c[] = 90.*pi/180.; // fixed contact angle...
        contact_angle = c;

        run();
}

/**
We set the boundary conditions, so to obtain a flow around a fixed cylinder. */
u.n[left]  = dirichlet(ramp(t,1.,1.));
p[left]    = neumann(0.);
pf[left]   = neumann(0.);

u.n[right] = neumann(0.);
p[right]   = dirichlet(0.);
pf[right]  = dirichlet(0.);

/**
Must impose no-slip on embedded boundaries!*/
u.n[embed] = dirichlet(0.);
u.t[embed] = dirichlet(0.);

event init (t = 0)
{
/**
Initialize solid geometry, with automatic refinement.*/
#if TREE
  astats s;
  int iter = 0;

  do {
    iter++;

    // Recompute embedded fractions on the current grid
                solid (cs, fs, cylinder_notch (x, y, XC, YC,R0,NS,THETA));

    // Force refinement where cs varies, i.e. near the embedded boundary
    s = adapt_wavelet ({cs},
                       (double[]){1.e-30},
                       maxlevel,
                       minlevel);

  } while ((s.nf || s.nc) && iter < 100);
#endif
/**
Initialize bubble location.*/

  fraction (f,notch_bubble (x, y,XC,YC,R0,NS*1.1,THETA));
}

event adapt (i++) {
  adapt_wavelet ({cs,u,f}, (double[]){1e-3,1e-2,1e-2,1e-3}, maxlevel, minlevel);
        solid (cs, fs, cylinder_notch (x, y, XC, YC,R0,NS,THETA));
        /**
        Unrefine all uninteresting areas...such as the inlet!*/
        unrefine (x < -L0/3.);

        /**
        Purge gas that is (somehow) trapped within the solid phase.*/
        foreach()
                if(cs[]<=0)
                        f[] = 0.;

}

event movie(t+=tMovie,t<=T){
  view(fov=0.4, tx = 0., ty = -R0/L0);
  draw_vof("cs", "fs",filled = -1);
        draw_vof ("f", lc = {1, 0, 0}, lw = 3);
//      squares ("u.x", linear = true);
        // Draw grid only on upper part of flow
        cells (lc = {0.7, 0.7, 0.7});

          vectors ("u", scale = 0.01);


  save("movie.mp4");
}



event compute_forces (i += 1, t<=T)
{
  coord Fp, Fmu;
  embed_force (p, u, mu, &Fp, &Fmu);
// Output forces
  static FILE * fp = NULL;
        if(pid()==0){
        if (!fp) {
          fp = fopen ("force_coeff.dat", "w");
          fprintf (fp, "# i t dt Cx Cy\n");
        }
                double Cx = (Fp.x+Fmu.x)*2.;
                double Cy = (Fp.y+Fmu.y)*2.;
        fprintf (fp, "%3.2e %06d %+6.5e %+6.5e %+6.5e %+6.5e\n",
                 Reynolds, i, t, dt, Cx, Cy);
        fflush (fp);
        }
}