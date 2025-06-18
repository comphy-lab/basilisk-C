/**
# Minimal example for Hybrid Lagrangian-Eulerian $\omega-\psi$ solver

![Works pretty good. Well done Lagrangian advection!](Lag_stream/mode3.mp4)
 */
#include "poisson.h"
#include "run.h"
#include "timestep.h"
#include "tracer-particles.h"
#include "view.h"
#include "scatter2.h"

scalar n[], omega[], psi[];

vector u[];

Particles p;

int maxlevel = 7;
//Function prototype for colorfull particles.
void scatter_color (Particles p, float s = 3, Colormap map = jet, double minv = -1, double maxv = 1);

int main() {
  L0 = 8;
  origin (-pi*4./3., -2*exp(1) + 1.2); // Not centered
  N = 1 << maxlevel;
  foreach_dimension() {
    psi[right] = dirichlet (0);
    psi[left] = dirichlet (0);
  }
  DT = 1e-1;
  run();
}

#define RAD (sqrt((sq(x) + sq(y))))
#define THETA(M) (M*asin(x/RAD))
#define RADP(P, M) ((1 + P*sin(THETA(M))))

double P = 1e-5, m = 3; // Perturbation and mode
double b = 1.45; //Vortex parameter

event init (t = 0) {
  CFL = 2;
  double betav = 1./(1 - sq(b));
  p = new_tracer_particles(0);
  scalar omega_uf[]; //unfiltered omega 
  foreach() {
    double rp = RAD*RADP(P,m), omg = 0;
    if (rp <= 1.)
      omg = 2;
    else if (rp > 1 && rp <= b) 
      omg = 2*betav;
    omega_uf[] = omg;
  }
  const face vector alphaf[] = {-sq(0.4/(2.*pi)), -sq(0.4/(2.*pi))};
  poisson (omega, omega_uf, alphaf, unity);  
  foreach() {
    double o = omega[];
    foreach_child() {
      particle pn;
      pn.x = x;
      pn.y = y;
      pn.z = o;
      add_particle(pn, p);
    }
  }
  P_RK2 = true; // There is no velocity update between step 3 and 1
}
/**
## The solver ...

... is impremented in this event.
 */
event velocity (i++) {
  dt = dtnext (timestep (u, DT));
  // Reset Eulerian fields
  foreach() {
    omega[] = n[] = 0;
  }
  // Diagnose "average" vorticity (poormans finite pointset method)
  foreach_particle_in(p) {
    foreach_point(x,y) {
      n[]++;
      omega[] += p().z;
    }
  }
  foreach() {
    if (n[] > 0)
      omega[] /= n[];
    else
      omega[] = HUGE;
  }
  // Apply a stupid filter to guess omega in "empty" cells 
  foreach() {
    if (!n[]) {
      int nn = 0;
      double on = 0;
      foreach_neighbor(1) {
	if (n[]) {
	  on += omega[];
	  nn++;
	}
      }
      if (nn)
	omega[] = on/nn;
      else
	omega[] = 0;
    }
  }
  // Find stream function and velocity field.
  poisson(psi, omega);
  foreach() {
    u.x[] = -(psi[0,1] - psi[0,-1])/(2*Delta);
    u.y[] = (psi[1] - psi[-1])/(2*Delta);
  }
}

event mov (t += 0.1, last) {
  view (width = 800, height = 800);
  scatter_color (p, s = 2, map = blue_white_red);
  scatter (p, s = 1);

  save ("mode3.mp4");

}

event stop (t = 20);

void scatter_color (Particles p, float s = 3, Colormap map = jet, double minv = -1, double maxv = 1) {
  double cmap[NCMAP][3];
  blue_white_red(cmap);
  bview * view = draw();
#if dimension == 2
  glTranslatef (0., 0., view->lc*view->fov/24.); //from draw_lines()
#endif
  glPointSize(s*view->samples);
  glBegin (GL_POINTS);
  foreach_particle_in(p) {
    Color b = colormap_color (cmap, p().z, minv, maxv);		
    glColor3f (b.r/255., b.g/255., b.b/255.);			
    
#if dimension == 2    
    glvertex2d (view, x, y);
#else // dimension == 3
    glvertex3d (view, x, y, z);
#endif
  }
  glEnd();
  view->ni++; 
}
