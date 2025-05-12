/**
# Sessile drop on an embedded wall

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

Here, a drop is initialized as a half-disk (i.e. the
initial contact angle is 90$^\circ$). The difference with the initial sessile.c 
test case is that we run the simulation of the sessile droplet on an embedded
horizontal solid and test if the contact angle are well imposed between the droplet
and the embedded solid.

*/
#include "grid/multigrid.h"
#include "embed.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "tension.h"
#include "tavares/contact-embed.h"

/**
Geometrical parameters of the embed wall. */
#if INCLINED
# define wall(x,y) (y - x - 0.97)
# define circle(x,y) (sq(x - 0) + sq(y - 1) - sq(0.25))
#else
#define wall(x,y) (y - 0.26)
# define circle(x,y) (sq(x) + sq(y - 0.26) - sq(0.25))
#endif

/** Boundary conditions and global parameters definition */

// u.n[embed] = dirichlet(0.);
u.t[embed] = dirichlet(0.);
u.t[bottom] = dirichlet(0.);
u.t[top]  = dirichlet(0.);

/**
Flow and contact angle parameters */
scalar tag[], kappa[], triple_cell[], f1[], Alpha, app_angle[];//,Re[];
vector ns[], nf[], nf1[]; 
int angle = 90;
int LEVEL = 6;
double volume_vof_init;

int main()
{
  L0=2;
  origin(-1.,0.);
  //size (2);
  init_grid (1 << LEVEL);
  
  rho1 = 1.;
  mu1 = 0.1;
  
  rho2 = 1.;
  mu2 = mu1;
  
  /**
  We use a constant viscosity. */
  
  
  /**
  We set the surface tension coefficient. */
  
  f.sigma = 1.;

  /**
  With bview we control interactively the contact angle and maximum level of refinement. */
  //display_control (angle, 0, 180, "Angle");
  //display_control (LEVEL, 5, 10, "Level");
  for (angle = 15; angle <= 165; angle += 15)
  run();
}

event init (t = 0)
{

  vertex scalar phi[];
  foreach_vertex() {
    phi[] = wall(x,y);
  }
  boundary ({phi});
  fractions (phi, cs, fs);
  fraction (f, - circle(x,y));
}

#if 0
event logfile (i++)
{
  fprintf (fout, "%g %g\n", t, normf(u.x).max);
}
event snapshots (t += 1)
{
  p.nodump = false;
  dump();
}
#endif

#define T 20.

#if 1
event contact(i++; t<=T){
  //scalar Alpha[];
  /**Identifier les cellules coupées et leurs voisins dans le solide*/
  triple_cell_dectection(angle, f, cs, ns, nf1, nf);
  /**Reconstruire la fraction volumique dans les ghost cell*/
  fraction_reconstruction(f, nf);
}
#endif

#if 0
event apparent_angle(i++;t<=T){
  foreach(){
    app_angle[] = 0.;
    if (triple_cell[]){
      //coord ntemp1 = interface_normal (point, f);
      foreach_dimension() app_angle[] += nf1.x[]*ns.x[] ;
      app_angle[]= (acos(app_angle[])*180)/pi;
    } 
  }
}

event effective_radius(i++,t<=T){
  double Re;
  foreach(){  
    if (triple_cell[]){
      if (interfacial(point, f)){
        Re[]=sqrt(sq(x)+sq(y-0.26));
      }  
    }
  }
}

#endif


#if 1
event curve(i++, i<=1){
  curvature (f, kappa);
}
#endif

#if 0
event logfile (t = 0) {
  stats s = statsf (f);
  volume_vof_init = s.sum;
}
#endif

#if 0
event volume (i++, i<=1){
  static FILE * fp = fopen ("volume","w");
  stats s = statsf (f);
  double erreur = ((volume_vof_init - s.sum)/volume_vof_init)*100;
  fprintf (fp, "%g %.5g\n", t, erreur); 
}
#endif

/**
At equilibrium (t = 15), we output the interface
shape and compute the (constant) curvature. */
#if 1
event end (t = T)
{
  output_facets (f, stdout);
  //curvature (f, kappa);
  stats s = statsf (kappa);
  double R = s.volume/s.sum, V = statsf(f).sum;
  fprintf (stderr, "%d %d %.5g %.3g\n", LEVEL, angle, R/sqrt(V/pi), s.stddev);
}
#endif

/**
We compare $R/R_0$ to the analytical expression, with $R_0=\sqrt{V/\pi}$.

~~~gnuplot
reset
set xlabel 'Contact angle (degrees)'
set ylabel 'R/R_0'
set arrow from 15,1 to 165,1 nohead dt 2
set xtics 15,15,165
plot 1./sqrt(x/180. - sin(x*pi/180.)*cos(x*pi/180.)/pi) t 'analytical', \
  'log' u 2:3 pt 7 t 'numerical'
~~~

*/
