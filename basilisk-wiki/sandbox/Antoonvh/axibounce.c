/**
![Astronaut Scott Kelly played ping pong in space. Video courtesy of
 [NASA Johnson](https://www.youtube.com/watch?v=TLbhrMCM4_0)
 ($\leftarrow$ click to watch the full
 movie)](https://lh3.googleusercontent.com/-HPsDGunAGPQ/VqPMV9unuEI/AAAAAAAB5zY/AZOVAQaaey8vT47GYfijW3scP_2FZQZmQCJoC/w500-h429-n/SpacePingPong.gif)

# A Bouncing Axisymetric Droplet

On this page a quasi three dimensional (3D) simulation of the bouncing
droplet is presented. The set-up is very similar to that of the planar
(2D) version [here](bounce.c), but now in the 3D-ish, axisymmetric
limit.
*/

#include "axi.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "tension.h"
#include "view.h"

/**
## The case set-up

Our "vertical" axis, if that concept is suitable to use inside a space
station, is now the swapped with the former horizontal axis. Hence,
what used to be associated with the `bottom` boundary in the 2D
example is now called `left`. Furthermore, the radial coordinate ($r$)
is associated with the $y$-direction in Basilisk's axisymmetric
formulation.
 */

f[left] = 0.;
u.t[left] = dirichlet(0.);
f[right] = 0.;
u.t[right] = dirichlet(0.);

/**
   The physical system consists of a liquid droplet with a radius $R$
   that travels trough the air towards the hydrophobic wall with a
   velocity $U$. Gauged from the movie, we approximate $R \approx 2
   \text{ cm}$ and $U\approx 3 \text{cm/s}. Together with the fluid
   properties of the air and water ($\rho_a,\rho_w,\mu_a,\mu_w,
   \sigma$) we can identify some dimensionless groups that describe
   the system.
   
   $$ Re = \frac{\rho_a UR}{\mu_a} \approx 30, $$ 
   
   $$ We = \frac{\rho_w U^2R}{\sigma} \approx 0.2, $$
   
   $$ \Pi_1 = \frac{\mu_w }{\mu_a} \approx 100, $$
   
   $$ \Pi_1 = \frac{\rho_w }{\rho_a} \approx 1000 \rightarrow 100, $$
   
   $$ D = 3 \rightarrow 2\text{; Axisymmetric} , $$
   
   where $D$ stands here for the number of spatial dimensions of the
   system and the '$\rightarrow$' symbols indicate dimensionless
   groups where we make a consession with respect to the reality
   onboard the space station to keep the numerical experiment feasable
   on a single-core system.

We set the model parameters such that we obtain the approximated
ratios, aproximately. It is done in such a way that we have a
normalized velocity scale ($U$), length scale ($R$) and gas phase
density ($\rho_a$) scale. Furthermore, we set the maximum grid
resolution to correspond to a $256^2$ equidistant grid. For some
minimally desired consistency, we take special care such that the
radial coordinate does not have any negative values.
*/

double R = 1.;
double U = 1.;
int maxlevel = 8;
double temp = 50.;

int main(){
  mu2 = 1./30.; //gas phase
  mu1 = 100./30.; //liquid phase
  rho2 = 1.; //gas phase
  rho1 = 100;//liquid phase
  f.sigma = 500.;
  init_grid(1 << 7);
  L0 = 10;
  X0 = -L0/2;
  Y0 = 0.;
  run();
}

/**
## initialization

After we have refined the grid with a ring of high resolution mesh,
the droplet is initialized and it is targeted at the left wall.
 */

event init(i = 0){
  refine(sq(x) + sq(y) < sq(R + 0.25) && sq(x) + sq(y) > sq(R - 0.25) && level < maxlevel);
  fraction(f, sq(R) - sq(x) - sq(y));
  foreach()
    u.x[] = -f[]*U;
}
/**
## Adaptation

Since the advective interface tracking and resolving for the surface
tencile waves only requires a high resolution mesh at the locations of
the interface, it seems smart to use grid adaptivity, we will check it
this was indeed the case.
*/
event adapt(i++)
  adapt_wavelet((scalar *){u,f}, (double []){0.02, 0.02, 0.001}, maxlevel);

/**
## Output

The general dynamics are visualized in a movie that shows the rebound
of the droplet. This requires to define some additional functions.
 */


static void glvertex3d (bview * view, double x, double y, double z) {
  if (view->map) {
    coord p = {x, y, z};
    view->map (&p);
    glVertex3d (p.x, p.y, p.z);
  }
  else
    glVertex3d (x, y, z);
}

struct _draw_vof_axi {
  char * c;
  char * s;
  bool edges;
  double larger;
  int filled;
  char * color;
  double min, max, spread;
  bool linear;
  colormap map;
  float fc[3], lc[3], lw;
};

trace
bool draw_vof_axi (struct _draw_vof_axi p)
{
  scalar c = lookup_field (p.c);
  if (c.i < 0) {
    fprintf (stderr, "draw_vof(): no field named '%s'\n", p.c);
    return false;
  }
  face vector s = lookup_vector (p.s);
  
  colorize_args (p);
  
  double cmin = 1e-3; // do not reconstruct fragments smaller than this

#if TREE
  // make sure we prolongate properly
  void (* prolongation) (Point, scalar) = c.prolongation;
  if (prolongation != fraction_refine) {
    c.prolongation = fraction_refine;
    boundary ({c});
  }
#endif // TREE
    
  bview * view = draw();
  colorize() {
    foreach()
      if (cfilter (point, c, cmin)) {
	coord n = facet_normal (point, c, s);
	double alpha = plane_alpha (c[], n);
	coord v[2];
	int m = facets (n, alpha, v);
	if (m == 2) {
	  coord g[4];
	  g[0].x = v[0].x; g[0].y = v[0].y; g[0].z = -y; 
	  g[1].x = v[0].x; g[1].y = v[0].y; g[1].z = y; 
	  g[2].x = v[1].x; g[2].y = v[1].y; g[2].z = y;
	  g[3].x = v[1].x; g[3].y = v[1].y; g[3].z = -y;
	  color_facet (p);
	  if (view->gfsview)
	    glNormal3d (- n.x, - n.y, - n.z);
	  else
	    glNormal3d (n.x, n.y, n.z);
	  for (double th=0; th < 2*pi ; th+=0.25){
	    glBegin (GL_POLYGON);
	    for (int i = 0; i < 4; i++) {
	      color_vertex (p, interp (point, g[i], col));
	      glvertex3d (view, x + g[i].x*Delta,
			  cos(th)*(y + g[i].y*Delta) + sin(th)*(z + g[i].z*Delta),
			  sin(th)*(y + g[i].y*Delta) + cos(th)*(z + g[i].z*Delta));
	    }
	    glEnd ();
	    view->ni++;
	  }
	}
      }
  }
#if TREE
  // revert prolongation
  if (prolongation != fraction_refine) {
    c.prolongation = prolongation;
    boundary ({c});
  }
#endif // TREE
  return true;
}

event viewer(t=0.05;t+=0.1;t<=temp){
  view(psi = -pi/2, theta = 0.1, phi = 0.2);
  box();
  draw_vof_axi("f");
  save("f.mp4");
}

/**
The simulation does display the ping pong dynamics;

![Well done draw_fov_axi!](axibounce/f.mp4)

The energy partitioning between kinetic and the surface free energy is
again calculated. The formulation as was used in the planar simulation
is modified to correct for the fact that the $y$ coordinate now
represents the radial coordinate.
 */

void find_facets(Point point, scalar f, double xy[4]){
  coord n;
  n = mycs (point, f);
  double alpha = plane_alpha (f[], n);
  coord segment[2];
  if (facets (n, alpha, segment) == 2){
    xy[0] = x + segment[0].x*Delta;
    xy[1] = y + segment[0].y*Delta;
    xy[2] = x + segment[1].x*Delta;
    xy[3] = y + segment[1].y*Delta;
  }else{
    printf("Warning:\nCould not find facets; expect unexpected behaviour.\n");
  }
}

event energy(i += 5){
  static FILE * fpe = fopen("axienergy","w");
  double e = 0;
  double a = 0;
  double xyf[4];
  foreach(reduction(+:e) reduction(+:a)){
    e += 0.5*(sq(Delta)*M_PI*2*y)*(sq(u.x[]) + sq(u.y[])) * ((rho1*f[]) + (rho2*(1-f[])));
    if (f[] > 0.00001 && f[] < 0.9999){
      find_facets(point, f, xyf);
      a += 2*M_PI*((xyf[1] + xyf[3])/2.)*pow(sq(xyf[0] - xyf[2]) + sq(xyf[1] - xyf[3]), 0.5);
    }
  }
  fprintf(fpe,"%g\t%d\t%g\t%g\t%g\t%g\n", t, i, e, (a-(4*M_PI*sq(R)))*f.sigma,a, ((a-(4*M_PI*sq(R)))*f.sigma) + e);
  printf("%g\t%d\t%g\t%g\t%g\t%g\n", t, i, e, (a - (4*M_PI*sq(R)))*f.sigma, a,((a - (4*M_PI*sq(R)))*f.sigma) + e);
}

/**
We can view the results;

~~~gnuplot
set xlabel 'time'
set ylabel 'Energy'
plot 'axienergy' u 1:3 w l lw 5 t 'Kinetic' ,\
     'axienergy' u 1:4 w l lw 3 t 'Potential' ,\
     'axienergy' u 1:6 w l lw 3 t 'Total'     
~~~

Overall, things seem familiar and consistent to what was learned from the planar case. Differences may be spotted when doing a more quantitative analysis of the energy partioning.

## The next step
In the movie, the droplet does not really appear to be axisymetric. It is generally quite challenging to make-up an a priori argument on how to dynamics would alter when the physical system is allotted with an additional spatial dimension. E.g. after having studied 2D turbulence (a la Kaichmann), would you be able to expect the main flow-behavioural characteristics of 3D turbulence (a la Kolmogorov)? Well, Not me! Therefore, it is required to break the axisymetry and study the flow in a fully-fleched 3D simulation. This study should also adresses the influence of the consession that was made with regards to the ratio of the phases' densities.
*/
