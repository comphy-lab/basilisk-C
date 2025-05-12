/**
# Sessile drop on an embedded boundary


The influence of gravity ($Eo=\rho_l g R0^2/\sigma > 0$) on the sessile droplet shape is presented


~~~gnuplot Equilibrium shapes for $0.5 \leq Eo \leq 50$
set term push
set term @SVG size 640,180
set size ratio -1
unset key
unset xtics
unset ytics
unset border
set xrange [-1.5:1.5]
set yrange [0:]
plot 'out' w l lt -1 lw 2 lc rgb "red" t '',\
     'out' u (-$1):2 w l lt -1 lw 2 lc rgb "red" t '' ,\
        0 lt -1 lw 3  t ''  
set term pop
~~~
*/
#include "grid/multigrid.h"
#include "embed.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "tension.h"
#include "../contact-embed.h"
#include "../geometrical_tools.h"
#include "view.h"

#define MINLEVEL 3
#define ADAPT 0
#define GRAVITY 1
#define R0 0.5
#define xc 0
#define yc 0 
# define wall(x,y) (y)
# define circle(x,y) (sq(x - xc) + sq(y - yc) - sq(R0))
# define T 10

double theta0=120, Eo=10;
double rf, hf; //final radius and height of the droplet
int LEVEL =7;

int main()
{
  size (2.);
  /**
  We shift the bottom boundary. */
  origin (0, -0.26);
  init_grid (1 << LEVEL);
  /**
  We use a variable viscosity and density. */
  mu1 = 0.1; mu2 = 0.1;
  rho1= 1; rho2 = 0.1;
  /**
  We set the surface tension coefficient. */
    f.sigma = 1;

  /** and vary the Eotvos number*/
  for (Eo=0.005; Eo<=50; Eo*=10) {
    const scalar c[] = theta0*pi/180.;
    contact_angle = c;
    run();
  }
}

event init (t = 0)
{
  /**
  We define the horizontal bottom wall and the initial (half)-circular
  interface. */
  solid (cs, fs, wall(x,y));
  fraction (f, - circle(x,y));
}

event acceleration (i++) {
  face vector av = a;
  foreach_face(y)
    av.y[] -= Eo*f.sigma/(rho1*sq(R0));
}

event logfile (i++; t <= T)
{
  /**
  If the curvature is almost constant, we stop the computation
  (convergence has been reached). */
  
  scalar kappa[];
  curvature (f, kappa);
  foreach()
    if (cs[] < 1.)
      kappa[] = nodata;
  if (statsf (kappa).stddev < 1e-6)
    return true;
}

/** Here we compute the final radius of the droplet (intersection between the fluid and the solid interface) */

#if 0
event final_radius (i++, t <= T){
  double dist=R0; 
  coord p;
  foreach() {
    if (interfacial(point, cs) && interfacial(point, f) && x>xc) { // potential triple cell
      coord n = interface_normal (point, f); 
      coord m = facet_normal (point, cs, fs); 
      coord o = {x, y, z};
      double alphan = plane_alpha (f[], n);
      double alpham = plane_alpha (cs[], m);
      double area = plane_area_center (m, alpham, &p);
      coord a[2]; coord b[2];
      if (facets (m, alpham, a) == 2 && facets (n, alphan, b) == 2) {
        foreach_dimension() {
          for (int i = 0; i <= 1; i+= 1) {
          a[i].x = o.x + a[i].x*Delta;
          b[i].x = o.x + b[i].x*Delta;
          }
        }
        coord It;
        foreach_dimension() It.x = 0.;
        if ( Intersect(a[0], a[1], b[0], b[1], &It) ){
          dist = sqrt(sq(It.x-xc)+ sq(It.y-yc));
          //ucl = embed_interpolate(point, unorm, p); // contact line velocity
        }
      }
    }
  }
  rf=dist; 
}
#endif

/** We can also compute the final height of the droplet considering the embedded boundary orientation */
#if 1
event final_height(i++, t <= T){
  coord n; 
  foreach(){
    if (cs[] < 1. && cs[] > 0) {
      n = facet_normal (point, cs, fs);
      normalize (&n);
      break;
    }
  }
  double detmin=10000, det;
  foreach(reduction(min:detmin)){
    if (f[] < 1.  && f[] > 0. && cs[]> 0 ){ 
      coord o = {x, y, z}; coord v = {o.x-xc, o.y-yc};
      normalize (&v);
      det=fabs(v.x*n.y-v.y*n.x);
      detmin=min(detmin,det);
      if (det==detmin) hf=sqrt(sq(o.x-xc)+sq(o.y-yc));
    }
  }  
}
#endif


#if ADAPT
event adapt (i++) {
  adapt_wavelet ({f,u}, (double[]){1e-4,1e-3,1e-3}, LEVEL, 3);
}
#endif


event end (t = end)
{
  /**
  At the end, we output the equilibrium shape. */
  output_facets (f, stdout);

   /**
   We also define the analytical radius and height for a given angle $\theta$ */
  double pitheta0 = theta0*pi/180.;
  double Rfa = R0*sqrt(pi/(2*(pitheta0-sin(pitheta0)*cos(pitheta0))));
  double hfa = Rfa*(1-cos(pitheta0));

  /**
   We write the analytical solutions in a file */
  double g = Eo*f.sigma/(rho1*sq(R0));
  double hfty = 2*sqrt(f.sigma/(rho1*g))*sin(pitheta0/2);
  fprintf (stderr, "%d %g %g %.5g %.5g\n", N, theta0, Eo, hf/hfa, hfty/hfa);
}

#if 0
event movie(i+=10,last){
  if (theta0 == 120 && Eo==5) {
    view(fov=20, tx = 0, ty = -0.5);
    draw_vof ("f", lw=2);
    draw_vof("cs", "fs",filled=-1);
    squares("f", linear = true, min = 0, max = 1);
    save("movie.mp4");
  }
}
#endif

/**
We compare $h_{\infty}/h_f$ to the analytical expression, with $h_{\infty}=2(\sqrt{\sigma/\rho g} )sin(\theta_s/2)$.

~~~gnuplot
reset
set logscale x
set format x "10^{%L}"
set xlabel 'Eo' 
set ylabel 'h_{/Symbol \245}/h_f' 
set xrange[0.0005:100]
set yrange[0:1.5]
set ytics 0,0.5,1.5
theta0=120*pi/180.

plot (2*sqrt(1/x)*sin(theta0/2.))/((1-cos(theta0))*sqrt(pi/(2*(theta0-sin(theta0)*cos(theta0))))) lt -1 dt 2 lc rgb "black" lw 2.5,\
     'file128_save.dat' u 3:4 w p pt 4 ps 2.5 lt -1 lw 2,\
     1 lt -1 dt 0 lw 2
~~~

![Relaxation toward a $120^\circ$ contact angle with gravity ($Eo=5$).](sessile-embed-gravity/movie.mp4)

*/
