/**
[Return to my homepage](http://basilisk.fr/sandbox/nlemoine/README)

# 2D test case for the Shallow Ice Approximation diffusive equation:

This is the 2D version of Halfar's similarity solution (see [halfar1D.c](halfar1D.c)) for a flat substratum ($z_b = 0$) and a zero mass balance. The solution has a radial symmetry:
$$\eta(r,t)=h(r,t)=(aV)^{1/3}f^2(t)\,g\!\left(\frac{r f(t)}{(aV)^{1/3}}\right)$$
*/

#include "grid/quadtree.h"
#include "run.h"
#include "diffusion.h"
#include "view.h"
#include "nlemoine/view-utils.h"

#define LEVEL         8
#define MINLEVEL      6
#define ETAE          20. // error on surface elevation
#define sec_per_year  31557600.
#define tfin          9000.
#define nGlen	      3.
#define AGlen	      5.0e-24
#define M_PI          acos(-1.0)

scalar zb[],eta[];
face vector D[];
double Volume,aV;
double t0,Gamma;

double dt, tini;
mgstats mgd;

double Beta(double xx, double yy){
 return exp(lgamma(xx)+lgamma(yy)-lgamma(xx+yy));
}

/**
## Diffusivity function */

int update_diffusivity(scalar zb, scalar eta, face vector D)
{
  //    /!\ reminder : D.x and D.y are not co-located
  foreach_face()
  {
    double hm = (eta[]+eta[-1,0]-zb[]-zb[-1,0])/2.; // ice thickness interpolated at face center
    hm = max(hm,0.);
    double sn = (eta[]-eta[-1,0])/Delta; // slope component along face normal
    double st = (eta[0,-1]+eta[-1,-1]-eta[0,1]-eta[-1,1])/4./Delta;	// slope component transverse to face normal
    D.x[] = sec_per_year*Gamma*pow(hm,5.)*(sq(sn)+sq(st));
  }

 return(0);
}

/**
## Halfar 2D similarity solution (Halfar, 1983)
In [Halfar's paper](https://doi.org/10.1029/JC088iC10p06043) the solution is scaled using total volume $V$ multiplied by a constant $a$. Noting $R(t)$ the radius at time $t$, the explanation is as follows:
$$\begin{aligned}
V(t) & = \int_{\theta=0}^{2\pi}\int_{r=0}^{R(t)}h(r,t)\,r\,dr\,d\theta =  2\pi\int_0^{R(t)=\frac{(aV)^{1/3}}{f(t)}}(aV)^{1/3}f^2(t)\,g\!\left(\frac{r}{R(t)}\right)r\,dr\\ 
& = 2\pi(aV)\int_0^1 g(\xi)\,\xi\,d\xi
\end{aligned}$$
Moreover,
$$\int_0^1 g(\xi)\,\xi\,d\xi = \int_0^1 \left(1-\xi^{\frac{n+1}{n}}\right)^{\frac{n}{2n+1}}\xi\,d\xi = \frac{n}{n+1}\int_0^1 \left(1-\xi'\right)^{\frac{n}{2n+1}}\xi'^{\frac{2n}{n+1}-1}d\xi'$$

The latter integral is the [Beta function](https://en.wikipedia.org/wiki/Beta_function) $\Beta\!\left(\frac{3n+1}{2n+1},\frac{2n}{n+1}\right)$, so we find the condition for $a$:
$$V=\frac{2\pi n}{n+1}(aV)\ \Beta\!\left(\frac{3n+1}{2n+1},\frac{2n}{n+1}\right)$$
hence
$$a=\frac{n+1}{2\pi n\ \Beta\!\left(\frac{3n+1}{2n+1},\frac{2n}{n+1}\right)}$$
*/
double Halfar2D(double xx, double yy, double tt)
{
  double rr = sqrt(sq(xx)+sq(yy));
  double ft = pow(t0/tt,1./(5.*nGlen+3.)); // eq. (5)
  double cbrtaV = cbrt(aV);
  double argg = rr*ft/cbrtaV;  // eq. (3)
  double res = argg < 1. ? cbrtaV*sq(ft)*pow( (1.0 - pow(argg,1.+1./nGlen)) , nGlen/(2.*nGlen+1.) ) : 0. ; // eq. (4)
  return res  ;
}

/**
## Main */

int main (int argc, char * argv[])
{
  N = 1 << LEVEL;
  L0 = 2.0e6;
  size (L0);
  origin (-L0/2.,-L0/2.);

  Volume = 2.0e15;  // Total volume (m3)

  double bet = Beta((3.*nGlen+1.)/(2.*nGlen+1.),2.*nGlen/(nGlen+1.));
  aV = Volume*(nGlen+1.)/(2.*M_PI*nGlen*bet);
  Gamma  = 2. * AGlen * pow(9.81*910.,nGlen)/(nGlen+2.) ;
  t0 = pow((2.*nGlen+1.)/(nGlen+1.),nGlen)*pow(aV,-nGlen/3.)/(5.*nGlen+3.); // eq. (6)
  t0 = t0/Gamma/sec_per_year;

  tini = 1.;

  run();
  return(0);
}

event init (i=0)
{
  foreach()
  {
    zb[] = 0.;
    eta[] = Halfar2D(x,y,tini);
  }
  boundary({zb,eta});
}

scalar l[];
scalar analytical[],numerical[];
scalar msk_ana[], msk_num[];

event movie (t=0.; t += 25.)
{
  foreach()
  {
    analytical[] = Halfar2D(x,y,t+tini);
    numerical[] = eta[];
    boundary({analytical,numerical});
    msk_ana[] = x > 0. ? -1 : 1.;
    msk_num[] = x > 0. ? 1. : -1;
  }

  float qview_x[4],qview_z[4],qview[4];
  (void) gl_axis_to_quat ((float[]){1,0,0}, 2.*M_PI/6., qview_x);
  (void) gl_axis_to_quat ((float[]){0,0,1}, M_PI/32., qview_z);
  gl_add_quats(qview_x, qview_z, qview);

  view (fov = 19., quat = {qview[0],qview[1],qview[2],qview[3]},
	tx = 0., ty = -0.05,
        sx = 1., sy = 1., sz = 200.,
        width = 1200, height = 768);
  char s[80];
  masked_squares ("analytical", linear = true, z = "analytical", 
                  min = 0., max = 4000., mask = msk_ana);
  masked_squares ("numerical", linear = true, z = "numerical",
                  min = 0., max = 4000., mask = msk_num);
  //surf_cells(analytical, mask = msk_ana);
  surf_cells(numerical, mask = msk_num);

  sprintf (s," left : analytical");
  draw_string (s, 0, size = 48, lw = 2);
  sprintf (s,"right : numerical ");
  draw_string (s, 3, size = 48, lw = 2);

  sprintf (s,"t = %3.1f ka ", t/1000.);
  draw_string (s, 2, size = 48, lw = 2);
  sprintf (s," Z-stretch x 200");
  draw_string (s, 1, size = 48, lw = 2);
  
  char fname[200];
  sprintf(fname,"frame-%4.4d.png",(int)t);  
  save (fname);
  
  stats s_ana = statsf (analytical);
  stats s_num = statsf (numerical);
  fprintf(stderr,"t = %.0f ka, Max. analytical solution: %.2lf m, Max. numerical solution: %.2lf m\n",t,s_ana.max,s_num.max);
}

event stop (t = tfin)
{  
  system ("for f in frame-*.png; do convert $f ppm:- && rm -f $f; done | "
	  "ppm2mp4 animate.mp4");
  fprintf(stderr,"Done.\n");
  return 1;
}

/**
## Time integration
The SIA is basically a diffusion equation, so we use the [Poisson solver](http://basilisk.fr/src/diffusion.h) for integration. Since this solver is time-implicit, there isn't any stability problem; however, since the diffusivity depends on the solution, setting a too large time step will cause spurious features to develop in the solution if the diffusivity field is not updated often enough. So we set the time step to a ''reasonable'' (up to 10) multiple of the explicit time step
$$\Delta t_\textrm{explicit} = \frac{1}{4}\frac{\Delta x^2}{D_\textrm{max}}$$
*/
event integration (i++)
{
  (void) update_diffusivity(zb,eta,D);

  stats s = statsf (D.x);
  double dtExplicit = 0.25*L0*L0/N/N/s.max;
  dt = dtnext(4.0*dtExplicit);

  mgd = diffusion(eta,dt,D);

  // preserve positivity of ice tickness
  foreach()
   eta[] = max(zb[],eta[]);
}

// Adaptivity

int adapt() {
#if TREE
  astats s = adapt_wavelet ({eta}, (double[]){ETAE},
			    LEVEL, MINLEVEL);
//  fprintf (ferr, "# refined %d cells, coarsened %d cells\n", s.nf, s.nc);
  return s.nf;
#else // Cartesian
  return 0;
#endif
}

// And finally we set the event that will apply our *adapt()* function at every timestep.

event do_adapt (i++) adapt();

/**
## Animation of the solution
![Animation of the solution. The left-hand part of the image displays the analytical solution, and the right-hand part displays the numerical solution and the adaptive grid.](./halfar2D/animate.mp4)(width=75% )
 */