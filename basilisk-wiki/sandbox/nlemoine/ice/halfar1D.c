/**
[Return to my homepage](http://basilisk.fr/sandbox/nlemoine/README)

# 1D test case for the Shallow Ice Approximation diffusive equation:

$$\frac{\partial\eta}{\partial t} + \nabla\cdot(h\mathbf{U}) = M$$
where $\eta$ is the ice surface elevation, $z_b$ the elevation of the substratum, $h = \eta - z_b$ the ice thickness, $\mathbf{U}$ the depth-average velocity and $M$ the local mass balance (difference between accumulation and ablation).

In the SIA, the flux density $\mathbf{q} = h\mathbf{U}$ boils down to a diffusive flux with a nonlinear diffusivity $D$ which depends on ice thickness and ice surface slope:
$$\mathbf{q} = h\mathbf{U} =  -\underbrace{\textstyle\frac{2A (\rho g)^n}{n+2} h^{n+2}\big\|\nabla \eta\big\|^{n-1}}_{D} \nabla\eta$$
$n$ is the exponent of Glen flow law, $n\approx 3$

Halfar 1D similarity solution ([Halfar, 1981](https://doi.org/10.1029/JC086iC11p11065)) holds for a flat substratum ($z_b = 0$) and a zero mass balance ($M=0$):
$$\eta(x,t)=h(x,t)=\sqrt{aV}f(t)\,g\!\left(\frac{x f(t)}{\sqrt{aV}}\right)$$
*/

//#include "grid/multigrid1D.h"
#include "grid/bitree.h"
#include "run.h"
#include "diffusion.h"

#define LEVEL         8
#define MINLEVEL      6
#define ETAE          20 // error on surface elevation
#define DE            1e7 // error on diffusivity (m2/yr)
#define sec_per_year  31557600.
#define tfin          100000.
#define nGlen	      3.
#define AGlen	      5.0e-24

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
  foreach_face()
  {
    double hm = (eta[]+eta[-1]-zb[]-zb[-1])/2.; // ice thickness interpolated at face center
    hm = max(hm,0.);
    double sn = (eta[]-eta[-1])/Delta; 				// slope component along face normal
    D.x[] = sec_per_year*Gamma*pow(hm,nGlen+2.)*pow(fabs(sn),nGlen-1.);
  }

 return(0);
}

/**
## Halfar 1D similarity solution (Halfar, 1981)
In [Halfar's paper](https://doi.org/10.1029/JC086iC11p11065) the solution is scaled using total volume $V$ multiplied by a constant $a$. Noting $W(t)$ the half-width at time $t$, the explanation is as follows:
$$\begin{aligned}
V(t) & = 2\int_0^{W(t)=\frac{\sqrt{aV}}{f(t)}}\sqrt{aV}f(t)\,g\!\left(\frac{x}{W(t)}\right)dx\\ 
& = 2(aV)\int_0^1 g(\xi)d\xi
\end{aligned}$$
Moreover,
$$\int_0^1 g(\xi)d\xi = \int_0^1 \left(1-\xi^{\frac{n+1}{n}}\right)^{\frac{n}{2n+1}}d\xi = \frac{n}{n+1}\int_0^1 \left(1-\xi'\right)^{\frac{n}{2n+1}}\xi'^{\frac{n}{n+1}-1}d\xi'$$
The latter integral is the [Beta function](https://en.wikipedia.org/wiki/Beta_function) $\Beta\!\left(\frac{3n+1}{2n+1},\frac{n}{n+1}\right)$, so we find the condition for $a$:
$$V=\frac{2n}{n+1}(aV)\ \Beta\!\left(\frac{3n+1}{2n+1},\frac{n}{n+1}\right)$$
hence
$$a=\frac{n+1}{2n\ \Beta\!\left(\frac{3n+1}{2n+1},\frac{n}{n+1}\right)}$$
*/
double Halfar1D(double xx, double tt)
{
    double ft = pow(t0/tt,1./(3.*nGlen+2.)); // eq. (12)
    double argg = fabs(xx*ft/sqrt(aV));	// eq. (10)
    double res = argg < 1. ? sqrt(aV)*ft*pow( (1.0 - pow(argg,1.+1./nGlen)) , nGlen/(2.*nGlen+1.) ) : 0. ; // eq. (11)
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

  Volume = 2.0e9;  // Total volume (m2)

  double bet = Beta((3.*nGlen+1.)/(2.*nGlen+1.),nGlen/(nGlen+1.));
  aV = Volume*(nGlen+1.)/(2.*nGlen*bet);
  Gamma  = 2. * AGlen * pow(9.81*910.,nGlen)/(nGlen+2.) ;
  t0 = pow((2.*nGlen+1.)/(nGlen+1.),nGlen)*pow(aV,-nGlen/2.)/(3.*nGlen+2.); // eq. (13)
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
    eta[] = Halfar1D(x,tini);
  }
  boundary({zb,eta});
}

event profile (t=0. ; t+=100.)
{
  FILE * fp;
  char name[80];
  double xp,eta_num,eta_ana;

  sprintf(name,"transect-%d.dat",(int)t);
  fp = fopen(name,"w");

  scalar Dc[];
  foreach()
     Dc[] = (D.x[]+D.x[1])/2.;

  for(int np = 0;np<N;np++)
  {
    xp = -L0/2.+(L0/N)*((float)np+0.5);
    eta_num = interpolate (eta, xp);
    eta_ana = Halfar1D(xp,t+tini);
    double Dp = interpolate (Dc, xp);
    Point point = locate(xp);
    fprintf (fp, "%g %g %g %g %d\n",xp/1000.,eta_num,eta_ana,Dp,point.level);
  }
  
  fclose(fp);
}

event stop (t = tfin);

/**
## Time integration
The SIA is a basically a diffusion equation, so we use the [Poisson solver](http://basilisk.fr/src/diffusion.h) for integration. Since this solver is time-implicit, there isn't any stability problem; however, since the diffusivity depends on the solution, setting a too large time step will cause spurious features to develop in the solution if the diffusivity field is not updated often enough. So we set the time step to a ''reasonable'' (up to 10) multiple of the explicit time step
$$\Delta t_\textrm{explicit} = \frac{1}{2}\frac{\Delta x^2}{D_\textrm{max}}$$
*/
event integration (i++)
{
  (void) update_diffusivity(zb,eta,D);

  stats s = statsf (D.x);
  double dtExplicit = 0.5*L0*L0/N/N/s.max;
  dt = dtnext(4.0*dtExplicit);

  mgd = diffusion(eta,dt,D);

  // preserve positivity of ice tickness
  foreach()
   eta[] = max(zb[],eta[]);
}

// Adaptivity

int adapt() {
#if TREE
  scalar Dc[];
  foreach()
     Dc[] = (D.x[]+D.x[1])/2.;

//  astats s = adapt_wavelet ({eta,Dc}, (double[]){ETAE,DE},
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
~~~gnuplot
   reset
   set xlabel "x (km)"
   set ylabel "Ice thickness (m)"
   ydim = 800 
   xdim = 480
   plot 'transect-0.dat' u 1:3 w lines lc rgb "black" lw 3 title 't = 0', \
        'transect-100.dat' u 1:3 w lines lc rgb "red" title 't = 0.1 ka (analytical)', \
        'transect-100.dat' u 1:2 w point pointtype 6 ps 1 lc rgb "red" title 't = 0.1 ka (numerical)', \
        'transect-1000.dat' u 1:3 w lines lc rgb "yellow" title 't = 1 ka (analytical)', \
        'transect-1000.dat' u 1:2 w point pointtype 6 ps 1 lc rgb "yellow" title 't = 1 ka (numerical)', \
        'transect-10000.dat' u 1:3 w lines lc rgb "green" title 't = 10 ka (analytical)', \
        'transect-10000.dat' u 1:2 w point pointtype 6 ps 1 lc rgb "green" title 't = 10 ka (numerical)', \
        'transect-100000.dat' u 1:3 w lines lc rgb "blue" title 't = 100 ka (analytical)', \
        'transect-100000.dat' u 1:2 w point pointtype 6 ps 1 lc rgb "blue" title 't = 100 ka (numerical)'   
~~~
*/

