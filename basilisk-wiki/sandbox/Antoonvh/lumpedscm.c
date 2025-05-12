/**
Note: This page is best viewed with the Firefox web browser.

# A Single-Column model for the Atmospheric Diurnal Cycle

We time integrate a 1D evolution equation for the vertical profiles of
the horizontal velocity components ($u,v$) and the buoyancy ($b$), 
using an adaptive bitree grid, a generic timeloop iterator and the
reaction-diffusion solver. 
*/

#include "grid/bitree.h"
#include "run.h"
#include "diffusion.h"
/**
For the turbulent mixing we use the following mixing functions
($f(Ri)$), for stable (`fris`) and unstable (`friu`) stratifications,
 */
//#define fris(Ri) (sq((1 - (Ri/0.20)))*(Ri < 0.20)) // Critical Ri
#define fris(x) (exp(-10.*x))                        // Exponential
//#define fris(x) (1./(1. + (10*x*(1.+8.*x))))       // Long Tail
#define friu(Ri) (sqrt(1. - (18.*Ri)))               // Holtslag en Boville 1992
// Note: no surface f(Ri)
/**
  The thermodynamic variable is forced at the surface with the surface
buoyancy flux ($B$) according to a a simple energy balance:

$$B = Q* + G,$$

where $Q*$ in a prescribed netto radiation and $G$ is a 'feed back
flux' that is dynamically evaluated as:

$$G = -\Lambda \left( b_{surf} - b_{ref} \right),$$

with $b_{ref}$ a constant reference temperature/buoyancy scale and
$\Lambda$ a coupling strength. The buoyancy at the surface ($b_{surf}$)
is evaluated using the lowest two cell values for $b$ ($b_1,
b_2$) and a sub-grid-scale model that assumes a logaritmic profile for
$b(z)$.

$$b_{surf} = \frac{b_2 - b_1*c}{1-c},$$

with, 

$$c = \frac{\mathrm{ln}\left(\frac{4\Delta}{z_{0,h}}\right)}{\mathrm{ln}\left(\frac{\Delta}{z_{0,h}}\right)}$$

that is the result from an integration exercise, assuming that the
grid-cell size ($\Delta$) is much larger than the roughness length
for heat ($z_{0,h}$); $\Delta \gg z_{0,h}$. The surface buoyancy can
be evaluated from the lowest to cells according to a lookup table that
stores the values for $c(\Delta)$. Furthermore, we `#define` macros to
compute $G$(`GFLX`) and $Q*$(`Qn`), such that we can readily evaluate
$B$.
 */
#define BSURF ((b[1] - b[]*c[level])/(1. - c[level]))
#define GFLX (-Lambda*(BSURF - bo))
/**
Prescribing,

$$Q*= \mathrm{max}\left[ B_1\mathrm{sin}\left(\frac{2 \pi t}{T} \right),
B_1 \right]$$.
 */

#define Qn (max(B0*sin(2.*M_PI*t/T), B1))
/**
We declare fields for $u, v$ and $b$ as `u`, `v` and `b`, respectively and
also set the value for $L_c = 1030$m and $T = 24\times3600$s.
*/		
scalar u[], v[], b[];
double Lc = 1030.;
double T = 24.*3600.;
/**
global parameters are declared for $\Lambda, B_0, B_1, N, U_{geo},
z_{0,m}$ and $z_{0,h}$.
 */
double Lambda, B0, B1, f, Nv, Ugeo, zom, zoh;
/**
The numerical values of which will be computed from the values of the
dimensionless groups:
 */
double Pi1 = -6;
double Pi2 = 2160.;
double Pi3 = 10.;
double Pi4 = 5366.;
double Pi5; // This one will be varied
double Pi6 = 5150.;
double Pi7 = 1.;

/**
We initialize/declare some usefull variables for the `deep soil`
buoyancy, the Von Karmann constant, a lookup table to help determine
$b_{surf}$ and the surface friction velocity ($u*$), the Mixed-layer
depth ($\mathcal{L}_c$), the number of averaging iterations (for the
inversion data), the maximum level of refinement, a structure for the
multigrid-solver details, the maximum windspeed, the inversion, $U$ at
$\approx40$m and a file pointer for the output of the inversion
strengths.
*/
double bo = 0., k = 0.4;
double c[20], lut[20], Lmix;

int ni, maxlevel = 9;
mgstats mgb;
double Umax, inv, u40;
FILE * fpi;
/**
We set boundary conditions where `left` and `right` are to be
interpreted as the surface and the domain top, respectively.
 */
b[left] = dirichlet (BSURF);
b[right] = dirichlet (bo + sq(Nv)*x);
u[left] = dirichlet (0.);
v[left] = dirichlet (0.);

int main(){
  /**
The file is opened and the parameter values are calculated from the
dimensionless groups. 
   */
  fpi = fopen ("Inversions", "w");
  Nv = Pi2/T;                    //Initial stratification strength
  f = Pi3/T;                     //Coriolis parameter
  B0 = sq(Lc*Nv)*M_PI/(2.*T);    //B_0
  B1 = B0/Pi1;                   //B_1
  Lambda = sqrt(B0*T)/Pi4;       //Lambda
  zom = Lc / Pi6;                //z_0,m
  zoh = zom * Pi7;               //z_0,h
  L0 = 3.*Lc;                    //Domain size
  /**
We vary the value of $\Pi_5$ from 1 to 8 in 57 runs and compute the
value of $U_{geo}$ accordingly.
   */
  for (Pi5 = 1.;  Pi5 <= 8; Pi5 += 0.125){
    Ugeo = Pi5*pow(B0*Lc, 1./3.);
    init_grid (128);
    run();
  }
}

/**
   This event initializes the setup.  
 */
event init (t = 0){
  /**
First, the inversion-strength data is reset.
   */
  ni = 0;
  inv = Umax = u40 = 0.;
  /**
For accurate time integration we *initially* set a small timestep (1
second) and only allow a small tolerance on the residual for the
Poisson problem. The timestep will be adaptive, the `TOLERANCE`
remains small.
   */
  DT = T/(3600.*24.);
  TOLERANCE = 1E-4;
  /**
Solution fields are initialized after the near-surface grid is refined
to the maximum resolution.
   */
  refine(x < (Lc/10.) && level < maxlevel);
  foreach(){
    u[] = Ugeo;
    v[] = 0.;
    b[] = bo + sq(Nv)*x;
  }
  /**
For the computations of $b_{surf}$(`c`) and the friction velocity $u*$
(`lut`), we create lookup tables for the various possible grid
resolutions thay *may* be used at the surface.
   */
  for (int j = 1; j <= maxlevel; j++){
    double d = (L0/((double)(1 << j)))/zoh;
    double d2 = (L0/((double)(1 << j)))/zom;
    c[j] = (log(4.*d) - 1.)/(log(d) - 1.);
    lut[j] = sq(k/(log(d2) - 1.));
  }
}
/**
## Time integration

We time integrate the evolution equation of the atmospheric profiles:
 */
event diff (i++){
  /**
We allocate temporary fields to store the tendency terms and the
diffusivity.
   */
  scalar rx[],ry[],rb[];
  face vector kh[];
  double B = 0;
  double ws = 0;
  /**
The momentum is forced with a height-and-time constant horizontal
pressure gradient ($-\nabla P$) and is affected by back ground
rotation according to the Coriolis Parameter ($f$). Introducing,

$$\overrightarrow{U_{geo}} = \frac{\overrightarrow{k}}{\rho f}\times
\nabla P$$,

that is known as the geostrophic wind.
   */
  foreach(){
    rx[] = f*v[];
    ry[] = f*(Ugeo - u[]);
    rb[] = 0.;
    /**
In the bottom cell, we add the tendency due to the buoyancy flux.
     */
    if (x < Delta){
      B = (Qn + GFLX);
      rb[] += B / Delta;
      /**
For the momentum flux at the surface, we use a robust law-off-the-wall in
the lowest grid cell ($z+z_0 < \Delta$).

$$u(z) = \frac{u*}{\kappa}\mathrm{ln}\left(\frac{z}{z_{0.m}}\right)$$

with $u*$ the friction velocity, $\kappa=0.4$ the Von Karmann
constant, $\mathrm{ln}(x)$ is natural logarithm of a dummy variable
$x$ and $z_{0,m}$ the roughness length for momentum. We can express a
tendency in the lowest cell $[u]$ due to the surface friction:

$$u*^2 =
\frac{u}{\|u\|} \left( \frac{uk}{\mathrm{ln}
\left(\frac{\Delta}{z_{0,m}} \right)} \right)^2.$$

The corresponding prefactors are already computed and stored in the
`lut` array.
       */
      rx[] += -sign(u[])*lut[level]*sq(u[])/Delta;
      ry[] += -sign(v[])*lut[level]*sq(v[])/Delta;
    }
  }
  /**
The local eddy diffusivity ($K$) is written as,

$$K = lV_*$$

where $l$ is the mixing length,

$$l = \mathrm{min} \left[ kz, 70m \right],$$

with $k= 0.4$ the Von Karman constant. $V_*$ is the mixing velocity
scale due to shear *and* convection, inspired by Troen en Mahrt (1986)
we write,

$$V_* = \sqrt{w_c^2 + \left(lSf(Ri))^2}$$

with w_c the vertical velocity variance for convection ($B>0$ and $z < \mathcal{L}_c$), 

$$w_c =3 w_d z/h \left( 1-z/\mathcal{L}_c \right) ^ 2,$$

with $w_d$ the Deardorf velocity scale,

$$w_d = (B\mathcal{L}_c)^{1/3},$$ 

and $\mathcal{L}_c$ the height of the well mixed layer,

$$\mathcal{L}_c =\sqrt{\frac{2\int b - N^2z \mathrm{d}z}{N^2}}$$

for which we followed the works of Van Heerwaarden, Mellado etc. 

Noting that $w_c = 0$ for the times with $Q*(t) < 0$ and heights $z > \mathrm{L}_c$. 

The term $lSf(Ri)$ is computed conforming to the definitions in Van
Hooft et al. (2018b).
   */
  Lmix = 0.;
  foreach()
    Lmix += (b[] - x*sq(Nv)) * Delta;
  if (Lmix > 0.)
    Lmix = sqrt (Lmix*2./sq(Nv));
  if (B > 0 && Lmix > 0)
    ws = pow(B*Lmix, 1./3.);
  foreach_face(){
    double sqd = (sq((u[] - u[-1])/(Delta)) + sq((v[] - v[-1])/(Delta)));
    double Ri = ((b[] - b[-1])/(Delta))/(sqd + 0.00001);
    double fRi;
    if (Ri < 0)
      fRi = friu(Ri);
    else
      fRi = fris(Ri);
    double l = min(k*x, (70./1030.)*Lc); // The mixing length
    double fraction = 0;
    if (Qn > 0)
      fraction = 3.*(x/Lmix * sq(1. - x/Lmix))*(x < Lmix);
    double Vs = sqrt(fraction*sq(ws) + sq(l)*sqd*sq(fRi));
    kh.x[] = l*Vs;
  }
  boundary(all);
  dt = dtnext(DT);
  /**
Now we have all the ingredients for the reaction-diffusion problem. 
   */
  int n = 0;
  mgb = diffusion (u, dt, kh, r = rx);
  n += mgb.i;
  mgb = diffusion (v, dt, kh, r = ry);
  n += mgb.i;
  mgb = diffusion (b, dt, kh, r = rb);
  n += mgb.i;
  /**
Based on the convergence properties we adapt the timestep. 
   */
  if (n > 10)  //Quickly reduce the timestep if things get rough
    DT = max (DT/(1+((double)n/10.)), T/(24.*3600.));
  if (n < 5)   //Slowly increase the timestep when time integration is easy.
    DT = min (DT*(1+((double)n/100.)), 20.*T/(24.*3600.));
}
/**
The grid is adapted based on the wavelet-estimated error for the
discretized representation of the solution field. The criteria are
chosen to be $\zeta_{u,v} = U_{geo}/20$ and $\zeta_b =
b_{c,\Lambda}/50$. These values work well enough. 
 */
event adapt (i++){
  double ue = Ugeo/20.;
  double be =  sq(Nv)*Lc/50.;
  adapt_wavelet({u,v,b}, (double[]){ue, ue, be}, maxlevel);
}

/**
## Output

Output consists of instantanious porofiles that are outputted 15 times
per phyisical hour.
*/

event profile(t += T/(24*15)){
  char fname[99];
  sprintf (fname, "ProfsPi5%g", Pi5);
  static FILE * fp = fopen (fname, "w");
  foreach(){
    if (x < Delta)
      fprintf (fp, "0 0 0 %g \n", BSURF);
    fprintf (fp, "%g %g %g %g \n", x, u[], v[], b[]);
  }
}
/**
We also monitor the occurance of a low level jet. The numerical values
are not very robust for the various mixing function ($f(\mathrm{Ri})$)
formulations, and hence not presented in the associated work.
 */
double xuvb[4][1000] = {0};
double trec = 0;
event llj (t = T/2; i += 5){
  bool new_record = false;
  foreach()
    if (sq(u[]) + sq(v[]) > sq(Umax))
      new_record = true;
  if (new_record){
    trec = t;
    int j = 0;
    foreach(){
      xuvb[0][j] = x;
      xuvb[1][j] = u[];
      xuvb[2][j] = v[];
      xuvb[3][j++] = b[];
    }
    xuvb[0][j] = -1.;
  }
}   
/**
Each timestep, we monitor the buoyancy difference (i.e. the inversion)
between $z = 0$ and $z = L_c/20$ and the windspeed at that
height. These are summed up over the 19th hour of simulation.
 */
event inversion (t = T*18./24.; t<= T*19./24.; i++){
  double zp1 = 0.;
  double zp2 = Lc / 20;
  inv += interpolate (b, zp2) - interpolate (b, zp1);
  u40 += sqrt(sq(interpolate (u, zp2)) + sq(interpolate (v, zp2)));
  ni++;
}
/**
Each ten timesteps we output some interesting solution diagnostics so
that we can see how the solution has evolved over time. We output
($Q*(t), G(t), B(t), b_{surf}, b_{mix}$ and $\mathcal{L}_c$ to a file
named: `timeseriesPi5` followed by a $\Pi_5$-value identifier.
 */
event timeseries (i += 10){
  char fname[99];
  sprintf (fname, "timeseriesPi5%g", Pi5);
  static FILE * fp = fopen (fname, "w");
  double B, G, bs, bmix;
  Lmix = 0;
  foreach()
    Lmix += (b[] - x*sq(Nv)) * Delta;
  if (Lmix > 0.)
    Lmix = sqrt (Lmix*2./sq(Nv));
bmix = Lmix * sq(Nv);
foreach_boundary (left){ 
    G = GFLX;
    B = Qn + G;
    bs = BSURF; 
  }
  if (i == 0)
    fprintf (fp, "t\tQn\tG\tB\tbs\tbmix\tLmix\n");
  fprintf (fp, "%g\t%g\t%g\t%g\t%g\t%g\t%g\n",
	  t, Qn, G, B, bs, bmix, Lmix); 
}

/**
## The last event

The `run()` stops once $t = T$ and then we output the (iteration)
averaged inversion data and the low-level jet data.
 */
event stop (t = T){
  int m = 0;
  char fname[99];
  sprintf (fname, "Profjet%g",Pi5);
  static FILE * fpj = fopen (fname, "w");
  while (xuvb[0][m] != -1.){
    fprintf(fpj, "%g\t%g\t%g\t%g\n",
	    xuvb[0][m], xuvb[1][m],xuvb[2][m], xuvb[3][m]);
    m++;
  }
  inv /= (double)ni;
  u40 /= (double)ni;
  fprintf (fpi, "%g\t%g\t%g\n", Pi5, inv, u40);
}

/**
## A movie

On systems with `gnuplot` and `ffmpeg` installed, we can choose to
generate animations of the simulation results:

![](lumpedscm/mov_scm.mp4)(width="800" height="400")
 */
#define GNUPLOT_AND_FFMPEG 1 //Switch for movie output
#if GNUPLOT_AND_FFMPEG
/**
We initialize a pipeline for the plots that will later be turned into a movie.
 */
static FILE * gnuplotPipe;
event init(t = 0){
  gnuplotPipe = popen ("gnuplot", "w");
  fprintf (gnuplotPipe,
           "set term pngcairo size 1200,600 enhanced font 'Times ,18'\n"
           "set yr [0: 1200]\n"
           "set ylabel 'height [m]'\n"
           "set grid\n"
           "set size square\n");
}
/**
We output 15 plots per physical hour. 
 */
int frame = 0;
event plot (t += T/(24*15)){
  if (fabs(Pi5 - round(Pi5)) < 0.001){
    fprintf (gnuplotPipe, "set output 'plot%d.png'\n", frame);
    fprintf (gnuplotPipe, "set multiplot layout 1,2 title 'Π_5 = %g,"
             "time %.02d:%.02d' font 'Times ,25'\n",
	    Pi5, (int)(t/(3600)), ((int)t%3600)/60);
    fprintf (gnuplotPipe,
             "set xr [-15: 20]\n"
             "set xlabel 'wind speed [m/s]'\n"
             "set key top left\n"
	    );
    fprintf (gnuplotPipe, "plot '-' w l lw 5 lt rgb'#11BB11' t 'u',"
             "'-' w l lw 5 lt rgb '#BB11BB' t 'v'\n");
    foreach()
      fprintf (gnuplotPipe, "%g %g\n",u[], x);
    fprintf (gnuplotPipe, "e\n");
    foreach()
      fprintf (gnuplotPipe, "%g %g\n",v[], x);
    fprintf (gnuplotPipe, "e\n");
    fprintf (gnuplotPipe,
             "set xr [260 : 290]\n"
             "set key off\n"
             "set xlabel 'θ [K]'\n"
            );
    fprintf (gnuplotPipe, "plot '-' w l lw 5 lt rgb '#CC1111'\n");
    /**
       Note that the movie displays the potential temperature ($\theta$)
       rather than the buoyancy, using $\theta_{ref} = 270K$ and $g = 10
ms^{-2}$. 
     */
    foreach()
      fprintf (gnuplotPipe, "%g %g\n", b[]*27 + 270, x);
    fprintf (gnuplotPipe, "e\n"
             "unset multiplot\n");
    fflush (gnuplotPipe);
    frame++;
  }
}
/**
Finally, we render a `.mp4` movie from these plots and then remove the plots
from the disk.
 */
event moviemaker(t = T){
  if (fabs(Pi5 - 8.) < 0.0001){ 
    system ("rm mov_scm.mp4");
    system ("ffmpeg -r 25 -f image2 -i plot%d.png -c:v libx264 -vf format=yuv420p -y mov_scm.mp4");
    system ("rm plot*");
  }
  return 1;
}
#endif