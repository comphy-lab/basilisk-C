/**
# Large-eddy simulation of a diunral cycle

This page discusses the setup and inplementation for the large-eddy
simulation (LES) of a idealized diunral cycle of the dry atmospheric
boundary layer (ABL). We use an adaptive octree grid, the
Navier-Stokes solver, a tracer and the subgrid-scale (SGS) model for
turbulent mixing of Vreman (2004). Furthermore we will `view` our
results.

This setup enherites quite some features of the single-column model
(SCM), that was run for the same case.
 */
#include "grid/octree.h" //<- Uncomment for *the* 3D simulation
#include "navier-stokes/centered.h"
#include "tracer.h"
#include "diffusion.h"
#include "SGS.h"
#include "view.h"
/**
   A macro is defined to find $b_{surf}$ (see SCM page)
 */
#define BSURF ((b[0,1] - b[]*c[level])/(1. - c[level]))
/**
The buoyancy field ($b$) is an *active* tracer field.
*/			    
scalar b[], * tracers = {b};
/**
To ensure no spurious fluxes at the top and bottom boundaries, we set
a boundary condition for the eddy viscosity.
 */
Evis[bottom] = dirichlet (0.);
Evis[top] = dirichlet (0.);
/**
The bottom buoyancy flux is computed as the sum of $G$ and $Q^*$ (see
the SCM page)
 */
#define Gbflx (Lambda*(beq - BSURF))
#define Qn (max (B0*sin(2.*M_PI*t/T), B1))
/**
A damping-layer tendency can be defined for a quantity `s`.  
 */
double S  = 0.1;
#define damping(s) (-S*min( sq(y - (L0*0.35)), 1.)*(y >(L0*0.35))*(s))
/** 
   The values for the *dimensionless* parameters are listed below;
*/
double Pi1 = -6.;
double Pi2 = 2160.;
double Pi3 = 10.;  
double Pi4 = 5366.;
double Pi5 = 0.5;  
double Pi6 = 5000.;
double Pi7 = 1.;  
double L0Lc = 3.;       // L0/Lc
/**
   Two dimensional parameters give rise to "lattice units" or
   simulation units. We choose a normalized length scale for $L_c$ and
   a cycle time-scale of 24*60=1440 units, i.e. inspired by the
   (dimensionless) number of minutes in a 24 hour period.
*/
double Lc  = 1.;      // A normalized length scale.
double T   = 1440.;   // = 24*60 time units in a cycle
double beq = 0.;      // Equal to $theta_{ref}$
double k = 0.4;       // Von Karmann constant
double lut[19], c[19];// Loop-up tables

/**
## Adaptive refinement parameters

These values should not affect the dynamics. However they do affect
the absolute values of all physical fields and how they are
represented in the simulation. Hence the refinement criteria should
be defined consistently.
*/ 
double be, ue;           // The criteria for buoyancy and the velocities
double uemin, bemin;     // Minimum for u and b
double fracu    = 2.;    // Fraction of the convective velocity scale
double fracb    = 1.;    // Fraction of the convective buoyancy scale
int    maxlevel = 8;    
/**
   These global variables are usefull and their values depend on the
   values provided above and are calculated in the "main" function.
*/

double B0, B1, Nb, Lambda, Ugeo, f, zom, zoh, zi;
/**
The buoyancy will be linked to the acceleration via the `av` face vector field.
*/
face vector av[];

int main(){
  /**
For the pressure field a third-order-accurate prolongation and
refinment attribute for grid-alligned variations are set. For the
momentum, we choose conservative refinement and the gradients in the
buoyancy field are computed with the generalized `mindmod2` slope
limiter, this helps with advection of sharp buoyancy fronts.
   */
  p.prolongation = p.refine = refine_linear;
  foreach_dimension()
    u.x.refine = refine_linear;
  b.gradient = minmod;
  /**
Next we compute the physical parameters from the dimensionless groups,
$L_c$ and $T$ (see the SCM page).
  */
  Nb = Pi2/T;
  f = Pi3/T;
  B0 = sq(Lc*Nb)*M_PI/(2.*T);
  B1 = B0/Pi1;
  Lambda = sqrt(B0*T)/Pi4;
  Ugeo = Pi5*pow(B0*Lc, 1./3.);
  zom = Lc / Pi6;
  zoh = zom * Pi7;
  L0 = L0Lc*Lc;
  /**
Aditionally, we express the minimum refinment criertia values as
$\zeta_{u,\mathrm{min}}= U_c/10$ and $\zeta_{b,\mathrm{min}}= b_c /
5$. The boundary conditions at the surface are consistent with
no-slip, and we enforce the initialized buoyancy value at the top.
   */
  uemin = pow(B0*Lc, 1./3.)/10.;
  bemin = pow(sq(B0)/Lc, 1./3.)/5.;
  ue = uemin;
  be = bemin;
  Y0 = 0;
  X0 = -L0/2;
  zi = Lc;
  periodic (left);
  a=av;
  u.t[bottom] = dirichlet (0.);
  b[top] = dirichlet (y*sq(Nb));
#if (dimension == 3)
  u.r[bottom] = dirichlet (0.);
  periodic (back);
  Z0 = X0;
#endif
  /**
For future reference, we save some parameter values of this `run()` to
a file called `systemparameters`.
   */
  FILE * fpstart = fopen ("systemparameters","w");
  fprintf (fpstart, "Pi1=%g\nPi2=%g \nPi4=%g \nPi6=%g \nPi7=%g\n "	\
	   "Nb=%g \nB0=%g \nB1=%g \nLambda=%g Ug =%g \nml = %d\n",	\
	   Pi1, Pi2, Pi4, Pi6, Pi7, Nb, B0, B1, Lambda, Ugeo, maxlevel);
  fprintf (fpstart,"zi = %g and %g \nbcq=%g\n bs = %g \nsqNb=%g\n" ,	\
	   Lc, pow((2.*B0*T/(M_PI*sq(Nb))), 0.5), pow((2*B0*T*sq(Nb)/M_PI), 0.5),
	   B1/Lambda, sq(Nb));
  fflush (fpstart);
  fclose (fpstart);
  N = 1 << 4;
  run();
}
/**
## Initialization

This event takes care of the initialization:
 */
event init(t=0){
  /**
Since we are not bothered to initialize a pressure, the solver will
have to find it in a short while. As such a small initial timestep
(`DT`) and tolerance on the residuals for the Poisson problems are
used.
   */
  TOLERANCE = 10E-5;
  DT = 10E-3;         //Small values
  /**
If possible, we restore for a dumpfile. The name of file can be
changed manually.
   */
  if (!restore("dumpt=390")){ //e.g. 
    /**
Else we refine the grid close to the surface and initialize mean
profiles.
     */
    refine (level < (maxlevel - 1) && y < (Lc/10.));
    refine (level < maxlevel && y < (Lc/20.));
    foreach(){
      b[] = sq(Nb)*y + (noise()*be/5.);
      u.y[] = noise()*ue/10.;
      u.z[] = Ugeo;
    }
  }
  boundary (all);
  lut[0] = 0.;
  /**
## Look-up tables

We also make a lookup table for the aerodynamic drag coefficient ($C_D$) at the wall, 

$$u_{*,i}^2 = C_D * u_i^2,$$

Using the law of the wall 

$$u = \frac{u_*}{k} \mathrm{log}\left( \frac{z}{z_0} \right) .$$

*And* a look-up table for the logaritmic extrapolation
coefficients of the buoyancy profile. (see the SCM page)
  */
  for (int j = 1; j <= maxlevel; j++){
    double d = (L0/((double)(1 << j)))/zoh;
    double d2 = (L0/((double)(1 << j)))/zom;
    c[j] = (log(4.*d) - 1.)/(log(d) - 1.);
    lut[j] = sq(k/log(d2));
  }
}

event tracer_diffusion (i++){
  /**
     The heatflux tendencies at the surface and within the buoyancy
     damping layer are time integrated using a rather unsophisticated
     forward-Euler method. Noting that since $\Pi_2 >> 1$, the
     evolution of the surface buoyancy flux is slow (time scale $T$)
     compared to the evolution of the flow at the discretization time
     scale (via CFL-criterion).
  */
  foreach(){
    b[] += dt*damping (b[] - sq(Nb)*y);
    if (y< Delta)
      b[] += dt*(Gbflx + Qn)/Delta;     // <-- This is the surface buoyancy flux
  }
}

/**
## Momentum forcing

Here we implement

* The acceleration of gravity is implemented with a buoyancy
   formulation.  
* The so-called "sponge layer".  
* The acceleration due to background rotation  
* The tendencies due to the surface shear stress, based on a law of the wall.  
*/

event acceleration (i++){
  /**
Notice how in the associated work, the coordinates referred to as
$x,y,z$ are substitued in the simulation by `z,x,y` and as such;
$\{u,v,w\} \rightarrow \{$`u.z ,u.x, u.z`$\}$.  respectively. This
curiosity is rooted in the fact that Basilisk's `y`-coordinate is
associated with the `bottom` and `top` boundary, and a corresponding
arbitrary choice.
   */
  foreach_face(x) //v
    av.x[] += f*(Ugeo - ((u.z[] + u.z[-1])/2.)) + damping((u.x[] + u.x[-1])/2.); 
  foreach_face(z) //u
    av.z[] += f*((u.x[] + u.x[0,0,-1])/2.) +  damping(((u.z[] + u.z[0,0,-1])/2. - Ugeo));
  foreach_face(y) //w
    av.y[]+=(gr.y*(b[]+b[0,-1])/2.) + damping((u.y[]+u.y[0,-1])/2.);
  foreach(){
    if (y + Y0 < Delta){
      u.x[] -= dt*sign(u.x[])*lut[level]*sq(u.x[])/Delta;
      u.z[] -= dt*sign(u.z[])*lut[level]*sq(u.z[])/Delta;
    }
  }
}
/**
   ## Adaptivity

   Due to the large degree of scale separation over the course of the
   simulation run, it seems attractive to use an adatpive grid
   approach. The grid adaptation strategy is based on a
   multi-resolution/wavelet analysis of the solution structure at the
   grid level. This may be viewed as a "feature detection algorithm". 
*/
double Lc, Uc;
event adapt (i++){
  /**
The small `TOLERANCE` is relaxed over time.
  */
  TOLERANCE = min (TOLERANCE*1.01, 1E-3);
  /**
     Since the simulation evolves over time, we use a dynamic meassure
     to express the refinement criteron (adaptive adaptation). For
     that purpose the convective velocity-fluctuation scale
     ($\mathcal{U}_c$) and the convective buoyancy-fluctuation scale:
     $b_c$ according to:

$$\mathcal{U}_c = \left( \mathrm{L}_c B \right) ^{1/3}, $$

$$b_c = \left(\frac{B}{L_c}\right)^{1/3},$$

with $\mathcal{L}_c$ the height of the well-mixed layer according to,

$$\mathcal{L}_c^2 = \frac{2}{N^2}\int \langle b\rangle -N^2z  \mathrm{d}z.$$

This assumes a positive $B$, for the stable part of the simulation we
default to the aforementioned minimal values for the refinement
criteria.
  */
  double uen    = 0.;
  double ben    = 0.;
  double H      = 0;
  double mxdbdz = 0;
  double bint   = 0;
  foreach (reduction(max:mxdbdz) reduction(+:H) reduction(+:bint)){
    bint += (b[] - (sq(Nb)*y))*dv();
    if ((b[0,1] - b[])/(Delta) > mxdbdz)
      mxdbdz = (b[0,1] - b[])/(Delta);
    if (y<Delta){
      double m = sq(Delta);
      H += (Qn + Gbflx) * m;
    }
  }
  DT = min (DT*1.02, 1./(Nb*((double)(1<<5))));
  Lc = 0.001;
  if (H > 0 && bint > 0){
    Lc = (1./(L0*pow(sq(Nb),0.5)))*pow(2*bint, 0.5);
    Uc = pow((H/sq(L0))*Lc, 1./3.);
    uen = max (Uc/fracu, uemin);
    double Bs = pow(sq(H/(sq(L0)))/Lc, 1./3.);
    ben = max (Bs/fracb, bemin);
  } else {
    uen = uemin;
    ben = bemin;
  }
  /**
To actual refinement criterion has a small history. 
   */
  ue = (ue + uen)/2.;
  be = (be + ben)/2.;
  fprintf (stderr, "ue = %g be = %g DT = %g, Lc = %g, 1/sqN*8 = %g\n",
	  ue, be, DT, Lc, 1./(sqrt(mxdbdz)*((double)(1<<3))));
  /**
An `if` statement ensures that the grid is not adapted to
be too coarse in the first few physical minutes before the atmospheric
dynamics occur...
   */
  if (t > 60./((double)(maxlevel - 4)))
    adapt_wavelet ((scalar*){b, u}, (double[]){be, ue, ue, ue}, maxlevel);
  fflush (stdout);
}
/**
## The last event

The simulation stops once $t = T$.
*/
event end (t = T);
/**
## Output

The main output features:

* Simulation dump files for post processing and simulation restarts.  
* Data from a virtual observation tower.  
* Vertical profiles (i.e. horizontal averages) of the solution fields.     
* A movie.  
* Time-serie data of domain-integrated and/or surface-integrated quantities.  

For a cleaner appearance of this page, it was chosen to separate the  
burden of defining the case setup from the output routines between  
different files. As such, the output routines are presented in a  
header file (click the link for details):  
 */
#include "diagnostics.h"
