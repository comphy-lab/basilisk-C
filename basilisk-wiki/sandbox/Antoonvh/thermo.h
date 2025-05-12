/**
# Moist thermodynamics

This page implements moist atmospheric thermodynamics, suitable for
turbulence-resolving boundarly-layer simulations. In order to simulate
clouds in the atmosphere we adopt a *model* for phase changes. Notice
this this happens entirely on the subgrid scale. We take inspirtation
from the formulations of
[DALES](https://github.com/dalesteam/)
and [Micro-HH](www.microhh.org).

A bunch of global constants are defined. Furthermore, fields for the
liquid water potential temperature ($\theta_v$, `thl` in Kelvin), the
total water specific humidity ($q_t$, `qt` in kg/kg) and an
approximate thermodynamic pressure ($p$, `pres` in Pascal) are
declared.
 */
double Rv = 461.5;        //in J/kg/K
double Rd = 287.;         //in J/kg/K
double L = 2500000.;      //in J/kg
double cpd = 1004.;       //in J/Kg/Kq
#define P00 100000.       //A constant reference pressure
double P0 = 101780., Pe = 7500;//Surface pressure and 1/e length scale
scalar thl[], qt[], pres[];
pres[bottom] = dirichlet (P0);
#ifndef RHO0            
#define RHO0 1.3 //Base density in kg/m^3  
#endif
double rho0 = RHO0;

#if dimension == 1
#define HEIGHT_ABOVE_SURF (x - X0)  
#else
#define HEIGHT_ABOVE_SURF (y - Y0)
#endif

/**
   A sequence of `DEFINE`s is used to compute the cloud liquid water
   content as a function of height and temperature. Note that there is
   no field dedicated to store the liquid water specific humidity.
*/
#define PRES  (P0*(exp(-HEIGHT_ABOVE_SURF/Pe)))   
#define EXNER (pow(pres[]/P00, Rd/cpd))
#define TK1   (thl[]*EXNER)
#define ES    (610.78*exp(17.27*(TK1 - 273.16)/(TK1 - 35.86)))
#define QSL   (Rd/Rv*(ES/(pres[] - (1. - Rd/Rv)*ES)))
#define QSAT  (QSL*((1. + ((sq(L)/(Rv*cpd*sq(TK1)))*qt[]))/(1 + ((sq(L)/(Rv*cpd*sq(TK1)))*QSL))))
#define QC    (max(qt[] - QSAT, 0))
/**
   The virtual potential temperature in a cell ($\theta_v$,`THV`) can
   be computed from $\theta_l$, $q_c$ and $q_t$.
 */
#define THETA (thl[] + L*QC/(cpd*EXNER)) 
#define THV (THETA*(1. - (1. - Rv/Rd)*qt[] - Rv/Rd*QC))
/**
   The virtual potential temperature is related to the buoyancy,
 */
double g_acc = 9.81, T_ref = 280.;
#define B_CENTERED (g_acc*(THV - T_ref)/T_ref)

/**
## The thermodynamic pressure

We need to have a quess at the thermodynamic pressure. It can
integrated from the surface pressure (`P0`) with height according to:

$$\frac{\partial \Pi}{\partial z} = -\frac{g}{c_p\theta_{v}}$$
with,
$$\Pi = \left(\frac{p}{p_0} \right)^{\frac{Rd}{c_p}}$$

Notice that this pressure definition is not modified by the flow. The
assumption being that $g$ is far greater than any acceleration in the
flow field. We use a multigrid solver for this purpose.
 */
#include "integrator.h"
struct STP {
  bool reset;   //Use a guess for the pressure
  double TH;    //The angle (in radians) of the vertical direction
};
int set_thermodynamic_pressure (struct STP stp) {
  if (!stp.TH)
    THETA_ANGLE = -pi/2.; //`Bottom` is the surface
  else
    THETA_ANGLE = stp.TH;
  scalar rhs[], _pres[], Pi[]; 
  double TOL = 0.1;  //Pa
  Pi[bottom] = dirichlet (pow(P0/P00, Rd/cpd));
  /**
     Notice that computing the rhs requires the pressure
itself. Therefore, we must call the iterative solver iteratively and
hope that the pressure converges. As a fist guess we use the current
values in `pres[]`, unless the argument `reset` is true, then we take
a typical profile. It is not wise to call this function very often.
  */
  foreach() {
    _pres[] = PRES;
    if (stp.reset)
      pres[]  = _pres[];
    Pi[] = EXNER;
  }
  int j = 0;
  do {
    foreach()
      rhs[] = -g_acc/(cpd*THV);
    integrate_dx (Pi, rhs);
    j++;
    foreach()
      pres[] = pow((Pi[]), cpd/Rd)*P00;
  } while (change (pres, _pres) > TOL && j < NITERMAX);
  if (j == NITERMAX)
    fprintf(ferr, "The thermodynamic pressure was not found within %d cycles\n", j);
  return j;
}
/**
   Alternatively, one may opt for the `set_pres ()` functionality
   (tip!). 
*/

#include "profile5c.h"
struct SP {
  int C;      // Levels
  bool guess; // Use a guess
};

void set_pres (struct SP sp) {
  int C = (1 << grid->maxdepth);
  if (sp.C)
    C = sp.C;
  double dz = 0.99999*L0/(double)(C - 1); // almost span the entire domain
  /**
     First we compute the $\theta_v$ field so we can find its
     horizontally averaged values.
  */
  scalar thv[];
  foreach() {
    if (sp.guess)
      pres[] = PRES;
    thv[] = THV;
  }
  boundary ({thv}); 
  /**
     A 1D look-up table is created for the pressure (`Y`), computed
     via $\theta_{v0}$, as described in van Heerwaarden et al. (2017).
  */
  double Y[C][3]; //y, thl, p
  memset (Y, 0, sizeof(Y));
  Y[0][0] = Y0 + SEPS;
  Y[0][2] = P0;
  average_over_yp ({thv}, &Y[0][1], Y[0][0]);
  for (int n = 1; n < C; n++) {
    Y[n][0] = Y[n - 1][0] + dz;
    average_over_yp ({thv}, &Y[n][1], Y[n][0]);
    if (pid() == 0)
      Y[n][2] = Y[n - 1][2]*exp(-g_acc*2.*dz/(Rd*pow(Y[n-1][2]/P00, Rd/cpd)*
					      (Y[n][1] + Y[n - 1][1])));
  }
  /**
In parallel, the root broadcasts its data.
   */
#if _MPI
  MPI_Bcast (&Y[0][0], 3*C, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif
  /**
     Each cell can find the corresponding pressure via interpolation
  */
  foreach() {
    double H = HEIGHT_ABOVE_SURF;
    int n =  (int)(H/dz) + 1;
    pres[] = ((H - Y[n - 1][0])*Y[n][2] +
	      (Y[n][0] - H)*Y[n - 1][2])/dz;
  }
  boundary ({pres});  
}

/**
## Moist dynamics
   
   Next, we couple the thermodynamics to the flow. That means $q_t$
   and $\theta_l$ become active tracers. It can optionally be swithed
   off by compiling with a `-DNO_DYNAMICS` tag.
*/
#ifndef NO_DYNAMICS
#include "tracer.h"
scalar * tracers = {qt, thl};

event defaults (i = 0) {
  if (is_constant(a.x)) {
    a = new face vector;
    foreach_face()
      a.x[] = 0;
    boundary ((scalar*){a});
  }
  pres.prolongation = refine_linear;
#if TREE
  pres.refine = refine_linear;          //conservative 3rd order vertical interp.
#endif
  thl.gradient = qt.gradient = minmod2; //A slope limiter
}
/**
   The effect of gravity on the flow is described via the buoyancy.
*/

event acceleration (i++) {
  face vector av = a;
  scalar buoy[];
  foreach()
    if (cm[] > 0)  //?
      buoy[] = B_CENTERED;
  boundary ({buoy});
  foreach_face(y) 
    av.y[] += fm.y[]*(buoy[0,-1] + buoy[])/2.;
}

/**
   In order to find the modified pressures without triggering
   spurrious currents, a small tolerance is used for a brief while.
*/
int back_it = 10;
double DEF_TOL = 1e-3, DEF_DT = 60.;

event init (t = 0) {
  DT = min(DT, 1.); 
  TOLERANCE = min(TOLERANCE, 1e-6);
}

event back_to_defaults (i = back_it) {
  TOLERANCE = DEF_TOL;
  DT = DEF_DT; 
}

#endif //NO_DYNAMICS
/**
## extras 

   A blue-to-white color bar for cloud visuals
 */
void cloud (double cmap[NCMAP][3]) {
  for (int i = 0; i < NCMAP; i++) {
    cmap[i][0]= 0.2 + 0.8 * (double)i/(double)NCMAP;
    cmap[i][1]= 0.2 + 0.8 * (double)i/(double)NCMAP;
    cmap[i][2]= 0.9 + 0.1 * (double)i/(double)NCMAP; 
  }
}
/**
## Usage

* [Shallow cumulus convection](bomex.c)  
* [A rising moist bubble](moist_bubble.c)

## Reference

Heerwaarden, C. C. V., Van Stratum, B. J., Heus, T., Gibbs, J. A.,
Fedorovich, E., & Mellado, J. P. (2017). MicroHH 1.0: a computational
fluid dynamics code for direct numerical simulation and large-eddy
simulation of atmospheric boundary layer flows. Geoscientific Model
Development, 10(8), 3145-3165.
 */
