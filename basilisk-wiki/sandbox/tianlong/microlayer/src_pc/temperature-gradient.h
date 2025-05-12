/**
# Temperature Gradient Phase Change Model

This phase change model is suitable for boiling conditions.
The vaporization rate is computed from an energy balance at
the interface (in every interfacial cell), assuming that
the interface is at saturation temperature and that the
heat conduction to the interface is used for the phase change
phenomena. Therefore, the vaporization rate is computed from
the following balance:

$$
\dot{m}_{Evap} \Delta h_{ev} =
    \mathbf{\dot{q}}_l \cdot \mathbf{n_\Gamma}
  + \mathbf{\dot{q}}_g \cdot \mathbf{n_\Gamma} =
  - \lambda_l \left.\dfrac{\partial T_l}{\partial \mathbf{n}_\Gamma}\right\vert_l
  - \lambda_g \left.\dfrac{\partial T_g}{\partial \mathbf{n}_\Gamma}\right\vert_g
$$

where $\lambda$ is the thermal conductivity, $\Delta h_{ev}$
is the enthalpy of evaporation and $\mathbf{n_\Gamma}$ is the
normal to the interface.
After the calculation of the vaporization rate, this model
solves the transport equation for the temperature field $T$
using a two-field formulation and including the presence of
the phase change, both in the diffusion equation and in the
transport of the temperature field which keeps into account
the Stefan convection.
*/

#include "intgrad.h"
#include "fracface.h"
#include "mydiffusion.h"
/**
## Memory Allocations

This phase change model defines a list of vaporization
rates *mEvapList* which is populated by a single scalar
*mEvap* because a single contribution to the vaporization
rate is considered.
*/

#ifndef PHASECHANGE
scalar mEvap[];
scalar * mEvapList = {mEvap};

/**
*fL* and *fG* store the value of volume fractions for the
calculation of the interface gradients, while the vectors
*fsL* and *fsG* contain the face fraction fields computed
using the [fracface.h](/sandbox/ecipriano/src/fracface.h)
module. */

scalar fL[], fG[], f0[];
face vector fsL[], fsG[];

face vector s[];
# define PHASECHANGE
#endif

/**
## Field Allocations

* *T* one-field temperature
* *TL* liquid-phase temperature field
* *TG* gas-phase temperature field
* *TInt* value of temperature at the interface
*/

scalar T[], TL[], TG[], TInt[];

//for small cut cells, the temperature needs to be fixed for numerical stablity
double fix_small_thd = 0.1;
void fixSmallCellsTemp(scalar TL, scalar TG);
void fixSmallCellsTempLevel(scalar TL, scalar TG, int l);

#if USE_CONJUGATE_HEAT
  scalar TS[];
#endif

/**
## User Data

Using this phase change model, the user should define the
following variables (SI units):

* *lambda1* Thermal conductivity in liquid phase
* *lambda2* Thermal conductivity in gas phase
* *dhev* Latent heat of vaporization
* *cp1* Specific heat capacity in liquid phase
* *cp2* Specific heat capacity in gas phase
* *TIntVal* Interface temperature value
*/

extern double lambda1, lambda2, dhev, cp1, cp2;
extern double TIntVal;
#if USE_CONJUGATE_HEAT
extern double lambdas, rhos, cps;
extern face vector is_solid_x_face;
extern face vector is_solid_y_face;
mgstats mgT;
#endif

#if (USE_CONJUGATE_HEAT)
extern scalar is_solid_y, is_solid_x, vof_heater;
extern double lambdas2, rhos2, cps2;
#endif

/**
Allocating diffusivity fields *lambdaf*, volume correction
*theta* and explicit source terms *sgT* and *sTS* for
the diffusion equation. */

face vector lambda1f[], lambda2f[];
scalar thetacorr1[], thetacorr2[];
scalar sgT[], slT[], sgTimp[], slTimp[];

#if USE_CONJUGATE_HEAT
face vector lambdasf[];
face vector thermocoefsf[];
scalar thetasolid[];
scalar betal[], betag[], betas[];
scalar ssT[];
extern double Rcc;
extern double (*delta_q)(double x, double t);
extern double (*delta_q_volume)(double x, double t);
void imposeTempBoundaryConjugate(scalar TL, scalar TG, scalar TS, bool with_Rcc);
extern const double tshift;
#endif

foreach_dimension()
static double mTempgradient_x (Point point, scalar c, scalar t)
{
  static const double cmin = 0.5;
  double cl = c[-1], cc = c[], cr = c[1];
  if (t.inverse)
    cl = 1. - cl, cc = 1. - cc, cr = 1. - cr;
  if (cc >= cmin && t.gradient != zero)
  {
    if (cr >= cmin)
    {
      if (cl >= cmin)
      {
        if (t.gradient)
          return t.gradient(t[-1] / cl, t[] / cc, t[1] / cr) / Delta;
        else
          return (t[1] / cr - t[-1] / cl) / (2. * Delta);
      }
      else
        return (t[1] / cr - t[] / cc) / Delta;
    }
    else if (cl >= cmin)
      return (t[] / cc - t[-1] / cl) / Delta;
  }
  return 0.;
}

static void myTemprefine (Point point, scalar s)
{
  scalar f = s.c;
  if (cm[] == 0. || (!s.inverse && f[] <= 0.) || (s.inverse && f[] >= 1.))
    foreach_child()
        s[] = 0.;
  else
  {
    coord g;
    foreach_dimension()
        g.x = Delta * mTempgradient_x(point, f, s);
    double sc = s.inverse ? s[] / (1. - f[]) : s[] / f[], cmc = 4. * cm[];
    foreach_child()
    {
      s[] = sc;
      foreach_dimension()
          s[] += child.x * g.x * cm[-child.x] / cmc;
      s[] *= s.inverse ? 1. - f[] : f[];
    }
  }
}

/**
## Defaults

In the defaults event we setup the tracer lists for the
advection of the temperature fields. */

event defaults (i = 0)
{
  TL.inverse = false;
  TG.inverse = true;

#ifdef CONSISTENTPHASE1
  fuext.tracers = list_append (fuext.tracers, TL);
#else
  fu.tracers = list_append (fu.tracers, TL);
#endif
#ifdef CONSISTENTPHASE2
  fuext.tracers = list_append (fuext.tracers, TG);
#else
  fu.tracers = list_append (fu.tracers, TG);
#endif

  /**
  On adaptive meshes, tracers need to use linear interpolation (rather
  than the default bilinear interpolation) to ensure conservation when
  refining cells. */

#if TREE
#if EMBED
      TL.refine = TL.prolongation = refine_embed_linear;
      TG.refine = TG.prolongation = refine_embed_linear;
#else
      TL.refine  = refine_linear;
      TG.refine  = refine_linear;
#endif
      TL.restriction = restriction_volume_average;
      TL.dirty = true; // boundary conditions need to be updated
      TG.restriction = restriction_volume_average;
      TG.dirty = true; // boundary conditions need to be updated
#endif

  scalar *myinterfaces = {fu, fuext};
  T.depends = list_add (T.depends, fu);
  T.depends = list_add (T.depends, fuext);
  for (scalar c in myinterfaces)
  {
    c.refine = c.prolongation = fraction_refine;
    c.dirty = true;
    scalar *tracers = c.tracers;
    for (scalar t in tracers)
    {
      t.depends = list_add (t.depends, T);
      t.restriction = restriction_volume_average;
      t.refine = t.prolongation = myTemprefine;
      t.dirty = true;
      t.c = c;
    }
  }
}

/**
## Init

In the init event, we avoid dumping all the fields that we
don't need to visualize. */

event init (i = 0)
{
  sgT.nodump = true;
  slT.nodump = true;
  fL.nodump  = true;
  fG.nodump  = true;
  thetacorr1.nodump = true;
  thetacorr2.nodump = true;
#if USE_CONJUGATE_HEAT
  thetasolid.nodump = true;
  ssT.nodump = true;
  betal.nodump = true;
  betag.nodump = true;
  betas.nodump = true;
#endif
  f0.nodump = true;
  slTimp.nodump = true;
  sgTimp.nodump = true;
  TInt.nodump = true;
  mEvap.nodump = true;
  stefanflow.nodump = true;
  stefanflowext.nodump = true;
}

/**
## Finalise

We deallocate the various lists from the memory. */

event cleanup (t = end)
{
  free (fu.tracers), fu.tracers = NULL;
  free (fuext.tracers), fuext.tracers = NULL;
}

/**
## Phase Change

In the *phasechange* event, the vaporization rate is computed
and the diffusion step for the mass fraction field (in liquid
and gas phase) is solved. */

event phasechange (i++)
{
  /**
  First we compute the value of the non volume-averaged
  temperature fields. This procedure allows a better
  calculation of the gradients close to the interface. */

  foreach() {
    f[] = clamp (f[], 0., 1.);
    f[] = (f[] > F_ERR) ? f[] : 0.;
    f0[] = f[];
    fL[] = f[]; fG[] = 1. - f[];
    TL[] = f[] > F_ERR ? TL[]/f[] : TIntVal;
    TG[] = ((1. - f[]) > F_ERR) ? TG[]/(1. - f[]) : TIntVal;
  }
  //boundary({fL,fG,TL,TG,f0});

  fixSmallCellsTemp(TL, TG);
#if USE_CONJUGATE_HEAT
  imposeTempBoundaryConjugate(TL, TG, TS, true);
#endif
  /**
  We compute the value of volume fraction *f* on the
  cell-faces using a geometric approach (necessary
  for interface gradients and diffusion equations). */

  face_fraction (fL, fsL);
  face_fraction (fG, fsG);

  /**
  We compute the vaporization rate from the interface
  jump condition, obtaining the vaporization rate per
  unit of interface surface, stored in *mEvap*. */

  foreach() {
    mEvap[] = 0.; TInt[] = 0.;
    if (f[] > F_ERR && f[] < 1.-F_ERR) {
      TInt[] = TIntVal;

      double ltrgrad = ebmgrad (point, TL, fL, fG, fsL, fsG, false, TIntVal, false);
      double gtrgrad = ebmgrad (point, TG, fL, fG, fsL, fsG, true,  TIntVal, false);

      mEvap[] += lambda1*ltrgrad/dhev;
      mEvap[] += lambda2*gtrgrad/dhev;

#ifdef SOLVE_LIQONLY
      mEvap[] = lambda1*ltrgrad/dhev;
#endif

#ifdef SOLVE_GASONLY
      mEvap[] = lambda2*gtrgrad/dhev;
#endif
    }
  }

  // if(grid->maxdepth > 10 && t < 20e-6)
  // {
  //   int total_num = 0;
  //   double total_mEvap = 0.0;
  //   foreach(reduction(+:total_num)reduction(+:total_mEvap))
  //   {
  //     if((int)is_solid[] == 1)
  //       continue;

  //     if(mEvap[] != 0.0)
  //     {
  //       total_num ++;
  //       total_mEvap += mEvap[];
  //     }
  //   }

  //   foreach ()
  //   {
  //     if ((int)is_solid[] == 1)
  //       continue;

  //     if (mEvap[] != 0.0)
  //     {
  //       mEvap[] = total_mEvap / (double)total_num;
  //     }
  //   }
  // }

  //boundary({mEvap});

  /**
  The calculation of the interface gradients is used
  also for the calculation of the source terms for the
  diffusion equation of the temperature fields. */

  //for me, this can be moved to trace_diffusion directly
//   foreach() {
//     sgT[] = 0., slT[] = 0.;
//     if (f[] > F_ERR && f[] < 1.-F_ERR) {
//       coord n = facet_normal (point, fL, fsL), p;
//       double alpha = plane_alpha (fL[], n);
//       double area = plane_area_center (n, alpha, &p);
//       normalize (&n);

//       double ltrgrad = ebmgrad (point, TL, fL, fG, fsL, fsG, false, TIntVal, false);
//       double gtrgrad = ebmgrad (point, TG, fL, fG, fsL, fsG, true, TIntVal, false);

//       double lheatflux = lambda1*ltrgrad;
//       double gheatflux = lambda2*gtrgrad;

// #ifdef AXI
//       slT[] = lheatflux/rho1/cp1*area*(y + p.y*Delta)/(Delta*y)*cm[];
//       sgT[] = gheatflux/rho2/cp2*area*(y + p.y*Delta)/(Delta*y)*cm[];
// #else
//       slT[] = lheatflux/rho1/cp1*area/Delta*cm[];
//       sgT[] = gheatflux/rho2/cp2*area/Delta*cm[];
// #endif
//     }
//   }

  /**
  We restore the tracer form of the liquid and gas-phase
  temperature fields. */

  foreach() {
    TL[] *= f[]*(f[] > F_ERR);
    TG[] *= (1. - f[])*((1. - f[]) > F_ERR);
    T[]  = TL[] + TG[];
  }
  //boundary({TL,TG,T});
}

/**
## Tracer Advection

We let the volume fractions *fu* and *fuext* to
advect the fields YL and YG, as implemented in
the tracer_advection event of [evaporation.h](evaporation.h)
*/

event tracer_advection (i++);

/**
## Tracer Diffusion

We solve the diffusion equations for the temperature fields
accounting for the phase change contributions. */

void fixSmallCellsTemp(scalar TL, scalar TG)
{
  extern scalar f;
  scalar Tfix[];
  //fix gas
  foreach()
  {
    Tfix[] = nodata;
#if USE_MY_SOLID
    if(is_solid[] == 1)
      continue;
#endif
    if(f[] >= (1.0 - fix_small_thd) && f[] < 1.0)
    {
      double total = 0.0;
      int num = 0;
      double frac = 0.0;
      foreach_neighbor(1)
      {
#if USE_MY_SOLID
        if (is_solid[] == 1)
          continue;
#endif
        if(f[] < 0.5)
        {
          total += (1.0 - f[]) * TG[];
          frac += 1.0 - f[];
          num ++;
        }
      }
      Tfix[] = frac == 0 ? TIntVal : total / frac;
      //printf("num = %d total = %g fixvalue = %g\n", num, total, Tfix[]);
    }
  }

  foreach()
  {
    TG[] = Tfix[] == nodata ? TG[] : Tfix[];
  }

  //fix liquid
  foreach()
  {
    Tfix[] = nodata;
#if USE_MY_SOLID
        if (is_solid[] == 1)
          continue;
#endif
    if(f[] > 0.0 && f[] < fix_small_thd)
    {
      double total = 0.0;
      int num = 0;
      double frac = 0.0;
      foreach_neighbor(1)
      {
#if USE_MY_SOLID
        if (is_solid[] == 1)
          continue;
#endif
        if(f[] > 0.5)
        {
          total += f[] * TL[];
          num ++;
          frac += f[];
        }
      }
      Tfix[] = frac == 0 ? TIntVal : total / frac;
      //printf("num = %d total = %g fixvalue = %g\n", num, total, Tfix[]);
    }
  }

  foreach()
  {
    TL[] = Tfix[] == nodata ? TL[] : Tfix[];
  }

  boundary({TL, TG});
}

void fixSmallCellsTempLevel(scalar TL, scalar TG, int l)
{
  extern scalar f;
  foreach_level(l)
  {
    if(f[] >= (1.0 - fix_small_thd) && f[] < 1.0)
    {
      TG[] = TIntVal;
    }
    else if(f[] > 0.0 && f[] < fix_small_thd)
    {
      TL[] = TIntVal;
    }
  }

  TL.dirty = TG.dirty = true;
  boundary_level({TL, TG}, l);
}

#if(USE_CONJUGATE_HEAT)
double getTbcS(double dx, double Tsol, double Tflu, double f1, double f2)
{
  double R_sol = dx / 2.0 / lambdas;
  double R_gas = dx / 2.0 / lambda2;
  double R_liq = dx / 2.0 / lambda1;

  double R_f = f1 * R_liq + f2 * R_gas;
  double R_tol = Rcc + R_sol + R_f;

  double Tfb = (R_f * Tsol + (Rcc + R_sol) * Tflu) / R_tol;
  double Tsb = ((Rcc + R_f) * Tsol + R_sol * Tflu) / R_tol;

  return Tsb;
}

double getTbcF(double dx, double Tsol, double Tflu, double f1, double f2)
{
  double R_sol = dx / 2.0 / lambdas;
  double R_gas = dx / 2.0 / lambda2;
  double R_liq = dx / 2.0 / lambda1;

  double R_f = f1 * R_liq + f2 * R_gas;
  double R_tol = Rcc + R_sol + R_f;

  double Tfb = (R_f * Tsol + (Rcc + R_sol) * Tflu) / R_tol;
  double Tsb = ((Rcc + R_f) * Tsol + R_sol * Tflu) / R_tol;

  return Tfb;
}

void imposeTempBoundaryConjugate(scalar TL, scalar TG, scalar TS, bool with_Rcc)
{
  extern scalar f;
  foreach()
  {
    if(is_solid[] == 1.0)
    {
      TL[] = TS[];
      TG[] = TS[];
      if(with_Rcc && is_solid[1] == 0.)
      {
        double T_sol = TS[];
        double T_liq = TL[1];
        double T_gas = TG[1];
        double tmp = getTbcF(Delta, T_sol, T_liq, 1.0, 0.0);
        double tmp2 = getTbcF(Delta, T_sol, T_gas, 0.0, 1.0);
        TL[] = 2.0 * getTbcF(Delta, T_sol, T_liq, 1.0, 0.0) - T_liq;
        TG[] = 2.0 * getTbcF(Delta, T_sol, T_gas, 0.0, 1.0) - T_gas;
      }
    }
    else
    {
      double T_liq = (fL[] > F_ERR) ? TL[] * (fL[]) : 0.;
      double T_gas = (fG[] > F_ERR) ? TG[] * (1.0 - fL[]): 0.;
      TS[] = T_liq + T_gas;
    }
  }

  TL.dirty = TG.dirty = TS.dirty = true;
  boundary({TL, TG, TS});
}
#endif

#if USE_CONJUGATE_HEAT

void coefDiffusionConjugate_1(double time, double dt)
{
  extern scalar is_solid_x, is_solid_y;
  face vector solid_face_b[];

  double dx = L0 / (1 << grid->maxdepth);

  double R_sol = dx / 2.0 / lambdas;
  double R_gas = dx / 2.0 / lambda2;
  double R_liq = dx / 2.0 / lambda1;

  double R_tol_l = Rcc + R_sol + R_liq;
  double R_tol_g = Rcc + R_sol + R_gas;

  foreach()
  {
    betal[] = 0.0;
    betag[] = 0.0;

    if(is_solid[] == 0. && is_solid_x[-1] == 1.)
    {
      double ratio_liq = R_liq / R_tol_l;
      double coefl = lambda1f.x[] * 2.0 * ratio_liq / Delta / Delta;
      betal[] = - coefl;
      slT[] += coefl * TS[-1];
      slT[] += coefl * R_sol * delta_q(y, time);

      double ratio_gas = R_gas / R_tol_g;
      double coefg = lambda2f.x[] * 2.0 * ratio_gas / Delta / Delta;
      betag[] = - coefg;
      sgT[] += coefg * TS[-1];
      sgT[] += coefg * R_sol * delta_q(y, time);
    }
  }

  foreach_face(x)
  {
    if(is_solid[] == 0. && is_solid_x[-1] == 1.)
    {
      lambda1f.x[] = 0.0;
      lambda2f.x[] = 0.0;
    }
  }

  foreach()
  {
    if(is_solid_x[] == 1. && is_solid[1] == 0.)
    {
      betal[] = 1.0 * dt;
      slT[] = -TS[] * dt;

      betag[] = 1.0 * dt;
      sgT[] = -TS[] * dt;
    }
  }

  //coefficients for heat transfer in solid
  //we use solid_face_b to store the coefficients
  foreach_face()
  {
    solid_face_b.x[] = is_solid_x_face.x[] * (1.0 - is_solid_y_face.x[]);
  }

  foreach_face(y)
  {
    if(is_solid_x[] == 1. && is_solid_y[] == 0. && is_solid_y[0,-1] == 1.)
    {
      solid_face_b.y[] = 1.0;
    }
  }

  foreach_face()
  {
    lambdasf.x[] = thermocoefsf.x[] * fm.x[] * solid_face_b.x[] * dt;
  }

  foreach()
  {
    betas[] = 0.0;
    ssT[] = 0.0;
    if(is_solid_x[] == 1. && is_solid[1] == 0.)
    {
      double ratio_liq = R_sol / R_tol_l;
      double ratio_gas = R_sol / R_tol_g;

      double coef_liq = lambdasf.x[1] * 2.0 * ratio_liq * fsL.x[1];
      double coef_gas = lambdasf.x[1] * 2.0 * ratio_gas * fsG.x[1];
      double coef_sol = (coef_liq + coef_gas) / Delta / Delta;
      betas[] = - coef_sol;
      ssT[] = coef_sol * TS[1];

      double coef_delta = (coef_liq * (Rcc + R_liq) + coef_gas * (Rcc + R_gas)) / Delta / Delta;
      ssT[] += coef_delta * delta_q(y, time);
    }
  }
  foreach()
  {
    //for solid cells
    if(vof_heater[] > 0.0)
    {
      double grad = delta_q_volume(y, time);
      double frac = vof_heater[];
      double rhom = rhos * frac + rhos2*(1.0-frac);
      double cpm = cps * frac + cps2 * (1.0-frac);
      extern double thickness_heater;
      grad = grad / thickness_heater * vof_heater[] * cm[] * dt;
      ssT[] += grad * (is_solid_x[]) * (1.0 - is_solid_y[]);
    }
  }
  foreach_face(x)
  {
    if(is_solid[] == 0. && is_solid_x[-1] == 1.)
    {
      lambdasf.x[] = 0.0;
    }
  }

  foreach()
  {
    if(is_solid[] == 0. && is_solid_x[-1] == 1.)
    {
      betas[] = 1.0 * dt;
      ssT[] = -TS[] * dt;
    }
  }

  foreach() {
    double frac = vof_heater[];
    double rhom = rhos * frac + (1.0 - frac) * rhos2;
    double cpm = cps * frac + (1.0 - frac) * cps2;
    double rho_cpm = rhos * cps * frac + (1.0 - frac) * rhos2 * cps2;
    //rho_cpm = rhom * cpm;
    //rho_cpm = 1.0 / rhos / cps * frac + (1.0 - frac) / rhos2 / cps2;
    //rho_cpm = 1.0 / rho_cpm;
    thetasolid[] = rho_cpm * cm[] * (is_solid_x[]) * (1.0 - is_solid_y[]);
  }
}


double R_liq_tol_l(double frac)
{
  double dx = L0 / (1 << grid->maxdepth);

  double lbd_liq = lambda1;
  double lbd_gas = lambda2;
  double lbd_sol = lambdas * frac + lambdas2 * (1.0 - frac);

  double R_sol = dx / 2.0 / lbd_sol;
  double R_gas = dx / 2.0 / lbd_gas;
  double R_liq = dx / 2.0 / lbd_liq;

  double R_tol_l = Rcc + R_sol + R_liq;
  double R_tol_g = Rcc + R_sol + R_gas;

  double ratio = R_liq / R_tol_l;
  return ratio;
}

double R_gas_tol_g(double frac)
{
  double dx = L0 / (1 << grid->maxdepth);

  double lbd_liq = lambda1;
  double lbd_gas = lambda2;
  double lbd_sol = lambdas * frac + lambdas2 * (1.0 - frac);

  double R_sol = dx / 2.0 / lbd_sol;
  double R_gas = dx / 2.0 / lbd_gas;
  double R_liq = dx / 2.0 / lbd_liq;

  double R_tol_l = Rcc + R_sol + R_liq;
  double R_tol_g = Rcc + R_sol + R_gas;

  double ratio = R_gas / R_tol_l;
  return ratio;
}

double R_sol_tol_l(double frac)
{
  double dx = L0 / (1 << grid->maxdepth);

  double lbd_liq = lambda1;
  double lbd_gas = lambda2;
  double lbd_sol = lambdas * frac + lambdas2 * (1.0 - frac);

  double R_sol = dx / 2.0 / lbd_sol;
  double R_gas = dx / 2.0 / lbd_gas;
  double R_liq = dx / 2.0 / lbd_liq;

  double R_tol_l = Rcc + R_sol + R_liq;
  double R_tol_g = Rcc + R_sol + R_gas;

  double ratio_liq = R_sol / R_tol_l;
  return ratio_liq;
}

double R_sol_tol_g(double frac)
{
  double dx = L0 / (1 << grid->maxdepth);

  double lbd_liq = lambda1;
  double lbd_gas = lambda2;
  double lbd_sol = lambdas * frac + lambdas2 * (1.0 - frac);

  double R_sol = dx / 2.0 / lbd_sol;
  double R_gas = dx / 2.0 / lbd_gas;
  double R_liq = dx / 2.0 / lbd_liq;

  double R_tol_l = Rcc + R_sol + R_liq;
  double R_tol_g = Rcc + R_sol + R_gas;

  double ratio_gas = R_sol / R_tol_g;
  return ratio_gas;
}

double ris_sol(double frac)
{
  double dx = L0 / (1 << grid->maxdepth);
  double lbd_sol = lambdas * frac + lambdas2 * (1.0 - frac);

  double R_sol = dx / 2.0 / lbd_sol;
  return R_sol;
}

void coefDiffusionConjugate_2(double time, double dt)
{
  extern scalar is_solid_x, is_solid_y;
  face vector solid_face_b[];

  double dx = L0 / (1 << grid->maxdepth);

  double R_gas = dx / 2.0 / lambda2;
  double R_liq = dx / 2.0 / lambda1;


  foreach ()
  {
    betal[] = 0.0;
    betag[] = 0.0;

    if(is_solid[] == 0. && is_solid_x[-1] == 1.)
    {
      double ratio_liq = R_liq_tol_l(1.0);
      double coefl = lambda1f.x[] * 2.0 * ratio_liq / Delta / Delta;
      slT[] += coefl * ris_sol(1.0) * delta_q(y, time);

      double ratio_gas = R_gas_tol_g(1.0);
      double coefg = lambda2f.x[] * 2.0 * ratio_gas / Delta / Delta;
      sgT[] += coefg * ris_sol(1.0) * delta_q(y, time);
    }
  }

  foreach_face(x)
  {
    if (is_solid[] == 0. && is_solid_x[-1] == 1.)
    {
      double ratio_liq = R_liq_tol_l(1.0);
      lambda1f.x[] = lambda1f.x[] * 2.0 * ratio_liq;

      double ratio_gas = R_gas_tol_g(1.0);
      lambda2f.x[] = lambda2f.x[] * 2.0 * ratio_gas;
    }
  }

  // coefficients for heat transfer in solid
  // we use solid_face_b to store the coefficients
  foreach_face()
  {
    solid_face_b.x[] = is_solid_x_face.x[] * (1.0 - is_solid_y_face.x[]);
  }

  foreach_face(y)
  {
    if (is_solid_x[] == 1. && is_solid_y[] == 0. && is_solid_y[0, -1] == 1.)
    {
      solid_face_b.y[] = 1.0;
    }
  }

  foreach_face()
  {
    lambdasf.x[] = thermocoefsf.x[] * fm.x[] * solid_face_b.x[] * dt;
  }

  foreach ()
  {
    betas[] = 0.0;
    ssT[] = 0.0;
    if(is_solid_x[] == 1. && is_solid[1] == 0.)
    {
      double ratio_liq = R_sol_tol_l(1.0);
      double ratio_gas = R_sol_tol_g(1.0);

      double coef_liq = lambdasf.x[1] * 2.0 * ratio_liq * fsL.x[1];
      double coef_gas = lambdasf.x[1] * 2.0 * ratio_gas * fsG.x[1];

      double coef_delta = (coef_liq * (Rcc + R_liq) + coef_gas * (Rcc + R_gas)) / Delta / Delta;
      ssT[] += coef_delta * delta_q(y, time);
    }
  }

  foreach()
  {
    //for solid cells
    if(vof_heater[] > 0.0)
    {
      double grad = delta_q_volume(y, time);
      extern double thickness_heater;
      grad = grad / thickness_heater * vof_heater[] * cm[] * dt;
      ssT[] += grad * (is_solid_x[]) * (1.0 - is_solid_y[]);
    }
  }

  foreach_face(x)
  {
    if (is_solid[] == 0. && is_solid_x[-1] == 1.)
    {
      double ratio_liq = R_sol_tol_l(1.0);
      double ratio_gas = R_sol_tol_g(1.0);

      double coef_liq = lambdasf.x[] * 2.0 * ratio_liq * fsL.x[];
      double coef_gas = lambdasf.x[] * 2.0 * ratio_gas * fsG.x[];
      
      double coef_sol = (coef_liq + coef_gas);

      lambdasf.x[] = coef_sol;
    }
  }

  foreach ()
  {
    double frac = vof_heater[];
    double rhom = rhos * frac + (1.0 - frac) * rhos2;
    double cpm = cps * frac + (1.0 - frac) * cps2;
    double rho_cpm = rhos * cps * frac + (1.0 - frac) * rhos2 * cps2;
    //rho_cpm = rhom * cpm;
    //rho_cpm = 1.0 / rhos / cps * frac + (1.0 - frac) / rhos2 / cps2;
    //rho_cpm = 1.0 / rho_cpm;
    thetasolid[] = rho_cpm * cm[] * (is_solid_x[]) * (1.0 - is_solid_y[]);
  }

  //the coefficients beyond the current region of interest will not influence the result
  //the bc will be enforced during each iteration
  //and those cells out of the region will not be taken into account for the 
  //computation of residual.
  //This will improve the convergence
  foreach_face(y)
  {
    if(is_solid[] == 0.)
    {
      lambdasf.y[] = lambda1f.y[] + lambda2f.y[];
    }

    if(is_solid_x[] == 1. && is_solid_y[] == 0.)
    {
      lambda1f.y[] = lambdasf.y[];
      lambda2f.y[] = lambdasf.y[];
    }
  }

  foreach_face(x)
  {
    if(is_solid_x[] == 1. && is_solid_y[] == 0.)
    {
      lambda1f.x[] = lambdasf.x[];
      lambda2f.x[] = lambdasf.x[];
    }

    if(is_solid[] == 0. && is_solid[-1] == 0.)
    {
      lambdasf.x[] = lambda1f.x[] + lambda2f.x[];
    }
  }
}

void setSolidThermoCoef(face vector coef)
{
#if(USE_CONJUGATE_HEAT)
  foreach_face()
  {
    coef.x[] = lambdas2;
  }

  foreach_face(x)
  {
    //we know for sure this is the heater solid face
    if(is_solid[] == 0.0 && is_solid[-1] == 1.0)
    {
      coef.x[] = lambdas;
    }

    if(vof_heater[] >= 1.)
    {
      coef.x[] = lambdas;
    }
  }

  foreach_face(y)
  {
    if(vof_heater[] > 0.0)
    {
      double frac = vof_heater[];

      double lbde = lambdas * frac + (1.0 - frac) * lambdas2;
      double rhose = rhos * frac + (1.0 - frac) * rhos2;
      double cpse = cps * frac + (1.0 - frac) * cps2;

      coef.y[] = lbde;
    }
    // if(is_solid[] == 1. && is_solid[1] == 0. && is_solid_y[] == 0.)
    // {
    //   double thick_titanium = 500e-9;
    //   double frac = thick_titanium / Delta;

    //   double lbde = lambdas * frac + (1.0 - frac) * lambdas2;
    //   double rhose = rhos * frac + (1.0 - frac) * rhos2;
    //   double cpse = cps * frac + (1.0 - frac) * cps2;

    //   coef.y[] = lbde / rhose / cpse;
    // }
  }
#else
  foreach_face()
  {
    coef.x[] = lambdas / rhos / cps;
  }
#endif
}

#endif //USE_CONJUGATE_HEAT

event tracer_diffusion (i++)
{
  /**
  We remove the fractions of f and mass fractions
  lower than F_ERR and we reconstruct the non-volume
  averaged form of the mass fraction fields, in order
  to improve the discretization of the face gradients
  in the diffusion equation. */

  foreach() {
    f[] = clamp (f[], 0., 1.);
    f[] = (f[] > F_ERR) ? f[] : 0.;
    //fL[] = f[]; fG[] = 1. - f[];
    //TL[] = fL[] > F_ERR ? TL[]/fL[] : 0.;
    //TG[] = ((1. - fL[]) > F_ERR) ? TG[]/(1. - fL[]) : 0.;

#ifdef CONSISTENTPHASE1
    TL[] = fuext[] > F_ERR ? TL[]/fuext[] : TIntVal;
#else
    TL[] = fu[] > F_ERR ? TL[]/fu[] : TIntVal;
#endif
#ifdef CONSISTENTPHASE2
    TG[] = ((1. - fuext[]) > F_ERR) ? TG[]/(1. - fuext[]) : TIntVal;
#else
    TG[] = ((1. - fu[]) > F_ERR) ? TG[]/(1. - fu[]) : TIntVal;
#endif
    fL[] = f[]; fG[] = 1. - f[];
  }
  //boundary({fL,fG,TL,TG});

  /**
  We compute the value of volume fraction *f* on the
  cell-faces using a geometric approach (necessary
  for interface gradients and diffusion equations). */

#if USE_CONJUGATE_HEAT
  imposeTempBoundaryConjugate(TL, TG, TS, true);
#endif

  face_fraction (fL, fsL);
  face_fraction (fG, fsG);

  foreach ()
  {
    sgT[] = 0., slT[] = 0.;
    if (f[] > F_ERR && f[] < 1. - F_ERR)
    {
      coord n = facet_normal(point, fL, fsL), p;
      double alpha = plane_alpha(fL[], n);
      double area = plane_area_center(n, alpha, &p);
      normalize(&n);

      double ltrgrad = ebmgrad(point, TL, fL, fG, fsL, fsG, false, TIntVal, false);
      double gtrgrad = ebmgrad(point, TG, fL, fG, fsL, fsG, true, TIntVal, false);

      double lheatflux = lambda1 * ltrgrad;
      double gheatflux = lambda2 * gtrgrad;

#ifdef AXI
      slT[] = lheatflux * area * (y + p.y * Delta) / (Delta * y) * cm[] * dt;
      sgT[] = gheatflux * area * (y + p.y * Delta) / (Delta * y) * cm[] * dt;
#else
      slT[] = lheatflux * area / Delta * cm[] * dt;
      sgT[] = gheatflux * area / Delta * cm[] * dt;
#endif
    }
  }

  /**
  We solve the diffusion equations, confined by means of
  the face fraction fields *fsL* and *fsG*. */

  foreach_face() {
    lambda1f.x[] = lambda1 * fsL.x[] * fm.x[] * dt;
    lambda2f.x[] = lambda2 * fsG.x[] * fm.x[] * dt;
  }
  //boundary((scalar *){lambda1f,lambda2f});

  foreach() {
    thetacorr1[] = cm[] * max(fL[], F_ERR) * rho1 * cp1;
    thetacorr2[] = cm[] * max(fG[], F_ERR) * rho2 * cp2;
  }
  //boundary({thetacorr1,thetacorr2});

#if USE_MY_SOLID
  //fixme: this should be modified when conjugate heat transfer is consisdered
  foreach()
  {
    sgT[] *= (1.0 - is_solid[]);
    slT[] *= (1.0 - is_solid[]);
    thetacorr1[] *= (1.0 - is_solid[]);
    thetacorr2[] *= (1.0 - is_solid[]);
  }

  extern scalar is_solid_x, is_solid_y;
  face vector solid_face_b[];
  foreach_face()
  {
    solid_face_b.x[] = is_solid_face.x[];
  }

  foreach_face(x)
  {
    if(is_solid[] == 0. && is_solid_x[-1] == 1.)
    {
      solid_face_b.x[] = 0.0;
    }
  }

  foreach_face()
  {
    lambda1f.x[] *= (1.0 - solid_face_b.x[]);
    lambda2f.x[] *= (1.0 - solid_face_b.x[]);
  }

#endif

#if USE_CONJUGATE_HEAT
  imposeTempBoundaryConjugate(TL, TG, TS, false);
  setSolidThermoCoef(thermocoefsf);
  scalar resL[], resG[], resS[];
#if (USE_CONJUGATE_HEAT == 1)
  coefDiffusionConjugate_1(t + dt + tshift, dt);
  mgT = diffusion(TS, 1.0, D = lambdasf, r = ssT, beta = betas, theta = thetasolid, res = {resS});
  diffusion (TG, 1.0, D=lambda2f, r=sgT, beta=betag, theta=thetacorr2, res = {resG});
  diffusion (TL, 1.0, D=lambda1f, r=slT, beta=betal, theta=thetacorr1, res = {resL});
  //mydiffusion(TL, TG, dt, D1 = lambda1f, D2 = lambda2f, r1 = slT, r2 = sgT, beta1 = betal, beta2 = betag,theta1 = thetacorr1, theta2 = thetacorr2);
#else //USE_CONJUGATE_HEAT == 1
  coefDiffusionConjugate_2(t + dt + tshift, dt);
  mgT = mydiffusionConjugate(TL, TG, TS, 1.0, Dl = lambda1f, Dg = lambda2f, Ds = lambdasf, rl = slT, rg = sgT, rs = ssT,
                                betal = betal, betag = betag, betas = betas, thetal = thetacorr1, thetag = thetacorr2, thetas = thetasolid,
                                resl = resL, resg = resG, ress = resS);
#endif //USE_CONJUGATE_HEAT == 1
  imposeTempBoundaryConjugate(TL, TG, TS, true);
#else //USE_CONJUGATE_HEAT
//#ifndef SOLVE_LIQONLY
  diffusion (TG, 1.0, D=lambda2f, r=sgT, theta=thetacorr2);
//#endif
//#ifndef SOLVE_GASONLY
  diffusion (TL, 1.0, D=lambda1f, r=slT, theta=thetacorr1);
//#endif
#endif


  //for the cut cells, the temperature of the phase with
  //smaller volume fraction will be more inaccurate
  fixSmallCellsTemp(TL, TG);

  foreach() {
    TL[] *= fL[];
    TG[] *= (1. - fL[]);
  }
  //boundary({TL,TG});

  /**
  We reconstruct the one-field temperature field summing
  the two fields $T_L$ and $T_G$ in tracer form. */

  foreach() {
    TL[] = (fL[] > F_ERR) ? TL[] : 0.;
    TG[] = (fG[] > F_ERR) ? TG[] : 0.;
    T[] = TL[] + TG[];
  }
  //boundary({TL,TG,T});


#if USE_MY_SOLID
  foreach()
    T[] = (int)is_solid[] == 1 ? TS[] : T[];
#endif
}

