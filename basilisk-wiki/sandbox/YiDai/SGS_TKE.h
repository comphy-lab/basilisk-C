/**
# The SGS TKE scheme
This header file implements SGS TKE scheme (Deardorff 1980; Heus et al 2010; template from [Antoon's SGS scheme](http://www.basilisk.fr/sandbox/Antoonvh/SGS.h).) We follow the parametrization of TKE equation from Heus et al 2010 (DALES code). 

$$\begin{aligned}
\frac{\partial e^{1 / 2}}{\partial t}= & -\tilde{u_j} \frac{\partial e^{1 / 2}}{\partial x_j}+\frac{1}{2 e^{1 / 2}} {\left[K_{\mathrm{m}}\left(\frac{\partial \tilde{u}_j}{\partial x_i}+\frac{\partial \tilde{u}_i}{\partial x_j}\right) \frac{\partial \tilde{u}_i}{\partial x_j}-K_{\mathrm{h}} \frac{g}{\theta_0} \frac{\partial\left(A \tilde{\theta}_{\mathrm{l}}+B \tilde{q}_{\mathrm{t}}\right)}{\partial z}\right] }+\frac{\partial}{\partial x_j}\left(2 K_{\mathrm{m}} \frac{\partial e^{1 / 2}}{\partial x_j}\right)-\frac{c_{\varepsilon} e}{2 \lambda}
\end{aligned}$$

We implictly solve the left first and fourth terms using tracer diffusion and directly add the shear and buoyancy production, dissipation terms as source terms in the solver.  
*/

#include "tracer.h"
#include "diffusion.h"

face vector Km[], Kh[], Ke[];
(const) face vector Pr;
scalar Evis[]; // Cell Centered diffusivity
// e120 should be declared together with buoyancy
// scalar e120[];
// scalar *tracers = {e120};
scalar e12p[];
double Cm, Cn;
double kappa;
double ce1, ce2;
double eps1 = 1E-10;

static inline void Evisprol(Point point, scalar s)
{
  foreach_child()
      Evis[] = bilinear(point, Evis) / 4.;
}

static inline void Evisres(Point point, scalar s)
{
  double sum = 0.;
  foreach_child()
      sum += s[];
  s[] = sum / 2.;
}

event defaults(i = 0)
{
  if (dimension != 3) // Allow to run, but give a warning
    fprintf(stdout, "Warning %dD grid. The used formulations only make sense for 3D turbulence simulations\n", dimension);
  mu = Km;
  Evis.prolongation = Evisprol;
  Evis.restriction = Evisres;
  Pr = unityf;
  kappa = 0.4;
  ce1 = 0.19;
  ce2 = 0.51;
  Cm = 0.12;
  Cn = 0.76;
  e120.refine = no_restriction;
  e120.coarsen = no_restriction;

#if TREE
  Evis.refine = no_restriction;
  Evis.coarsen = no_restriction;
  foreach_dimension()
  {
    Kh.x.refine = no_restriction;
    Km.x.coarsen = no_restriction;
  }
#endif
}

//Here first diagnose Km using mixing length scale
event Km_tke(i++)
{
  scalar lambda[];
  // stable or unstable condition !! this requries buoyancy field
  // how to transfer face gradient to cell center gradient and calculate the mixing length
  // - here we get face vector of buoyancy gradient
  // vector db[];
  // gradients({b}, {db});
  scalar dbdz[];
  foreach ()
  {
    dbdz[] = (b[0, 1] - b[0, -1]) / (Delta * 2);
    // limiting the temp gradient value
    if (dbdz[] < 0 && dbdz[] >= -eps1)
      dbdz[] = -eps1;
    if (dbdz[] >= 0 && dbdz[] <= eps1)
      dbdz[] = eps1;

    if (dbdz[] <= 0)
    {
      lambda[] = Delta;
    }
    else
    {
      e120[] = (e120[] <= Emin) ? Emin : e120[];
      lambda[] = pow(pow(kappa * y, -1) + pow(Cn * e120[] / sqrt(fabs(dbdz[])), -1), -1);
    }
    Evis[] = Cm * lambda[] * e120[];
  }
  boundary({Evis});

  foreach_face()
  {
    Km.x[] = (Evis[] + Evis[-1]) / 2;
    Kh.x[] = Km.x[] / Pr.x[];
    Ke.x[] = Km.x[] * 2.;
  }
  boundary({Km, Kh, Ke});

  /**
  TKE prognostic -- 1. shear production term sbshr
    $$\left(\frac{\partial \tilde{u}_j}{\partial x_i}+\frac{\partial \tilde{u}_i}{\partial x_j}\right) \frac{\partial \tilde{u}_i}{\partial x_j}$$
  
  */
  scalar tdef2[]; // cell center S^2/2
  foreach(){
    // the diagonal part
    tdef2[] = 2 / sq(Delta) * (sq(u.x[1, 0, 0] - u.x[]) + sq(u.y[0, 1, 0] - u.y[]) + sq(u.z[0, 0, 1] - u.z[]));
    // other parts
    tdef2[] += 0.25 * (sq((u.y[0, 0, 1] - u.y[-1, 0, 1]) / Delta + (u.x[0, 0, 1] - u.x[0, 0, 0]) / Delta) +
                       sq((u.y[0, 0, 0] - u.y[-1, 0, 0]) / Delta + (u.x[0, 0, 0] - u.x[0, 0, -1]) / Delta) +
                       sq((u.y[1, 0, 0] - u.y[0, 0, 0]) /  Delta + (u.x[1, 0, 0] - u.x[1, 0, -1]) / Delta) +
                       sq((u.y[1, 0, 1] - u.y[0, 0, 1]) /  Delta + (u.x[1, 0, 1] - u.x[1, 0, 0]) / Delta));
    tdef2[] += 0.25 * (sq((u.x[0, 1, 0] - u.x[0, 0, 0]) /  Delta + (u.z[0, 1, 0] - u.z[-1, 1, 0]) / Delta) +
                       sq((u.x[0, 0, 0] - u.x[0, -1, 0]) / Delta + (u.z[0, 0, 0] - u.z[-1, 0, 0]) / Delta) +
                       sq((u.x[1, 0, 0] - u.x[1, -1, 0]) / Delta + (u.z[1, 0, 0] - u.z[0, 0, 0]) / Delta) +
                       sq((u.x[1, 1, 0] - u.x[1, 0, 0]) /  Delta + (u.z[1, 1, 0] - u.z[0, 1, 0]) / Delta));
    tdef2[] += 0.25 * (sq((u.z[0, 0, 1] - u.z[0, 0, 0]) /  Delta + (u.y[0, 0, 1] - u.y[0, -1, 1]) / Delta) +
                       sq((u.z[0, 0, 0] - u.z[0, 0, -1]) / Delta + (u.y[0, 0, 0] - u.y[0, -1, 0]) / Delta) +
                       sq((u.z[0, 1, 0] - u.z[0, 1, -1]) / Delta + (u.y[0, 1, 0] - u.y[0, 0, 0]) / Delta) +
                       sq((u.z[0, 1, 1] - u.z[0, 1, 0]) /  Delta + (u.y[0, 1, 1] - u.y[0, 0, 1]) / Delta));
  }

  // here we cancel out the e120 in the equations to avoid division
  foreach ()
    e12p[] = (Cm * lambda[] * tdef2[] / 2 - 3 * Cm * lambda[] * dbdz[] / 2)   // shear and buoyancy
             - (ce1 + ce2 * lambda[] / Delta) * sq(e120[]) / (2 * lambda[]); // dissipation
}

// diffusiona and advection terms happen here
mgstats mgb;
event tracer_diffusion(i++)
{
  mgb = diffusion(e120, dt, Ke, r = e12p);
}
/**
## Test
* [Convective boundary layer case from Antoon](vhm_tke.c)
* There should be more tests

##Reference
Heus, T., van Heerwaarden, C.C., Jonker, H.J., Pier Siebesma, A., Axelsen, S., Van Den Dries, K., Geoffroy, O., Moene, A.F., Pino, D., De Roode, S.R. and Vilà-Guerau de Arellano, J., 2010. Formulation of the Dutch Atmospheric Large-Eddy Simulation (DALES) and overview of its applications. Geoscientific Model Development, 3(2), pp.415-444.
*/