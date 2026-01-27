/**
# Conjugate Heat Transfer

This module implements Conjugate Heat Transfer (CHT) in a naive manner, which
does not require embedded boundaries. One of the boundaries of the (square)
domain is considered as the solid-fluid interface; the solid is resolved on the
same grid as the fluid, and the boundary conditions on the solid-fluid interface
are computed and set in order to achieve the heat exchange between the different
phases.

The main advantage of this method is its simplicity, and the fact that it can
exploit the default ghost cells method to impose the boundary conditions.
Therefore, it does not require dealing with interfacial gradients as tipically
done with the embedded boundaries. However, the cost is the use of simple
geometries, and the introduction of additional cells into the solution of the
fluid phase. Conversely from methods using a master-slave coupling approach,
this approach automatically guarantees MPI parallelization, and it does not
introduce the overhead due to the interpolation of values from another domain,
potentially in another processor. In many systems (e.g., droplet combustion or
nucleate boiling problems) the solid heater is simple and can be effectively
described by this approach.

## Usage

Include this module after the Navier--Stokes and two-phase solver, there is not
strict rule about it. Set the global variables which store the solid material
properties.  Set additional heat sources such as the joule effect. Other source
terms can be introduced by directly modifying the field `QSexp` introduced at
the next section.  Finally, decide the shape of the solid by setting the volume
and surface fractions `cw` and `fw`, this can be done using the `solid()`
function but consider that the solid-fluid boundary must coicinde with one of the
domain boundaries. The following lines assume that the fluid-solid boundary is
the left one.

~~~literatec
#include "conjugate.h"

T[left] = dirichlet (TSB[]);
TS[left] = dirichlet (TFB[]);

int main() {
  ...
  rho3 = 8000., cp3 = 1080., lambda3 = 0.2;
  qjoule = 1e4;
  ...
}

event init (i = 0) {
  solid (cw, fw, x - LW);
}
~~~

## Improvements

1. Consider using a better adaptation method for `TS`.
1. Improve handling of Dirichlet boundary condition on the solid (embedded)
boundary which is not in contact with a domain boundary. */

/**
## Variables

We introduce the *wall* metric factors `cw` and `fw` which are used to correct
the integration of the solid heat equation. */

scalar cw[];
face vector fw[];

/**
The solid temperature `TS` is resolved in this file with boundary conditions
provided by the user. The fields `TSB` and `TFB` are the temperature boundary
fields on the solid, and on the fluid side, respectively. The heat source/sink
terms can be introduced using the field `QSexp`. If there is no thermal heat
transfer resistance, the temperature is continuous at the solid-fluid
boundaries. */

scalar TS[], TSB[], TFB[], TLB[], TGB[], QSexp[];

/**
### Physical properties

We define the density, heat capacity, and thermal conductivity of the solid. The
thermal power of the Joule effect can be set using the respective variable,
initialized as null. The initial solid temperature is set through `TS0`. */

double rho3 = 1., cp3 = 1., lambda3 = 1., qjoule = 0., TS0 = 300.;

/**
### Heat transfer resistance

We introduce the values of heat transfer resistance (from Hertz-Knudsen theory) and
the accommodation coefficient. The former is a function of the material
properties, while the latter is a case-dependent factor. */

double RInt = 0., acc = 1.;

/**
## Initialization

We specify an initial solid temperature and the function used to refine or
coarsen the solid temperature and the solid volume fraction. */


event defaults (i = 0) {
  foreach()
    TS[] = TS0;

#if TREE
  cw.refine = cw.prolongation = fraction_refine;
  for (scalar t in {TS}) {
    // fixme: solid tracers deserve a better refinement method
    t.refine = t.prolongation = refine_injection;
    t.dirty = true;
    //t.c = cw;
    //t.depends = list_add (t.depends, cw);
  }
#endif
}

event reset_sources (i++) {
  foreach()
    QSexp[] = 0.;
}

/**
## Solid-Fluid coupling

Before solving the diffusion equations on the solid and on the fluid domains, we
compute the temperatures at both sides of the interface from the energy balance.
These temperatures are imposed as Dirichlet boundary conditions. */

double conjugate_dirichlet = 0.;

static void conjugate (void) {
  scalar TL = liq->T, TG = gas->T;
  scalar lambdal = liq->lambda, lambdag = gas->lambda;
  foreach() {
    double RL = Delta/(2.*lambdal[] + 1e-10);
    double RG = Delta/(2.*lambdag[] + 1e-10);
    double RS = Delta/(2.*lambda3 + 1e-10);

    double lambdaf = f[]*lambda1 + (1. - f[])*lambda2;
    double RF = Delta/(2.*lambdaf);

    double TLL = (f[] > 1e-10) ? TL[]/f[] : 0.;
    double TGG = (1. - f[] > 1e-10) ? TG[]/(1. - f[]) : 0.;

    T[] = TL[] + TG[];

    TSB[] = ((acc*RInt + RF)*TS[] + RS*T[])/(acc*RInt + RS + RF);
    TLB[] = ((acc*RInt + RS)*TLL + RL*TS[])/(acc*RInt + RS + RL);
    TGB[] = ((acc*RInt + RS)*TGG + RG*TS[])/(acc*RInt + RS + RG);
    TFB[] = ((acc*RInt + RS)*T[] + RF*TS[])/(acc*RInt + RS + RF);
  }

  /**
  We use the interface gradients calculation in order to impose Dirichlet BCs on
  the solid heater. */

  if (conjugate_dirichlet) {
    scalar fl[], fg[];
    foreach() {
      fl[] = cw[];
      fg[] = 1. - cw[];
    }

    face vector fsl[], fsg[];
    face_fraction (fl, fsl);
    foreach_face()
      fsg.x[] = fsl.x[];

    foreach_interfacial_plic (cw, F_ERR) {
      double bc = conjugate_dirichlet;
      double intgrad = ebmgrad (point, TS, fl, fg, fsl, fsg, false, bc, false);
      QSexp[] += lambda3*intgrad*dirac;
    }
  }
}

/**
## Heat equation

We solve the heat equation for the solid phase:

$$
\rho_s cp_s \dfrac{\partial T_s}{\partial t} =
\nabla\cdot\left(\lambda_s\nabla T_s\right) + \dot{Q}_s
$$

considering metric factors and heat source/sink terms. Explicit (`e`) and
implicit (`i`) source terms to the solid diffusion equation can be added
exploiting the function pointer `energy_sources` defined below. */

void no_sources (scalar i, scalar e) {
  return;
}

void (* energy_sources) (scalar i, scalar e) = no_sources;

event tracer_diffusion (i++) {
  conjugate();
  face_fraction (cw, fw);

  scalar cmm[];
  foreach()
#if AXI
    cmm[] = max (fabs (Y0 - (y - Y0)), 1./HUGE)*cw[];
#else
    cmm[] = cm[]*cw[];
#endif

  face vector fmm[];
  foreach_face()
#if AXI
    fmm.x[] = max (fabs (Y0 - (y - Y0)), 1./HUGE)*fw.x[];
#else
    fmm.x[] = fm.x[]*fw.x[];
#endif

  scalar theta[];
  foreach()
    theta[] = max (cmm[]*rho3*cp3, 1e-10);

  face vector D[];
  foreach_face()
    D.x[] = fmm.x[]*lambda3;

  energy_sources (zeroc, QSexp);

  diffusion (TS, dt, D = D, theta = theta, r = QSexp);
}

/**
## Solid properties

The materials of the solid heaters/suspenders are pretty much the same in all
experimental works. Therefore, for common materials we provide a functions that
helps initializing solid properties. */

enum solid_props {
  QUARTZ, SiC, SAPPHIRE
};

void solid_properties (enum solid_props s) {
  switch (s) {
    case QUARTZ:
      rho3 = 2220, cp3 = 760., lambda3 = 1.5;
      break;
    case SiC:
      rho3 = 2740., cp3 = 670., lambda3 = 5.2;
      break;
    case SAPPHIRE:
      rho3 = 4890., cp3 = 410., lambda3 = 10.9;
      break;
  }
}

