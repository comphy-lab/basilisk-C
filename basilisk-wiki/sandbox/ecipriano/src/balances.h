/**
# Mass Conservation Balances

This module implements mass consevation tests for the phase
change models. The operations implemented here are too boring
and long to be written in every simulation. Therefore, an external
header that automatically computes the balances is useful. The
implementation is based on the nomenclature and on the data
structures used in [multicomponent.h](multicomponent.h). However, it
stores everything in a separate structure, making this module general
also for other phase change models.

The mass conservation tests can be defined in different ways. If a
liquid phase is placed in a closed domain, and no phase change
occurs, the mass balance reads:
$$
  \int_V\rho_l dV = m_L^0
$$
which means that the total liquid mass (assumed to be a function of
the volume fraction only) must be constant and equal to the initial
amount of mass in liquid phase $m_L^0$.
If phase change is present, the variation of the liquid mass depends
on the total amount of evaporated liquid:
$$
  \int_V \rho_l dV - \int_0^{\tau}dt\int_V\dot{m}\delta_\Gamma dV
  = m_L^0
$$
where $\dot{m} < 0$ if the liquid evaporates. If the domain is not
closed and the liquid (or the gas-phase) is allowed to leave the
domain, the mass balance must include also the contribution of the
fluxes across the boundaries:
$$
  \int_V \rho_l dV +
  \int_0^{\tau}dt\oint_{S_b} \rho_l \mathbf{u}\cdot\mathbf{n} dS_b -
  \int_0^{\tau}dt\int_V\dot{m}\delta_\Gamma dV =
  m_L^0
$$
Exactly the same concept applies for the gas phase. If the same
balances are written for every chemical species, the equation must
account for the mass fraction of the chemical species and the
diffusive flux:
$$
  \int_V \rho_l\omega_{i,l} dV +
  \int_0^{\tau}dt\oint_{S_b} \rho_l\omega_{i,j}
  \mathbf{u}\cdot\mathbf{n} dS_b -
  \int_0^{\tau}dt\oint_{S_b}
  \rho_l\mathcal{D}_{i,l}\nabla\omega_{i,j}\cdot\mathbf{n}dS_b -
  \int_0^{\tau}dt\int_V\dot{m}\delta_\Gamma dV =
  m_{L,i}^0
$$
*/

#define BALANCES

/**
We define a structure that gathers data that are used or
modified in this module. */

struct MassBalance {
  // Data
  scalar * YLList;
  scalar * YGList;
  scalar * mEvapList;
  double * liq_start;
  double * gas_start;
  char ** liq_species;
  char ** gas_species;
  double rho1;
  double rho2;
  double * inDmix1;
  double * inDmix2;
  bool boundaries;
  int maxlevel;
  FILE * fb;
  // Results
  double totmass1;
  double totmass2;
  double totmass1init;
  double totmass2init;
  double totmass1bd;
  double totmass2bd;
  double totmass1bdnow;
  double totmass2bdnow;
  double evapmass;
  double evapmassnow;
  double * liq_mass;
  double * gas_mass;
  double * liq_massinit;
  double * gas_massinit;
  double * liq_massbd;
  double * gas_massbd;
  double * liq_massbdnow;
  double * gas_massbdnow;
  double * tot_evap;
  double * tot_evapnow;
  double * liq_relerr;
  double * gas_relerr;
};

struct MassBalance mb = {0};

/**
## Helper Functions

We define functions that facilitate loops over all the domain
boundaries, adjusting the stencil index according to the boundary ID.
*/

@if _MPI
# define allreduce(vector,size) {                                   \
    MPI_Allreduce (MPI_IN_PLACE, vector, size, MPI_DOUBLE, MPI_SUM, \
        MPI_COMM_WORLD);                                            \
}
@endif

/**
### *face_gradient_bid()*: Return the face gradient at a given boundary.
*/

static double face_gradient_bid (Point point, scalar Y, int bid) {
  double grad = 0.;
  switch (bid) {
    case 0: grad = (Y[1,0]  - Y[])/Delta;  break;  // right
    case 1: grad = (Y[-1,0] - Y[])/Delta;  break;  // left
    case 2: grad = (Y[0,1]  - Y[])/Delta;  break;  // top
    case 3: grad = (Y[0,-1] - Y[])/Delta;  break;  // bottom
  }
  return grad;
}

/**
### *face_value_bid()*: Return the face value at a given boundary.
*/

static double face_value_bid (Point point, scalar Y, int bid) {
  double val = 0.;
  switch (bid) {
    case 0: val = 0.5*(Y[1,0]  + Y[]); break;  // right
    case 1: val = 0.5*(Y[-1,0] + Y[]); break;  // left
    case 2: val = 0.5*(Y[0,1]  + Y[]); break;  // top
    case 3: val = 0.5*(Y[0,-1] + Y[]); break;  // bottom
  }
  return val;
}

/**
### *face_flux_bid()*: Return the face flux at a given boundary.

We must choose a convention: the flux is positive if leaving the
domain from the right and top boundaries, and negative if it enters
the fomain from left and bottom. */

static double face_flux_bid (Point point, scalar Y, int bid) {
  double flux = 0., faceval = face_value_bid (point, Y, bid);
  switch (bid) {
    case 0: flux = +faceval*uf.x[]*Delta; break;  // right
    case 1: flux = -faceval*uf.x[]*Delta; break;  // left
    case 2: flux = +faceval*uf.y[]*Delta; break;  // top
    case 3: flux = -faceval*uf.y[]*Delta; break;  // bottom
  }
  return flux;
}

/**
### *diffusion_boundary()*: Sum the diffusive fluxes across a given boundary.
*/

static void diffusion_boundary (Point point, int bid) {
  foreach_elem (YGList, jj) {
    scalar YG = YGList[jj];
    double gradYG = face_gradient_bid (point, YG, bid);
    mb.gas_massbdnow[jj] -= mb.rho2*mb.inDmix2[jj]*gradYG*Delta*dt;
  }
  foreach_elem (YLList, jj) {
    scalar YL = YLList[jj];
    double gradYL = face_gradient_bid (point, YL, bid);
    mb.liq_massbdnow[jj] -= mb.rho1*mb.inDmix1[jj]*gradYL*Delta*dt;
  }
}

/**
### *advection_boundary()*: Sum the convective flux across a given boundary.
*/

scalar U[];

static void advection_boundary (Point point, int bid) {
  foreach_elem (YGList, jj) {
    scalar YG = YGList[jj];
    double fluxYG = face_flux_bid (point, YG, bid);
    mb.gas_massbdnow[jj] += mb.rho2*fluxYG*dt;
  }
  foreach_elem (YLList, jj) {
    scalar YL = YLList[jj];
    double fluxYL = face_flux_bid (point, YL, bid);
    mb.liq_massbdnow[jj] += mb.rho1*fluxYL*dt;
  }

  /**
  Linear interpolation approximation for the calculation of the face
  volume fraction. */

  double ff = face_value_bid (point, f, bid);
  mb.totmass1bdnow += mb.rho1*face_flux_bid (point, U, bid)*ff*dt;
  mb.totmass2bdnow += mb.rho2*face_flux_bid (point, U, bid)*(1. - ff)*dt;
}

/**
### *write_balances()*: write the results in the correct order.
*/

static void write_balances (void) {
  fprintf (mb.fb, "%g %g %g %g ",
      t, mb.totmass1, mb.totmass2, mb.evapmass);
  foreach_elem (mb.YLList, jj)
    fprintf (mb.fb, "%g ", mb.liq_mass[jj]);
  foreach_elem (mb.YGList, jj)
    fprintf (mb.fb, "%g ", mb.gas_mass[jj]);
  foreach_elem (mb.mEvapList, jj)
    fprintf (mb.fb, "%g ", mb.tot_evap[jj]);
  foreach_elem (mb.YLList, jj)
    fprintf (mb.fb, "%g ", mb.liq_relerr[jj]);
  foreach_elem (mb.YGList, jj)
    fprintf (mb.fb, "%g ", mb.gas_relerr[jj]);
  fprintf (mb.fb, "\n");
}

/**
### *compute_balances()*: compute the volume and surface integrals.
*/

static void compute_balances (void) {
  /**
  The total liquid and gas phase masses for each species are set to
  zero before summing up over every cell and boundary. */

  foreach_elem (YGList, jj) {
    mb.gas_mass[jj] = 0.;
    mb.gas_massbdnow[jj] = 0.;
    mb.tot_evapnow[jj] = 0.;
  }
  foreach_elem (YLList, jj) {
    mb.liq_mass[jj] = 0.;
    mb.liq_massbdnow[jj] = 0.;
  }
  mb.totmass1 = 0.;
  mb.totmass2 = 0.;
  mb.evapmassnow = 0.;
  mb.totmass1bdnow = 0.;
  mb.totmass2bdnow = 0.;

  /**
  We compute the total mass in liquid and gas phase for every
  chemical species of the domain. The total evaporated mass is also
  computed for each species. */

  foreach(serial) {
    vofrecon vr = vof_reconstruction (point, f);
    foreach_elem (mb.mEvapList, jj) {
      scalar mEvap = mEvapList[jj];
      scalar YG    = YGList[jj];

      mb.tot_evapnow[jj] += mEvap[]*vr.dirac*dv()*dt;
      mb.gas_mass[jj] += mb.rho2*YG[]*dv();
    }
    foreach_elem (mb.YLList, jj) {
      scalar YL = mb.YLList[jj];
      mb.liq_mass[jj] += mb.rho1*YL[]*dv();
    }
    mb.totmass1 += mb.rho1*f[]*dv();
    mb.totmass2 += mb.rho2*(1. - f[])*dv();
    mb.evapmassnow += mEvapTot[]*vr.dirac*dv()*dt;
  }

  /**
  If *boundaries* are included (not by default), the total masses are
  corrected according to the convective and diffusive fluxes through
  the boundaries. */

  if (mb.boundaries) {
    for (int b = 0; b < nboundary; b++) {
      foreach_boundary(b) {
        diffusion_boundary (point, b);
        advection_boundary (point, b);
      }
    }
@if _MPI
    allreduce (&mb.totmass1bdnow, 1);
    allreduce (&mb.totmass2bdnow, 1);
    allreduce (mb.liq_massbdnow, NLS);
    allreduce (mb.gas_massbdnow, NGS);
@endif

    foreach_elem (YLList, jj)
      mb.liq_massbd[jj] += mb.liq_massbdnow[jj];
    foreach_elem (YGList, jj)
      mb.gas_massbd[jj] += mb.gas_massbdnow[jj];
    mb.totmass1bd += mb.totmass1bdnow;
    mb.totmass2bd += mb.totmass2bdnow;
  }

  /**
  We sum the masses contributions from every processor. */

@if _MPI
  allreduce (&mb.totmass1, 1);
  allreduce (&mb.totmass2, 1);
  allreduce (&mb.evapmassnow, 1);
  allreduce (mb.liq_mass, NLS);
  allreduce (mb.gas_mass, NGS);
  allreduce (mb.tot_evapnow, NGS);
@endif

  /**
  We add the boundaries contributions, and we compute the differences
  between the initial and the current mass. */

  foreach_elem (YLList, jj)
    mb.liq_mass[jj] = mb.liq_massinit[jj] - (mb.liq_mass[jj] + mb.liq_massbd[jj]);
  foreach_elem (YGList, jj) {
    mb.tot_evap[jj] += mb.tot_evapnow[jj];
    mb.gas_mass[jj] = (mb.gas_mass[jj] + mb.gas_massbd[jj]) - mb.gas_massinit[jj];
  }
  mb.totmass1 = mb.totmass1init - (mb.totmass1 + mb.totmass1bd);
  mb.totmass2 = (mb.totmass2 + mb.totmass2bd) - mb.totmass2init;
  mb.evapmass += mb.evapmassnow;
}

/**
## Defaults

We initialize the varibales of the mass balances struct. */

event defaults (i = 0) {
  mb.YLList = NULL;
  mb.YGList = NULL;
  mb.mEvapList = NULL;
  mb.liq_start = NULL;
  mb.gas_start = NULL;
  mb.totmass1 = 0.;
  mb.totmass2 = 0.;
  mb.totmass1init = 0.;
  mb.totmass2init = 0.;
  mb.totmass1bd = 0.;
  mb.totmass2bd = 0.;
  mb.totmass1bdnow = 0.;
  mb.totmass2bdnow = 0.;
  mb.evapmass = 0.;
  mb.liq_mass = NULL;
  mb.gas_mass = NULL;
  mb.liq_massinit = NULL;
  mb.gas_massinit = NULL;
  mb.liq_massbd = NULL;
  mb.gas_massbd = NULL;
  mb.liq_massbdnow = NULL;
  mb.gas_massbdnow = NULL;
  mb.tot_evap = NULL;
  mb.tot_evapnow = NULL;
  mb.liq_relerr = NULL;
  mb.gas_relerr = NULL;
  mb.boundaries = false;
  mb.maxlevel = 5;

  /**
  We allocate the fields that store the resuls of the mass
  balance calculations. */

  mb.liq_mass      = (double *)malloc(NLS*sizeof(double));
  mb.gas_mass      = (double *)malloc(NGS*sizeof(double));
  mb.tot_evap      = (double *)malloc(NGS*sizeof(double));
  mb.tot_evapnow   = (double *)malloc(NGS*sizeof(double));
  mb.liq_relerr    = (double *)malloc(NLS*sizeof(double));
  mb.gas_relerr    = (double *)malloc(NGS*sizeof(double));
  mb.liq_massbd    = (double *)malloc(NLS*sizeof(double));
  mb.gas_massbd    = (double *)malloc(NGS*sizeof(double));
  mb.liq_massbdnow = (double *)malloc(NLS*sizeof(double));
  mb.gas_massbdnow = (double *)malloc(NGS*sizeof(double));
  mb.liq_massinit  = (double *)malloc(NLS*sizeof(double));
  mb.gas_massinit  = (double *)malloc(NGS*sizeof(double));

  /**
  We reset the initial values of mass in liquid and gas phase and the
  total vaporizing mass. */

  foreach_elem (mb.YLList, jj) {
    mb.liq_mass[jj] = 0.;
    mb.liq_massbd[jj] = 0.;
    mb.liq_massbdnow[jj] = 0.;
    mb.liq_massinit[jj] = 0.;
    mb.liq_relerr[jj] = 0.;
  }

  foreach_elem (mb.YGList, jj) {
    mb.gas_mass[jj] = 0.;
    mb.gas_massbd[jj] = 0.;
    mb.gas_massbdnow[jj] = 0.;
    mb.gas_massinit[jj] = 0.;
    mb.tot_evap[jj] = 0.;
    mb.gas_relerr[jj] = 0.;
  }

  /**
  We set a unit helper field, used for the boundary fluxes. */

  foreach()
    U[] = 1.;
}

/**
## Initialization

In the initialization event, we write the header of the
output file, in order to facilitate the identification
of the correct column to be plotted. */

/**
Write the *balances* file header and compute the intial mass of each
phase and species. */

event init (i = 0) {
  char name[80];
  sprintf (name, "balances-%d", mb.maxlevel);
  mb.fb = fopen (name, "w");
  fprintf (mb.fb, "time[s](1) totmass1(2) totmass2(3) evapmass(4) ");
  int counter = 5;
  foreach_elem (mb.YLList, jj)
    fprintf (mb.fb, "mL_%s(%d) ", mb.liq_species[jj], counter++);
  foreach_elem (mb.YGList, jj)
    fprintf (mb.fb, "mG_%s(%d) ", mb.gas_species[jj], counter++);
  foreach_elem (mb.mEvapList, jj)
    fprintf (mb.fb, "mE_%s(%d) ", mb.gas_species[jj], counter++);
  foreach_elem (mb.YLList, jj)
    fprintf (mb.fb, "liqerr_%s(%d) ", mb.liq_species[jj], counter++);
  foreach_elem (mb.YGList, jj)
    fprintf (mb.fb, "gaserr_%s(%d) ", mb.gas_species[jj], counter++);
  fprintf (mb.fb, "\n");

  /**
  We compute the initial mass in gas and liquid phase. */

  foreach(serial) {
    mb.totmass1init += mb.rho1*f[]*dv();
    mb.totmass2init += mb.rho2*(1. - f[])*dv();
  }

@if _MPI
  allreduce (&mb.totmass1init, 1);
  allreduce (&mb.totmass2init, 1);
@endif

  foreach_elem (YLList, jj)
    mb.liq_massinit[jj] = mb.totmass1init*mb.liq_start[jj];
  foreach_elem (YGList, jj)
    mb.gas_massinit[jj] = mb.totmass2init*mb.gas_start[jj];
}

/**
## Cleanup

We de-allocate the dynamic arrays from memory. */

event cleanup (t = end) {
  fclose (mb.fb);
  free (mb.liq_mass);
  free (mb.gas_mass);
  free (mb.tot_evap);
  free (mb.tot_evapnow);
  free (mb.liq_relerr);
  free (mb.gas_relerr);
  free (mb.liq_massbd);
  free (mb.gas_massbd);
  free (mb.liq_massbdnow);
  free (mb.gas_massbdnow);
  free (mb.liq_massinit);
  free (mb.gas_massinit);
}

/**
## Compute Balances

First, we write a function that is responsible for the
writing operations, based on the data in the structrue
*mb*. */

event balances (i++,last) {
  compute_balances();
  if (pid()==0)
    write_balances();
}

/**
## TODO

1. Add compatibility with metrics
2. Check the bid numbering system, what happens if new bid are
defined? And in 3D?
3. For some reasons *foreach_boundary()* does not work correctly when
a grid/multigrid.h is used.
*/

