/**
# Multicomponent Phase Change Model with Variable Properties

This is the variable-properties extension of the [multicomponent.h](/sandbox/ecipriano/src/multicomponent.h)
phase change model. Eventually, the two models will be merged in a
single header.

This phase change model describes the evaporation of non-isothermal
systems with multiple chemical species. The interface jump condition
is computed coupling the species and temperature equation in every
interfacial cell, assuming that thermodynamic equilibrium conditions
occur at the gas-liquid interface. The following balance is solved:

$$
  \begin{cases}
    \dot{m}_i = 
    \dot{m}\hat \omega_{i,l} - \rho_l \mathcal{D}_{i,l} \left.\dfrac{\partial
    \omega_{i,l}}{\partial\mathbf{n}_\Gamma}\right\vert_l 
    =
    \dot{m}\hat\omega_{i,g} - \rho_g \mathcal{D}_{i,g} \left.\dfrac{\partial
    \omega_{i,g}}{\partial\mathbf{n}_\Gamma}\right\vert_g
    \\
    \sum_{j=0}^{NS} \dot{m}_i \Delta h_{ev,i} = \lambda_l \left.\dfrac{
    \partial T_l}{\partial \mathbf{n}_\Gamma}\right\vert_l + \lambda_g
    \left.\dfrac{\partial T_g}{\partial \mathbf{n}_\Gamma}\right\vert_g
    \\
    \hat x_{i,g} = k_{eq,i}(\hat T) \hat x_{i,l}
    \\
  \end{cases}
$$

obtaining the vaporization rate $\dot{m}_i$ for each chemical species,
the liquid interface mass fractions $\hat{\omega}_{i,l}$, the gas
interface mass fractions $\hat{\omega}_{i,g}$, and the interface temperature
$\hat{T}$.
model
*/

#define PCM_MULTICOMPONENT

/**
We include the definitions of the interface gradients for the jump condition,
the face fraction for the diffusion equation, and thermodynamic functions
for the equilibrium at the gas-liquid interface.
*/

#include "intgrad.h"
#include "fracface.h"
#include "diffusion.h"
#include "thermodynamics.h"

/**
## Model Options

The following compiler directives can be defined to modify some options
of this code:

* *NGS*: number of chemical species in gas-phase
* *NLS*: number of chemical species in liquid-phase
* *SOLVE_TEMPERATURE*: non-isothermal system: solve the temperature equation
* *USE_CLAPEYRON*: Clausius-Clapeyron equation for the thermodynamic equilibrium
* *USE_ANTOINE*: Antoine equation for the thermodynamic equilibrium (see [thermodynamics.h](thermodynamics.h))

If *USE_GSL* is defined, we can solve the jump condition as a non-linear
system of equations. The calculation of the interface temperature from the
vaporization rate also relies of this keyword:

* *undefined*: temperature not solved, only decoupled solution of jump conditions
* *USE_GSL = 0*: decoupled solution of jump condition for the chemical species.
The interface temperature is obtained from the non-linear algebraic equations solver.
Required when *SOLVE_TEMPERATURE* is defined.
* *USE_GSL = 1*: decoupled solution of jump condition and interface temperature
to obtain fist guess values which are then refined by the non-linear algebraic equations
solver.
*/

#ifdef USE_GSL
# include "ijc-coupled-gsl.h"
#endif

/**
## User Data

The following data should be provided by the user:

* *gas_species*: list of gas phase species
* *liq_species*: list of liquid phase species
* *inert_species*: name of the inert species
* *gas_start*: initial mass fractions of gas species
* *liq_start*: initial mass fractions of liquid species
* *inDmix1*: diffusivities of species in liquid phase
* *inDmix2*: diffusivities of species in gas phase
* *inKeq*: thermodynamic equilibrium constant
* *Tboil*: boiling temperature of liquid species
(only if *USE_CLAPEYRON*)
*/

extern char* gas_species[NGS];
extern char* liq_species[NLS];
extern char* inert_species[1];

extern double gas_start[NGS];
extern double liq_start[NLS];

extern double inDmix1[NLS];
extern double inDmix2[NGS];

extern double inKeq[NLS];
#ifdef USE_CLAPEYRON
extern double Tboil[NLS];
#endif

/**
By default, the molecular weights are set to 1 for each
species. This can be modified in an init event. The default
thermodynamic reference pressure of the system is 1atm.
*/

double inMW[NGS];
double Pref = 101325.;

/**
If *SOLVE_TEMPERATURE* is defined, additional data must be
defined:

* *lambda1*: thermal conductivity in liquid phase
* *lambda2*: thermal conductivity in gas phase
* *dhev*: enthalpy of evaporation
* *cp1*: specific heat for the liquid phase
* *cp2*: specific heat for the gas phase
* *TG0, TL0*: inital temperatures in gas and liquid phase
*/

#ifdef SOLVE_TEMPERATURE
extern double lambda1, lambda2, dhev, cp1, cp2;
extern double TL0, TG0;

scalar T[], TInt[];
scalar TL, TG;
scalar sgT[], slT[], sgTimp[], slTimp[];
face vector lambda1f[], lambda2f[];

mgstats mgTL, mgTG;
#endif

/**
## Fields

Since we don't know a-priori how many species
we want to simulate, we initialize empty lists
containing the vaporization rate for every species,
the mass fractions, the interface properties and
phase change source terms for the transport equations.
*/

scalar * mEvapList = NULL;    // [NGS]
scalar * YList     = NULL;    // [NGS]
scalar * YLList    = NULL;    // [NLS]
scalar * YGList    = NULL;    // [NGS]
scalar * YLIntList = NULL;    // [NLS]
scalar * YGIntList = NULL;    // [NGS]
scalar * slexpList = NULL;    // [NLS]
scalar * slimpList = NULL;    // [NLS]
scalar * sgexpList = NULL;    // [NGS]
scalar * sgimpList = NULL;    // [NGS]
#ifdef VARPROP
scalar * Dmix1List = NULL;    // [NLS]
scalar * Dmix2List = NULL;    // [NGS]
scalar * dhevList  = NULL;    // [NLS]
scalar * Cp1List   = NULL;    // [NLS]
scalar * Cp2List   = NULL;    // [NGS]
#endif
#ifdef MOLAR_DIFFUSION
scalar * XLList    = NULL;    // [NLS]
scalar * XGList    = NULL;    // [NGS]
scalar * XLIntList = NULL;    // [NLS]
scalar * XGIntList = NULL;    // [NGS]
#endif

/**
We declare useful fields used for loops over chemical
species:

* *LSI*: vector with Liquid Species Indices [NLS]
* *GOSI*: vector with Gas-Only Species Indices [NGOS]
* *NGOS*: number of Gas-Only Species Indices
* *inertIndex*: index/position of the inert species.
*/

int * LSI = NULL;
int * GOSI = NULL;
int NGOS;
int inertIndex;

/**
We initilize other useful fields. */

bool success;
bool init_fields = true;

scalar fG[], fL[], fuT[];
face vector fsL[], fsG[];
scalar f0[];
scalar MW1mix[], MW2mix[];

scalar divu[], fold[];

/**
Variable properties stuff. */

#ifdef VARPROP
scalar dummy[];
scalar rho1v0[], rho2v0[];
#endif

double TEMPERATURE_TOLERANCE = 1.e-3 [*];

/**
## Defaults

In the defaults event, the lists containing the species
are initialized from the lists of species names in gas
and liquid phase. The tracers are assigned in this event.
*/

event defaults (i = 0)
{
  /**
  Fill Liquid Species Indices *LSI* array. */

  Array * arrLSI = array_new();
  for (int ii=0; ii<NGS; ii++) {
    for (int jj=0; jj<NLS; jj++) {
      if (strcmp(gas_species[ii], liq_species[jj]) == 0) {
        int idx = ii; array_append (arrLSI, &idx, sizeof(int));
      }
    }
  }
  LSI = (int *) array_shrink (arrLSI);

  /**
  Fill Gas Only Species Indices *GOSI* array. */

  Array * arrGOSI = array_new();
  for (int ii=0; ii<NGS; ii++) {
    bool thisSpeciesIsAlsoLiquid = false;
    for (int jj=0; jj<NLS; jj++) {
      if (strcmp(liq_species[jj], gas_species[ii]) == 0) {
        thisSpeciesIsAlsoLiquid = true;
      }
    }
    if (!thisSpeciesIsAlsoLiquid) {
      int idx = ii; array_append (arrGOSI, &idx, sizeof(int));
    }
  }
  NGOS = arrGOSI->len/sizeof(int);
  GOSI = (int *) array_shrink (arrGOSI);
  assert (NGOS == (NGS - NLS));

  /**
  Fill index of inert species among the gas species. */

  inertIndex = -1;
  for (int ii=0; ii<NGS; ii++) {
    if (strcmp(gas_species[ii], inert_species[0]) == 0)
      inertIndex = ii;
  }

  /**
  We create species fields from the list of species names. */

  YList     = NULL;
  YLList    = NULL;
  YGList    = NULL;
  YGIntList = NULL;
  YLIntList = NULL;
  mEvapList = NULL;
  slexpList = NULL;
  slimpList = NULL;
  sgexpList = NULL;
  sgimpList = NULL;
#ifdef VARPROP
  Dmix1List = NULL;
  Dmix2List = NULL;
  dhevList  = NULL;
  Cp1List   = NULL;
  Cp2List   = NULL;
#endif
#ifdef MOLAR_DIFFUSION
  XLList    = NULL;
  XGList    = NULL;
  XLIntList = NULL;
  XGIntList = NULL;
#endif

  for (int jj=0; jj<NLS; jj++) {
    scalar a = new scalar;
    free (a.name);
    char name[20];
    sprintf (name, "%s_L", liq_species[jj]);
    a.name = strdup (name);
    YLList = list_append (YLList, a);
  }
  reset (YLList, 0.);

  for (int jj=0; jj<NGS; jj++) {
    scalar a = new scalar;
    free (a.name);
    char name[20];
    sprintf (name, "%s_G", gas_species[jj]);
    a.name = strdup (name);
    YGList = list_append (YGList, a);
  }
  reset (YGList, 0.);

  for (int jj=0; jj<NLS; jj++) {
    scalar a = new scalar;
    free (a.name);
    char name[20];
    sprintf (name, "%s_LInt", liq_species[jj]);
    a.name = strdup (name);
    YLIntList = list_append (YLIntList, a);
  }
  reset (YLIntList, 0.);

  for (int jj=0; jj<NGS; jj++) {
    scalar a = new scalar;
    free (a.name);
    char name[20];
    sprintf (name, "%s_GInt", gas_species[jj]);
    a.name = strdup (name);
    YGIntList = list_append (YGIntList, a);
  }
  reset (YGIntList, 0.);

  for (int jj=0; jj<NGS; jj++) {
    scalar a = new scalar;
    free (a.name);
    char name[20];
    sprintf (name, "%s", gas_species[jj]);
    a.name = strdup (name);
    YList = list_append (YList, a);
  }
  reset (YList, 0.);

  for (int jj=0; jj<NGS; jj++) {
    scalar a = new scalar;
    free (a.name);
    char name[20];
    sprintf (name, "mEvap_%s", gas_species[jj]);
    a.name = strdup (name);
    mEvapList = list_append (mEvapList, a);
  }
  reset (mEvapList, 0.);

  for (int jj=0; jj<NLS; jj++) {
    scalar a = new scalar;
    scalar b = new scalar;
    free (a.name);
    free (b.name);
    char aname[20];
    char bname[20];
    sprintf (aname, "slexp_%s", liq_species[jj]);
    sprintf (bname, "slimp_%s", liq_species[jj]);
    a.name = strdup (aname);
    b.name = strdup (bname);
    a.nodump = true;
    b.nodump = true;
    slexpList = list_append (slexpList, a);
    slimpList = list_append (slimpList, b);
  }
  reset (slexpList, 0.);
  reset (slimpList, 0.);

  for (int jj=0; jj<NGS; jj++) {
    scalar a = new scalar;
    scalar b = new scalar;
    free (a.name);
    free (b.name);
    char aname[20];
    char bname[20];
    sprintf (aname, "sgexp_%s", gas_species[jj]);
    sprintf (bname, "sgimp_%s", gas_species[jj]);
    a.name = strdup (aname);
    b.name = strdup (bname);
    a.nodump = true;
    b.nodump = true;
    sgexpList = list_append (sgexpList, a);
    sgimpList = list_append (sgimpList, b);
  }
  reset (sgexpList, 0.);
  reset (sgimpList, 0.);

#ifdef VARPROP
  for (int jj=0; jj<NLS; jj++) {
    scalar a = new scalar;
    free (a.name);
    char name[20];
    sprintf (name, "Dmix1_%s", liq_species[jj]);
    a.name = strdup (name);
    //a.nodump = true;
    Dmix1List = list_append (Dmix1List, a);
  }
  for (int jj=0; jj<NGS; jj++) {
    scalar a = new scalar;
    free (a.name);
    char name[20];
    sprintf (name, "Dmix2_%s", gas_species[jj]);
    a.name = strdup (name);
    //a.nodump = true;
    Dmix2List = list_append (Dmix2List, a);
  }
  reset (Dmix1List, 0.);
  reset (Dmix2List, 0.);

  for (int jj=0; jj<NLS; jj++) {
    scalar a = new scalar;
    free (a.name);
    char name[20];
    sprintf (name, "dhev_%s", liq_species[jj]);
    a.name = strdup (name);
    //a.nodump = true;
    dhevList = list_append (dhevList, a);
  }
  reset (dhevList, 0.);

  for (int jj=0; jj<NLS; jj++) {
    scalar a = new scalar;
    free (a.name);
    char name[20];
    sprintf (name, "Cp1_%s", liq_species[jj]);
    a.name = strdup (name);
    a.nodump = true;
    Cp1List = list_append (Cp1List, a);
  }
  for (int jj=0; jj<NGS; jj++) {
    scalar a = new scalar;
    free (a.name);
    char name[20];
    sprintf (name, "Cp2_%s", gas_species[jj]);
    a.name = strdup (name);
    a.nodump = true;
    Cp2List = list_append (Cp2List, a);
  }
  reset (Cp1List, 0.);
  reset (Cp2List, 0.);
#endif
#ifdef MOLAR_DIFFUSION
  for (int jj=0; jj<NLS; jj++) {
    scalar a = new scalar;
    free (a.name);
    char name[20];
    sprintf (name, "XL_%s", liq_species[jj]);
    a.name = strdup (name);
    a.nodump = true;
    //a.dirty = false;
    a.dirty = true;
    XLList = list_append (XLList, a);
  }
  reset (XLList, 0.);

  for (int jj=0; jj<NGS; jj++) {
    scalar a = new scalar;
    free (a.name);
    char name[20];
    sprintf (name, "XG_%s", gas_species[jj]);
    a.name = strdup (name);
    a.nodump = true;
    //a.dirty = false;
    a.dirty = true;
    XGList = list_append (XGList, a);
  }
  reset (XGList, 0.);

  for (int jj=0; jj<NLS; jj++) {
    scalar a = new scalar;
    free (a.name);
    char name[20];
    sprintf (name, "XLInt_%s", liq_species[jj]);
    a.name = strdup (name);
    a.nodump = true;
    XLIntList = list_append (XLIntList, a);
  }
  reset (XLIntList, 0.);

  for (int jj=0; jj<NGS; jj++) {
    scalar a = new scalar;
    free (a.name);
    char name[20];
    sprintf (name, "XGInt_%s", gas_species[jj]);
    a.name = strdup (name);
    a.nodump = true;
    XGIntList = list_append (XGIntList, a);
  }
  reset (XGIntList, 0.);
#endif

  fL.nodump = true;
  fG.nodump = true;

  for (scalar s in YGList)
    s.inverse = true;

  for (scalar s in YLList)
    s.inverse = false;

#ifdef MOLAR_DIFFUSION
  for (scalar s in XGList)
    s.inverse = true;

  for (scalar s in XLList)
    s.inverse = false;
#endif

  /**
  We restore to NULL the tracers associated with fu
  and fuext in order to re-run simulations correctly. */

  fuext.tracers = NULL;
  fu.tracers    = NULL;

#ifdef CONSISTENTPHASE1
  fuext.tracers = list_concat (fuext.tracers, YLList);
#else
  fu.tracers = list_concat (fu.tracers, YLList);
#endif
#ifdef CONSISTENTPHASE2
  fuext.tracers = list_concat (fuext.tracers, YGList);
#else
  fu.tracers = list_concat (fu.tracers, YGList);
#endif

  /**
  On adaptive meshes, tracers need to use linear interpolation (rather
  than the default bilinear interpolation) to ensure conservation when
  refining cells. */

#if TREE
    for (scalar s in YGList) {
#if EMBED
      s.refine = s.prolongation = refine_embed_linear;
#else
      s.refine  = refine_linear;
#endif
      s.restriction = restriction_volume_average;
      s.dirty = true; // boundary conditions need to be updated
    }
#endif

#if TREE
    for (scalar s in YLList) {
#if EMBED
      s.refine = s.prolongation = refine_embed_linear;
#else
      s.refine  = refine_linear;
#endif
      s.restriction = restriction_volume_average;
      s.dirty = true; // boundary conditions need to be updated
    }
#endif

#ifdef SOLVE_TEMPERATURE

  TL = new scalar;
  TG = new scalar;

  TL.inverse = false;
  TG.inverse = true;

  sgT.nodump    = true;
  slT.nodump    = true;
  sgTimp.nodump = true;
  slTimp.nodump = true;

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

# if TREE
  TL.refine  = refine_linear;
  TL.restriction = restriction_volume_average;
  TL.dirty = true; // boundary conditions need to be updated
  TG.refine  = refine_linear;
  TG.restriction = restriction_volume_average;
  TG.dirty = true; // boundary conditions need to be updated
# endif
#endif

  /**
  Set default multicomponent properties. */

  for (int jj=0; jj<NGS; jj++)
    inMW[jj] = 1.;

#ifdef VARPROP
  if (is_constant (rho2v)) {
    scalar * l = list_copy (all);
    rho2v = new scalar;
    free (all);
    all = list_concat ({rho2v}, l);
    free (l);
  }
  if (is_constant (mu2v)) {
    scalar * l = list_copy (all);
    mu2v = new scalar;
    free (all);
    all = list_concat ({mu2v}, l);
    free (l);
  }
  if (is_constant (cp2v)) {
    scalar * l = list_copy (all);
    cp2v = new scalar;
    free (all);
    all = list_concat ({cp2v}, l);
    free (l);
  }
  if (is_constant (lambda2v)) {
    scalar * l = list_copy (all);
    lambda2v = new scalar;
    free (all);
    all = list_concat ({lambda2v}, l);
    free (l);
  }
  if (is_constant (rho1v)) {
    scalar * l = list_copy (all);
    rho1v = new scalar;
    free (all);
    all = list_concat ({rho1v}, l);
    free (l);
  }
  if (is_constant (mu1v)) {
    scalar * l = list_copy (all);
    mu1v = new scalar;
    free (all);
    all = list_concat ({mu1v}, l);
    free (l);
  }
  if (is_constant (cp1v)) {
    scalar * l = list_copy (all);
    cp1v = new scalar;
    free (all);
    all = list_concat ({cp1v}, l);
    free (l);
  }
  if (is_constant (lambda1v)) {
    scalar * l = list_copy (all);
    lambda1v = new scalar;
    free (all);
    all = list_concat ({lambda1v}, l);
    free (l);
  }
#endif
}

/**
## Init

In the init event, We initialize the temperature and
chemical species and set to zero additional fields. */

event init (i = 0)
{
  foreach() {
    for (int jj=0; jj<NLS; jj++) {
      scalar s = YLList[jj];
      if (init_fields)
        s[] = f[]*liq_start[jj];
      else
        s[] = 0.;
    }
    for (int jj=0; jj<NGS; jj++) {
      scalar s = YGList[jj];
      if (init_fields)
        s[] = (1. - f[])*gas_start[jj];
      else
        s[] = 0.;
    }
    for (int jj=0; jj<NLS; jj++) {
      scalar s = YLIntList[jj];
      s[] = 0.;
    }
    for (int jj=0; jj<NGS; jj++) {
      scalar s = YGIntList[jj];
      s[] = 0.;
    }
    for (int jj=0; jj<NGS; jj++) {
      scalar s = mEvapList[jj];
      s[] = 0.;
    }
    for (int jj=0; jj<NGOS; jj++) {
      scalar s  = YList[GOSI[jj]];
      scalar sg = YGList[GOSI[jj]];
      s[] = sg[];
    }
    for (int jj=0; jj<NLS; jj++) {
      scalar s  = YList[LSI[jj]];
      scalar sg = YGList[LSI[jj]];
      scalar sl = YLList[jj];
      s[] = sg[] + sl[];
    }
    for (int jj=0; jj<NLS; jj++) {
      scalar slexp = slexpList[jj];
      scalar slimp = slimpList[jj];
      slexp[] = 0.;
      slimp[] = 0.;
    }
    for (int jj=0; jj<NGS; jj++) {
      scalar sgexp = sgexpList[jj];
      scalar sgimp = sgimpList[jj];
      sgexp[] = 0.;
      sgimp[] = 0.;
    }
    for (int jj=0; jj<NGS; jj++) {
      scalar mEvap = mEvapList[jj];
      mEvap[] = 0.;
    }
  }

#ifdef SOLVE_TEMPERATURE
  if (init_fields) {
    foreach() {
      TL[] = TL0*f[];
      TG[] = TG0*(1. - f[]);
      T[]  = TL[] + TG[];
    }
  }
#endif
}

/**
## Cleanup

We deallocate the various lists from the memory. */

event cleanup (t = end)
{
  delete (mEvapList), free (mEvapList), mEvapList = NULL;
  delete (YList), free (YList), YList = NULL;
  delete (YLList), free (YLList), YLList = NULL;
  delete (YGList), free (YGList), YGList = NULL;
  delete (YLIntList), free (YLIntList), YLIntList = NULL;
  delete (YGIntList), free (YGIntList), YGIntList = NULL;
  delete (slexpList), free (slexpList), slexpList = NULL;
  delete (slimpList), free (slimpList), slimpList = NULL;
  delete (sgexpList), free (sgexpList), sgexpList = NULL;
  delete (sgimpList), free (sgimpList), sgimpList = NULL;
#ifdef VARPROP
  delete (Dmix1List), free (Dmix1List), Dmix1List = NULL;
  delete (Dmix2List), free (Dmix2List), Dmix2List = NULL;
  delete (dhevList), free (dhevList), dhevList = NULL;
  delete (Cp1List), free (Cp1List), Cp1List = NULL;
  delete (Cp2List), free (Cp2List), Cp2List = NULL;
#endif
#ifdef MOLAR_DIFFUSION
  delete (XLList), free (XLList), XLList = NULL;
  delete (XGList), free (XGList), XGList = NULL;
  delete (XLIntList), free (XLIntList), XLIntList = NULL;
  delete (XGIntList), free (XGIntList), XGIntList = NULL;
#endif
  delete (fu.tracers), free (fu.tracers), fu.tracers = NULL;
  delete (fuext.tracers), free (fuext.tracers), fuext.tracers = NULL;
  free (LSI);
  free (GOSI);
#ifdef SOLVE_TEMPERATURE
  delete ({TL,TG});
#endif
}

/**
## Reset Source Terms

We set to zero phase change source terms, in order to add other
source term contributions from outside this module. */

event reset_sources (i++)
{
  foreach() {
#ifdef SOLVE_TEMPERATURE
    sgT[] = 0.;
    slT[] = 0.;
    sgTimp[] = 0.;
    slTimp[] = 0.;
#endif
    for (int jj=0; jj<NLS; jj++) {
      scalar slexp = slexpList[jj];
      scalar slimp = slimpList[jj];
      slexp[] = 0.;
      slimp[] = 0.;
    }
    for (int jj=0; jj<NGS; jj++) {
      scalar sgexp = sgexpList[jj];
      scalar sgimp = sgimpList[jj];
      sgexp[] = 0.;
      sgimp[] = 0.;
    }
  }
}

/**
## Helper functions

We define a function which is used to update the mixture molecular
weight and the mole fractions (if needed). */

#include "multicomponent-properties.h"

void update_mw_moles (void) {
  double MW1[NLS], MW2[NGS];
  for (int jj=0; jj<NLS; jj++)
    MW1[jj] = inMW[LSI[jj]];
  for (int jj=0; jj<NGS; jj++)
    MW2[jj] = inMW[jj];

//#ifdef MOLAR_DIFFUSION
//  for (scalar s in XLList)
//    s.dirty = false;
//  for (scalar s in XGList)
//    s.dirty = false;
//#endif

  foreach() {
    double x1[NLS], y1[NLS];
    for (int jj=0; jj<NLS; jj++) {
      scalar YL = YLList[jj];
      y1[jj] = (NLS == 1.) ? 1. : YL[];
    }
    correctfrac (y1, NLS);
    mass2molefrac (x1, y1, MW1, NLS);
    MW1mix[] = mass2mw (y1, MW1, NLS);
#ifdef MOLAR_DIFFUSION
    for (int jj=0; jj<NLS; jj++) {
      scalar XL = XLList[jj];
      XL[] = x1[jj];
    }
#endif

    double x2[NGS], y2[NGS];
    for (int jj=0; jj<NGS; jj++) {
      scalar YG = YGList[jj];
      y2[jj] = (NGS == 1.) ? 1. : YG[];
    }
    correctfrac (y2, NGS);
    mass2molefrac (x2, y2, MW2, NGS);
    MW2mix[] = mass2mw (y2, MW2, NGS);
#ifdef MOLAR_DIFFUSION
    for (int jj=0; jj<NGS; jj++) {
      scalar XG = XGList[jj];
      XG[] = x2[jj];
    }
#endif
  }

  double MW_TOL = T_PROP;

  foreach() {
    if (f[] <= MW_TOL) {
      double MW1mixvgh = 0.;

      int counter = 0;
      foreach_neighbor(1) {
        if (f[] > MW_TOL) {
          counter++;
          MW1mixvgh += MW1mix[];
        }
      }
      MW1mix[] = (counter != 0.) ? MW1mixvgh/counter : 0.;
    }
  }

  foreach() {
    if ((1. - f[]) <= MW_TOL) {
      double MW2mixvgh = 0.;

      int counter = 0;
      foreach_neighbor(1) {
        if ((1. - f[]) > MW_TOL) {
          counter++;
          MW2mixvgh += MW2mix[];
        }
      }
      MW2mix[] = (counter != 0.) ? MW2mixvgh/counter : 0.;
    }
  }

  foreach() {
    if (f[] > F_ERR) {
      double x1[NLS], y1[NLS];
      for (int jj=0; jj<NLS; jj++) {
        scalar YL = YLList[jj];
        y1[jj] = (NLS == 1) ? 1. : YL[];
      }
      correctfrac (y1, NLS);
      mass2molefrac (x1, y1, MW1, NLS);
      MW1mix[] = mass2mw (y1, MW1, NLS);
#ifdef MOLAR_DIFFUSION
      for (int jj=0; jj<NLS; jj++) {
        scalar XL = XLList[jj];
        XL[] = x1[jj];
      }
#endif
    }

    if ((1. - f[]) > F_ERR) {
      double x2[NGS], y2[NGS];
      for (int jj=0; jj<NGS; jj++) {
        scalar YG = YGList[jj];
        y2[jj] = YG[];
      }
      correctfrac (y2, NGS);
      mass2molefrac (x2, y2, MW2, NGS);
      MW2mix[] = mass2mw (y2, MW2, NGS);
#ifdef MOLAR_DIFFUSION
      for (int jj=0; jj<NGS; jj++) {
        scalar XG = XGList[jj];
        XG[] = x2[jj];
      }
#endif
    }
  }

#ifdef MOLAR_DIFFUSION
  boundary (XGList);
  boundary (XLList);
#endif
  boundary ({MW1mix,MW2mix});

//  for (int b = 0; b < nboundary; b++) {
//    foreach_boundary (b) {
//
//      double x1[NLS], y1[NLS];
//      for (int jj=0; jj<NLS; jj++) {
//        scalar YL = YLList[jj];
//        y1[jj] = 0.5*(YL[] + get_ghost (point, YL, b));
//      }
//      correctfrac (y1, NLS);
//      mass2molefrac (x1, y1, MW1, NLS);
//      double MW1face = mass2mw (y1, MW1, NLS);
//      set_ghost (point, MW1mix, b, MW1face);
//#ifdef MOLAR_DIFFUSION
//      for (int jj=0; jj<NLS; jj++) {
//        scalar XL = XLList[jj];
//        set_ghost (point, XL, b, x1[jj]);
//      }
//#endif
//
//      double x2[NGS], y2[NGS];
//      for (int jj=0; jj<NGS; jj++) {
//        scalar YG = YGList[jj];
//        y2[jj] = 0.5*(YG[] + get_ghost (point, YG, b));
//      }
//      correctfrac (y2, NGS);
//      mass2molefrac (x2, y2, MW2, NGS);
//      double MW2face = mass2mw (y2, MW2, NGS);
//      set_ghost (point, MW2mix, b, MW2face);
//#ifdef MOLAR_DIFFUSION
//      for (int jj=0; jj<NGS; jj++) {
//        scalar XG = XGList[jj];
//        set_ghost (point, XG, b, x2[jj]);
//      }
//#endif
//    }
//  }
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

#ifdef SOLVE_TEMPERATURE
    TL[] = f[] > F_ERR ? TL[]/f[] : 0.;
    TG[] = ((1. - f[]) > F_ERR) ? TG[]/(1. - f[]) : 0.;
#endif

    /**
    We recover value of mass fractions from their tracer form. */

    for (int jj=0; jj<NLS; jj++) {
      scalar YL = YLList[jj];
      YL[] = f[] > F_ERR ? YL[]/f[] : 0.;
    }

    for (int jj=0; jj<NGS; jj++) {
      scalar YG = YGList[jj];
      YG[] = ((1. - f[]) > F_ERR) ? YG[]/(1. - f[]) : 0.;
    }
  }

  /**
  We compute the mole fraction fields if the diffusivity is
  molar-based. */

#ifdef MOLAR_DIFFUSION
  update_mw_moles();
#endif

  /**
  The thermodynamic and transport properties are updated at the
  beginning of the each time-step. */

#if defined (VARPROP) && !defined (NO_UPDATE_PROPERTIES)
  update_properties();
#endif

  /**
  We compute the value of volume fraction *f* on the
  cell-faces using a geometric approach (necessary
  for interface gradients and diffusion equations). */

  face_fraction (fL, fsL);
  face_fraction (fG, fsG);

  /**
  We assign the interface temperature value. */

#ifdef SOLVE_TEMPERATURE
  foreach() {
    TInt[] = 0.;
    if (f[] > F_ERR && f[] < 1.-F_ERR)
#ifdef FIXED_INTERFACE_TEMPERATURE
      TInt[] = FIXED_INTERFACE_TEMPERATURE;
#else
      TInt[] = avg_neighbor (point, TL, f);
#endif
  }
#endif

  /**
  We compute the molecular weight of the gas-only
  species mixture. */

  scalar MWGmix[];
  foreach() {
    MWGmix[] = 0.;
    if (f[] > F_ERR && f[] < 1.-F_ERR) {
      double yGinorm[NGOS];
      double sumYGi = 0.;
      for (int jj=0; jj<NGOS; jj++) {
        scalar YG = YGList[GOSI[jj]];
        yGinorm[jj] = YG[];
        sumYGi += yGinorm[jj];
      }
      for (int jj=0; jj<NGOS; jj++) {
        yGinorm[jj] /= (sumYGi + 1.e-10);
      }
      for (int jj=0; jj<NGOS; jj++) {
        MWGmix[] += yGinorm[jj] / inMW[GOSI[jj]];
      }
      MWGmix[] = 1. / (MWGmix[] + 1.e-10);
    }
  }

  /**
  We compute total vaporization flowrate. */

  foreach() {
    mEvapTot[] = 0.;

    /**
    We reset to zero mEvap for every species, and we set to zero
    the interface mass fraction fields. */

    for (int jj=0; jj<NGS; jj++) {
      scalar mEvap = mEvapList[jj];
      scalar YGInt = YGIntList[jj];
      mEvap[] = 0.;
      YGInt[] = 0.;
    }
    for (int jj=0; jj<NLS; jj++) {
      scalar YLInt = YLIntList[jj];
      YLInt[] = 0.;
    }

    if (f[] > F_ERR && f[] < 1.-F_ERR) {

      /**
      We create fields to store the local mole fractions, the
      conversion is required by the thermodynamic equilibrium. */

      double XGIntConv[NLS+1], YGIntConv[NLS+1];
      double XLIntConv[NLS], YLIntConv[NLS];
      double inMWG[NLS+1], inMWL[NLS];

      for (int jj=0; jj<NLS; jj++) {
        inMWL[jj] = inMW[LSI[jj]];
        inMWG[jj] = inMW[LSI[jj]];
        XGIntConv[jj] = 0.;
        YGIntConv[jj] = 0.;
        XLIntConv[jj] = 0.;
        YLIntConv[jj] = 0.;
      }

      /**
      We convert *YLInt* to mole fractions *XLInt*. */

      for (int jj=0; jj<NLS; jj++) {
        scalar YL    = YLList[jj];
        scalar YLInt = YLIntList[jj];
        YLIntConv[jj] = avg_neighbor (point, YL, f);
        YLInt[] = YLIntConv[jj];
      }
      mass2molefrac (XLIntConv, YLIntConv, inMWL, NLS);
#ifdef MOLAR_DIFFUSION
      for (int jj=0; jj<NLS; jj++) {
        scalar XLInt = XLIntList[jj];
        XLInt[] = XLIntConv[jj];
      }
#endif

      /**
      We compute *XGInt* from the thermodynamic equilibrium.
      Different equilibrium options are available: constant
      thermodynamic equilibrium constant, Clausius-Clapeyron
      relation, and Antoine equation. */

      double sumXGi = 0.;
      for (int jj=0; jj<NLS; jj++) {
        XGIntConv[jj] = inKeq[jj]*XLIntConv[jj];
#ifdef USE_CLAPEYRON
        double dhevvh = dhev;
# ifdef VARPROP
        scalar dhevjj = dhevList[jj];
        dhevvh = dhevjj[];
# endif
        XGIntConv[jj] = clapeyron (min (TInt[], Tboil[jj]-1.),
            Tboil[jj], dhevvh, inMW[jj])*XLIntConv[jj];
#endif
#ifdef USE_ANTOINE
        scalar YL = YLList[jj];
        XGIntConv[jj] = min (YL.antoine (TInt[], Pref), 0.98)*XLIntConv[jj];
#endif
#ifdef USE_ANTOINE_OPENSMOKE
        XGIntConv[jj] = min (opensmoke_antoine (TInt[], Pref, jj), 0.98)*XLIntConv[jj];
#endif
        sumXGi += XGIntConv[jj];
      }
      XGIntConv[NLS] = 1. - sumXGi;
      inMWG[NLS] = MWGmix[];
      mole2massfrac (YGIntConv, XGIntConv, inMWG, NLS+1);

      /**
      We set the gas phase interface mass fraction values using
      the converted fractions. */

      for (int jj=0; jj<NLS; jj++) {
        scalar YGInt = YGIntList[LSI[jj]];
        YGInt[] = YGIntConv[jj];
#ifdef MOLAR_DIFFUSION
        scalar XGInt = XGIntList[LSI[jj]];
        XGInt[] = XGIntConv[jj];
#endif
      }

      /**
      We adjust the interface mass fractions of the gas-only
      species in the system. */

      double yGinorm[NGOS];
      double sumYGi = 0., sumYGl = 0.;
      if (f[] > F_ERR && f[] < 1.-F_ERR) {
        for (int jj=0; jj<NLS; jj++) {
          sumYGl += YGIntConv[jj];
        }
        for (int jj=0; jj<NGOS; jj++) {
          scalar YG = YGList[GOSI[jj]];
          yGinorm[jj] = YG[];
          sumYGi += yGinorm[jj];
        }
        for (int jj=0; jj<NGOS; jj++) {
          scalar YGInt = YGIntList[GOSI[jj]];
          YGInt[] = (1. - sumYGl)*yGinorm[jj]/(sumYGi + 1.e-10);
        }
      }

      /**
      We repeat the same operations to close the molar fractions of the
      gas-only species in the system. */

#ifdef MOLAR_DIFFUSION
      double xGinorm[NGOS];
      double sumXGigo = 0., sumXGl = 0.;
      if (f[] > F_ERR && f[] < 1.-F_ERR) {
        for (int jj=0; jj<NLS; jj++) {
          sumXGl += XGIntConv[jj];
        }
        for (int jj=0; jj<NGOS; jj++) {
          scalar XG = XGList[GOSI[jj]];
          xGinorm[jj] = XG[];
          sumXGigo += xGinorm[jj];
        }
        for (int jj=0; jj<NGOS; jj++) {
          scalar XGInt = XGIntList[GOSI[jj]];
          XGInt[] = (1. - sumXGl)*xGinorm[jj]/(sumXGigo + 1.e-10);
        }
      }
#endif

      /**
      We compute the sum of the diffusive fluxes in gas phase, to be
      used for the Fick corrected approach. */

      double jGtot = 0.;
#ifdef FICK_CORRECTED
      for (int jj=0; jj<NGS; jj++) {
        double rho2vh = rho2;
        double inDmix2vh = inDmix2[jj];
# ifdef VARPROP
        rho2vh = rho2v[];
        scalar Dmix2v = Dmix2List[jj];
        inDmix2vh = Dmix2v[];
# endif
# ifdef MOLAR_DIFFUSION
        scalar XGInt = XGIntList[jj];
        scalar XG    = XGList[jj];
        double gtrgrad = ebmgrad (point, XG, fL, fG, fsL, fsG, true, XGInt[], &success);
        gtrgrad *= (MW2mix[] > 0.) ? inMW[jj]/MW2mix[] : 0.;
# else
        scalar YGInt = YGIntList[jj];
        scalar YG    = YGList[jj];
        double gtrgrad = ebmgrad (point, YG, fL, fG, fsL, fsG, true, YGInt[], &success);
# endif
        jGtot += -rho2vh*inDmix2vh*gtrgrad;
      }
#endif

      /**
      We compute the total vaporization rate from the sum
      of the interface mass balances over all the chemical
      species in liquid phase. */

      double sum_jG = 0., sum_YGInt = 0.;

      for (int jj=0; jj<NLS; jj++) {
        //scalar YGInt = YGIntList[LSI[jj]];
        //scalar YLInt = YLIntList[jj];
        //scalar YL    = YLList[jj];
        //scalar YG    = YGList[LSI[jj]];

        //YLInt[] = avg_neighbor (point, YL, f);
        //YGInt[] = YGIntConv[jj];

        double rho2vh = rho2;
        double inDmix2vh = inDmix2[LSI[jj]];
#ifdef VARPROP
        rho2vh = rho2v[];
        scalar Dmix2v = Dmix2List[LSI[jj]];
        inDmix2vh = Dmix2v[];
#endif

        scalar YGInt = YGIntList[LSI[jj]];
#ifdef MOLAR_DIFFUSION
        scalar XGInt = XGIntList[LSI[jj]];
        scalar XG    = XGList[LSI[jj]];
        double gtrgrad = ebmgrad (point, XG, fL, fG, fsL, fsG, true, XGInt[], &success);
        gtrgrad *= (MW2mix[] > 0.) ? inMW[LSI[jj]]/MW2mix[] : 0.;
#else
        scalar YG    = YGList[LSI[jj]];
        double gtrgrad = ebmgrad (point, YG, fL, fG, fsL, fsG, true, YGInt[], &success);
#endif
        sum_jG += -rho2vh*inDmix2vh*gtrgrad - YGInt[]*jGtot;
        sum_YGInt += YGInt[];
      }
#ifdef DIFFUSIVE
      mEvapTot[] = sum_jG;
#else
      mEvapTot[] = sum_jG/min(1. - sum_YGInt, 0.99);
#endif
    }
  }

  /**
  From the knowledge of the total vaporization rate, we
  compute the vaporization rate for each species. */

  foreach() {
    if (f[] > F_ERR && f[] < 1.-F_ERR) {

      double jGtot = 0.;
#ifdef FICK_CORRECTED
      for (int jj=0; jj<NGS; jj++) {
        double rho2vh = rho2;
        double inDmix2vh = inDmix2[jj];
# ifdef VARPROP
        rho2vh = rho2v[];
        scalar Dmix2v = Dmix2List[jj];
        inDmix2vh = Dmix2v[];
# endif

# ifdef MOLAR_DIFFUSION
        scalar XGInt = XGIntList[jj];
        scalar XG    = XGList[jj];
        double gtrgrad = ebmgrad (point, XG, fL, fG, fsL, fsG, true, XGInt[], &success);
        gtrgrad *= (MW2mix[] > 0.) ? inMW[jj]/MW2mix[] : 0.;
# else
        scalar YGInt = YGIntList[jj];
        scalar YG    = YGList[jj];
        double gtrgrad = ebmgrad (point, YG, fL, fG, fsL, fsG, true, YGInt[], &success);
# endif
        jGtot += -rho2vh*inDmix2vh*gtrgrad;
      }
#endif

      for (int jj=0; jj<NLS; jj++) {
        scalar mEvap = mEvapList[LSI[jj]];
        scalar YGInt = YGIntList[LSI[jj]];

        double rho2vh = rho2;
        double inDmix2vh = inDmix2[LSI[jj]];
#ifdef VARPROP
        rho2vh = rho2v[];
        scalar Dmix2v = Dmix2List[LSI[jj]];
        inDmix2vh = Dmix2v[];
#endif

#ifdef MOLAR_DIFFUSION
        scalar XGInt = XGIntList[LSI[jj]];
        scalar XG    = XGList[LSI[jj]];
        double gtrgrad = ebmgrad (point, XG, fL, fG, fsL, fsG, true, XGInt[], &success);
        gtrgrad *= (MW2mix[] > 0.) ? inMW[LSI[jj]]/MW2mix[] : 0.;
#else
        scalar YG    = YGList[LSI[jj]];
        double gtrgrad = ebmgrad (point, YG, fL, fG, fsL, fsG, true, YGInt[], &success);
#endif

#ifdef DIFFUSIVE
        mEvap[] = - rho2vh*inDmix2vh*gtrgrad - YGInt[]*jGtot;
#else
        mEvap[] = mEvapTot[]*YGInt[] - rho2vh*inDmix2vh*gtrgrad - YGInt[]*jGtot;
#endif
      }
    }
  }

  /**
  If *SOLVE_TEMPERATURE*, the value of interface temperature
  is obtained from the total vaporization rate just computed
  and the interface temperature gradients. This result in an
  algebraic equation which is numerically zero-ed. */

#ifdef USE_GSL

# ifdef SOLVE_TEMPERATURE
  ijc_CoupledTemperature();
# endif

  /**
  Finally, we can refine the first guess values from the decoupled
  solution from the fully-coupled interface jump condition. Solved
  as a non-linear system of equations. */

# if USE_GSL > 0
  ijc_CoupledNls();
# endif

#endif

  /**
  The source terms for the diffusion equation of the species
  mass fractions in gas an liquid phase are computed here. */

  foreach() {

    if (f[] > F_ERR && f[] < 1.-F_ERR) {
      coord n = facet_normal (point, fL, fsL), p;
      double alpha = plane_alpha (fL[], n);
      double area = plane_area_center (n, alpha, &p);
      normalize (&n);

      for (int jj=0; jj<NGS; jj++) {
        scalar sgexp = sgexpList[jj];
        scalar sgimp = sgimpList[jj];
        scalar mEvap = mEvapList[jj];

        //scalar YGInt = YGIntList[jj];
#ifdef AXI
        //sgexp[] += -(mEvap[] - mEvapTot[]*YGInt[])*area*(y + p.y*Delta)/(Delta*y)*cm[];
        //sgimp[] += 0.;
        sgexp[] += -mEvap[]*area*(y + p.y*Delta)/(Delta*y)*cm[];
        sgimp[] += +mEvapTot[]*area*(y + p.y*Delta)/(Delta*y)*cm[];
#else
        //sgexp[] += -(mEvap[] - mEvapTot[]*YGInt[])*area/Delta*cm[];
        //sgimp[] += 0.;
        sgexp[] += -mEvap[]*area/Delta*cm[];
        sgimp[] += +mEvapTot[]*area/Delta*cm[];
#endif
      }

      for (int jj=0; jj<NLS; jj++) {
        scalar slexp = slexpList[jj];
        scalar slimp = slimpList[jj];
        scalar mEvap = mEvapList[LSI[jj]];

        //scalar YLInt = YLIntList[jj];
#ifdef AXI
        //slexp[] += +(mEvap[] - mEvapTot[]*YLInt[])*area*(y + p.y*Delta)/(Delta*y)*cm[];
        //slimp[] += 0.;
        slexp[] += +mEvap[]*area*(y + p.y*Delta)/(Delta*y)*cm[];
        slimp[] += -mEvapTot[]*area*(y + p.y*Delta)/(Delta*y)*cm[];
#else
        //slexp[] += +(mEvap[] - mEvapTot[]*YLInt[])*area/Delta*cm[];
        //slimp[] += 0.;
        slexp[] += +mEvap[]*area/Delta*cm[];
        slimp[] += -mEvapTot[]*area/Delta*cm[];
#endif
      }

#ifdef SOLVE_TEMPERATURE
      double bc = TInt[];
      double gtrgrad = ebmgrad (point, TG, fL, fG, fsL, fsG, true, bc, &success);
      double ltrgrad = ebmgrad (point, TL, fL, fG, fsL, fsG, false, bc, &success);

      double lambda1vh = lambda1;
      double lambda2vh = lambda2;
# ifdef VARPROP
      lambda1vh = lambda1v[];
      lambda2vh = lambda2v[];
# endif
      double lheatflux = lambda1vh*ltrgrad;
      double gheatflux = lambda2vh*gtrgrad;

# ifndef VARPROP  // Without VARPROP this form converges better
      lheatflux /= (rho1*cp1);
      gheatflux /= (rho2*cp2);
# endif

# ifdef AXI
      slT[] += lheatflux*area*(y + p.y*Delta)/(Delta*y)*cm[];
      sgT[] += gheatflux*area*(y + p.y*Delta)/(Delta*y)*cm[];
# else
      slT[] += lheatflux*area/Delta*cm[];
      sgT[] += gheatflux*area/Delta*cm[];
# endif
#endif
    }
  }

  /**
  We include the velocity correction for the diffusive fluxes.
  Using the FICK_CORRECTED approach this method enforces mass
  conservation in multicomponent diffusion.
  Using MOLAR_DIFFUSION, this procedure includes a missing term
  in the implicit diffusion step, which is a function of the
  mixture molecular weight. */

#ifdef MASS_DIFFUSION_ENTHALPY
  foreach() {
    double mde1 = 0., mde2 = 0.;
    coord gTL = {0., 0., 0.,}, gTG = {0., 0., 0.};
    coord gYLj = {0., 0., 0.}, gYGj = {0., 0., 0.};
    coord gYLsum = {0., 0., 0.}, gYGsum = {0., 0., 0.};

    foreach_dimension() {
      gTG.x = (TG[1] - TG[-1])/(2.*Delta);
      gTL.x = (TL[1] - TL[-1])/(2.*Delta);
    }

    foreach_dimension() {
      for (int jj=0; jj<NGS; jj++) {
        scalar Dmix2v = Dmix2List[jj];
# ifdef MOLAR_DIFFUSION
        scalar XG = XGList[jj];
        gYGsum.x -= (MW2mix[] > 0.) ?
          rho2v[]*Dmix2v[]*inMW[jj]/MW2mix[]*(XG[1] - XG[-1])/(2.*Delta) : 0.;
# else
        scalar YG = YGList[jj];
        gYGsum.x -= rho2v[]*Dmix2v[]*(YG[1] - YG[-1])/(2.*Delta);
# endif
      }

      for (int jj=0; jj<NLS; jj++) {
        scalar Dmix1v = Dmix1List[jj];
# ifdef MOLAR_DIFFUSION
        scalar XL = XLList[jj];
        gYLsum.x -= (MW1mix[] > 0.) ?
          rho1v[]*Dmix1v[]*inMW[LSI[jj]]/MW1mix[]*(XL[1] - XL[-1])/(2.*Delta) : 0.;
# else
        scalar YL = YLList[jj];
        gYLsum.x -= rho1v[]*Dmix1v[]*(YL[1] - YL[-1])/(2.*Delta);
# endif
      }

      for (int jj=0; jj<NGS; jj++) {
        scalar YG = YGList[jj];
        scalar Cp2v = Cp2List[jj];
        scalar Dmix2v = Dmix2List[jj];
# ifdef MOLAR_DIFFUSION
        scalar XG = XGList[jj];
        gYGj.x = (MW2mix[] > 0.) ?
          -rho2v[]*Dmix2v[]*inMW[jj]/MW2mix[]*(XG[1] - XG[-1])/(2.*Delta) : 0.;
# else
        gYGj.x = -rho2v[]*Dmix2v[]*(YG[1] - YG[-1])/(2.*Delta);
# endif
        mde2 += Cp2v[]*(gYGj.x - YG[]*gYGsum.x)*gTG.x;
      }

      for (int jj=0; jj<NLS; jj++) {
        scalar YL = YLList[jj];
        scalar Cp1v = Cp1List[jj];
        scalar Dmix1v = Dmix1List[jj];
# ifdef MOLAR_DIFFUSION
        scalar XL = XLList[jj];
        gYLj.x = (MW1mix[] > 0.) ?
          -rho1v[]*Dmix1v[]*inMW[LSI[jj]]/MW1mix[]*(XL[1] - XL[-1])/(2.*Delta) : 0.;
# else
        gYLj.x = -rho1v[]*Dmix1v[]*(YL[1] - YL[-1])/(2.*Delta);
# endif
        mde1 += Cp1v[]*(gYLj.x - YL[]*gYLsum.x)*gTL.x;
      }
    }
    slT[] -= mde1*cm[]*(f[] > (1. - F_ERR));
    sgT[] -= mde2*cm[]*(f[] < F_ERR);
  }
#endif

  /**
  The divergence of the velocity field is computed after the calculation
  of source terms for species and temperature. */

#if defined (VARPROP) && !defined (NO_DIVERGENCE)
  update_divergence();
  //update_divergence_density();
#endif

  /**
  We restore the tracer form of the liquid and gas-phase
  mass fraction fields. */

  foreach() {
    for (scalar YL in YLList)
      YL[] *= f[]*(f[] > F_ERR);
    for (scalar YG in YGList) 
      YG[] *= (1. - f[])*((1. - f[]) > F_ERR);
#ifdef SOLVE_TEMPERATURE
    TL[] *= f[]*(f[] > F_ERR);
    TG[] *= (1. - f[])*((1. - f[]) > F_ERR);
#endif
  }
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

We solve the diffusion equations for species and temperature
accounting for the phase change contributions. */

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

    for (scalar YL in YLList)
      YL[] = (fuext[] > F_ERR) ? YL[]/fuext[] : 0.;
    for (scalar YG in YGList)
      YG[] = ((1. - fu[]) > F_ERR) ? YG[]/(1. - fu[]) : 0.;

    fL[] = f[]; fG[] = 1. - f[];

#ifdef SOLVE_TEMPERATURE
    TL[] = (fuext[] > F_ERR) ? TL[]/fuext[] : 0.;
    TG[] = ((1. - fu[]) > F_ERR) ? TG[]/(1. - fu[]) : 0.;
#endif
  }

#ifdef MOLAR_DIFFUSION
  update_mw_moles();
#endif

  /**
  We compute the value of volume fraction *f* on the
  cell-faces using a geometric approach (necessary
  for interface gradients and diffusion equations). */

  face_fraction (fL, fsL);
  face_fraction (fG, fsG);

  /**
  We include the velocity correction for the diffusive fluxes.
  Using the FICK_CORRECTED approach this method enforces mass
  conservation in multicomponent diffusion.
  Using MOLAR_DIFFUSION, this procedure includes a missing term
  in the implicit diffusion step, which is a function of the
  mixture molecular weight. */

  face vector phicGtot[];
  foreach_face() {
    phicGtot.x[] = 0.;
    for (int jj=0; jj<NGS; jj++) {
#ifdef FICK_CORRECTED
      double rho2f = 0.5*(rho2v[] + rho2v[-1]);
# ifdef VARPROP
      scalar Dmix2 = Dmix2List[jj];
      double Dmix2f = 0.5*(Dmix2[] + Dmix2[-1]);
# else
      double Dmix2f = inDmix2[jj];
# endif
# ifdef MOLAR_DIFFUSION
      scalar XG = XGList[jj];
      double MW2mixf = 0.5*(MW2mix[] + MW2mix[-1]);
      phicGtot.x[] += (MW2mixf > 0.) ?
        rho2f*Dmix2f*inMW[jj]/MW2mixf*face_gradient_x (XG, 0) : 0.;
# else
      scalar YG = YGList[jj];
      phicGtot.x[] += rho2f*Dmix2f*face_gradient_x (YG, 0);
# endif // MOLAR_DIFFUSION
#else
      phicGtot.x[] = 0.;
#endif  // FICK_CORRECTED
    }
    phicGtot.x[] *= fsG.x[]*fm.x[];
  }

  face vector phicLtot[];
  foreach_face() {
    phicLtot.x[] = 0.;
    for (int jj=0; jj<NLS; jj++) {
#ifdef FICK_CORRECTED
      double rho1f = 0.5*(rho1v[] + rho1v[-1]);
# ifdef VARPROP
      scalar Dmix1 = Dmix1List[jj];
      double Dmix1f = 0.5*(Dmix1[] + Dmix1[-1]);
# else
      double Dmix1f = inDmix1[jj];
# endif
# ifdef MOLAR_DIFFUSION
      scalar XL = XLList[jj];
      double MW1mixf = 0.5*(MW1mix[] + MW1mix[-1]);
      phicLtot.x[] += (MW1mixf > 0.) ?
        rho1f*Dmix1f*inMW[LSI[jj]]/MW1mixf*face_gradient_x (XL, 0) : 0.;
# else
      scalar YL = YLList[jj];
      phicLtot.x[] += rho1f*Dmix1f*face_gradient_x (YL, 0);
# endif // MOLAR_DIFFUSION
#else
      phicLtot.x[] = 0.;
#endif  // FICK_CORRECTED
    }
    phicLtot.x[] *= fsL.x[]*fm.x[];
  }

  for (int jj=0; jj<NGS; jj++) {
    face vector phicjj[];
    foreach_face() {
      phicjj.x[] = phicGtot.x[];
#ifdef MOLAR_DIFFUSION
      scalar Dmix2 = Dmix2List[jj];
      double Dmix2f = 0.5*(Dmix2[] + Dmix2[-1]);
      double rho2f = 0.5*(rho2v[] + rho2v[-1]);
      double MW2mixf = 0.5*(MW2mix[] + MW2mix[-1]);

      phicjj.x[] -= (MW2mixf > 0.) ?
        rho2f*Dmix2f/MW2mixf*face_gradient_x (MW2mix, 0)*fsG.x[]*fm.x[] : 0.;
#endif  // MOLAR_DIFFUSION
    }
    scalar YG = YGList[jj];
    double (* gradient_backup)(double,double,double) = YG.gradient;
    YG.gradient = zero;
    face vector flux[];
    tracer_fluxes (YG, phicjj, flux, dt, zeroc);
    YG.gradient = gradient_backup;

    foreach()
      foreach_dimension()
        YG[] += (rho2v[] > 0.) ? dt/(rho2v[])*(flux.x[] - flux.x[1])/(Delta*cm[]) : 0.;
  }

  for (int jj=0; jj<NLS; jj++) {
    face vector phicjj[];
    foreach_face() {
      phicjj.x[] = phicLtot.x[];
#ifdef MOLAR_DIFFUSION
      scalar Dmix1 = Dmix1List[jj];
      double Dmix1f = 0.5*(Dmix1[] + Dmix1[-1]);
      double rho1f = 0.5*(rho1v[] + rho1v[-1]);
      double MW1mixf = 0.5*(MW1mix[] + MW1mix[-1]);

      phicjj.x[] -= (MW1mixf > 0.) ?
        rho1f*Dmix1f/MW1mixf*face_gradient_x (MW1mix, 0)*fsL.x[]*fm.x[] : 0.;
#endif  // MOLAR_DIFFUSION
    }
    scalar YL = YLList[jj];
    double (* gradient_backup)(double,double,double) = YL.gradient;
    YL.gradient = zero;
    face vector flux[];
    tracer_fluxes (YL, phicjj, flux, dt, zeroc);
    YL.gradient = gradient_backup;

    foreach()
      foreach_dimension()
        YL[] += (rho1v[] > 0.) ? dt/(rho1v[])*(flux.x[] - flux.x[1])/(Delta*cm[]) : 0.;
  }

  /**
  We solve the diffusion equations, confined by means of
  the face fraction fields *fsL* and *fsG*. */

  scalar theta1[], theta2[];

#if TREE
  theta1.refine = theta1.prolongation = fraction_refine;
  theta2.refine = theta2.prolongation = fraction_refine;
  theta1.dirty = true;
  theta2.dirty = true;
#endif

  for (int jj=0; jj<NLS; jj++) {
    face vector Dmix1f[];
    foreach_face() {
#ifdef VARPROP
      scalar Dmix1jj = Dmix1List[jj];
      double Dmix1vh = 0.5*(Dmix1jj[] + Dmix1jj[-1]);
      double rho1vh = 0.5*(rho1v[] + rho1v[-1]);
      Dmix1f.x[] = rho1vh*Dmix1vh*fsL.x[]*fm.x[];
#else
      Dmix1f.x[] = rho1*inDmix1[jj]*fsL.x[]*fm.x[];
      //Dmix1f.x[] = inDmix1[jj]*fsL.x[]*fm.x[];
#endif
    }

    foreach()
#ifdef VARPROP
      theta1[] = cm[]*max(rho1v[]*fL[], F_ERR);
#else
      theta1[] = cm[]*max(rho1*fL[], F_ERR);
#endif

    scalar YL = YLList[jj];
    scalar slexp = slexpList[jj];
    scalar slimp = slimpList[jj];

    //foreach() {
    //  slexp[] = (f[] > F_ERR) ? slexp[] : 0.;
    //  slimp[] = (f[] > F_ERR) ? slimp[] : 0.;
    //}

    diffusion (YL, dt, D=Dmix1f, r=slexp, beta=slimp, theta=theta1);
  }

  for (int jj=0; jj<NGS; jj++) {

    face vector Dmix2f[];
    foreach_face() {
#ifdef VARPROP
      scalar Dmix2v = Dmix2List[jj];
      double Dmix2vh = 0.5*(Dmix2v[] + Dmix2v[-1]);
      double rho2vh = 0.5*(rho2v[] + rho2v[-1]);
      Dmix2f.x[] = rho2vh*Dmix2vh*fsG.x[]*fm.x[];
#else
      Dmix2f.x[] = rho2*inDmix2[jj]*fsG.x[]*fm.x[];
#endif
    }

    foreach()
#ifdef VARPROP
      theta2[] = cm[]*max(rho2v[]*fG[], F_ERR);
#else
      theta2[] = cm[]*max(rho2*fG[], F_ERR);
#endif

    scalar YG = YGList[jj];
    scalar sgexp = sgexpList[jj];
    scalar sgimp = sgimpList[jj];

    //foreach() {
    //  sgexp[] = (f[] > F_ERR) ? sgexp[] : 0.;
    //  sgimp[] = (f[] > F_ERR) ? sgimp[] : 0.;
    //}

    diffusion (YG, dt, D=Dmix2f, r=sgexp, beta=sgimp, theta=theta2);
  }

#ifdef SOLVE_TEMPERATURE

  foreach_face() {
#ifdef VARPROP
    //double alpha1l = (rho1v[] != 0.) ? lambda1v[]/rho1v[]/cp1v[] : 0.;
    //double alpha2l = (rho2v[] != 0.) ? lambda2v[]/rho2v[]/cp2v[] : 0.;
    //double alpha1r = (rho1v[-1] != 0.) ? lambda1v[-1]/rho1v[-1]/cp1v[-1] : 0.;
    //double alpha2r = (rho2v[-1] != 0.) ? lambda2v[-1]/rho2v[-1]/cp2v[-1] : 0.;
    //lambda1f.x[] = 0.5*(alpha1r + alpha1l)*fsL.x[]*fm.x[];
    //lambda2f.x[] = 0.5*(alpha2r + alpha2l)*fsG.x[]*fm.x[];
    lambda1f.x[] = 0.5*(lambda1v[] + lambda1v[-1])*fsL.x[]*fm.x[];
    lambda2f.x[] = 0.5*(lambda2v[] + lambda2v[-1])*fsG.x[]*fm.x[];
#else
    lambda1f.x[] = lambda1/rho1/cp1*fsL.x[]*fm.x[];
    lambda2f.x[] = lambda2/rho2/cp2*fsG.x[]*fm.x[];
#endif
  }

  foreach() {
#ifdef VARPROP
    theta1[] = cm[]*max(fL[]*rho1v[]*cp1v[], F_ERR);
    theta2[] = cm[]*max(fG[]*rho2v[]*cp2v[], F_ERR);
#else
    theta1[] = cm[]*max(fL[], F_ERR);
    theta2[] = cm[]*max(fG[], F_ERR);
#endif
  }

  /**
  Solve diffusion equations for temperature. */

  double CACHE_TOLERANCE = TOLERANCE;
  TOLERANCE = TEMPERATURE_TOLERANCE;
  mgTL = diffusion (TL, dt, D=lambda1f, r=slT, theta=theta1);
  mgTG = diffusion (TG, dt, D=lambda2f, r=sgT, theta=theta2);
  TOLERANCE = CACHE_TOLERANCE;

# if PRINT_ITERS
  fprintf (stderr, "TL: iter = %d - resa = %g - nrelax = %d\n", mgTL.i, mgTL.resa, mgTL.nrelax);
  fprintf (stderr, "TG: iter = %d - resa = %g - nrelax = %d\n", mgTG.i, mgTG.resa, mgTG.nrelax);
# endif

#endif

  foreach() {
    /**
    We reconstruct the mass fractions from the total mass
    in every cell. */

    double totmassliq = 0.;
    for (scalar YL in YLList)
      totmassliq += YL[];
    for (scalar YL in YLList)
      YL[] = (totmassliq > 0.) ? YL[]/totmassliq : 0.;

    double totmassgas = 0.;
    for (scalar YG in YGList)
      totmassgas += YG[];
    for (scalar YG in YGList)
      YG[] = (totmassgas > 0.) ? YG[]/totmassgas : 0.;

    for (scalar YL in YLList)
      YL[] *= f[];
    for (scalar YG in YGList)
      YG[] *= (1. - f[]);
#ifdef SOLVE_TEMPERATURE
    TL[] *= f[];
    TG[] *= (1. - f[]);
#endif
  }

  /**
  We reconstruct the volume-averaged mass fractions and
  the one-field mass fraction field. */

  foreach() {
    for (scalar YL in YLList) {
      YL[] = clamp (YL[], 0., 1.);
      YL[] = (YL[] > F_ERR) ? YL[] : 0.;
    }
    for (scalar YG in YGList) {
      YG[] = clamp (YG[], 0., 1.);
      YG[] = (YG[] > F_ERR) ? YG[] : 0.;
    }
    for (int jj=0; jj<NGOS; jj++) {
      scalar Y  = YList[GOSI[jj]];
      scalar YG = YGList[GOSI[jj]];

      Y[] = YG[];
    }
    for (int jj=0; jj<NLS; jj++) {
      scalar Y  = YList[LSI[jj]];
      scalar YG = YGList[LSI[jj]];
      scalar YL = YLList[jj];

      Y[] = YL[] + YG[];
    }
  }

#ifdef SOLVE_TEMPERATURE
  foreach()
    T[] = TL[] + TG[];
#endif
}

