/**
# Memory allocation for variable properties
Allocation of memory for species mass fractions, diffusion coefficients, temperatures, etc.
List are used to dynamically allocate the scalars according to the number of species in the kinetic scheme.
*/

#include "opensmoke.h"

unsigned int NGS, NSS;
scalar omega[];

scalar* YGList_G    = NULL; // gas species in the gas phase
scalar* YGList_S    = NULL; // gas species in the porous phase
scalar* YGList_Int  = NULL; // gas species at the interface
scalar* sSexpList   = NULL; // source terms for gas species in the porous phase
scalar* sGexpList   = NULL; // source terms for gas species in the gas phase
scalar* YSList      = NULL; // solid species in the porous phase
scalar* DmixGList_G = NULL; // diffusion coefficients for gas species in the gas phase
scalar* DmixGList_S = NULL; // diffusion coefficients for gas species in the porous phase
#ifdef MOLAR_DIFFUSION
scalar* XGList_G    = NULL; // mole fractions for gas species in the gas phase
scalar* XGList_S    = NULL; // mole fractions for gas species in the porous phase
scalar* XGList_Int  = NULL; // mole fractions for gas species at the interface
#endif
#ifdef MASS_DIFFUSION_ENTHALPY
scalar* cpGList_G = NULL;   // specific heats for gas species in the gas phase
scalar* cpGList_S = NULL;   // specific heats for gas species in the porous phase
#endif

double* gas_start; // starting mass fractions for gas species
double* sol_start; // starting mass fractions for solid species
double* gas_MWs; // molecular weights for gas species
double* sol_MWs; // molecular weights for solid species
double Pref = 101325.; // reference pressure [Pa]

#ifdef SOLVE_TEMPERATURE
scalar T[]; // temperature field
double lambdaS = 1.; double lambdaG = 1.; // thermal conductivities [W/m/K]
double cpS = 1.; double cpG = 1.; // specific heats [J/kg/K]

double TS0 = 300.; double TG0 = 300.; // initial temperatures [K]

scalar TInt[]; // interface temperature
scalar TS, TG; // solid and gas temperatures
scalar sST[], sGT[]; // source terms for solid and gas temperatures
face vector lambda1f[], lambda2f[]; // face vector thermal conductivities for porous and gas phases
vector lambda1v[], lambda2v[]; // thermal conductivities for porous and gas phases

/**
We want to consider anisotropic thermal conductivity for the solid phase.
Thus, we set the boundary conditions for the vector field lambda1v (solid phase) accordingly.
For simplicity, we do the same for the gas phase (lambda2v), although it is not strictly necessary.
These boundary conditions are set to Neumann (zero gradient) on all sides are necessary for the correct
calculation of the fluxes at the boundaries.
*/

lambda2v.n[bottom] = neumann(0.); 
lambda2v.t[bottom] = neumann(0.);
lambda2v.n[left] = neumann(0.);
lambda2v.t[left] = neumann(0.);
lambda2v.n[right] = neumann(0.);
lambda2v.t[right] = neumann(0.);
lambda2v.n[top] = neumann(0.);
lambda2v.t[top] = neumann(0.);

lambda1v.n[bottom] = neumann(0.);
lambda1v.t[bottom] = neumann(0.);
lambda1v.n[left] = neumann(0.);
lambda1v.t[left] = neumann(0.);
lambda1v.n[right] = neumann(0.);
lambda1v.t[right] = neumann(0.);
lambda1v.n[top] = neumann(0.);
lambda1v.t[top] = neumann(0.);
#endif

bool success;
scalar MWmixG_G[], MWmixG_S[]; // mixture molecular weights
scalar fG[], fS[]; // phase fractions
face vector fsS[], fsG[]; // face fractions

event defaults (i = 0) {

  /*
  We load the kinetic scheme using OpenSMOKE++ functions.
  */

  char kinfolder_root[120];
  sprintf (kinfolder_root, "%s/kinetics/%s/kinetics",
      getenv ("OPENSMOKE_INTERFACE"), kinfolder);

  OpenSMOKE_Init();
  OpenSMOKE_ReadKinetics (kinfolder_root);
  OpenSMOKE_ReadSolidKinetics (kinfolder_root);
  NGS = OpenSMOKE_NumberOfSpecies(); // Number of gas species
  NSS = OpenSMOKE_NumberOfSolidSpecies(); // Number of solid species
  
  YGList_G    = NULL;
  YGList_S    = NULL;
  YGList_Int  = NULL;
  sSexpList   = NULL;
  sGexpList   = NULL;
  YSList      = NULL;
  DmixGList_G = NULL;
  DmixGList_S = NULL;
  #ifdef MOLAR_DIFFUSION
  XGList_G    = NULL;
  XGList_S    = NULL;
  XGList_Int  = NULL;
  #endif
  #ifdef MASS_DIFFUSION_ENTHALPY
  cpGList_G = NULL;
  cpGList_S = NULL;
  #endif

  //Allocate gas species fields
  for (int jj = 0; jj<NGS; jj++) {
    scalar a = new scalar;
    free (a.name);
    char name[20];
    sprintf (name, "%s_S",OpenSMOKE_NamesOfSpecies(jj));
    a.name = strdup (name);
    YGList_S = list_append (YGList_S, a);
  }
  reset (YGList_S, 0.);

  for (int jj = 0; jj<NGS; jj++) {
    scalar a = new scalar;
    free (a.name);
    char name[20];
    sprintf (name, "%s_G",OpenSMOKE_NamesOfSpecies(jj));
    a.name = strdup (name);
    YGList_G = list_append (YGList_G, a);
  } 
  reset (YGList_G, 0.);

  for (int jj = 0; jj<NGS; jj++) {
    scalar a = new scalar;
    free (a.name);
    char name[20];
    sprintf (name, "%s_Int",OpenSMOKE_NamesOfSpecies(jj));
    a.name = strdup (name);
    YGList_Int = list_append (YGList_Int, a);
  }
  reset (YGList_Int, 0.);

  //Allocate solid species fields
  for (int jj = 0; jj<NSS; jj++) {
    scalar a = new scalar;
    free (a.name);
    char name[20];
    sprintf (name, "%s",OpenSMOKE_NamesOfSolidSpecies(jj));
    a.name = strdup (name);
    YSList = list_append (YSList, a);
  }
  reset (YSList, 0.);

  //Allocate diff coeff fields
  for (int jj = 0; jj<NGS; jj++) {
    scalar s = new scalar;
    free (s.name);
    char name[20];
    sprintf (name, "D_%s_S",OpenSMOKE_NamesOfSpecies(jj));
    s.name = strdup (name);
    DmixGList_S = list_append (DmixGList_S, s);
  }
  reset (DmixGList_S, 0.);

  for (int jj = 0; jj<NGS; jj++) {
    scalar s = new scalar;
    free (s.name);
    char name[20];
    sprintf (name, "D_%s_G",OpenSMOKE_NamesOfSpecies(jj));
    s.name = strdup (name);
    DmixGList_G = list_append (DmixGList_G, s);
  }
  reset (DmixGList_G, 0.);

// fields for the source therms
for (int jj=0; jj<NGS; jj++) {
    scalar a = new scalar;
    scalar b = new scalar;
    free (a.name);
    free (b.name);
    char aname[20];
    char bname[20];
    sprintf (aname, "sSexp_%s", OpenSMOKE_NamesOfSpecies(jj));
    sprintf (bname, "sGexp_%s", OpenSMOKE_NamesOfSpecies(jj));
    a.name = strdup (aname);
    b.name = strdup (bname);
    a.nodump = true;
    b.nodump = true;
    sSexpList = list_append (sSexpList, a);
    sGexpList = list_append (sGexpList, b);
  }
  reset (sSexpList, 0.);
  reset (sGexpList, 0.);

  #ifdef MOLAR_DIFFUSION
  //Allocate gas molar fraction fields
  for (int jj = 0; jj<NGS; jj++) {
    scalar a = new scalar;
    free (a.name);
    char name[20];
    sprintf (name, "XG_%s_S",OpenSMOKE_NamesOfSpecies(jj));
    a.name = strdup (name);
    XGList_S = list_append (XGList_S, a);
  }
  reset (XGList_S, 0.);

  for (int jj = 0; jj<NGS; jj++) {
    scalar a = new scalar;
    free (a.name);
    char name[20];
    sprintf (name, "XG_%s_G",OpenSMOKE_NamesOfSpecies(jj));
    a.name = strdup (name);
    XGList_G = list_append (XGList_G, a);
  }
  reset (XGList_G, 0.);

  for (int jj = 0; jj<NGS; jj++) {
    scalar a = new scalar;
    free (a.name);
    char name[20];
    sprintf (name, "XG_%s_Int",OpenSMOKE_NamesOfSpecies(jj));
    a.name = strdup (name);
    XGList_Int = list_append (XGList_Int, a);
  }
  reset (XGList_Int, 0.);
  #endif
  #ifdef MASS_DIFFUSION_ENTHALPY
  for (int jj = 0; jj<NGS; jj++) {
    scalar a = new scalar;
    free (a.name);
    char name[20];
    sprintf (name, "cpG_%s_G",OpenSMOKE_NamesOfSpecies(jj));
    a.name = strdup (name);
    cpGList_G = list_append (cpGList_G, a);
  }
  reset (cpGList_G, 0.);
  for (int jj = 0; jj<NGS; jj++) {
    scalar a = new scalar;
    free (a.name);
    char name[20];
    sprintf (name, "cpG_%s_S",OpenSMOKE_NamesOfSpecies(jj));
    a.name = strdup (name);
    cpGList_S = list_append (cpGList_S, a);
  }
  reset (cpGList_S, 0.);
  #endif

  //initialize vector with initial values
  gas_start = (double *)malloc(NGS * sizeof(double));
  sol_start = (double *)malloc(NSS * sizeof(double));
  gas_MWs = (double *)malloc(NGS * sizeof(double));
  sol_MWs = (double *)malloc(NSS * sizeof(double));

  for (int jj=0; jj<NGS; jj++) {
    gas_start[jj] = 0.;
    gas_MWs[jj] = OpenSMOKE_MW(jj);
  }

  for (int jj=0; jj<NSS; jj++) {
    sol_start[jj] = 0.;
    sol_MWs[jj] = OpenSMOKE_MW_Solid(jj);
  }

  // Set correct side for inverse fields
  for (scalar s in YGList_S)
    s.inverse = false;

  for (scalar s in YGList_G)
    s.inverse = true;

  for (scalar s in YSList)
    s.inverse = false;

#ifdef MOLAR_DIFFUSION
  for (scalar s in XGList_S)
    s.inverse = false;

  for (scalar s in XGList_G)
    s.inverse = true;
#endif

  // Add species to f tracer list
  scalar *temp;
  temp = list_concat(f.tracers, YGList_S);
  if (f.tracers)
    free(f.tracers);
  f.tracers = temp;

  temp = list_concat(f.tracers, YGList_G);
  if (f.tracers)
    free(f.tracers);
  f.tracers = temp;

  temp = list_concat(f.tracers, YSList);
  if (f.tracers)
    free(f.tracers);
  f.tracers = temp;

  fS.nodump = true;
  fG.nodump =true;

#ifdef SOLVE_TEMPERATURE
  TS = new scalar;
  TG = new scalar;

  TS.inverse = false;
  TG.inverse = true;

  sST.nodump = true;
  sGT.nodump = true;

  f.tracers = list_append (f.tracers, TS);
  f.tracers = list_append (f.tracers, TG);
#endif

#ifdef VARPROP
  if (is_constant(rhoGv_G)) {
    scalar * l = list_copy(all);
    rhoGv_G = new scalar;
    free(all);
    all = list_concat({rhoGv_G}, l);
    free(l);
  }
  if (is_constant(rhoGv_S)) {
    scalar * l = list_copy(all);
    rhoGv_S = new scalar;
    free(all);
    all = list_concat({rhoGv_S}, l);
    free(l);
  }
  if (is_constant(rhoSv)) {
    scalar * l = list_copy(all);
    rhoSv = new scalar;
    free(all);
    all = list_concat({rhoSv}, l);
    free(l);
  }
  if (is_constant(muGv_G)) {
    scalar * l = list_copy(all);
    muGv_G = new scalar;
    free(all);
    all = list_concat({muGv_G}, l);
    free(l);
  }
  if (is_constant(muGv_S)) {
    scalar * l = list_copy(all);
    muGv_S = new scalar;
    free(all);
    all = list_concat({muGv_S}, l);
    free(l);
  }
  if (is_constant(lambdaGv_G)) {
    scalar * l = list_copy(all);
    lambdaGv_G = new scalar;
    free(all);
    all = list_concat({lambdaGv_G}, l);
    free(l);
  }
  if (is_constant(lambdaGv_S)) {
    scalar * l = list_copy(all);
    lambdaGv_S = new scalar;
    free(all);
    all = list_concat({lambdaGv_S}, l);
    free(l);
  }
  if (is_constant(lambdaSv)) {
    scalar * l = list_copy(all);
    lambdaSv = new scalar;
    free(all);
    all = list_concat({lambdaSv}, l);
    free(l);
  }
  if (is_constant(cpGv_G)) {
    scalar * l = list_copy(all);
    cpGv_G = new scalar;
    free(all);
    all = list_concat({cpGv_G}, l);
    free(l);
  }
  if (is_constant(cpGv_S)) {
    scalar * l = list_copy(all);
    cpGv_S = new scalar;
    free(all);
    all = list_concat({cpGv_S}, l);
    free(l);
  }
  if (is_constant(cpSv)) {
    scalar * l = list_copy(all);
    cpSv = new scalar;
    free(all);
    all = list_concat({cpSv}, l);
    free(l);
  }
  #endif
}

event init (i = 0) {
  //initialize gas fractions fields
  // check if the sum of the fractions is 1
  double sum = 0.;
  for (int jj=0; jj<NGS; jj++)
    sum += gas_start[jj];
  if (fabs(sum-1.) > 1e-10) {
    fprintf(stderr, "Sum of gas fractions is not 1. Exiting...\n");
    exit(1);
  }

  foreach() {
    for (int jj=0; jj<NGS; jj++) {
      scalar YG = YGList_G[jj];
      YG[] = gas_start[jj]*(1-f[]);
    }

    for (int jj=0; jj<NGS; jj++) {
      scalar YG = YGList_S[jj];
      YG[] = gas_start[jj]*f[];
    }
  }

  foreach() {
    for (int jj=0; jj<NGS; jj++) {
      scalar YGInt = YGList_Int[jj];
      YGInt[] = gas_start[jj];
    }
  }

  //initialize solid fractions fields
  // check if the sum of the fractions is 1
  sum = 0.;
  for (int jj=0; jj<NSS; jj++)
    sum += sol_start[jj];
  if (fabs(sum-1.) > 1e-6) {
    fprintf(stderr, "Sum of solid fractions is not 1. Exiting...\n");
    exit(1);
  }

  foreach() {
    for (int jj=0; jj<NSS; jj++) {
      scalar YS = YSList[jj];
      YS[] = sol_start[jj]*f[];
    }
  }

  foreach() {
    for (int jj=0; jj<NGS; jj++) {
      scalar sSexp = sSexpList[jj];
      scalar sGexp = sGexpList[jj];
      sSexp[] = 0.;
      sGexp[] = 0.;
    }
  }

#ifdef SOLVE_TEMPERATURE
  foreach() {
    TS[] = TS0*f[];
    TG[] = TG0*(1. - f[]);
    T[]  = TS[] + TG[];
    TInt[] = f[]<1-F_ERR && f[] > F_ERR ? TS0 : 0.;
  }
#endif
}

event cleanup (t = end) {
  delete (YSList),      free (YSList),      YSList = NULL;
  delete (YGList_S),    free (YGList_S),    YGList_S = NULL;
  delete (YGList_G),    free (YGList_G),    YGList_G = NULL;
  delete (YGList_Int),  free (YGList_Int),  YGList_Int = NULL;
  delete (sSexpList),   free (sSexpList),   sSexpList = NULL;
  delete (sGexpList),   free (sGexpList),   sGexpList = NULL;
  delete (DmixGList_S), free (DmixGList_S), DmixGList_S = NULL;
  delete (DmixGList_G), free (DmixGList_G), DmixGList_G = NULL;
  free(gas_start),  gas_start = NULL;
  free(sol_start),  sol_start = NULL;
  free(gas_MWs),    gas_MWs = NULL;
  free(sol_MWs),    sol_MWs = NULL;
  delete(f.tracers),  free(f.tracers),   f.tracers = NULL;

#ifdef MOLAR_DIFFUSION
  delete (XGList_S), free (XGList_S), XGList_S = NULL;
  delete (XGList_G), free (XGList_G), XGList_G = NULL;
  delete (XGList_Int), free (XGList_Int), XGList_Int = NULL;
#endif

#ifdef MASS_DIFFUSION_ENTHALPY
  delete (cpGList_G), free (cpGList_G), cpGList_G = NULL;
  delete (cpGList_S), free (cpGList_S), cpGList_S = NULL;
#endif

#ifdef SOLVE_TEMPERATURE
  delete ({TS,TG});
#endif

#ifdef TEMPERATURE_PROFILE
  TemperatureProfile_Free();
#endif
}