/**
# Batch Reactors

Collection of ideal reactor functions for the operator splitting
solution of reactive systems ([Cuoci et al. 2013](#cuoci2013numerical)).
The temperature and chemical species concentrations are considered
perfectly mixed inside each cell. Neglecting the convective and diffusive
transport in the conservation equations, the equations reduce to the ODE
system of a batch reactor:
$$
\begin{cases}
  \dfrac{d\omega_i}{dt} = \dfrac{1}{\rho} \sum \limits_{j=1}^{NR} \nu_{i,j} r_j MW_i \\
  \dfrac{d T}{dt} = \dfrac{1}{\rho c_p} \sum \limits_{j=1}^{NR} r_j \Delta H_{R_j}
\end{cases}
$$
where $\nu$ is the matrix of the stoichiometric coefficients, while
$r_j$ is the rate of the $j$-th reaction ([Fogler 2020](#fogler2020elements)).
The term $\Delta H_R$ is the enthalpy change associated with the
chemical reactions. The reaction rate $r_j$ is usually computed as:
$$
  r_j = k_j \prod \limits_{i=1}^{NS} c_i^{\nu_{i,j}}
$$
and the kinetic constant $k_j$ is computed from a modified Arrhenius equation:
$$
  k_j = AT^n \exp^{-E_a/(RT)}
$$
This procedure is managed by the OpenSMOKE++ library ([opensmoke.h](opensmoke.h)).
*/

#include "radiation.h"

/**
## *batch_isothermal_constantpressure()*: solve chemical species reactions

* *y*: vector (length = NS) containing the initial values of the system
* *dt*: simulation time step
* *dy*: right-hand-side of the batch system of equations
* *args*: structure with additional arguments to be passed to the system
*/

void batch_isothermal_constantpressure (const double * y, const double dt, double * dy, void * args) {

  /**
  Unpack data for the ODE system. */

  UserDataODE * data = (UserDataODE *)args;
  int NS = data->NS;
  double rho = data->rho;
  double Pressure = data->P;
  double Temperature = data->T;
  double * sources = data->sources;

  /**
  Set temperature and pressure of the system,
  for the calculation of the kinetic constants. */

  OpenSMOKE_GasProp_SetTemperature (Temperature);
  OpenSMOKE_GasProp_SetPressure (Pressure);

  /**
  We create a vector with the mass fractions and we
  convert it to mole fractions. */

  double massfracs[NS], molefracs[NS];
  for (int jj=0; jj<NS; jj++) {
    massfracs[jj] = y[jj];
  }

  double MWMix = 0.;
  OpenSMOKE_MoleFractions_From_MassFractions (molefracs, &MWMix, massfracs);

  double ctot = Pressure/(R_GAS*1000.)/Temperature;
  double ci[NS], ri[NS];
  for (int jj=0; jj<NS; jj++) {
    ci[jj] = ctot*molefracs[jj];
    ri[jj] = 0.;
  }
  rho = ctot*MWMix;

  /**
  Compute the kinetic constant, reaction rates,
  and formation rates. */

  OpenSMOKE_GasProp_KineticConstants ();
  OpenSMOKE_GasProp_ReactionRates (ci);
  OpenSMOKE_GasProp_FormationRates (ri);

  /**
  Equation for the chemical species. */

  for (int jj=0; jj<OpenSMOKE_NumberOfSpecies(); jj++) {
    dy[jj] = OpenSMOKE_MW(jj)*ri[jj]/rho;
    sources[jj] = dy[jj]*rho;
  }
}

/**
## *batch_nonisothermal_constantpressure()*: solve chemical species reactions

* *y*: vector (length = NS+1) containing the initial values of the system
* *dt*: simulation time step
* *dy*: right-hand-side of the batch system of equations
* *args*: structure with additional arguments to be passed to the system
*/

void batch_nonisothermal_constantpressure (const double * y, const double dt, double * dy, void * args) {

  unsigned int NS = OpenSMOKE_NumberOfSpecies();

  /**
  Unpack data for the ODE system. */

  UserDataODE * data = (UserDataODE *)args;
  double rho = data->rho;
  double cp = data->cp;
  double Pressure = data->P;
  double Temperature = y[OpenSMOKE_NumberOfSpecies()];
  double * sources = data->sources;

  /**
  Set temperature and pressure of the system,
  for the calculation of the kinetic constants. */

  OpenSMOKE_GasProp_SetTemperature (Temperature);
  OpenSMOKE_GasProp_SetPressure (Pressure);

  /**
  We create a vector with the mass fractions and we
  convert it to mole fractions. */

  double massfracs[NS], molefracs[NS], MWs[NS];
  for (int i = 0; i < NS; i++) {
    massfracs[i] = y[i];
    MWs[i] = OpenSMOKE_MW (i);
  }
  correctfrac (massfracs, NS);
  double MWMix = mass2mw (massfracs, MWs, NS);
  mass2molefrac (molefracs, massfracs, MWs, NS);

  //double MWMix = 0.;
  //OpenSMOKE_MoleFractions_From_MassFractions (molefracs, &MWMix, massfracs);

//#ifdef VARPROP

  /**
  We compute the density and the heat capacity as
  a function of the thermodynamic state. */

  double ctot = Pressure/(R_GAS*1000.)/Temperature;  // [kmol/m3]
  rho = ctot*MWMix;  // [kg/m3]

  cp = OpenSMOKE_GasProp_HeatCapacity (molefracs);

  double ci[NS], ri[NS];
  for (int i = 0; i < NS; i++)
    ci[i] = ctot*molefracs[i];

  data->rho = rho;
  data->cp = cp;

//#endif

  //double ci[NS], ri[NS];
  //for (int jj=0; jj<NS; jj++) {
  //  //ci[jj] = massfracs[jj]*rho/OpenSMOKE_MW(jj);
  //  //ci[jj] = ctot*molefracs[jj];
  //  //ci[jj] = (ci[jj] < 0.) ? 0. : ci[jj];
  //  ri[jj] = 0.;
  //}

  //OpticallyThinProperties otp;
  //otp.T = Temperature;
  //otp.P = Pressure;
  //otp.xH2O = molefracs[otm.indexH2O];
  //otp.xCO2 = molefracs[otm.indexCO2];
  //otp.xCO  = molefracs[otm.indexCO];
  //otp.xCH4 = molefracs[otm.indexCH4];

  ThermoState ts;
  ts.T = Temperature, ts.P = Pressure, ts.x = molefracs;
  double divq_rad = opensmoke_optically_thin (&ts);

  /**
  Compute the kinetic constant, reaction rates,
  and formation rates. */

  OpenSMOKE_GasProp_KineticConstants ();
  OpenSMOKE_GasProp_ReactionRates (ci);     // [kmol/m3]
  OpenSMOKE_GasProp_FormationRates (ri);    // [kmol/m3/s]

  /**
  Equation for the chemical species. */

  for (int jj=0; jj<OpenSMOKE_NumberOfSpecies(); jj++) {
    dy[jj] = OpenSMOKE_MW(jj)*ri[jj]/rho;
    sources[jj] = dy[jj]*rho;
  }

  /**
  Get the heat of reaction and compute the equation for the
  temperature. We add non-linear contributions such as the
  heat dissipated for radiation */

  double QR = OpenSMOKE_GasProp_HeatRelease (ri);
  //dy[OpenSMOKE_NumberOfSpecies()] = (QR + divq_rad (&otp))/rho/cp;
  dy[OpenSMOKE_NumberOfSpecies()] = (QR + divq_rad)/rho/cp;
  sources[OpenSMOKE_NumberOfSpecies()] = dy[OpenSMOKE_NumberOfSpecies()]*rho*cp;
}

/**
# References

~~~bib
@article{cuoci2013numerical,
  title={Numerical modeling of laminar flames with detailed kinetics based on the operator-splitting method},
  author={Cuoci, Alberto and Frassoldati, Alessio and Faravelli, Tiziano and Ranzi, Eliseo},
  journal={Energy \& fuels},
  volume={27},
  number={12},
  pages={7730--7753},
  year={2013},
  publisher={ACS Publications}
}

@book{fogler2020elements,
  title={Elements of chemical reaction engineering},
  author={Fogler, H Scott},
  year={2020},
  publisher={Pearson Boston}
}
~~~
*/

