/**
# Selective Frequency Damping method in the stream-vorticity form

The low-pass filtered version of the initial flow is calculated at each iteration using an Euler explicite scheme ($\Delta$ a control parameter of the method) :
$$
    \mathbf{\bar{q}}^{n+1}_i \simeq
    \left( \mathbf{q}^n_i-\mathbf{\bar{q}}^n_i \right)
    \frac{\delta t}{\Delta}
    +\mathbf{\bar{q}}^n_i
$$

The flow is then forced by adding to the acceleration term ($\chi$ the second control parameter) : 
$$
\chi(\mathbf{q}^n_i-\mathbf{\bar{q}}^n_i)
$$

Because it has to be integrated at each time step :
$$
- \delta t\cdot\chi(\mathbf{q}^n_i-\mathbf{\bar{q}}^n_i)
$$

# Code */

scalar omegabar[];

/**
SFD's control parameters :*/
double SFD_Delta;
double SFD_chi;

/**
Variables that should be defined in the user's main .c code : */

double freq_SFD; // frequency of the instability to kill (approximative)
bool SFD_toggle; // conditions for activation of the SFD

/**
Calculation of the SFD's control parameters. */

event init (i = 0, last) {
 
  SFD_Delta = 1/(M_PI*freq_SFD);
  SFD_chi = 1/SFD_Delta;

/**
We set the initial conditions for the filtered flow, they should be the same as the calculated flow.*/
  foreach()
    omegabar[] = omega[];
}

/**
Calculation of the correction term and of the next step's filtered field. */

event SFD (i++) {
  foreach() {
    omega[] += SFD_toggle ? - dt * SFD_chi * (omega[] - omegabar[])/2. : 0.;
    omegabar[] = (omega[] - omegabar[]) * (dt/SFD_Delta) + omegabar[];
  }
}