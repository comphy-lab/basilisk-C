/**
# Selective Frequency Damping method

The low-pass filtered version of the initial flow is calculated at each iteration using an Euler explicite scheme ($\Delta$ a control parameter of the method) :
$$
    \mathbf{\bar{q}}^{n+1}_i \simeq
    \left( \mathbf{q}^n_i-\mathbf{\bar{q}}^n_i \right)
    \frac{\delta t}{\Delta}
    +\mathbf{\bar{q}}^n_i
$$

The flow is then forced by adding to the acceleration term ($\chi$ the second control parameter) : 
$$
-\chi(\mathbf{q}^n_i-\mathbf{\bar{q}}^n_i)
$$

# Code

Variables that should be defined in the user's main .c code : */

double freq_SFD; // frequency of the instability to kill (approximative)
bool SFD_toggle; // conditions for activation of the SFD

/**
SFD's control parameters*/

double SFD_Delta;
double SFD_chi;

vector ubar[];
scalar pbar[];

/**
av[] the correction term is imposed as the acceleration. */

face vector av[];
event defaults (i = 0) {
  a = av;
}

/**
We set the initial conditions for the filtered flow, they should be the same as the calculated flow.*/

event init (i = 0, last) {
  
  SFD_Delta = 1/(M_PI*freq_SFD);
  SFD_chi = 1/SFD_Delta;
  
  foreach() {
    foreach_dimension()
      ubar.x[] = u.x[];
    pbar[] = p[];
  }
}

/**
We add the SFD correction term to the acceleration. */

event acceleration (i++) {
  foreach_face()
    av.x[] = SFD_toggle ? - SFD_chi * ((u.x[] - ubar.x[]) + (u.x[-1] - ubar.x[-1]))/2. : 0.;
}

/**
We calculate the next step’s filtered field. */

event bar (i++) {
  foreach() {
    foreach_dimension()
      ubar.x[] = (u.x[] - ubar.x[]) * (dt/SFD_Delta) + ubar.x[];
    pbar[] = (p[] - pbar[]) * (dt/SFD_Delta) + pbar[];
  }
}
