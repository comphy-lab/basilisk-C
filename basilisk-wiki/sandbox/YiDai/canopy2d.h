/**
## Canopy model
We implement canopy model based on Patton et al 2016; Boekee et al 2023. This parameterization should incorporate with large eddy simulation. This header file only solve energy balance for the leaf.
To solve the flow part, the convective heat exchange (source term for temp) and wind drag should be including in other files (!!!)
*/
// #include "runge-kutta.h"                 // for now we use forward euler, later on maybe we will just use time integgrator from basilisk
// define where is the canopy
#define Zh 20. // canopy top height [m]   mind the float
// #define CANOPY (y <= Zh)? 1:0              // create a mask for canopy layer
#define PAD(s) ((s <= Zh) ? 2 * sin(s / Zh * M_PI) : 0)
// wind drag
#define Cd 0.15 // dimensionless drag coefficient (Shaw and Schumann 1992)

// default parameters
#define boltz 5.67E-8
#define vis 1.37E-5 // kinematic viscosity of the air [m^2 s^-1]
#define Gconst 9.81
#define T_ref 273.
#define Kd 0.024 // thermal conductivity of air [W m-1 K-1]
// other parameters
#define L_l 4E-2 // representative length scale [m]
#define VF_s 0.1
#define VF_g VF_s
#define VF_l (1. - VF_s) // view factor of sky, ground, leaf
#define eps_s 0.8
#define eps_g 0.98
#define eps_l 0.96 // emissivity of sky, ground, leaf
#define T_s 260.
#define T_g 280.   // temperature of sky and surface in [kelvin]
#define Cp_l 2.0E6 // leaf heat capacity [J m^-3 K^-1]
#define Cp_a 1005. // dry air heat capacity [J kg^-1 K^-1]
#define rho_a 1.   // air density at near 10 C [kg m^-3]
// geometry of leaf
// #define R_l 0.00657                              // radius of leaf surface [m]
#define d_l 2.0E-4                // thickness of the leaf  [m]
#define A_l (2. * M_PI * sq(L_l)) // surface area leaf      [m^2]
#define V_l (A_l / 2. * d_l)      // volume leaf            [m^3]
// leaf temperature
scalar TV[]; // leaf temp [K]
scalar H[];  // convective heat exchange [W m^-2]

event leaf_flow(i++)
{
    // refine the canopy layer
    refine(y <= Zh * 1.5 && (level <= maxlevel));
    // radiation - Stefan–Boltzmann law
    // double Lwin, Lwout;
    scalar Lwnet[];
    foreach ()
    {
        if (y <= Zh)
        {
            double Lwin = 0.5 * VF_s * eps_s * boltz * pow(T_s, 4) +
                          0.5 * VF_g * eps_g * boltz * pow(T_g, 4) +
                          1. * VF_l * eps_l * boltz * pow(TV[], 4);
            double Lwout = eps_l * boltz * pow(TV[], 4);
            Lwnet[] = Lwin - Lwout;
        }
    }
    // Convective -- resistance of turbulent heat exchange -- Boekee et al 2023
    foreach ()
    {
        if (y <= Zh)
        {
            double T_a = b[] * T_ref / Gconst + T_ref;
            double gstar = Gconst * (TV[] - T_a) / T_a;

            // watch out for the sign here for gstar
            double M = sqrt(sq(u.x[]) + sq(u.y[]) + sq(u.z[]) + fabs(2. * L_l * gstar));
            double Re = fabs(M) * L_l / vis;
            double Nu = (Re > 2E4) ? 0.032 * pow(Re, 0.8) : 0.6 * pow(Re, 0.5);
            double rH = L_l / Nu / Kd * Cp_a * rho_a;

            H[] = Cp_a * rho_a / rH * (TV[] - T_a);
            // diagnose leaf temp with forward euler
            // for vegetation is negative, in the source term of b, it should be positive
            TV[] += dt * (Lwnet[] - H[]) * A_l / (Cp_l * V_l);
        }
    }
}

/**
## Test
* [Inversion layer with caonpy case in 2D](canopy_test.c)
* There should be more tests

##Reference
Patton, E.G., Sullivan, P.P., Shaw, R.H., Finnigan, J.J. and Weil, J.C., 2016. Atmospheric stability influences on coupled boundary layer and canopy turbulence. Journal of the Atmospheric Sciences, 73(4), pp.1621-1647.

Boekee, J., Dai, Y., Schilperoort, B., van de Wiel, B.J. and ten Veldhuis, M.C., 2023. Plant–atmosphere heat exchange during wind machine operation for frost protection. Agricultural and Forest Meteorology, 330, p.109312.
*/
