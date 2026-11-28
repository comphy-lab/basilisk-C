/**

# Convection with the multilayer

This code aims to represent a oceanic convection case. The water is stratified
with $N^2=\sqrt{g\frac{\partial \rho}{\partial z}}$. 

The equation of state is linear, following $\rho = \rho_0(1-\alpha(T-T_0))$ or
in Basilisk's dr.h notation $\Delta \rho = -\alpha(T-T_0)$.

The experiment depends on the Richardson number $Ri=N^2.H^2/w_*^2$ and Prandtl
number which is the ratio of viscosity to temperature diffusion $Pr=\nu/diff_T$.
We use the convective velocity $w_*=(-B0.h)^{1/3}$, $B_0$ the buoyancy flux at the
surface $B_0=\frac{g\alpha}{\rho.cp}Q$ with $Q$ the heat flux in $W.m^{-2}$, $H$
the depth of the mixed layer. An analytical expression for the mixed layer height
evolution [Van Roeckel, Souza] can be derived:

$$
\begin{aligned}
  h(t) = \sqrt{a \frac{-B_0}{N^2}t }
\end{aligned}
$$

where $a$ is a prefactor. [Souza] found that their prefactor was $a=3.0$ when
they take into account entrainment of cold water into the mixed layer, while
[Roeckel] has $a=2.8$ without entrainment.

## Warning

Please note that this code uses the nh-staggered.h code from
sandbox/bderembl/libs/nh-staggered.h. Copy it to src/layered to be able to use
it easily with ml_convection.c

*/

#include "grid/multigrid1D.h"
#include "layered/hydro.h"
//#include "layered/nh.h"
#include "layered/nh-staggered.h"
/**
We use remapping, but this might not be necessary
*/

#include "layered/remap.h"
#include "layered/perfs.h"
#include "bderembl/libs/extra.h"      // parameters from namlist
#include "bderembl/libs/netcdf_bas.h" // read/write netcdf files

#include "hugoj/lib/diffusionH.h" // Neuman bott and top
//#define g_ 9.81


/*
DEFAULT PARAMETERS

These parameters are changed by the values in the namelist

Dimensions : [Length, Time, Temperature, Energy, mass]
*/
char namlist[80] = "namelist.toml";    // file name of namlist
char file_out[20] = "out.nc";          // file name of output
// -> Initial conditions
double strat = 0.000002 [0,-2];       // [s-2] N^2 stratification
double Ts = 20. [0,0,1];              // [K] Surface temperature (arbitrary)
// -> Forcing
double qt = 100. [-2,-1,0,1];         // [W.m-2] Heat flux
// -> Domain definition
int N_grid = 64 [];                       // 2^N_grid : number of x and y gridpoints
double L = 200.0 [1];                 // domain size
int N_layer = 5 [];                      // number of layers
double h0 = 1.0 [1];                  // depth of water
// -> Runtime parameters
int restart = 0 [];                  // 1: restart, 0: no restart
double tend = 2.0 [0,1];              // end time of simulation
// -> saving outputs
double dtout = 2.0 [0,1];             // dt for output in netcdf
double smalltime = 1e-10 [0,1];       // (s) small time increment
int pad = 4 [0];                      // number of 0-padding for ouput files
int nout = 1 [0];                     // number of the outfile
char fileout[100];                    // name of outfile
// -> physical properties
double nu0 = 0.00025 [2,-1];          // Viscosity for vertical diffusion
double thetaH = 0.5 [0];              // theta_h for dumping fast barotropic modes
// -> stratification related
double rho0 = 1025. [-3,0,0,0,1];     // [kg.m-3] reference density
double cp = 4.2e3 [0,0,-1,1,-1];      // [J.kg-1.K-1] heat capacity water
double alphaT = 2e-4 [0,0,-1];        // [K-1] Thermal expansion coeff for water
double Pr = 1. [0,0,0,0,0];                       // Prandtl number = nu/diff_T
double diff_T = 1.5e-5 [2,-1];         // [m2.s-1] Scalar vertical diffusion coeff
double T0 = 20. [0,0,1];              // [°C] Reference temperature
double Trand = 0.0001 [0,0,1];          // [°C] Random temperature perturbution
const double g_ = 9.81 [1,-2];        // [m.s-2] Gravity
double mld_i = 10.;                   // initial MLD
#define drho(T) (alphaT*(T0-T))       // Linear equation of state (Vallis 2.4)

#include "layered/dr.h"

// diag
double *T_profile;
double *W_profile;
double *Tflx_profile;
double dt_mean = 1.;

/** Initial profile of temperature */
double T_ini(double z, double mld_i, double strat, double g_, double alphaT, double Ts){
  if (z>=-mld_i){
    return Ts;
  }
  else{
    return strat/(g_*alphaT)*(z + mld_i) + Ts;
  }
}

/** main program*/
int main(int argc, char *argv[])  
{
 
  // Building a 'params' array with all parameters from the namlist
  params = array_new();
  add_param("N_grid", &N_grid, "int");
  add_param("L", &L, "double");
  add_param("N_layer", &N_layer, "int");
  add_param("h0", &h0, "double");
  add_param("tend", &tend, "double");
  add_param("nu0", &nu0, "double");
  add_param("thetaH", &thetaH, "double");
  add_param("restart", &restart, "int");
  add_param("dtout", &dtout, "double");
  add_param("strat", &strat, "double");
  add_param("Ts", &Ts, "double");
  add_param("rho0", &rho0, "double");
  add_param("cp", &cp, "double");
  add_param("alphaT", &alphaT, "double");
  add_param("qt", &qt, "double");
  add_param("dt_mean", &dt_mean, "double");
  add_param("Pr", &Pr, "double");
  add_param("mld_i", &mld_i, "double");

  // Search for the configuration file with a given path or read params.in
  if (argc == 2)
    strcpy(file_param, argv[1]);
  else
    strcpy(file_param, namlist);
  read_params(file_param);
  
  /** Settings solver values from namlist values */
  L0 = L;
  nu = nu0;
  diff_T = nu0/Pr;
  /** If the horizontal grid is too coarse, a wave develops at the free surface.
      In 1D and for N>=128 the surface stays static */
  N = N_grid; 
  nl = N_layer;
  G = g_;
  theta_H = thetaH;
  CFL_H = 1.;
  CFL = 0.5;
  TOLERANCE = 1e-2;

  origin (-L0/2., -L0/2.);
  
  //periodic (right);
  #if dimension==2
    periodic (top);
  #endif
  periodic (left);

  // T profile
  T_profile = (double*)calloc(nl, sizeof(double));
  W_profile = (double*)calloc(nl, sizeof(double));
  Tflx_profile = (double*)calloc(nl, sizeof(double));

  fprintf (stderr, "Read in parameters!\n");
  fprintf (stderr, "CFL_H=%g CFL=%g\n", CFL_H,CFL);

  run();
  
}

event init(i =  0) {
  // reset({T},0.);
  geometric_beta (1., true); // if !=0, varying layer thickness 1./3.
  // step 1: set eta and h
  foreach(cpu) {
    zb[] = -h0;
    eta[] = 0.;
    double H = - zb[];
    foreach_layer() 
      h[] = H*beta[point.l];
  }

  reset({T},0.);
  // step 2: set currents
  foreach(cpu) {
    double z = zb[];
    //dimensional (T[]==Trand);
    foreach_layer() {
      z += h[]/2.;
      foreach_dimension()
        u.x[] = 0.;
      w[] = 0.;

/**
Initial temperature profile is linear. We add a white noise in the upper layers
to trigger the instability.
*/
      T[] = T_ini(z,mld_i,strat,g_,alphaT,Ts) + Trand * noise() * exp(-z/(10.0*h0));
      z += h[]/2.;
    }
  }
  
  // initializing diag arrays
  // T_profile[0] = Trand; // <- this passes the dimensional analysis
  // dimensional (T_profile[0] == Trand); // <- this doesnt
   for (int i=0; i<nl; ++i) {
     T_profile[i] = Trand*0.;
     Tflx_profile[i] = Trand*L0/DT*0.;
     W_profile[i] = L0/DT*0.;
   }
  
  show_dimension (T);
  
  // original version, bug for now (currents writes = 0)
  create_nc({zb, h, u, w, eta, T}, file_out);

  fprintf (stderr,"Done initialization!\n");
}

/** Vertical diffusion on T
This implements a heat flux at the surface */
event viscous_term (i++)
{
  foreach()
    vertical_diffusion2 (point, h, T, dt, diff_T, qt/(diff_T*rho0*cp), 0.); // strat/(g_*alphaT)
  horizontal_diffusion ({T}, diff_T, dt);
  horizontal_diffusion ({u, w}, nu, dt);
}

/** Writing a netcdf file*/
event output(t = 0.; t<= tend+smalltime; t+=dtout){
  write_nc();
}

event end(t=tend+smalltime){
  dump();
}


double* h_avg(scalar var, double* profile){
  /*
  This function computes the horizontal average of var.
  
  INPUTS:
    var: scalar (C Basilisk), the variable to average.
    profile: array of double of length nl
  OUTPUTS:
    an array of double with the average values inside
   */

  int NB;
  #if dimension==1
    NB=N;
  #else
    NB=N*N;
  #endif

  foreach_layer() {
    double sum = 0.;
    foreach (reduction(+:sum))
      sum += var[];
    profile[_layer] += sum / NB;
  }

  return profile;
}

int write_profile(char* name, double* profile){
  // main worker is writing the file
  if (pid()==0) {
    for (int i=0; i<nl; ++i) {
      fprintf (stdout, "%f %d %g\n", t, i, profile[i]);
    }
    fprintf(stdout,"\n");
  }
  return 0;
}


// This event compute layer average of T, w
event compute_horizontal_avg (t+=dt_mean; t<=tend+smalltime){
  T_profile = h_avg(T, T_profile);
  // show_dimension (T_profile[0]);
  // show_dimension (W_profile[0]);
  //W_profile = h_avg(w, W_profile);
}

// event compute_flx_wT(t+=dt_mean, t<=tend+1e-10){
//   foreach(reduction(+:Tflx_profile[:nl]))
//     foreach_layer(){
//     #if dimension==1
//       Tflx_profile[point.l] += ((w[]-W_profile[point.l])*(T[]-T_profile[point.l])) / N;
//     #else
//       Tflx_profile[point.l] += ((w[]-W_profile[point.l])*(T[]-T_profile[point.l])) / (N*N);
//     #endif
//   }
// }


// This even writes to a file the layer average
event write_diag(t=0., t+=dt_mean){
    //write_profile("T_profile.dat", T_profile, fp1);
    if (pid()==0){
      for (int ik=0; ik<nl; ++ik) {
          fprintf (stdout, "%f %d %g\n", t, ik, T_profile[ik]);
        }
      fprintf(stdout,"\n");
    }
    //write_profile("Tflx_profile.dat", Tflx_profile, fp2);

    // Reset the profile for all workers
    for (int i=0; i<nl; ++i) {
      T_profile[i] = 0.0;
      //Tflx_profile[i] = 0.0;
    }
}


event cleanup(t=end){
  free(T_profile);
  //free(Tflx_profile);
}

/*

## TO DO : 
- add dimensions
*/


/**
Results: plots
~~~pythonplot

import numpy as np
import matplotlib.pyplot as plt

# Load data: columns = time, layer, temperature
with open("out", "r") as f:
    lines = []
    for line in f:
        if line.lstrip().startswith("#"):
            break
        lines.append(line)

data = np.loadtxt(lines)
# data = np.loadtxt("out", usecols=(0,1,2))

times = np.unique(data[:, 0])

# Plasma colormap
cmap = plt.cm.plasma
norm = plt.Normalize(times.min(), times.max())

fig, ax = plt.subplots(figsize=(10, 7))

# Plot one temperature profile per time instant
for t in times:
    
    profile = data[data[:, 0] == t]

    layer = profile[:, 1]
    temperature = profile[:, 2]

    ax.plot(
        temperature,
        layer,
        color=cmap(norm(t)),
        linewidth=1.5
    )

ax.set_xlabel("Temperature (°C)")
ax.set_ylabel("Layer")

ax.set_xlim(19.90, 20.01)
ax.set_ylim(0, 127)

# No legend or time labels
ax.grid(True, alpha=0.25)

plt.tight_layout()
plt.savefig("T_profiles.png", dpi=300)
plt.show()


~~~

**/
