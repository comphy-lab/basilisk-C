/**

# How to use diffusion.h and dr.h together

This 1D case is an example on how to set boundary condition for a scalar T.
Temperature is initialized as affine function of $z$, at the surface a (heat) flux is
imposed and at the bottom we impose a that the initial gradient is conserved at
all $t$.

Note: imposing flux = 0 at the bottom could be done using $\lambda_b$ = HUGE but
its not clean. A proper Neumann condition would need a modification of
diffusion.h. I started that but adding an argument with default value (allowing
existing code to continue to run without modification) to the function
'vertical_diffusion' raise some error. Working on this modified diffusion.h
requires to copy/link nh.h, implicit.h, hydro.h in the directory of diffusion.h.
There might be a better way to do this...

*/
#include "grid/multigrid.h"
//#include "grid/multigrid1D.h"
#include "layered/hydro.h"
#include "layered/nh.h"
#include "layered/remap.h"
#include "layered/perfs.h"

//#include "bderembl/libs/netcdf_bas.h" // read/write netcdf files

double tend = 3600.0;
double smalltime = 1e-10;

// -> stratification related
double strat = 0.000002;       // [s-2] N^2 stratification
double Ts = 20.;              // [K] Surface temperature (arbitrary)
double qt = -500.;         // [W.m-2] Heat flux
double rho0 = 1025.;     // [kg.m-3] reference density
double cp = 4.2e3;      // [J.kg-1.K-1] heat capacity water
double alphaT = 2e-4;        // [K-1] Thermal expansion coeff for water
double Pr = 1.;                       // Prandtl number = nu/kappa
double kappa = 1.5e-5;         // [m2.s-1] Scalar vertical diffusion coeff
double T0 = 20.;              // [°C] Reference temperature
double Trand = 0.001;           // [°C] Random temperature perturbution
const double g_ = 9.81;        // [m.s-2] Gravity
#define drho(T) (alphaT*(T0-T))       // Linear equation of state (Vallis 2.4)
#define Tini(z) strat/(g_*alphaT)*z + Ts
#include "layered/dr.h"

static FILE * fp;

int main(int argc, char *argv[])  
{
  L0 = 50.;
  nu = 0.01;
  kappa = nu;
  N = 1; 
  nl = 30;
  G = 9.81;
  theta_H = 0.51;
  CFL_H = 8.;
  CFL = 0.8;
  
  origin (-L0/2., -L0/2.);
  
  #if dimension==2
    periodic (top);
  #endif
  periodic (left);

  run();
}

event init(i =  0) {

  foreach() {
    zb[] = -100.;
    eta[] = 0.;
    double H = - zb[];
    double z = zb[];
    foreach_layer() {
      h[] = H/nl;
      z += h[]/2.;
      foreach_dimension()
        u.x[] = 0.;
      w[] = 0.;
      T[] = Tini(z) + Trand * noise()*exp(z/100.) ;
      z += h[]/2.;
    }
  }

  fp  = fopen("T_profile.dat","w"); // reset file
  fclose(fp);

}


event viscous_term (i++)
{
  foreach() {
    vertical_diffusion (point, h, T, dt, kappa, qt/(kappa*rho0*cp), T[0,0,0]-strat/(g_*alphaT) , 1.); // 
  }
}

event log (i++){
  fp  = fopen("T_profile.dat","a");
  if (fp == NULL){
    fprintf(stderr, "Error opening file T_profile.dat");
    return 2;
  }
  foreach() {
    foreach_layer()
      fprintf (fp, "%f %d %d %g\n", t, i, point.l, T[]);
  }
  fprintf(fp,"\n\n");
  fclose(fp);
}

event stop (t = tend);


/**

~~~pythonplot
import numpy as np
import matplotlib.pyplot as plt
# --- load data ---
data = []
with open("T_profile.dat") as f:
    for line in f:
        line = line.strip()
        if line.startswith("#") or line == "":
            continue
        t, it, layer, T = line.split()
        data.append((float(t), int(it), int(layer), float(T)))

data = np.array(data)

# --- group by iteration ---
iterations = np.unique(data[:, 1].astype(int))

fig, ax = plt.subplots(figsize=(8, 6))
cmap = plt.get_cmap("viridis", len(iterations))

for idx, it in enumerate(iterations):
    mask = data[:, 1].astype(int) == it
    block = data[mask]
    layer = block[:, 2]
    T     = block[:, 3]
    t_val = block[0, 0]

    ax.plot(T, layer, color=cmap(idx), marker="+", linestyle="-", label=f"t = {t_val:.1f}")

ax.set_xlabel("T")
ax.set_ylabel("Layer")
ax.set_title("Temperature profiles")
ax.set_xlim([19.875,20.025])
#ax.legend(loc="upper left", bbox_to_anchor=(1, 1))
plt.tight_layout()
plt.savefig("T_profiles.png", dpi=150)
plt.show()
~~~


*/
