/**

# How to use diffusion.h and dr.h together

This 1D case is an example on how to set boundary condition for a scalar T.
Temperature is initialized as an affine function of $z$, at the surface a (heat) flux is
imposed and at the bottom we impose a that the initial gradient is conserved at
all $t$.

If using  this with a machine with 2 gpus, make sure to select the right one.
For my nvidia gpu, I do
__NV_PRIME_RENDER_OFFLOAD=1 __GLX_VENDOR_LIBRARY_NAME=nvidia make

*/
#include "grid/multigrid.h"
#include "layered/hydro.h"
#include "layered/nh.h"
#include "layered/remap.h"
#include "layered/perfs.h"
#include "hugoj/lib/diffusionH.h"
#include "bderembl/libs/netcdf_bas.h"

double tend = 3600.0;
double smalltime = 1e-10;
double H0 = 100.;
// -> stratification related
double strat = 0.000002;       // [s-2] N^2 stratification
double Ts = 20.;              // [K] Surface temperature (arbitrary)
double qt = -500.;         // [W.m-2] Heat flux
double rho0 = 1025.;     // [kg.m-3] reference density
double cp = 4.2e3;      // [J.kg-1.K-1] heat capacity water
double alphaT = 2e-4;        // [K-1] Thermal expansion coeff for water
double Pr = 1.;                       // Prandtl number = nu/diffT
double diffT = 1.5e-5;         // [m2.s-1] Scalar vertical diffusion coeff
double T0 = 20.;              // [°C] Reference temperature
double Trand = 0.001;           // [°C] Random temperature perturbution
const double g_ = 9.81;        // [m.s-2] Gravity
#define drho(T) (alphaT*(T0-T))       // Linear equation of state (Vallis 2.4)
#define Tini(z) strat/(g_*alphaT)*z + Ts
#include "layered/dr.h"
double fluxbot,fluxtop;

double* Temp;

static FILE * fp;

int main(int argc, char *argv[])  
{
  L0 = 50.;
  nu = 0.01;
  diffT = nu/Pr;
  N = 1; 
  nl = 30;
  G = 9.81;
  theta_H = 0.51;
  CFL_H = 8.;
  CFL = 0.8;
  
  #if HDIFF || T_CST
  N = 8;
  #endif

  Temp = (double *)calloc(nl, sizeof(double));

  origin (-L0/2., -L0/2.);
  
  /**

  Boundary condition for temperature. We test two cases: i) imposing the initial
  stratification at top and bottom, ii) impose a destabilising flux at top and
  initial stratification at bottom. 

  */
  fluxbot = strat/(g_*alphaT);
#if HEATING
  fluxtop = qt/(diffT*rho0*cp);
#else
  fluxtop = strat/(g_*alphaT); 
#endif // HEATING

#if NEUMANN0
  fluxbot = 0.;
  fluxtop = 0.;
#endif // NEUMANN0

#if NOT_PERIODIC
#else
#if dimension==2
    periodic (top);
 #endif
  periodic (left);
#endif // NOT_PERIODIC
  run();
}

event init(i =  0) {

  foreach() {
    zb[] = -H0;
    eta[] = 0.;
    double H = - zb[];
    double z = zb[];
    foreach_layer() {
      h[] = H/nl;
      z += h[]/2.;
      foreach_dimension()
        u.x[] = 0.;
      w[] = 0.;
    #if T_CST
      T[] = Ts;
    #else
      T[] = Tini(z); // + Trand * noise()*exp(z/100.) ;
    #endif // T_CST
      z += h[]/2.;
    }
  }

  fp  = fopen("T_profile.dat","w"); // reset file
  fclose(fp);
  create_nc({zb, h, u, w, eta, T}, "out.nc");
}

/**
We impose the stratification at the top and the bottom to be the initial
stratification. We should not see any change in the profile.
*/

event viscous_term (i++)
{
  foreach() {
    vertical_diffusion2 (point,   // point
                        h,        // h
                        T,        // scalar
                        dt,       // dt
                        diffT,    // D
                        fluxtop,  // dst
                        fluxbot); // dsb
  }
  #if HDIFF
  horizontal_diffusion ({T}, diffT, dt);
  horizontal_diffusion ({u, w}, nu, dt);
  #endif // HDIFF
}

event log (i++){
  fp  = fopen("T_profile.dat","a");
  if (fp == NULL){
    fprintf(stderr, "Error opening file T_profile.dat");
    return 2;
  }

  #if HDIFF || T_CST
  foreach (reduction(+:Temp[:nl])){
    foreach_layer (){
      Temp[point.l] += T[];
    }
  }
  for (int kl=0; kl<nl; kl++){
    Temp[kl] /= N*N;
    fprintf (fp, "%f %d %d %g\n", t, i, kl, Temp[kl]);
    Temp[kl] = 0.;
  }
  #else
  foreach() {
    foreach_layer()
      fprintf (fp, "%f %d %d %g\n", t, i, point.l, T[]);
  }
  #endif 
  fprintf(fp,"\n\n");
  fclose(fp);
  write_nc();
}

event stop (t = tend){
  free(Temp);
}


/**
## ----------------------
## HOLD STRATIF

~~~pythonplot Profile of T with imposed gradient at the top and the bottom of the domain imposed to be the initial stratification: we expect that nothing moves
import numpy as np
import matplotlib.pyplot as plt
data = np.loadtxt("T_profile.dat")
nl=30
nt = data.shape[0]//nl
fig, ax = plt.subplots(figsize=(8, 6))
cmap = plt.get_cmap("viridis", nt)
for t in range(nt):
    layer=data[t*nl:(t+1)*nl,2]
    T = data[t*nl:(t+1)*nl,3]
    ax.plot(T,layer,color=cmap(t), marker="+", linestyle="-")
ax.set_xlabel("T")
ax.set_ylabel("Layer")
ax.set_title("Temperature profiles (hold stratif)")
ax.set_xlim([19.875,20.025])
#ax.legend(loc="upper left", bbox_to_anchor=(1, 1))
plt.tight_layout()
plt.savefig("T_profiles.png", dpi=150)
~~~


## ----------------------
## HEATING

~~~pythonplot Profile of T with surface cooling
data = np.loadtxt("../diff_test_H/T_profile.dat")
nl=30
nt = data.shape[0]//nl
fig, ax = plt.subplots(figsize=(8, 6))
cmap = plt.get_cmap("viridis", nt)
for t in range(nt):
    layer=data[t*nl:(t+1)*nl,2]
    T = data[t*nl:(t+1)*nl,3]
    ax.plot(T,layer,color=cmap(t), marker="+", linestyle="-")
ax.set_xlabel("T")
ax.set_ylabel("Layer")
ax.set_title("Temperature profiles (Q<0)")
ax.set_xlim([19.875,20.025])
#ax.legend(loc="upper left", bbox_to_anchor=(1, 1))
plt.tight_layout()
plt.savefig("T_profiles_H.png", dpi=150)
~~~

## ----------------------
## NEUMANN0

~~~pythonplot Profile of T Neumann boundary conditions
data = np.loadtxt("../neumann0/T_profile.dat")
nl=30
nt = data.shape[0]//nl
fig, ax = plt.subplots(figsize=(8, 6))
cmap = plt.get_cmap("viridis", nt)
for t in range(nt):
    layer=data[t*nl:(t+1)*nl,2]
    T = data[t*nl:(t+1)*nl,3]
    ax.plot(T,layer,color=cmap(t), marker="+", linestyle="-")
ax.set_xlabel("T")
ax.set_ylabel("Layer")
ax.set_title("Temperature profiles (dTdz=0 top, bot)")
ax.set_xlim([19.875,20.025])
#ax.legend(loc="upper left", bbox_to_anchor=(1, 1))
plt.tight_layout()
plt.savefig("T_profiles_neumann0.png", dpi=150)
~~~

## ----------------------
## NEUMANN0, HDIFF

~~~pythonplot Profile of T Neumann boundary conditions and horizontal diffusion
data = np.loadtxt("../with_diff/T_profile.dat")
nl=30
nt = data.shape[0]//nl
fig, ax = plt.subplots(figsize=(8, 6))
cmap = plt.get_cmap("viridis", nt)
for t in range(nt):
    layer=data[t*nl:(t+1)*nl,2]
    T = data[t*nl:(t+1)*nl,3]
    ax.plot(T,layer,color=cmap(t), marker="+", linestyle="-")
ax.set_xlabel("T")
ax.set_ylabel("Layer")
ax.set_title("Temperature profiles (w/ hdiff, dTdz=0 bot and top)")
ax.set_xlim([19.875,20.025])
#ax.legend(loc="upper left", bbox_to_anchor=(1, 1))
plt.tight_layout()
plt.savefig("T_profiles_with_diff.png", dpi=150)
~~~

## ----------------------
## NEUMANN0, HDIFF, CUDA
~~~pythonplot Profile of T Neumann boundary conditions and horizontal diffusion (cuda)
data = np.loadtxt("../with_diff.cuda/T_profile.dat")
nl=30
nt = data.shape[0]//nl
fig, ax = plt.subplots(figsize=(8, 6))
cmap = plt.get_cmap("viridis", nt)
for t in range(nt):
    layer=data[t*nl:(t+1)*nl,2]
    T = data[t*nl:(t+1)*nl,3]
    ax.plot(T,layer,color=cmap(t), marker="+", linestyle="-")
ax.set_xlabel("T")
ax.set_ylabel("Layer")
ax.set_title("Temperature profiles (w/ hdiff, dTdz=0 bot and top, cuda)")
ax.set_xlim([19.875,20.025])
#ax.legend(loc="upper left", bbox_to_anchor=(1, 1))
plt.tight_layout()
plt.savefig("T_profiles_with_diff_cuda.png", dpi=150)
~~~

## ----------------------
## NEUMANN0, T_CST

~~~pythonplot Profile of T, with Neumann boundary conditions and constant initial temperature
data = np.loadtxt("../T_cst/T_profile.dat")
nl=30
nt = data.shape[0]//nl
fig, ax = plt.subplots(figsize=(8, 6))
cmap = plt.get_cmap("viridis", nt)
for t in range(nt):
    layer=data[t*nl:(t+1)*nl,2]
    T = data[t*nl:(t+1)*nl,3]
    ax.plot(T,layer,color=cmap(t), marker="+", linestyle="-")
ax.set_xlabel("T")
ax.set_ylabel("Layer")
ax.set_title("Temperature profiles (T=cst, dT/dz=0 bot and top)")
ax.set_xlim([19.999,20.001])
#ax.legend(loc="upper left", bbox_to_anchor=(1, 1))
plt.tight_layout()
plt.savefig("T_profiles_T_cst.png", dpi=150)
~~~

## ----------------------
## NEUMANN0, T_CST, CUDA

~~~pythonplot Profile of T, with Neumann boundary conditions and constant initial temperature (cuda)
data = np.loadtxt("../T_cst.cuda/T_profile.dat")
nl=30
nt = data.shape[0]//nl
fig, ax = plt.subplots(figsize=(8, 6))
cmap = plt.get_cmap("viridis", nt)
for t in range(nt):
    layer=data[t*nl:(t+1)*nl,2]
    T = data[t*nl:(t+1)*nl,3]
    ax.plot(T,layer,color=cmap(t), marker="+", linestyle="-")
ax.set_xlabel("T")
ax.set_ylabel("Layer")
ax.set_title("Temperature profiles (T=cst, dT/dz=0 bot and top, cuda)")
ax.set_xlim([19.999,20.001])
#ax.legend(loc="upper left", bbox_to_anchor=(1, 1))
plt.tight_layout()
plt.savefig("T_profiles_T_cst_cuda.png", dpi=150)
~~~

## ----------------------
## NEUMANN0, T_CST, NOT_PERIODIC

~~~pythonplot Profile of T, with Neumann boundary conditions and constant initial temperature (cuda)
data = np.loadtxt("../T_cst_notperiodic/T_profile.dat")
nl=30
nt = data.shape[0]//nl
fig, ax = plt.subplots(figsize=(8, 6))
cmap = plt.get_cmap("viridis", nt)
for t in range(nt):
    layer=data[t*nl:(t+1)*nl,2]
    T = data[t*nl:(t+1)*nl,3]
    ax.plot(T,layer,color=cmap(t), marker="+", linestyle="-")
ax.set_xlabel("T")
ax.set_ylabel("Layer")
ax.set_title("Temperature profiles (T=cst, dT/dz=0 bot and top, notperiodic)")
ax.set_xlim([19.999,20.001])
#ax.legend(loc="upper left", bbox_to_anchor=(1, 1))
plt.tight_layout()
plt.savefig("T_profiles_T_cst_not_periodic.png", dpi=150)
~~~

## ----------------------
## NEUMANN0, T_CST, NOT_PERIODIC, CUDA

~~~pythonplot Profile of T, with Neumann boundary conditions and constant initial temperature (cuda)
data = np.loadtxt("../T_cst_notperiodic.cuda/T_profile.dat")
nl=30
nt = data.shape[0]//nl
fig, ax = plt.subplots(figsize=(8, 6))
cmap = plt.get_cmap("viridis", nt)
for t in range(nt):
    layer=data[t*nl:(t+1)*nl,2]
    T = data[t*nl:(t+1)*nl,3]
    ax.plot(T,layer,color=cmap(t), marker="+", linestyle="-")
ax.set_xlabel("T")
ax.set_ylabel("Layer")
ax.set_title("Temperature profiles (T=cst, dT/dz=0 bot and top, notperiodic,cuda)")
ax.set_xlim([19.999,20.001])
#ax.legend(loc="upper left", bbox_to_anchor=(1, 1))
plt.tight_layout()
plt.savefig("T_profiles_T_cst_not_periodic_cuda.png", dpi=150)
plt.show()
~~~

*/
