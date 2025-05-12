/**
# The Immiscible Rayleigh-Taylor Instability

The incompressible, immiscible Rayleigh-Taylor instability occurs when a lighter
fluid pushes against a heavier one, causing finger-like structures to form at
their interface. This leads to characteristic mushroom-shaped patterns as the
fluids mix, driven by gravitational or accelerational forces. A typical
simulation looks something like this:

![The movie shows the interface between two immiscible fluids corresponding to a
case with $\lbrace\mathcal{A}=0.05$, $Re=12.5$, $\mathsf{B}=0.14$,
$\mathsf{S}=3.00$, $\mathsf{D}=0.00\rbrace$. Simulation uses a maximum
refinement level 10 and ran using 4096 MPI cores at the CCRT-TGCC. The complete
parameters can be found [here](default_parameters_rti.json)](rti3D.mp4)

## 1. General equations

### 1.1. Dimensional equations

We consider a system of two immiscible, incompressible fluids, each with density
$\rho_i$, where the subscript $i \in \lbrace 1; 2 \rbrace$ identifies the
respective fluid. Here, $\rho_1$ denotes the density of the lighter fluid and
$\rho_2 > \rho_1$ represents the denser fluid. The fluids are confined within a
cuboidal container. In a Cartesian coordinate frame $\vec{x} = (x, y, z)$ moving
with the container, the system is subjected to the ambient gravitational field
$g$, which acts in the vertical $z$-direction.

This system is governed by the incompressible Navier-Stokes equations, which
describe the dynamics of the velocity field $\vec{U}_i(\vec{x},t)$ and the
pressure $P(\vec{x},t)$. For each phase, this yields:

$$
\begin{aligned}
  \partial_t\vec{U}_{i} + \vec{U}_{i}\cdot(\nabla\vec{U}_{i}) &= -\frac{1}{\rho_i}\nabla\cdot\vec{\tau}_i - \vec{a}(t)
  \\
  \partial_t\rho_{i} + \vec{U}_{i}\cdot(\nabla\rho_{i}) &= 0
  \\
  \nabla\cdot\vec{U}_{i} &= 0 
\end{aligned}
$$

which express the conservation of momentum, the conservation of mass and the
incompressibility condition, respectively. The external acceleration is
$\vec{a}(t) = g\vec{e}_z$, and $\vec{\tau}_i$ is the stress tensor:

$$
\begin{aligned}
  \vec{\tau}_i = -P_i\vec{I} + \mu_i (\nabla\vec{U}_i+\nabla\vec{U}_i^T) + \gamma\kappa\delta_s(\vec{n}\otimes\vec{n})
\end{aligned}
$$

where $\mu_i$, $\gamma$, $\kappa$, $\delta_s$, and $\vec{n}$ represent the
viscosity, surface tension coefficient, local curvature, a surface Dirac delta
function, and the unit normal vector of the free surface, respectively. We may
also write the viscous stress tensor $\Pi_{\nu_i} = (\nabla \vec{U}_i +
{\nabla \vec{U}_i}^T)$ and the effect of surface tension $\vec{f}_\gamma =
\kappa \delta_s \vec{n}$ for compactness.

### 1.2. Dimensionless equations

These variables, along with the spatial coordinate $\vec{x}$ and time $t$, are
nondimensionalized using the characteristic length $1/k_0$, linear growth rate
$\sigma=\sqrt{\mathcal{A} g k_0}$, and a reference density $\rho_0 = (\rho_1 +
\rho_2 )/2$. Let us introduce the following dimensionless variables:

$$
\begin{aligned}
  \vec{x}^* = k_0 x, \quad 
  t^* = \sqrt{ \mathcal{A} g k_0 }t, \quad 
  \vec{U}_i^* = \frac{\vec{U}_i}{\sqrt{{\mathcal{A} g}/{k_0}}}, \quad 
  P_i^* = \frac{P_i}{\rho_0 {\mathcal{A} g}/{k_0}}, \quad 
  \rho_i^* = \frac{\rho_i}{\rho_0}, \quad 
\end{aligned}
$$

The nondimensionalized governing equations become:

$$
\begin{aligned}
    \partial_{t}\vec{U}_{1} + \vec{U}_{1}\cdot(\nabla\vec{U}_{1}) &=     
    -\frac{1}{\rho_1}\nabla P_1     
    + \frac{\Upsilon}{\mathrm{Re}} \nabla \cdot \Pi_{\nu_1}
    + \frac{2\mathcal{A}}{\mathrm{We}} \frac{\vec{f}_\gamma}{\rho_1}
    - \vec{e}_z,
    \\
    \partial_{t}\vec{U}_{2} + \vec{U}_{2}\cdot(\nabla\vec{U}_{2}) &=     
    -\frac{1}{\rho_2}\nabla P_2    
    +\frac{1}{\mathrm{Re}}\nabla\cdot\Pi_{\nu_2}
    +\frac{2\mathcal{A}}{\mathrm{We}} \frac{\vec{f}_\gamma}{\rho_2}
    - \vec{e}_z,
    \\
    \partial_{t}\rho_{i} + \vec{U}_{i}\cdot(\nabla\rho_{i}) &= 0
    \\
    \nabla\cdot\vec{U}_{i} &= 0,
\end{aligned}
$$

where the $^*$ symbols are omitted for clarity. The key nondimensional
parameters are:

1. The Reynolds number: $\mathrm{Re} = \sqrt{\mathcal{A} g}/\nu_2 k_0^{3/2}$
2. The Weber number: $\mathrm{We} = (\rho_2-\rho_1) \mathcal{A} g/ k_0^2 \gamma$
3. The Atwood number: $\mathcal{A} = (\rho_2-\rho_1)/(\rho_2+\rho_1)$
4. The viscosity ratio: $\Upsilon = \mu_1/\mu_2$

### 1.3. One-fluid formulation

In the one-fluid formulation, we may introduce a volume fraction field
$\phi_i(\vec{x},t)$, representing the spatial distributions of each fluid:

$$
\begin{aligned}
  \phi_1(\vec{x},t) = \frac{V_1}{V_1+V_2}, \quad
  \phi_2(\vec{x},t) = \frac{V_2}{V_1+V_2}, \quad 
  \phi_1(\vec{x},t) + \phi_2(\vec{x},t) = 1
\end{aligned}
$$

In practice, we may use only the volume fraction of heavy fluid and drop the
subscript. Since fluids are non-miscible, $\phi(\vec{x},t)$ is either 0 or 1
everywhere but on the interface.

The renormalized density $\rho$ is related to the volume fraction of heavy fluid
$\phi$ 

$$
\begin{aligned}
  \rho(\phi) = (1 + \mathcal{A})\phi   + (1 - \mathcal{A})(1-\phi).
\end{aligned}
$$

Additionally, we introduce average velocity and pressure fields,
$\vec{U}=\phi\vec{U}_2 + (1-\phi)\vec{U}_1$ and $P=\phi P_2 + (1-\phi)P_1$.
Also, we introduce a viscous stress tensor
$\Pi_{\nu}=\phi\Pi_{\nu_2}+(1-\phi)\Upsilon\Pi_{\nu_1}$. This yields

$$
\begin{aligned}
  \partial_{t}\vec{U} + \vec{U}\cdot(\nabla\vec{U}) &= 
  -\frac{1}{\rho(\phi)}\nabla P    
  + \frac{1}{\mathrm{Re}} \nabla \cdot \Pi_{\nu}
  + \frac{2\mathcal{A}}{\mathrm{We}} \frac{\vec{f}_\gamma}{\rho(\phi)}
  - \vec{e}_z,
  \\
  \partial_{t}\phi + \vec{U}\cdot(\nabla \phi) &= 0,
  \\
  \nabla\cdot\vec{U} &= 0,
\end{aligned}
$$

In the nondimensionalized system, the viscous effects are characterized by the
Reynolds number, $\mathrm{Re} = \sqrt{\mathcal{A} g}/\nu_2 k_0^{3/2}$, where
$\nu_2$ is the kinematic viscosity of the denser fluid. The viscous stress
tensor $\Pi_\nu$ differs between the two Newtonian fluids due to their distinct
viscosities: In fluid 2 (the denser fluid), the viscous stress tensor is given
by $\Pi_\nu = (\nabla\vec{U}+\nabla\vec{U}^T)$. In fluid 1 (the lighter fluid),
the viscous stress tensor is modified by the viscosity ratio $\Upsilon =
\mu_1/\mu_2$, and is expressed as $\Pi_\nu =
\Upsilon(\nabla\vec{U}+\nabla\vec{U}^T)$. The viscosity ratio $\Upsilon$
introduces asymmetry in the viscous contributions between the two fluids,
influencing the dynamics at the interface and throughout the flow.  The effect of
surface tension $\gamma$ between the two fluids is incorporated into the system
through the Weber number, $\mathrm{We} = (\rho_2-\rho_1) \mathcal{A} g/ k_0^2
\gamma$, which quantifies the relative importance of inertial forces to surface
tension forces at the interface. A higher Weber number indicates that inertia
dominates, while a lower Weber number suggests that surface tension is more
significant in stabilizing the interface. The influence of surface tension is
expressed through a renormalized force, $\vec{f}_\gamma$, which is localized at
the fluid interface. This force depends on the renormalized curvature of the
interface $\kappa$ and acts to minimize the interfacial area, thereby
stabilizing distortions caused by the flow. In the nondimensionalized form, the
force expressing surface tension can be written as $\vec{f}_\gamma =
\kappa\delta_s\vec{n}$ where $\vec{n}$ is the unit normal vector to the
interface and $\delta_s$ the surface Dirac function localized at the interface
([Tryggvason et al. (2011)](@Tryggvason2011)).

To facilitate the implementation in Basilisk, we may write the equations in 
the following form:

$$
\begin{aligned}
  \partial_{t}\vec{U} + \vec{U}\cdot(\nabla\vec{U}) &= 
  \frac{1}{\rho(\phi)} \left[
  - \nabla P
  + \mu(\phi) \nabla \cdot (\nabla\vec{U} + \nabla\vec{U}^T)
  + \sigma\vec{f}_\sigma
  \right]
  + \vec{a},
  \\
  \partial_{t}\phi + \vec{U}\cdot(\nabla \phi) &= 0,
  \\
  \nabla\cdot\vec{U} &= 0,
\end{aligned}
$$

which can be solved using existing Basilisk code.

### 1.4. Initial conditions

The instability growth rate depends on the parameters $\lbrace\mathrm{Re},
\mathrm{We}, \mathcal{A}, \Upsilon\rbrace$, but also on the initial shape of the
interface. We consider a flat interface with an annular spectrum for the
interface perturbation as in [Dimonte et al. (2004)](@dimonte2004) characterized
by:

- a mean wavenumber $k_0$
- a perturbation bandwidth $\Delta k$
- and a rms amplitude $\eta_0$.

As in [Thévenin et al (2024)](@thevenin2024), the initial perturbation is defined 
by two dimensionless numbers
$$
  \mathsf{B} \equiv \frac{\Delta k}{k_0}, \quad
  \mathsf{S} \equiv k_0 \eta_0
$$
while the renormalized diffusion thickness of the interface $\mathrm{D}$ tends
to zero. Here, $\mathsf{B}$ indicates a narrow or large band multi-mode
perturbation, while $\mathsf{S}$ measures the steepness of the interface
perturbation. 

## 2. Implementation in Basilisk
### 2.1. Include solver blocks
We use a combination of the two-phase incompressible solver with embedded
boundaries and periodic boundaries.

*/

#define MAXLEVEL 8
#define FILTERED
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "navier-stokes/conserving.h"
#include "tension.h"
#ifdef REDUCED
#include "reduced.h"
#endif

/**
### 2.2. Set problem parameters

In order to comply with the nondimensionalization proposed, the parameters are
listed as follows:

| Parameter    | Value                                   |
|--------------|-----------------------------------------|
| Width        | $W k_0$                                 |
| Height       | $H k_0$                                 |
| Depth        | $D k_0$                                 |
| Density 1    | $(1 + \mathcal{A})$                     |
| Density 2    | $(1 - \mathcal{A})$                     |
| Viscosity 1  | $\Upsilon/\mathrm{Re}$                  |
| Viscosity 2  | $1/\mathrm{Re}$                         |
| S. tension   | $2\mathcal{A}/\mathrm{We}$              |
| Gravity      | $1$                                     |

We'll read the problem parameters from an external file using the
[`.json`](https://en.wikipedia.org/wiki/JSON) file format and the `cJSON`
standard library. To look at the full list of parameters we can try
```
>> jq . default_parameters_rti2d.json
{
  "control": {
    "R": 1.0,
    "B": 1.0,
    "S": 1.0,
    "D": 0
  },
  "fluids": {
    "atwood": -0.5,
    "density1": 1.0,
    "density2": 3.0,
    "viscosity1": 0.00075,
    "viscosity2": 0.00225,
    "tension": 0.020
  },
  ...
}
```

*/
#include "rti_parameters.h"
double default_values[26]; // Array to hold parsed values

vector h[];

int main(){

  /** First, we parse the `.json` file and store the results inside
  `default_values` */
  #if dimension == 2
    const char *filename = "parameters_rti2d.json";
  #else
    const char *filename = "parameters_rti3d.json";
  #endif  
  if (file_exists(filename)){
    read_and_parse_json_mpi(filename, default_values);
  }
  else {
    read_and_parse_json_mpi("../default_parameters_rti.json", default_values);
  }

  /** Then, assign the fluid properties, */ 
  ref.At = default_values[0];
  rho1 = default_values[1];            // for f = 1
  rho2 = default_values[2];            // for f = 0
  mu1 = default_values[3];             // for f = 1
  mu2 = default_values[4];             // for f = 0
  f.sigma = default_values[5];
  f.height = h;
  f.refine = f.prolongation = fraction_refine;
  p.nodump = false;

  /** the geometric parameters, */ 
  L0 = default_values[7] * default_values[9] / ref.Lx;
  X0 = -L0 / 2.;
  Y0 = -L0 / 2.;
  #if dimension == 3
    Z0 = -L0 / 2.;
  #endif

  /** the forcing parameters (if any), */ 
  force.G0 = default_values[10];
  force.V0 = default_values[11];
  force.freq0 = default_values[12];
  force.V0prev = default_values[13];
  force.omega0 = force.freq0 * (2 * pi);
  force.period0 = (2 * pi) / force.omega0;
  force.F0 = force.V0 * force.omega0 / force.G0;
  force.F0prev = force.V0prev * force.omega0 / force.G0;
  #ifdef REDUCED
    #if dimension == 2
      G.y = -force.G0;
    #else
      G.z = -force.G0;
    #endif    
  #endif

  /** and the initial perturbation */ 
  defect.eta0 = default_values[19] / ref.Lx;
  defect.k0 = default_values[20];
  defect.dk = default_values[21];
  defect.kmin = defect.k0 - defect.dk / 2.;
  defect.kmax = defect.k0 + defect.dk / 2.;

  /** Now, we set-up the numerical parameters */ 
  N = 1 << 4;
  CFL = default_values[21];
  DT = default_values[22] * force.period0;
  TOLERANCE = default_values[23];
  NITERMIN = default_values[24];

  /**   
  and compute the dimensionless numbers
  */ 
  sim.At = (rho2-rho1)/(rho2+rho1);                         // Atwood number
  sim.Re = (rho2/mu2)*(force.omega0);                       // Reynolds number
  sim.Fr = sq(force.omega0)/(force.G0);                     // Froude number
  sim.We = sq(force.omega0) * (rho2-rho1)/f.sigma;          // Weber number
  sim.Bo = sim.Fr / sim.We;                                 // Bond number
  sim.Oh = sqrt(sim.We) / sim.Re;                           // Ohnesorge number
  sim.B = defect.dk / defect.k0;
  sim.S = defect.k0 * defect.eta0;
  sim.D = 0.;
  if (pid() == 0){
    fputs("RTI\n", stderr);
    print_parameters(ref, force, defect, sim);
  }

  periodic(left);
  #if dimension == 3
    periodic(top);
  #endif
  run();
}

/** 
### 2.3. Set grid refinement
*/

#if TREE
event adapt(i++){
  #if dimension == 2
    adapt_wavelet({f, rhov, u}, (double[]){default_values[23], default_values[23], default_values[23], default_values[23]}, maxlevel = MAXLEVEL);
  #else
    adapt_wavelet({f, rhov, u}, (double[]){default_values[23], default_values[23], default_values[23], default_values[23], default_values[23]}, maxlevel = MAXLEVEL);
  #endif  
}
#endif

/** 
### 2.4. Remove tiny droplets
Every now and then, we remove the droplets that are too small to be properly 
resolved.
*/

#include "tag.h"
event drop_remove (i += 20) {
  remove_droplets (f, 1, 0);
}

/** 
### 2.5. Set initial conditions
We consider an annular spectrum for the interface perturbation as in [Dimonte et
al. (2004)](https://doi.org/10.1063/1.1688328). Additionally, the velocity field
is set to zero.
*/

#define _W0 (default_values[9] / ref.Lx)
#define _mindel (L0 / (1 << MAXLEVEL))
#if dimension == 2
  #define rectanglebox(extra) intersection((_W0 / 2 + extra - x), (_W0 / 2 + extra + x))
#else
  #define rectanglebox(extra) intersection((L0 / 8 + extra - y), (L0 / 8 + extra + y))
#endif

#include "view.h"
#if dimension == 2
  #include "../input_fields/initial_conditions_dimonte_fft1.h"
#else
  #include "../input_fields/initial_conditions_dimonte_fft2.h"
#endif
event init(i = 0){
  if (!restore(file = "./backup")){

    #if TREE
      #if dimension == 2
        refine(((y < L0/16) && (y > -L0/16)) && level < MAXLEVEL);
      #else
        refine(((z < L0/16) && (z > -L0/16)) && level < MAXLEVEL);
      #endif      
    #endif

    #if dimension == 2
    {
      vertex scalar phi[];
      initial_condition_dimonte_fft(phi, amplitude=defect.eta0, NX=(1 << MAXLEVEL), kmin = defect.kmin, kmax = defect.kmax);
      fractions(phi, f);
    }
    #else
    {
      vertex scalar phi[];
      initial_condition_dimonte_fft2(phi, amplitude=defect.eta0, NX=(1 << MAXLEVEL), NY=(1 << MAXLEVEL), kmin = defect.kmin, kmax = defect.kmax, isvof=1);
      fractions(phi, f);
    }
    #endif

    foreach(){
      foreach_dimension(){
        u.x[] = 0.;
      }
    }
  }

/** Then, visualize the initial conditions just to make sure  */
  {
    #if dimension == 3
    view(camera="bottom");
    #endif
    draw_vof("f");
    squares("f", linear = false, n = {0, 0, 1}, alpha = 0.);
    save("init_f.png");

    box();
    cells(n = {0, 0, 1}, alpha = 0.);
    save("init_grid.png");
  }

  #if dimension == 3
  {
    view(camera="iso");
    box();
    squares("f", linear = false, n = {1, 0, 0}, alpha = 0.);
    squares("f", linear = false, n = {0, 1, 0}, alpha = 0.);
    squares("f", linear = false, n = {0, 0, 1}, alpha = 0.);
    save("init2_f.png");

    view(camera="iso");
    box();
    cells(n = {1, 0, 0}, alpha = 0.);
    cells(n = {0, 1, 0}, alpha = 0.);
    cells(n = {0, 0, 1}, alpha = 0.);
    save("init2_grid.png");
  }
  #endif
}

/** 
### 2.6. Set boundary conditions
We consider periodic boundaries.

### 2.7. Apply the external forcing

We have the option of using an horizontal acceleration to stabilize the
interface,
$$
\vec{a} = b \omega^2 \sin(\omega t) \hat{e}_x - g \hat{e}_y
$$
where $b$ is the oscillation amplitude and $\omega$ the forcing frequency.
To prevent sloshing, we tipically use a ramping function.

#### 2.7.1. Functions to calculate a ramp value based on input parameters
We have the option of using a sigmoid function
*/
double sigmoid(double t1, double k){
  return 1.0 / (1.0 + exp(-k * (t - t1)));
}

/**
or a soft-ramp function 
*/
double soft_ramp(double t1, double p0, double p1, double k, double m){
  double ramp;
  double var1 = log(exp((p1 - p0) * k) - 1) / (k * m);
  if (t < t1 + var1){
    ramp = (1 / k) * log(1 + exp(k * m * (t - t1))) + p0;
  }
  else{
    ramp = p1;
  }
  return ramp;
}
/**
which multiply a periodic acceleration in the horizontal direction

#### 2.7.2. Add the external forcing to the acceleration term
*/
event acceleration(i++)
{
  double p0 = 0.0, p1 = force.F0 / force.F0prev;
  double t0 = default_values[25], t1 = default_values[15] / ref.T;
  double m = default_values[16] * ref.T;
  double k = default_values[17];
  if (default_values[14] == 0){
    force.ramp = soft_ramp(t1 - t0, p0, p1, k, m);
    force.Gn = (force.G0) * (force.F0prev * sin(force.omega0 * t));
  }
  else{
    force.ramp = sigmoid(t1 - t0, k);
    force.Gn = (force.G0) * (force.F0 * sin(force.omega0 * t));
  }

  face vector av = a;
  foreach_face(x){
      av.x[] -= force.Gn * force.ramp;
  }

  force.Gnm1 = force.Gn;
  #ifndef REDUCED
    #if dimension == 2
      foreach_face(y){
        av.y[] -= force.G0;
      }
    #else
      foreach_face(z){
        av.z[] -= force.G0;
      }
    #endif
  #endif
}

/** #### 2.7.3. We also tweak the CFL condition to take into account the acceleration */

#include "acceleration_cfl.h"



/** 
## 3. Outputs

To examine the dynamics of the Rayleigh-Taylor Instability, we focus on the
zero- dimensional (0D) turbulent quantities derived from averaging across the
entire width of the mixing layer, which are commonly utilized in mixing models.
Given that the RTI is statistically homogeneous in the horizontal plane, it is
practical to introduce the horizontal average $\bar{Q}$ and the fluctuation $Q'$
of a quantity $Q$, as determined by the Reynolds decomposition $Q = \bar{Q} +
Q'$. 

Within the context of the Boussinesq approximation, the RTI exhibits a zero mean
velocity, $\bar{\vec{u}} = 0$. Moreover, the volume fraction field can be
expressed as $\mathsf{f}(\vec{x}, t) = \bar{\mathsf{f}}(z, t) +
\mathsf{f}'(\vec{x}, t)$, allowing to define a mixing zone size as
$$
L(t) = 6 \int \bar{\mathsf{f}}(1 - \bar{\mathsf{f}}) dz
$$
derived from a piecewise mean concentration profile representing a turbulent
Rayleigh-Taylor mixing zone ([Andrews & Spalding,
1990](https://www.doi.org/10.1063/1.857652)).

For now, we'll extract the horizontally averaged profiles for different
quantities, such that 0D quantities may be obtained in post-process. 
*/


event logfile(i++; t <= default_values[18]){
  fprintf(stdout, " res: %.5f \t %g %d %d \t %g %d %d \n", t, mgp.resa, mgp.i, mgp.nrelax, mgu.resa, mgu.i, mgu.nrelax);
}

#include "../output_fields/xdmf/output_xdmf.h"
#include "../output_fields/vtu/output_vtu.h"
#include "../output_fields/profiles_foreach_region.h"

#ifndef _NOOUTPUTS
  #if dimension ==2 
    #include "rti2D_outputs.h"
  #else
    #include "rti3D_outputs.h"
  #endif
#endif
event finalize(t = end){
  dump("backup");
  squares("f", linear = false);
  save("final_f.png");

  cells();
  save("final_grid.png");
}

event backups(t += 0.05){
  dump("backup");

  squares("f", linear = false);
  save("backup_f.png");

  cells();
  save("backup_grid.png");
}



/**
# References

~~~bib

@article{grea2019,
  title={Frozen waves in turbulent mixing layers},
  author={Gr{\'e}a, Beno{\^\i}t-Joseph and Briard, Antoine},
  journal={Physical Review Fluids},
  volume={4},
  number={6},
  pages={064608},
  year={2019},
  publisher={APS}
}

@book{tryggvason2011,
  title={Direct numerical simulations of gas--liquid multiphase flows},
  author={Tryggvason, Gr{\'e}tar and Scardovelli, Ruben and Zaleski, St{\'e}phane},
  year={2011},
  publisher={Cambridge university press}
}

@article{bunner2002,
  title={Dynamics of homogeneous bubbly flows Part 1. Rise velocity and microstructure of the bubbles},
  author={Bunner, Bernard and Tryggvason, Gr{\'e}tar},
  journal={Journal of Fluid Mechanics},
  volume={466},
  pages={17--52},
  year={2002},
  publisher={Cambridge University Press}
}

@article{dimonte2004,
  title={Dependence of turbulent Rayleigh-Taylor instability on initial perturbations},
  author={Dimonte, Guy},
  journal={Physical Review E—Statistical, Nonlinear, and Soft Matter Physics},
  volume={69},
  number={5},
  pages={056305},
  year={2004},
  publisher={APS}
}

@article{thevenin2024,
  title={The memory of Rayleigh-Taylor turbulence},
  author={Th{\'e}venin, S{\'e}bastien and Gr{\'e}a, B-J and Kluth, Gilles and Nadiga, Balu},
  journal={arXiv preprint arXiv:2403.17832},
  year={2024}
}

~~~
*/