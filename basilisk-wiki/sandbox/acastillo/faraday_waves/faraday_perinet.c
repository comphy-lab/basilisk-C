/** 
# Super-Square Patterns in Faraday Waves

This section describes our numerical simulation of super-square patterns in
Faraday waves, specifically replicating the findings of [Kahouadji et al. (2015)
](#kahouadji2015). We use a two-fluid system (air and silicon oil) subjected to
vertical oscillation, which is governed by the incompressible Navier-Stokes
equations.

## 1. General equations

### 1.1. Dimensional equations

We model two immiscible, incompressible fluids (air, $\rho_1$; silicon oil,
$\rho_2 > \rho_1$) in a cuboidal container. The container oscillates vertically
at frequency $\omega$ and amplitude $a$, adding to the ambient gravitational
acceleration $g$.

The system's dynamics are described by the **incompressible Navier-Stokes
equations**:

$$
\begin{aligned}
  \partial_t\vec{U}_{i} + \vec{U}_{i}\cdot(\nabla\vec{U}_{i}) &= 
  \frac{1}{\rho_i}\nabla\cdot\vec{\tau}_i - \vec{a}(t)
  \\
  \nabla\cdot\vec{U}_{i} &= 0 
\end{aligned}
$$

where $\vec{U}_i$ is the velocity, $\rho_i$ is the density, $\vec{a}(t) =
(g+a\omega^2\cos(\omega t))\vec{e}_z$ is the external acceleration, and
$\vec{\tau}_i$ is the stress tensor, incorporating pressure, viscous forces
($\mu_i$), and surface tension ($\gamma$). Fluid properties are listed below:

| Fluid       | Density ($\rho$) | Viscosity ($\mu$) | Thickness ($h$) |
| :---------- | :--------------- | :---------------- | :-------------- |
| Air         | 1.205 kg/m³      | 1.82e-5 Pa·s      | 29.5 mm         |
| Silicon Oil | 965 kg/m³        | 0.02 Pa·s         | 14.5 mm         |

Surface tension ($\gamma$) at the oil-air interface is **0.02 N/m**. The authors
investigate three oscillation frequencies:

| $\omega/2\pi$ | Critical Wavenumber ($k_c$) | Critical Wavelength ($\lambda_c$) | Instability Threshold |
| :------------ | :-------------------------- | :-------------------------------- | :-------------------- |
| 30 Hz         | 551.36 1/m                  | 11.40 mm                          | 0.733g                |
| 60 Hz         | 1065.59 1/m                 | 5.90 mm                           | 2.650g                |
| 90 Hz         | 1467.84 1/m                 | 4.28 mm                           | 5.320g                |

But we are going to focus on the 30 Hz case. For this case, the authors
used a $132 \text{ mm} \times 132 \text{ mm} \times 44 \text{ mm}$ domain with
$512 \times 512 \times 256$ grid points (44 grid points per wavelength).

### 1.2. Dimensionless equations

To generalize, we non-dimensionalize the equations using characteristic length
$[L]=k^{-1}$, time $[T]=\omega^{-1}$, and reference density $\rho_0=\rho_2$.
This yields key dimensionless parameters:

| Parameter             | Symbol                 | Definition                                   | Value                |
| :-------------------- | :--------------------- | :------------------------------------------- | :------------------- |
| Reynolds number       | $\mathrm{Re}$          | $\rho_0 \omega / (k^2\mu_2)$                 | $29.918$             |
| Weber number          | $\mathrm{We}$          | $\rho_0 \omega^2 / (\gamma k^3)$             | $10.228$             |
| Froude number         | $\mathrm{Fr}$          | $\omega^2 / (gk)$                            | $6.569$              |
| Forcing term          | $\mathrm{F}$           | $a\omega^2/g$                                | $> 0.733$            |
| Viscosity ratio       | $\Upsilon$             | $\mu_1/\mu_2$                                | $9.1 \times 10^{-4}$ |
| Density ratio         | $\hat{\rho}$           | $\rho_1/\rho_2$                              | $0.00125$            |

The dimensionless governing equations become:

$$
\begin{aligned}
  \partial_{t}\vec{U}_{1} + \vec{U}_{1}\cdot(\nabla\vec{U}_{1}) &= 
  \frac{1}{\rho_1}\left[
  -\nabla P_1     
  + \frac{\Upsilon}{\mathrm{Re}} \nabla \cdot \Pi_1
  + \frac{1}{\mathrm{We}} \vec{f}_\gamma
  \right]
  - \frac{1}{\mathrm{Fr}} \left(1 + \mathrm{F}\cos(t)\right)\vec{e}_z,
  \\
  \partial_{t}\vec{U}_{2} + \vec{U}_{2}\cdot(\nabla\vec{U}_{2}) &=   
  \frac{1}{\rho_2}\left[
  - \nabla P_2    
  + \frac{1}{\mathrm{Re}} \nabla \cdot \Pi_2
  + \frac{1}{\mathrm{We}} \vec{f}_\gamma
  \right]
  - \frac{1}{\mathrm{Fr}} \left(1 + \mathrm{F}\cos(t)\right)\vec{e}_z,
  \\
  \nabla\cdot\vec{U}_{i} &= 0,
\end{aligned}
$$

### 1.3. One-fluid formulation

In practice, we use a **one-fluid formulation** with a volume fraction field
$\phi(\vec{x},t)$ for the heavy fluid. This allows us to define effective
density $\rho(\phi)$ and viscosity $\mu(\phi)$ based on the fluid distribution:

$$
\begin{aligned}
  \rho(\phi) =& \phi + \hat{\rho}(1-\phi)
  \\
  \mu(\phi) =& \frac{1}{\mathrm{Re}}\phi + \frac{\Upsilon}{\mathrm{Re}} (1-\phi)
\end{aligned}
$$

The averaged velocity $\vec{U}$ and pressure $P$ then allow for a single set of
coupled equations:

$$
\boxed{
\begin{aligned}
  \partial_{t}\vec{U} + \vec{U}\cdot(\nabla\vec{U}) &= 
  \frac{1}{\rho(\phi)}
  \left[
  -\nabla P    
  + \mu(\phi) \nabla \cdot \Pi
  + \frac{1}{\mathrm{We}} \vec{f}_\gamma
  \right]
  - \frac{1}{\mathrm{Fr}} \left(1 + \mathrm{F}\cos(t)\right)\vec{e}_z,
  \\
  \partial_{t}\phi + \vec{U}\cdot(\nabla \phi) &= 0,
  \\
  \nabla\cdot\vec{U} &= 0,
\end{aligned}
}
$$
which can be solved using existing Basilisk code.

## 2. Implementation in Basilisk

### 2.1. Include solver blocks

Our simulation utilizes Basilisk's two-phase incompressible solver with embedded
boundaries, along with modifications from
[contact-embed.h](http://basilisk.fr/sandbox/tavares/contact-embed.h) to achieve
correct aspect ratios. More details are available in [Tavares et al.
(2024)](#tavares2024).

*/

#define FILTERED
#include "embed.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "navier-stokes/conserving.h"
#include "tension.h"
#include "./my_contact-embed.h"


/**
### 2.2. Setting Problem Parameters

To align with the dimensionless equations, we've set the simulation parameters
as follows:

| Parameter         | Symbol              |
| :---------------- | :------------------ |
| Width             | $W k$               |
| Height            | $H k$               |
| Depth             | $D k$               |
| Density 1         | $\hat{\rho}$        |
| Density 2         | $1$                 |
| Viscosity 1       | $\Upsilon/\mathrm{Re}$ |
| Viscosity 2       | $1/\mathrm{Re}$     |
| Surface Tension   | $1/\mathrm{We}$     |
| Gravity           | $1/\mathrm{Fr}$     |
| Forcing Amplitude | $\mathrm{F}$        |
| Forcing Frequency | $1$                 |

For practical reasons, we are going to define some structures to hold 
parameters related to the initial perturbation and the periodic forcing.

*/
struct InitialPerturbation
{
  double a0;
  double m;
  double B;
  double kmin;
  double kmax;
};
struct InitialPerturbation perturb;

struct PeriodicForcing
{
  double G0;
  double F0;
  double F0prev;
  double freq0;
  double omega0;
  double period0;
  double Gn;
  double ramp;
};
struct PeriodicForcing force;

/** 
 
For this test we used 512 points along the horizontal direction, which is
comparable to the spatial resolution used in the referenced article. We can 
lower this resolution and the simulation lenght for the website.

*/

int maxlevel=7;     // Maximum level of refinement 
int minlevel=4;     // Minimum level of refinement 
int mesh_tol=1e-1;  // Refinement mesh criteria
double tend=1.0;  // Time end

/** First, we'll create some macros to deal with the geometry, */ 
#define _H0 ((29.5/11.40)*2*pi)
#define _mindel (L0 / (1 << maxlevel))
#if dimension == 2
  #define _D0 (0.)
  #define rectanglebox(extra) intersection((_H0 + extra - y), (_H0 + extra + y))
#else
  #define _D0 (L0)
  #define rectanglebox(extra)                                     \  
    intersection(                                                 \
    intersection((_H0 + extra - z), (_H0 + extra + z)), \
    intersection((_D0 / 2. + extra - y), (_D0 / 2. + extra + y)))
#endif

vector h[];
int main(){

  /** then, we assign the fluid properties, */ 
  rho1 = 0.0012487;
  rho2 = 1.0;
  mu1  = 3.0416796e-05;
  mu2  = 0.0334250;
  f.sigma = 0.097770;
  f.height = h;
  f.refine = f.prolongation = fraction_refine;
  p.nodump = false;

  /** the geometric parameters, */ 
  L0 = 24*pi; 
  X0 = -L0/2.;
  Y0 = -(14.5/11.40)*2*pi;
  #if dimension == 3
    Y0 = -L0 / 2.;
    Z0 = -(14.5/11.40)*2*pi;
  #endif

  /** the forcing parameters */ 
  force.G0 = 0.152230;
  force.F0 = 1.0;
  force.F0prev = 1.0;
  force.omega0 = 1.0;  
  force.period0 = (2 * pi) / force.omega0; 

  /** and the initial perturbation */ 
  perturb.a0 = 0.05; // Perturbation amplitude
  perturb.m  = 4;    // Mode index
  perturb.B  = 2;    // Bandwidth
  perturb.kmin = pi * (perturb.m - perturb.B/2) / L0;
  perturb.kmax = pi * (perturb.m + perturb.B/2) / L0;

  /** We'll also set-up the numerical parameters */ 
  N = 1 << minlevel;
  CFL = 0.50;
  DT = 0.50*force.period0;
  TOLERANCE = 1e-3;
  NITERMIN = 4;
  
  /** and display everything, */ 
  if (pid() == 0){
    fputs("SIMULATION PARAMETERS\n", stderr);
    fprintf(stderr, " Mode       : %g \n", perturb.m);
    fprintf(stderr, " Length     : %g \n", L0);
    fprintf(stderr, " Wavenumber : %g \n", perturb.kmin);
    fprintf(stderr, " Wavenumber : %g \n", perturb.kmax);
    fprintf(stderr, " Gravity    : %g \n", force.G0);
    fprintf(stderr, " Forcing    : %g \n", force.F0);
    fprintf(stderr, " Frequency  : %g \n", force.omega0);
    fprintf(stderr, " Period     : %g \n", force.period0);
    fprintf(stderr, " rho1       : %g \n", rho1);
    fprintf(stderr, " rho2       : %g \n", rho2);
    fprintf(stderr, " rho1/rho2  : %g \n", rho1 / rho2);
    fprintf(stderr, " mu1        : %g \n", mu1);
    fprintf(stderr, " mu2        : %g \n", mu2);
    fprintf(stderr, " mu1/mu2    : %g \n", mu1 / mu2);
    fprintf(stderr, " nu1        : %g \n", mu1 / rho1);
    fprintf(stderr, " nu2        : %g \n", mu2 / rho2);
    fprintf(stderr, " nu1/nu2    : %g \n", (mu1 / rho1) / (mu2 / rho2));
    fprintf(stderr, " sigma      : %g \n\n", f.sigma);
  }
  run();
}

/** 
### 2.3. Set grid refinement
*/

event adapt(i++){
  #if dimension == 2
    adapt_wavelet({f,cs,u}, (double[]){mesh_tol, mesh_tol, mesh_tol, mesh_tol}, maxlevel = maxlevel, minlevel = minlevel);
  #else
    adapt_wavelet({f,cs,u}, (double[]){mesh_tol, mesh_tol, mesh_tol, mesh_tol, mesh_tol}, maxlevel = maxlevel, minlevel = minlevel);
  #endif
}


/** 
### 2.4. Set initial conditions

We start from a flat surface perturbed by annular spectrum for the interface
perturbation as in [Dimonte et al. (2004)](#dimonte2004). Additionally, the
velocity field is set to zero.

*/

#include "view.h"
#if dimension == 2
  #include "../input_fields/initial_conditions_dimonte_fft1.h"
#else  
  #include "../input_fields/initial_conditions_dimonte_fft2.h"  
#endif

event init(i = 0){
  if (!restore(file = "./backup")){
    
    for (int l = 5; l <= maxlevel; ++l){
      refine( rectanglebox( 4.*L0/(1 << l) ) > 0 && level < l);
    }

    #if dimension == 2
    {
      vertex scalar phi[];
      initial_condition_dimonte_fft(phi, amplitude=perturb.a0, kmin = perturb.kmin, kmax = perturb.kmax, NX=(1<<maxlevel));
      foreach_vertex(){
        phi[] *= -1;
      }
      fractions(phi, f);
    }
    #else
    {
      vertex scalar phi[];
      initial_condition_dimonte_fft2(phi, amplitude=perturb.a0, kmin = perturb.kmin, kmax = perturb.kmax, NX=(1<<maxlevel), NY=(1<<maxlevel), isvof=1);
      foreach_vertex(){
        phi[] *= -1;
      }
      fractions(phi, f);      
    }
    #endif

    /** and zero velocity field  */
    foreach(){
      foreach_dimension(){
        u.x[] = 0.;
      }
    }    
  }

  /** If there are walls, we ensure the fields are zero insider the walls */
  foreach(){
    f[]*= (rectanglebox(1.*_mindel) > 0.);
    p[]*= (rectanglebox(1.*_mindel) > 0.);
    foreach_dimension(){
      u.x[]*= (rectanglebox(1.*_mindel) > 0.);
    }
  }
  N = 1 << maxlevel;
  solid(cs, fs, rectanglebox(0.));
  fractions_cleanup(cs, fs, smin=1e-3, opposite=1);
  restriction({fs, cs, f, u});
}

/** 
### 2.5. Set boundary conditions
We considered no-slip conditions on all boundaries, but stree-free conditions 
should give similar results.
*/

#if dimension == 2
  u.n[left] = dirichlet(0);
  u.t[left] = dirichlet(0);
  u.n[right] = dirichlet(0);
  u.t[right] = dirichlet(0);
  u.n[embed] = dirichlet(0);
  u.t[embed] = dirichlet(0);
#else
  u.n[left] = dirichlet(0);
  u.t[left] = dirichlet(0);
  u.r[left] = dirichlet(0);
  u.n[right] = dirichlet(0);
  u.t[right] = dirichlet(0);
  u.r[right] = dirichlet(0);
  u.n[embed] = dirichlet(0);
  u.t[embed] = dirichlet(0);
  u.r[embed] = dirichlet(0);
#endif

/** 
### 2.6. Apply the external forcing

#### 2.6.1. Functions to calculate a ramp value based on input parameters
We may use a sigmoid function
*/
double sigmoid(double t1, double k){
  return 1.0 / (1.0 + exp(-k * (t - t1)));
}

/**
which multiplies a periodic acceleration in the vertical direction

#### 2.6.2. Add the external forcing to the acceleration term
*/

event acceleration(i++){
  
  double t0 = 0.; 
  double t1 = 1.;
  double k = 4.0;
  
  force.ramp = sigmoid(t1-t0, k);
  force.Gn = (force.G0) * (force.F0 * sin(force.omega0 * t));
  
  face vector av = a;
  #if dimension == 2
    foreach_face(y)
      av.y[] -= force.G0 + force.Gn*force.ramp;
  #else
    foreach_face(z)
      av.z[] -= force.G0 + force.Gn*force.ramp;
  #endif 
}

/** 
## 3. Outputs

### 3.1. Logfile
We follow the residuals from the iterative solvers,

*/

event logfile(i++; t <= tend){
  fprintf(stderr, " res: %.5f \t %g %d %d \t %g %d %d \n", t, mgp.resa, mgp.i, mgp.nrelax, mgu.resa, mgu.i, mgu.nrelax);
}

/** 

### 3.2. Animation of the surface

To get the supersquare patterns, it is essential to utilize a large domain for
the simulations. Due to their substantial size, running these simulations on the
Basilisk server is impractical.

For example, a specific simulation was conducted on the Ruche cluster at the
Moulon Mesocenter, taking approximately 60 hours to complete. Detailed
information about the cluster can be found
[here](https://mesocentre.pages.centralesupelec.fr/).

The simulation, which spans 250 dimensionless time units (about 40 forcing
periods) is characterized by the following statistics:

```
# Octree, 5881 steps, 225809 CPU, 2.276e+05 real, 1.13e+06 points.step/s, 62 var
# 40 procs, MPI: min 2.8e+04 (12%) avg 3.6e+04 (16%) max 4.4e+04 (19%)
```

The development of the Faraday waves can be seen in the video below:

<iframe width="640" height="640" src="https://www.youtube.com/embed/4HFcfVq1Quk" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>

*/

event movie(t += pi/16){
#if dimension == 2
  view(camera = "front", ty = -0.1, tx=0.05);
  box();
  draw_vof("f", lc = {1, 0, 0}, lw = 2.);
  squares("cs[] > 0. ? f[] : nodata", linear = false);
#elif dimension == 3
  view(camera = "iso");
  scalar Z[];
  foreach()
    Z[] = z < ((14.5/11.40)*2*pi) ? f[] : 1.0;
  view(camera="iso");
  draw_vof("Z");
#endif
  save("front.mp4");
}

/**
### 3.3. Evolution of the kinetic energy

Our simulation tracks kinetic energy over time, revealing the initial
perturbation, exponential growth, and saturation to a finite amplitude.

~~~pythonplot Evolution kinetic energy
import numpy as np
import matplotlib.pyplot as plt
	
n = 64
fig, (ax1, ax2) = plt.subplots(nrows=2, sharex=True)
data = np.loadtxt('globals_from_ruche.asc', delimiter=' ', usecols=[0,1,2])
ax1.plot(data[:,0], data[:,1])
ax2.plot(data[:,0], data[:,2])
ax1.set_ylabel(r'External acceleration')
ax2.set_ylabel(r'Kinetic energy')
ax2.set_xlabel(r'Time')
ax1.set_xlim([0,250])
plt.tight_layout()
plt.savefig('plot_ek_vs_t.svg')
~~~ 
*/

event time_series(t += pi/16){
  double ekin = 0;
  foreach (reduction(+ : ekin)){
    foreach_dimension()	{
      ekin += dv() * 0.5 * rhov[] * sq(u.x[]);
    }
  }
  if (pid() == 0){
    FILE *fp = fopen("globals.asc", "a");
    fprintf(fp, "%f %.9g %.9g %.9g \n", t, force.Gn, ekin, force.ramp);
    fclose(fp);
  }
}

#ifdef _OUTPUTS
#include "../output_fields/vtu/output_vtu.h"
#include "curvature.h"
event save_facets(t += pi/32){
  scalar kappa[];
  curvature(f, kappa);
  char fname[99];
  sprintf(fname, "interface_t%.5f", t);
  output_facets_vtu(f, kappa, fname);
}
#endif

/**
# References

~~~bib

@article{tavares2024,
  title={A coupled VOF/embedded boundary method to model two-phase flows on arbitrary solid surfaces},
  author={Tavares, Mathilde and Josserand, Christophe and Limare, Alexandre and Lopez-Herrera, Jos{\'e} Ma and Popinet, St{\'e}phane},
  journal={Computers \& Fluids},
  volume={278},
  pages={106317},
  year={2024},
  doi={10.1016/j.compfluid.2024.106317},
  publisher={Elsevier}
}

@article{kahouadji2015,
  title={Numerical simulation of supersquare patterns in Faraday waves},
  author={Kahouadji, Lyes and P{\'e}rinet, Nicolas and Tuckerman, Laurette S and Shin, Seungwon and Chergui, Jalel and Juric, Damir},
  journal={Journal of Fluid Mechanics},
  volume={772},
  pages={R2},
  year={2015},
  publisher={Cambridge University Press},
  doi={10.1017/jfm.2015.213}
}

@article{dimonte2004,
  author = {Dimonte, Guy and Youngs, D. L. and Dimits, A. and Weber, S. and Marinak, M. and Wunsch, S. and Garasi, C. and Robinson, A. and Andrews, M. J. and Ramaprabhu, P. and Calder, A. C. and Fryxell, B. and Biello, J. and Dursi, L. and MacNeice, P. and Olson, K. and Ricker, P. and Rosner, R. and Timmes, F. and Tufo, H. and Young, Y.-N. and Zingale, M.},
  title = {A comparative study of the turbulent Rayleigh–Taylor instability using high-resolution three-dimensional numerical simulations: The Alpha-Group collaboration},
  journal = {Physics of Fluids},
  volume = {16},
  number = {5},
  pages = {1668-1693},
  year = {2004},
  month = {05},
  issn = {1070-6631},
  doi = {10.1063/1.1688328},
  url = {https://doi.org/10.1063/1.1688328},
  eprint = {https://pubs.aip.org/aip/pof/article-pdf/16/5/1668/19152852/1668\_1\_online.pdf},
}


~~~
*/
