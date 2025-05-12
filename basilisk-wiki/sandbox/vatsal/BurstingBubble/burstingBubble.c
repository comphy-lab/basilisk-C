/**
This repository contains the codes used for simulating the cases discussed in the manuscript: Bursting bubble in a viscoplastic medium. The results presented here are currently under review in Journal of Fluid Mechanics. The preprint of the article is available [here](https://arxiv.org/abs/2101.07744).

The supplementary videos are available [here](https://youtube.com/playlist?list=PLf5C5HCrvhLFETl6iaRr21pzr5Xab1OCM)

# Introduction:
We investigate the classical problem of bubble bursting at a liquid-gas interface, but now in the presence of a viscoplastic liquid medium. Here are the schematics of the problem. This simualtion will start from Figure 1(c). 

<p align="center">
  <img src="schematic.png" width="30%">
  <caption><p align="center">Schematics for the process of a bursting bubble: (a) A gas bubble in bulk. (b) The bubble approaches the free surface forming a liquid film (thickness $\delta$) between itself and the free surface. (c) A bubble cavity forms when the thin liquid film disappears.</caption>
</p>

## Some typical simulations
<div id="wrapper">
<iframe width="560" height="315" src="https://www.youtube-nocookie.com/embed/DY6TpfzZ2fM" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>
<iframe width="560" height="315" src="https://www.youtube-nocookie.com/embed/vYxNehfcKac" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>
<iframe width="560" height="315" src="https://www.youtube-nocookie.com/embed/PJbKsEjIbjQ" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>
<iframe width="560" height="315" src="https://www.youtube-nocookie.com/embed/aJbnD8wAPmQ" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>
<div class="center-block"></div> 
</div>

Bursting bubble dynamics for different capillary-Bingham numbers. (a) $\mathcal{J} = 0.0$: A typical case with a Newtonian liquid medium, (b) $\mathcal{J} =0.1$: A weakly viscoplastic liquid medium in which the process still shows all the major characteristics of the Newtonian liquid, (c) $\mathcal{J} = 0.5$: A case of moderate yield stress whereby the jetting is suppressed, nonetheless the entire cavity still yields, and (d) $\mathcal{J} = 1.0$: A highly viscoplastic liquid medium whereby a part of the cavity never yields. The left part of each video shows the magnitude of the velocity field, and the right part shows the magnitude of the deformation tensor on a $\log_{10}$ scale. The transition to the black region (low strain rates) marks the yield-surface location in the present study. For all the cases in this figure, $\mathcal{O}h = 10^{-2}$.

# Numerical code
Id 1 is for the Viscoplastic liquid pool, and Id 2 is Newtonian gas.
*/

#include "axi.h"
#include "navier-stokes/centered.h"
#define FILTERED // Smear density and viscosity jumps
/**
To model Viscoplastic liquids, we use a modified version of [two-phase.h](http://basilisk.fr/src/two-phase.h). [two-phaseVP.h](two-phaseVP.h) contains these modifications.
*/
#include "two-phaseAxiVP.h"
/**
 You can use: conserving.h as well. Even without it, I was still able to conserve the total energy (also momentum?) of the system if I adapt based on curvature and vorticity/deformation tensor norm (see the adapt even). I had to smear the density and viscosity anyhow because of the sharp ratios in liquid (Bingham) and the gas.
*/
#include "navier-stokes/conserving.h"
#include "tension.h"
#include "reduced.h"
#include "distance.h"
/**
We use a modified adapt-wavelet algorithm available [(here)](http://basilisk.fr/sandbox/pairetti/bag_mode/adapt_wavelet_limited.h). It is written by *César Pairetti* (Thanks :)). We use to ensure that refinement is higher near the substrate.
*/
#include "pairetti/bag_mode/adapt_wavelet_limited.h"

#define tsnap (0.001)
// Error tolerancs
#define fErr (1e-3)                                 // error tolerance in f1 VOF
#define KErr (1e-4)                                 // error tolerance in f2 VOF
#define VelErr (1e-2)                               // error tolerances in velocity
#define OmegaErr (1e-3)                             // error tolerances in vorticity

// Numbers!
#define RHO21 (1e-3)
#define MU21 (2e-2)
#define Ldomain 8

// boundary conditions
u.n[right] = neumann(0.);
p[right] = dirichlet(0.);

int MAXlevel;
double Oh, Bond, tmax;
char nameOut[80], dumpFile[80];

int  main(int argc, char const *argv[]) {
  L0 = Ldomain;
  origin (-L0/2., 0.);
  init_grid (1 << 6);
  // Values taken from the terminal
  MAXlevel = atoi(argv[1]);
  tauy = atof(argv[2]);
  Bond = atof(argv[3]);
  Oh = atof(argv[4]);
  tmax = atof(argv[5]);

  // Ensure that all the variables were transferred properly from the terminal or job script.
  if (argc < 6){
    fprintf(ferr, "Lack of command line arguments. Check! Need %d more arguments\n",6-argc);
    return 1;
  }
  fprintf(ferr, "Level %d, Oh %2.1e, Tauy %4.3f, Bo %4.3f\n", MAXlevel, Oh, tauy, Bond);

  // Create a folder named intermediate where all the simulation snapshots are stored.
  char comm[80];
  sprintf (comm, "mkdir -p intermediate");
  system(comm);
  // Name of the restart file. See writingFiles event.
  sprintf (dumpFile, "dump");

  /**
  We consider the burst of a small axisymmetric bubble at a surface of an incompressible Bingham fluid. To nondimensionalise the governing equations, we use the initial bubble radius $R_0$, and inertia-capillary velocity $V_\gamma = \sqrt{\gamma/(\rho_lR_0)}$, respectively. Pressure and stresses are scaled with the characteristic capillary pressure, $\gamma/R_0$. The dimensionless equations for mass and momentum conservation, for the liquid phase, then read

  $$
  \nabla \cdot \boldsymbol{u} = 0
  $$
  $$
  \frac{\partial\boldsymbol{u}}{\partial t} + \nabla\boldsymbol{\cdot}\left(\boldsymbol{uu}\right) = -\nabla p + \nabla\boldsymbol{\cdot}\boldsymbol{\tau} - \mathcal{B}o\,\hat{\boldsymbol{e}}_{\boldsymbol{\mathcal{Z}}},
  $$
  where $\boldsymbol{u}$ is the velocity vector, $t$ is time, $p$ is the pressure and $\boldsymbol{\tau}$ represents the deviatoric stress tensor. We use the regularized Bingham model with
  $$
  \boldsymbol{\tau} = 2\,\text{min}\left(\frac{\mathcal{J}}{2\|\boldsymbol{\mathcal{D}}\|} + \mathcal{O}h, \mathcal{O}h_\text{max}\right)\boldsymbol{\mathcal{D}}
  $$
  where $\|\boldsymbol{\mathcal{D}}\|$ is the second invariant of the deformation rate tensor, $\boldsymbol{\mathcal{D}}$, and $\mathcal{O}h_{max}$ is the viscous regularisation parameter. The three dimensionless numbers controlling the equations above are the capillary-Bingham number $\left(\mathcal{J}\right)$, which accounts for the competition between the capillary and yield stresses, the Ohnesorge number $\left(\mathcal{O}h\right)$ that compares the inertial-capillary to inertial-viscous time scales, and the Bond number $\left(\mathcal{B}o\right)$, which compares gravity and surface tension forces:
  $$
  \mathcal{J} = \frac{\tau_yR_0}{\gamma},\,\,\mathcal{O}h = \frac{\mu_l}{\sqrt{\rho_l\gamma R_0}},\,\,\mathcal{B}o = \frac{\rho_l gR_o^2}{\gamma}.
  $$

  Here, $\gamma$ is the liquid-gas surface tension coefficient, and $\tau_y$ and $\rho_l$ are the liquid's yield stress and density, respectively. Next, $\mu_l$ is the constant viscosity in the Bingham model. Note that in our simulations, we also solve the fluid's motion in the gas phase, using a similar set of equations (Newtonian). Hence, the further relevant non-dimensional groups in addition to those above are the ratios of density $\left(\rho_r = \rho_g/\rho_l\right)$ and viscosity $\left(\mu_r = \mu_g/\mu_l\right)$. In the present study, these ratios are kept fixed at $10^{-3}$ and $2 \times 10^{-2}$, respectively (see above). 
  */

  mumax = 1e8*Oh;  // The regularisation value of viscosity
  rho1 = 1., rho2 = RHO21;
  mu1 = Oh, mu2 = MU21*Oh;
  f.sigma = 1.0;
  G.x = -Bond;
  run();
}

/**
This event is specific to César's adapt_wavelet_limited.
*/
int refRegion(double x, double y, double z){
  return (y < 1.28 ? MAXlevel+2 : y < 2.56 ? MAXlevel+1 : y < 5.12 ? MAXlevel : MAXlevel-1);
}

/**
## Initial Condition

The initial shape of the bubble at the liquid-gas interface can be calculated by solving the Young-Laplace equations ([Lhuissier & Villermaux, 2012](https://doi.org/10.1017/jfm.2011.418​)) and it depends on the $\mathcal{B}o$ number. 
Resources: 

* [Alex Berny's Sandbox](http://www.basilisk.fr/sandbox/aberny/bubble/bubble.c)
* [My Matlab code:](https://github.com/VatsalSy/Bursting-Bubble-In-a-Viscoplastic-Medium/blob/main/InitialCondition.m) Also see the results for different $\mathcal{B}o$ number [here](https://youtu.be/Z_vdsOW5fsg).

<p align="center">
  <img src="VillermauxComparision.png" width="25%">
  <caption><p align="center">Comparision of the initial shape calculated using the Young-Laplace equations. Ofcourse, for this study, we only need: $\mathcal{B}o = 10^{-3}$. In the figure, $a$ is the capillary length, $a = \sqrt{\gamma/(\rho_lg)}$</caption>
</p>

Since we do not focus on the influence of $\mathcal{B}o$ number in the present study, I am not elaborating on it here. For all the simulations, I use the interfacial shape calculated for $\mathcal{B}o = 10^{-3}$. The resultant data file is available [here](https://raw.githubusercontent.com/VatsalSy/Bursting-Bubble-In-a-Viscoplastic-Medium/main/Bo0.0010.dat).

**Note:** The curvature diverges at the cavity-free surface intersection. We fillet this corner to circumvent this singularity, introducing a rim with a finite curvature that connects the bubble to the free surface. We ensured that the curvature of the rim is high enough such that the subsequent dynamics are independent of its finite value.
*/
event init (t = 0) {
  if (!restore (file = dumpFile)){

    char filename[60];
    sprintf(filename,"Bo%5.4f.dat",Bond);
    FILE * fp = fopen(filename,"rb");
    if (fp == NULL){
      fprintf(ferr, "There is no file named %s\n", filename);
      return 1;
    }
    coord* InitialShape;
    InitialShape = input_xy(fp);
    fclose (fp);
    scalar d[];
    distance (d, InitialShape);
    while (adapt_wavelet_limited ((scalar *){f, d}, (double[]){1e-8, 1e-8}, refRegion).nf);
    /**
    The distance function is defined at the center of each cell, we have
    to calculate the value of this function at each vertex. */
    vertex scalar phi[];
    foreach_vertex(){
      phi[] = -(d[] + d[-1] + d[0,-1] + d[-1,-1])/4.;
    }
    /**
    We can now initialize the volume fraction of the domain. */
    fractions (phi, f);
  }
}

/**
## Adaptive Mesh Refinement
*/
event adapt(i++){
  /**
  We adapt based on curvature, $\kappa$ and vorticity $\omega$. 
  Adaptation based on $\kappa$ ensures a constant grid resolution across the interface. See [this](http://basilisk.fr/sandbox/Antoonvh/rc.c) for further reading. 

  We also adapt based on vorticity in the liquid domain. I have noticed that this refinement helps resolve the fake-yield surface accurately (see the black regions in the videos below). 
  */
  scalar KAPPA[], omega[];
  curvature(f, KAPPA);
  vorticity (u, omega);
  foreach(){
    omega[] *= f[];
  }
  boundary ((scalar *){KAPPA, omega});
  adapt_wavelet_limited ((scalar *){f, u.x, u.y, KAPPA, omega},
     (double[]){fErr, VelErr, VelErr, KErr, OmegaErr},
     refRegion);
  /**
  ## Alternatively
  At higher $\mathcal{O}h$ and $\mathcal{J}$ numbers, vorticities in the liquid cease to be interesting. In that case, one might want to adapt based on the norm of deformation tensor, $\mathbf{\mathcal{D}}$. I already calculate $\|\mathbf{\mathcal{D}}\|$ in [two-phaseAxiVP.h](two-phaseAxiVP.h).

  **Note:** $\mathbf{\mathcal{D}}$ based refinement is way more expensive than $\omega$ based refinement.
  */
  // adapt_wavelet_limited ((scalar *){f, u.x, u.y, KAPPA, D2},
  //    (double[]){fErr, VelErr, VelErr, KErr, 1e-3},
  //    refRegion);
}

/**
## Dumping snapshots
*/
event writingFiles (t = 0; t += tsnap; t <= tmax) {
  dump (file = dumpFile);
  sprintf (nameOut, "intermediate/snapshot-%5.4f", t);
  dump(file=nameOut);
}

/**
## Ending Simulation
*/
event end (t = end) {
  fprintf(ferr, "Done: Level %d, Oh %2.1e, Tauy %4.3f, Bo %4.3f\n", MAXlevel, Oh, tauy, Bond);
}

/**
## Log writing
*/
event logWriting (i+=100) {
  double ke = 0.;
  foreach (reduction(+:ke)){
    ke += (2*pi*y)*(0.5*(f[])*(sq(u.x[]) + sq(u.y[])))*sq(Delta);
  }
  static FILE * fp;
  if (i == 0) {
    fprintf (ferr, "i dt t ke\n");
    fp = fopen ("log", "w");
    fprintf (fp, "i dt t ke\n");
    fprintf (fp, "%d %g %g %g\n", i, dt, t, ke);
    fclose(fp);
  } else {
    fp = fopen ("log", "a");
    fprintf (fp, "%d %g %g %g\n", i, dt, t, ke);
    fclose(fp);
  }
  fprintf (ferr, "%d %g %g %g\n", i, dt, t, ke);
  if (ke > 1e3 || ke < 1e-6){
    if (i > 1e2){
      return 1;
    }
  }
}

/**
## Running the code
~~~bash
#!/bin/bash
qcc -fopenmp -Wall -O2 burstingBubble.c -o burstingBubble -lm
export OMP_NUM_THREADS=8
./burstingBubble 10 0.25 1e-3 1e-2 5.0
~~~

# Output and Results
The post-processing codes and simulation data are available at: [PostProcess](https://github.com/VatsalSy/Bursting-Bubble-In-a-Viscoplastic-Medium/tree/main/PostProcess)
*/
