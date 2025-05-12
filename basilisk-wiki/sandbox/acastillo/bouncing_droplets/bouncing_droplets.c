/** 
# 1. Bouncing droplets impacting a fluid bath

This code is adapted from the set-up used in "*Inertio-capillary rebound of a
droplet impacting a fluid bath*" by [Alventosa, Cimpeanu & Harris
(2023)](#Alventosa2023). Here, the original implementation of
[two-phaseDOD.h](https://github.com/VatsalSy/Lifting-a-sessile-drop/blob/master/CaseI/two-phaseDOD.h)
is replaced in favor of using [no-coalescence.h](src/no-coalescence.h).
The original code posted by [R.
Cimpeanu](http://basilisk.fr/sandbox/rcimpeanu/README) can be found in their
[github.](https://github.com/rcsc-group/BouncingDroplets/) 

<center>
<table>
<tr>
<td><center>![](Alventosa2023a.png){ width="60%" }</center></td>
<td><center>![](Alventosa2023b.png){ width="100%" }</center></td>
</tr>
</table>
<td><center>A small water droplet (R ≈ 0.4 mm) rebounds from a bath of the
same fluid and a schematic of the problem. Both images are taken from 
[Alventosa, Cimpeanu & Harris (2023)](#Alventosa2023)</center></td>
</center>
<br>

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
which can be solved using (mostly) existing Basilisk code.

*/

/**
# 2. Implementation in Basilisk

## 2.1. Include solver blocks
*/

#include "axi.h"                     // Axisymmetric geometry
#include "navier-stokes/centered.h"  // Solve Navier-Stokes equations
#define FILTERED                     // Smear density and viscosity jumps
#define mu(f)  (1./(clamp(f,0,1)*(1./mu1 - 1./mu2) + 1./mu2))
#include "two-phase.h"               // Two-phase VoF approach
#include "no-coalescence.h"          // Prevent numerical coalescence
#include "tension.h"                 // Include surface tension between phases
                    
/** 
## 2.2. Set problem parameters
*/

/**   
We setup a structure to hold dimensional problem parameters. The key dimensional
parameters are:

| Parameter                   | Symbol           | Value       |
|-----------------------------|------------------|-------------|
| Impact speed                | $V_0$            | 20-100 cm/s |
| Droplet radius              | $R_d$            | 0.035 cm    |
| Density (water)             | $\rho_\ell$      | 0.998 g/cm³ |
| Kinematic viscosity (water) | $\nu_\ell$       | 0.978 cSt   |
| Surface tension (water)     | $\sigma$         | 72.2 mN/m   |
| Gravitational acceleration  | $g$              | 9.81 m/s²   |

*/
struct DProblemParameters
{
  double rho_l;   // Density, fluid 1 (kg/m^3)
  double rho_g;   // Density, fluid 2 (kg/m^3)
  double mu_l;    // Dynamic viscosity, fluid 1 (Pa.s)
  double mu_g;    // Dynamic viscosity, fluid 2 (Pa.s)
  double sigma;   // Surface tension (N/m)
  double g;       // Gravitational acceleration (m/s^2)
  double Rd;      // Drop radius (m)
  double Vo;      // Initial drop velocity (m/s)
};
struct DProblemParameters Dim = {0};

/**   
and a structure to hold non-dimensional problem parameters. The key
nondimensional parameters are:

| Parameter  | Symbol           | Definition                              |
|------------|------------------|-----------------------------------------|
| Weber      | $\mathrm{We}$    | $\rho_\ell V_0^2 R_d/\sigma$            |
| Bond       | $\mathrm{Bo}$    | $\rho_\ell g R_d^2/\sigma$              |
| Ohnesorge  | $\mathrm{Oh}$    | $\mu_\ell /\sqrt{\sigma R_d \rho_\ell}$ |
| Reynolds   | $\mathrm{Re}$    | $\rho_\ell V_0 R_d/\mu_\ell$            |
| Density ratio   | $\rho^*$    | $\rho_g/\rho_\ell$                      |
| Viscosity ratio | $\mu^*$     | $\mu_g/\mu_\ell$                        |

*/
struct NDProblemParameters
{
  double We;    // Weber number
  double Re;    // Reynolds number
  double Fr;    // Froude number
  double Bo;    // Bond number 
  double Oh;    // Ohnesorge
  double rho;   // Density ratio
  double mu;    // Viscosity ratio
  double H0;    // Pool height (in radii)
  double h0;    // Film height (in radii)
  double L0;    // Computational box size (in radii)
};
struct NDProblemParameters ND = {0};

/**   

*/

int minLevel = 5;
int maxLevel; // = 11;
double tend;

int main() {

  /**   
  First, we get the dimensional parameters
  */ 
  Dim.rho_l = 998.0;      // Density, fluid 1 (kg/m^3)
  Dim.rho_g = 1.21;       // Density, fluid 2 (kg/m^3)
  Dim.mu_l = 0.998e-3;    // Dynamic viscosity, fluid 1 (Pa.s)
  Dim.mu_g = 1.81e-5;     // Dynamic viscosity, fluid 2 (Pa.s)
  Dim.sigma = 0.0722;     // Surface tension (N/m)
  Dim.g = 9.81;           // Gravitational acceleration (m/s^2)
  Dim.Rd = 0.35e-3;       // Drop radius (m)
  Dim.Vo = 0.20;          // Initial drop velocity (m/s)
  tend = 6.0;             // prescribed simulation end time
  maxLevel = 9;           // prescribed maximum resolution level
  
  /**   
  and compute the dimensionless numbers
  */ 
  ND.We = (Dim.rho_l * pow(Dim.Vo, 2.0) * Dim.Rd) / Dim.sigma;
  ND.Re = (Dim.rho_l * Dim.Vo * Dim.Rd) / Dim.mu_l;
  ND.Fr = Dim.Vo / pow(Dim.Rd * Dim.g, 0.5);
  ND.Bo = Dim.rho_l * Dim.g * pow(Dim.Rd, 2.0) / Dim.sigma;
  ND.Oh = Dim.mu_l / pow(Dim.rho_l * Dim.sigma * Dim.Rd, 0.5);
  ND.rho = Dim.rho_g / Dim.rho_l;
  ND.mu = Dim.mu_g / Dim.mu_l;

  ND.H0 = 4.0;
  ND.L0 = 8.0;
  ND.h0 = -ND.L0/2. + ND.H0;

  fputs("SIMULATION PARAMETERS\n", stdout);
  fprintf(stdout, " Reynolds number   : %0.6f \n", ND.Re); 
  fprintf(stdout, " Weber number      : %0.6f \n", ND.We); 
  fprintf(stdout, " Froude number     : %0.6f \n", ND.Fr); 
  fprintf(stdout, " Bond number       : %0.6f \n", ND.Bo); 
  fprintf(stdout, " Ohnesorge number  : %0.6f \n", ND.Oh); 
  fprintf(stdout, " Density ratio     : %0.6f \n", ND.rho); 
  fprintf(stdout, " Viscosity ratio   : %0.6f \n", ND.mu); 
  fprintf(stdout, " Pool height       : %0.6f \n", ND.H0); 
  fprintf(stdout, " Film height       : %0.6f \n", ND.h0); 
  fprintf(stdout, " Domain size       : %0.6f \n", ND.L0); 
  
  /** Then, assign the fluid properties, */ 
  rho1 = 1.;          // for f = 1
  rho2 = ND.rho;      // for f = 0  
  mu1 = 1./ND.Re;     // for f = 1
  mu2 = ND.mu/ND.Re;  // for f = 0  
  f.sigma = 1./ND.We;

  /** the geometric parameters, */ 
  init_grid(1 << 8);
  size(ND.L0);                     
  origin(-0.5*ND.L0, 0.0);

  /** Now, we set-up the numerical parameters and run */ 
  DT = 1e-2;
  NITERMIN = 1; // default 1
  NITERMAX = 200; // default 100
  TOLERANCE = 1e-6; // default 1e-3
  run();
}

/** 
## 2.3. Set grid refinement
*/
event adapt(i++){
  adapt_wavelet ((scalar *){f, u}, (double[]){1e-5, 1e-3, 1e-3}, maxLevel, minLevel);
}

/** 
## 2.4. Set initial conditions
*/
#define CIRCLE(H,R,a,b) (sq(R*(1.0+a)) - sq(x - (H+R*(1.0+b))) - sq(y))
event init (t = 0.0) {

  // Strong refinement around the interfacial regions
  refine( 
    union(
      intersection(
        (CIRCLE( ND.h0, 1.0, 0.05, 0.5) > 0), 
        (CIRCLE( ND.h0, 1.0,-0.05, 0.5) < 0)
      ), 
      intersection(
        (-x + ND.h0 + 0.005 > 0.0), 
        (-x + ND.h0 - 0.005 < 0.0)
      )
    ) && level < maxLevel);
  
  scalar f1[], f2[];
  // Create active liquid phase as union between drop and film
  fraction (f1, CIRCLE( ND.h0, 1.0, 0.0, 0.5) );
  fraction (f2, - x + ND.h0);

  foreach()
    f[] = f1[] + f2[];
  
  // Initialise uniform velocity field inside droplet
  foreach(){
    u.x[] = -1.0*f1[];
    u.y[] = 0.0;
    p[] = 0.0;
  }  
}
#undef CIRCLE

/** 
## 2.5. Set boundary conditions

When using axisymmetry, `x` indicates the axial direction, while `y` indicates
the radial direction. We set the boundary conditions accordingly.

- `left` $\rightarrow$ Pool bottom : no-slip, no permeability
- `top` $\rightarrow$ Pool side : no-slip, no permeability
- `right` $\rightarrow$ Above the pool : outflow condition
- `bottom`$\rightarrow$  Rotation axis : symmetry

*/

u.n[left] = dirichlet(0.);
u.t[left] = dirichlet(0.);
p[left] = neumann(0.);
pf[left] = neumann(0.);

u.n[top] = dirichlet(0.); 
u.t[top] = neumann(0.); 

u.n[right] = neumann(0.);
p[right] = dirichlet(0.);
pf[right] = dirichlet(0.);

/** 
## 2.6. Apply the external forcing
*/
event acceleration (i++) {
  face vector av = a;
  foreach_face(x)  
    av.x[] -= 1./pow(ND.Fr,2.0);
  foreach_face(y)  
    av.y[] += 0.0;
}

/** 
# 3. Outputs
## 3.1. Logging residuals and number of iterations from each solver
*/

event logfile(i++; t <= tend){
  fprintf(stderr, " residuals: %.5f \t %g %d %d \t %g %d %d \n", t, mgp.resa, mgp.i, mgp.nrelax, mgu.resa, mgu.i, mgu.nrelax);
}

/** 
## 3.2. Saving slices in VTKHDF format
*/

#include "view.h"
#include "../output_fields/vtkhdf/output_vtkhdf.h"
event save_snapshots(t += 0.1) {
  char filename[99];	
  sprintf(filename, "snapshot_t%0.1f.vtkhdf", t);
  output_vtkhdf({f, rhov, p}, {u}, filename);
}

/** 
## 3.3. Saving interfaces in VTU format
*/

#include "../output_fields/vtu/output_vtu.h"
#include "curvature.h"
event save_facets(t += 0.1){
  scalar kappa[];
  char fname[99];

  for (scalar s in interfaces){
    curvature(s, kappa);
    sprintf(fname, "interface_%s_t%0.1f", s.name, t);
    output_facets_vtu(s, kappa, fname);
  }
}

/** 
## 3.4. Counting droplets and storing position
*/

#include "tag.h"
event droplets(t += 0.01){
  scalar m[];
  for (scalar s in interfaces){
    foreach ()
      m[] = s[] > 1e-3;
    int n = tag(m);

    double v[n], ymin[n], ymax[n];
    coord dpos[n];
    coord dvel[n];
    for (int j = 0; j < n; j++){
      v[j] = 0.;
      ymin[j] =  1e100;
      ymax[j] = -1e100;
      dpos[j].x = dpos[j].y = dpos[j].z = 0.;
      dvel[j].x = dvel[j].y = dvel[j].z = 0.;
    }
    foreach (serial){
      if (m[] > 0){
        int j = m[] - 1;
        v[j] += dv() * s[];
        coord p = {x, y, z};
        if (p.y > ymax[j]) ymax[j] = p.y;
        if (p.y < ymin[j]) ymin[j] = p.y;
        foreach_dimension(){
          dpos[j].x += dv() * s[] * p.x;
          dvel[j].x += dv() * s[] * u.x[];
        }
      }
    }

#if _MPI
    MPI_Allreduce(MPI_IN_PLACE, v, n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, ymax, n, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, ymin, n, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, dpos, 3 * n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, dvel, 3 * n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif

    if (pid() == 0){
      FILE *fp = fopen("droplets.asc", "a");
      for (int j = 0; j < n; j++)
        fprintf(fp, "%s %d %g %d %g %g %g %g %g %g %g\n", s.name, i, t, j, v[j],
                  dpos[j].x / v[j], dpos[j].y / v[j],
                  dvel[j].x / v[j], dvel[j].y / v[j],
                  ymax[j], ymin[j]);
      fclose(fp);
      fflush(fp);
    }
  }
}


/** 
## 3.5. Saving animations

![Interfaces/Grid refinement](bouncing_droplets/movie_summary.mp4)(width="75%")

*/

#include "draw.h"
event movies (t += 0.05){

  // Define auxiliar fields
  char timestring[100];
  scalar omega[], l[];
    foreach(){
    omega[] = (u.y[1,0] - u.y[-1,0])/(2.*Delta) - (u.x[0,1] - u.x[0,-1])/(2.*Delta);    
    l[] = level;
  }

  // Movie 1
  view(width=1900, height=1050, fov=7.0, ty = 0.0, quat = { 0, 0, -0.707, 0.707 });    
  {
    {
      int index = 0;
      float * colors[] = {(float[]){1,0,0},(float[]){0,0,1}};
      for (scalar s in interfaces)
        draw_vof (s.name, fc = colors[index++], filled = 1, min = 0, max = 3, lw = 2.);
    }
    mirror({0,1}) {
      for (scalar s in interfaces)
        draw_vof(s.name, lw=2);
      cells(lw=0.5);
      squares("l", map = cool_warm, min = minLevel, max = maxLevel);
    }
    sprintf(timestring, "t=%2.03f",t);
    draw_string(timestring, pos=1, lc= { 0, 0, 0 }, lw=2);  
    save ("movie_summary.mp4");
  }

  // Movie 2
  {
    for (scalar s in interfaces)
      draw_vof(s.name, lw=2);
    squares("u.x", map = cool_warm, min = -1., max = 0.5,);    
    mirror({0,1}) {
      for (scalar s in interfaces)
        draw_vof(s.name, lw=2);
      squares("u.y", map = cool_warm, min = -0.5, max = 2.);
    }
    sprintf(timestring, "t=%2.03f",t);
    draw_string(timestring, pos=1, lc= { 0, 0, 0 }, lw=2);
    save ("movie_velocities.mp4");
  } 

  // Movie 3
  { 
    for (scalar s in interfaces)
      draw_vof(s.name, lw=2);
    squares("omega", map = cool_warm, min = -3., max = 3.);
    mirror({0,1}) {
      for (scalar s in interfaces)
        draw_vof(s.name, lw=2);
      squares("p", map = cool_warm, min = -0.25, max = 4.0);
    } 
    sprintf(timestring, "t=%2.03f",t);
    draw_string(timestring, pos=1, lc= { 0, 0, 0 }, lw=2);    
    save ("movie_vorticity.mp4");
  }
}

event finalize(t = end)
  dump("backup");

event backups(t += 1.0)
  dump("backup");

/**
# References

~~~bib

@article{Alventosa2023, 
  title={Inertio-capillary rebound of a droplet impacting a fluid bath}, 
  volume={958}, 
  DOI={10.1017/jfm.2023.88},
  journal={Journal of Fluid Mechanics}, 
  author={Alventosa, Luke F.L. and Cimpeanu, Radu and Harris, Daniel M.}, 
  year={2023}, 
  pages={A24}
}

~~~
*/
