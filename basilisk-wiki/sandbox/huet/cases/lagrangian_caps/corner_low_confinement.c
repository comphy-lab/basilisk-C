/**
# Motion of a train of capsules through a corner at $Re = 25$

![Train of initially spherical capsules flowing through a corner channel](https://damienhuet.github.io/images/basilisk_sandbox_files/corner_post_processing/un.mp4)(width="30%")

## Relevant solvers and macros
We first define some macros relevant to the physical parameters:

* the aspect ratio $\beta = w/2a$ of the capsule diameter $2a$ over the channel width $b$
* the channel width $w$
* the initial radius $a$ of the capsule, set to unity
* the length of the domain $L_0$
* the characteristic velocity $U_{avg}$ of the flow, which we will set as the inlet velocity
* the density $\rho$, set to unity
* the channel Reynolds number $Re_w = \rho U_{avg} w/\mu$, set to 25 in this simulation
* the viscosity \mu
* the Capillary number $Ca = \mu U_{avg} a/E_s$, set to $\approx 0.12$
* the elastic modulus $E_s$ of the membrane
* the reduced bending coefficient $E_b^\star = E_b/(E_s a^2)$, set to 0.005
* the bending coefficient $E_b$
* we set the ```STOKES``` boolean to ```false``` since we consider $Re=25$
* we tell the capsule solver that up to 100 capsules will exist in this simulation
* we set the non-dimensional initial gap between the two capsules, ```ND_INITIAL_GAP```, which will be multiplied by the capsules' initial
diameter.
*/

#define BETA_RATIO 6.
#define WIDTH 6.
#define RADIUS (WIDTH/BETA_RATIO)
#define L0 (WIDTH*20 + 1.e-10)
#define U_AVG 1.
#define RHO 1.
#define REw 25.
#define MU (RHO*U_AVG*WIDTH/REw)
#define CA (0.35/3.)
#define E_S (MU*U_AVG/CA)
#define ND_EB 0.005
#define E_B (ND_EB*E_S*sq(RADIUS))
#define TAU_VISCOUS (WIDTH*REw/(4*U_AVG)/10.)
#define STOKES false
#define NCAPS 100
#define INITIAL_POSITION (-30*RADIUS - WIDTH/2.)
#define ND_INITIAL_GAP (.25)


/**
Then, some macros relevant to the solvers:

* the simulation time $t_{end}$
* the maximum allowed time step $\Delta t_{max}$
* the tolerance for the wavelet adaptivity solver
* the minimum and maximum levels of the Eulerian mesh
* the level of the Lagrangian triangulation (```LAGLEVEL``` 0 is an
icosahedron, and the triangulation approaches a sphere as ```LAGLEVEL```
increases)
* the tolerance for the Poisson solver
* the Jacobian preconditioner for the viscous solver
* the frequency at which pictures and restart files are generated
* a boolean indicating whether the capsule configuration is outputted in
paraview in addition to bview pictures.
*/

#ifndef T_END
  #define T_END (TAU_VISCOUS + 80.)
#endif
#define DT_MAX (7.5e-4)
#define U_TOL (0.05*U_AVG)
#define MINLEVEL 2
#define MAXLEVEL 11
#define LAGLEVEL 4
#define MY_TOLERANCE (1.e-3*U_AVG)
#define JACOBI 1
#define OUTPUT_FREQ ( t += 0.1 )
#define DUMP_OUTPUT_FREQ ( t += 1. )
#define PARAVIEW_CAPSULE 1
#define PARAVIEW_FLOW_FIELD 0
#define RESTART_CASE 0
#define RESTORE_FRAME -1


/**
Before introducing the capsule, we let the initial flow field develop, as
indicated by the macro ```INITIAL_FLOW_FIELD```. We redefine the end time, the
time step, the output frequency and we force the number of capsules to be zero.
*/
#ifndef INITIAL_FLOW_FIELD
  #define INITIAL_FLOW_FIELD 0
#endif
#if INITIAL_FLOW_FIELD
  #undef T_END
  #define T_END (TAU_VISCOUS)
  #undef DT_MAX
  #define DT_MAX (T_END/1000.)
  #undef OUTPUT_FREQ
  #define OUTPUT_FREQ ( i += 10 )
  #define NCAPS 0
#endif
#define STR(s) #s
#define STRING(s) STR(s)


/**
We call the solvers and various header files needed to perform this simulation:

* the octree grid
* the embedded boundaries
* the centered Navier-Stokes solver
* the front-tracking solver for immersed capsules
* the finite element solver for elastic membranes (using the neo-Hookean law)
* the paraboloid-fitting solver for bending stresses on membranes
* a header file containing routines to define common initial capsule shapes
* the viewing functions in their front-tracking flavour
*/

#include "grid/octree.h"
#include "embed.h"
#include "navier-stokes/centered.h"
#include "lagrangian_caps/capsule-ft.h"
#include "lagrangian_caps/neo-hookean-ft.h"
#include "lagrangian_caps/bending-ft.h"
#include "lagrangian_caps/common-shapes-ft.h"
#include "lagrangian_caps/view-ft.h"
#include "lagrangian_caps/vtu_output.h"

/**
## Simulation setup

Because embedded boundaries are used, we have to multiply the density,
viscosity and 1/density fields by the fluid volume or face fraction.

We also declare file pointers to output performance and flow convergence data.
*/
scalar rhov[];
face vector muv[];
face vector alphav[];

FILE* fperf = NULL;
FILE* iff_stats = NULL;
int nb_pic;
/**
### Generalities
In the main function, we initialize the grid and set basic parameters to their
values defined above, before starting the time iterations.
*/
int main(int argc, char* argv[]) {
  origin(-.25*L0 + 1.e-11, -.75*L0 + 1.e-11, -.5*L0 + 1.e-11);
  N = 1 << MINLEVEL;
  init_grid(N);
  stokes = STOKES;
  DT = DT_MAX;
  TOLERANCE = MY_TOLERANCE;
  TOLERANCE_MU = MY_TOLERANCE;
  rho = rhov;
  mu = muv;
  alpha = alphav;
  run();
}

/**
### Boundary conditions
We set no slip boundary conditions everywhere except at the inlet and outlet,
where in the former case a constant (normal) velocity profile is imposed, and
in the latter case Neumann boundary conditions are set.
*/

u.n[left] = dirichlet(0.);
u.t[left] = dirichlet(0.);
u.r[left] = dirichlet(0.);
u.n[right] = neumann(0.);
u.t[right] = dirichlet(0.);
u.r[right] = dirichlet(0.);
u.n[top] = dirichlet(0.);
u.n[bottom] = dirichlet(((fabs(x) < WIDTH/2 && fabs(z) < WIDTH/2) ? U_AVG : 0.));
u.t[top] = dirichlet(0.);
u.t[bottom] = dirichlet(0.);
u.r[top] = dirichlet(0.);
u.r[bottom] = dirichlet(0.);
u.n[front] = dirichlet(0.);
u.n[back] = dirichlet(0.);
u.t[front] = dirichlet(0.);
u.t[back] = dirichlet(0.);
u.r[front] = dirichlet(0.);
u.r[back] = dirichlet(0.);
u.n[embed] = dirichlet(0.);
u.t[embed] = dirichlet(0.);
u.r[embed] = dirichlet(0.);
p[bottom] = neumann(0.);
pf[bottom] = neumann(0.);
p[right] = dirichlet(0.);
pf[right] = dirichlet(0.);

/**
### Initialization
For performance monitoring, we output some quantities related to the
number of grid cells and the Poisson solver. Following Arthur Ghigo's sandbox,
we choose a third order interpolation stencil for the velocity field.


We then define a spherical membrane. Note that the centroid of the capsule needs
to be the at origin at the end of the ```init``` event, in order to define
outward-pointing normal vectors. Moreover, the elastic stresses are assumed to
be zero at the end of this event. Any initial stress or initial shift needs to
take place in the following ```init_adapt``` event.
*/

#define vertical_channel (intersection(intersection(WIDTH/2 - fabs(x),\
  WIDTH/2 - fabs(z)), WIDTH/2 - y))
#define horizontal_channel (intersection(intersection(WIDTH/2 - fabs(y),\
  WIDTH/2 - fabs(z)), x - WIDTH/2))

scalar mcs[];
vector prev_u[];
int last_introduced_caps;
event init (i = 0) {
  fperf = fopen("perf.csv", "w");
  fprintf(fperf, "total_nb_cells, nb_iter_viscous, resb_viscous, resa_viscous, nb_relax_viscous, nb_iter_pressure, resb_pressure, resa_pressure, nb_relax_pressure\n");

  for (scalar s in {u})
    s.third = true;

  #if INITIAL_FLOW_FIELD
    iff_stats = fopen("iff_stats.csv", "w");
    fprintf(iff_stats, "time iteration number_of_cells avg_velocity_change max_velocity_change max_x_velocity outlet_flow_rate avg_output_vel\n");

/**
Then the geometry of the embedded boundaries
is generated. The considered geometry is a square duct which centerline forms
a corner at the origin. This geometry is the union of a vertical and a horizontal channels as defined below. Initially mesh is adapted around the
channel walls. When we generate the initial flow field, we obviously don't
create and initialize the capsule.
*/
    astats ss;
    int ic = 0;
    do {
      ic++;
      solid (cs, fs, union(vertical_channel, horizontal_channel));
      foreach() mcs[] = cs[];
      ss = adapt_wavelet ({mcs}, (double[]) {1.e-30}, maxlevel = MAXLEVEL,
        minlevel = MINLEVEL);
    } while ((ss.nf || ss.nc) && ic < 100);
  #else
    char dump_name[64];
    #if RESTART_CASE
      sprintf(dump_name, "flow_%s.dump", STRING(RESTORE_FRAME));
      restore(dump_name);
      sprintf(dump_name, "caps_%s.dump", STRING(RESTORE_FRAME));
      #if RESTORE_OLD_DUMP
        tweaked_restore_capsules(dump_name);
      #else
        restore_capsules(dump_name);
      #endif
      last_introduced_caps = NCAPS - 1;
    #else
      sprintf(dump_name, "initial_flow_field_re%s_lvl%s.dump", STRING(REw),
        STRING(MAXLEVEL));
      restore(dump_name);
      activate_spherical_capsule(&CAPS(0), level = LAGLEVEL, radius = RADIUS);
      last_introduced_caps = 0;
      for(int k=0; k<CAPS(0).nln; k++) CAPS(0).nodes[k].pos.y += INITIAL_POSITION;
      generate_lag_stencils(no_warning = true);
      comp_centroid(&CAPS(0));
    #endif
      /**
  With the membrane properly created to an initially strain-free spherical shape,
  we now proceed to refining the Eulerian mesh around the membrane.
      */
    astats ss;
    int ic = 0;
    do {
      ic++;
      tag_ibm_stencils();
      solid (cs, fs, union(vertical_channel, horizontal_channel));
      foreach() mcs[] = cs[];
      ss = adapt_wavelet ({mcs, stencils, u}, (double[]) {1.e-30, 1.e-30,
        U_TOL, U_TOL, U_TOL}, maxlevel = MAXLEVEL, minlevel = MINLEVEL);
      generate_lag_stencils(no_warning = true);
    } while ((ss.nf || ss.nc) && ic < 100);

  /**
  We also erase the content of the files where the membrane position and triangle
  areas will be stored, and store the initial configuration.
  */
    if (pid() == 0) {
      for(int k=0; k<NCAPS; k++) {
        char fposname[64];
        char ftriname[64];
        sprintf(fposname, "mb_%d_pos.csv", k);
        sprintf(ftriname, "mb_%d_tri.csv", k);
        #if !RESTART_CASE
          FILE* file = fopen(fposname, "w");
          fclose(file);
          file = fopen(ftriname, "w");
          fclose(file);
        #endif
        if (CAPS(k).isactive) {
          dump_plain_nodes_pos(&CAPS(k), fposname);
          dump_plain_triangles(&CAPS(k), ftriname);
        }
      }
    }
    #if PARAVIEW_CAPSULE
      pv_output_ascii();
    #endif
    #if PARAVIEW_FLOW_FIELD
      scalar * list = {p};
      vector * vlist = {u};
      save_data(list, vlist, t);
    #endif
  #endif
  nb_pic = RESTORE_FRAME + 1;
}

/**
In this event, we introduce a new capsule in the simulation as soon as the
previous capsule has advanced by $(1 + G^\star)D$, where $D$ is the initial
diameter of the capsules, and $G^\star > 0$ is the non-dimensional initial gap
between the new capsule and the previous one.
*/
event manage_capsules_active_status (i++) {
  if (last_introduced_caps + 1 < NCAPS) {
    if (CAPS(last_introduced_caps + 1).isactive == false &&
      (CAPS(last_introduced_caps).centroid.y > INITIAL_POSITION +
      2*RADIUS*(1 + ND_INITIAL_GAP))) {
      activate_spherical_capsule(&CAPS(last_introduced_caps + 1),
        level = LAGLEVEL, radius = RADIUS);
      for(int k=0; k<CAPS(last_introduced_caps + 1).nln; k++)
        CAPS(last_introduced_caps + 1).nodes[k].pos.y += INITIAL_POSITION;
      comp_centroid(&CAPS(last_introduced_caps + 1));
      generate_lag_stencils(no_warning = true);
      astats ss;
      int ic = 0;
      do {
        ic++;
        tag_ibm_stencils();
        solid (cs, fs, union(vertical_channel, horizontal_channel));
        foreach() mcs[] = cs[];
        ss = adapt_wavelet ({mcs, stencils, u}, (double[]) {1.e-30, 1.e-30,
          U_TOL, U_TOL, U_TOL}, maxlevel = MAXLEVEL, minlevel = MINLEVEL);
        generate_lag_stencils(no_warning = true);
      } while ((ss.nf || ss.nc) && ic < 100);
      last_introduced_caps++;
    }
  }
  for(int i=0; i<NCAPS; i++) {
    if (CAPS(i).isactive) {
      double max_x = -HUGE;
      for(int j=0; j<CAPS(i).nln; j++)
        if (CAPS(i).nodes[j].pos.x > max_x) max_x = CAPS(i).nodes[j].pos.x;
      if (max_x >= .75*L0 - RADIUS) CAPS(i).isactive = false;
    }
  }
}

/**
### Updating properties and adapting the mesh
We ensure that the following density, viscosity and 1/density fields are
multiplied by the face or volume fractions.
*/
event properties (i++) {
  foreach_face() {
    muv.x[] = MU*fm.x[];
    alpha.x[] = 1./RHO*fm.x[];
  }
  foreach() rhov[] = RHO*cm[];
}

/**
Throughout the simulation, the mesh is adapted according to the velocity field,
and a maximum level is enforced around the channel walls and around the
membrane (if applicable).
*/
event adapt (i++) {
  #if INITIAL_FLOW_FIELD
    adapt_wavelet ({mcs, u}, (double[]) {1.e-30, U_TOL,  U_TOL, U_TOL},
      maxlevel = MAXLEVEL, minlevel = MINLEVEL);
  #else
    tag_ibm_stencils();
    adapt_wavelet ({mcs, stencils, u}, (double[]) {1.e-30, 1.e-30,
      U_TOL,  U_TOL, U_TOL}, maxlevel = MAXLEVEL, minlevel = MINLEVEL);
    generate_lag_stencils();
  #endif
}

/**
### Data output
If there is no capsule and we are generating the initial flow field, we monitor
the flow convergence by outputting the changes in the velocity field as well
as the flow rate and average velocity at the outflow.

If a capsule is considered, we dump the position of every Lagrangian node in
order to analyze the membrane position and shape later with a script language
(e.g. Julia).
*/
#if INITIAL_FLOW_FIELD
  event output (i ++) {
#else
  event output (OUTPUT_FREQ) {
#endif
  #if INITIAL_FLOW_FIELD
    double max_deltau = -HUGE;
    double avg_deltau = 0;
    double flow_rate = 0.;
    int nbc = 0;
    foreach(reduction(max:max_deltau)
      reduction(+:avg_deltau) reduction(+:nbc) reduction(+:flow_rate)) {
      if (cm[] > 1.e-10) {
        nbc++;
        double nu = sqrt(sq(u.x[] - prev_u.x[]) + sq(u.y[] - prev_u.y[])
          + sq(u.z[] - prev_u.z[]));
        avg_deltau += nu;
        if (nu > max_deltau) max_deltau = nu;
        foreach_dimension() prev_u.x[] = u.x[];

        /* We also output the flow rate at the boundary of the domain */
        if (fabs(x - 3.*L0/4) < Delta/2.)
          flow_rate += u.x[]*sq(Delta)*fs.x[];
      }
    }
    if (nbc > 0) avg_deltau /= nbc;
    double avg_output_vel = flow_rate/sq(WIDTH);
    fprintf(iff_stats, "%g %d %d %g %g %g %g %g\n", t, i, nbc, avg_deltau,
      max_deltau, normf(u.x).max, flow_rate, avg_output_vel);
    fflush(iff_stats);
  #else
    for(int k=0; k<NCAPS; k++) {
      char fposname[64];
      char ftriname[64];
      sprintf(fposname, "mb_%d_pos.csv", k);
      sprintf(ftriname, "mb_%d_tri.csv", k);
      dump_plain_nodes_pos(&CAPS(k), fposname);
      dump_plain_triangles(&CAPS(k), ftriname);
    }
    #if PARAVIEW_CAPSULE
      pv_output_ascii();
    #endif
    #if PARAVIEW_FLOW_FIELD
      scalar * list = {p};
      vector * vlist = {u};
      save_data(list, vlist, t);
    #endif
  #endif

  /**
Some figures about perfromance are always insightful:
  */
  fprintf(fperf, "%ld %d %g %g %d %d %g %g %d\n", grid->tn, mgu.i, mgu.resb,
    mgu.resa, mgu.nrelax, mgp.i, mgp.resb, mgp.resa, mgp.nrelax);
  fflush(fperf);
}

/**
Monitor progress
*/
event progress ( i++ ) {
  fprintf(stderr, "t = %g, i = %d\n", t, i);
}

/**
The capsule and the flow field are saved every ```10*OUTPUT_FREQ``` iterations.
*/
#if NCAPS > 0
  event save (DUMP_OUTPUT_FREQ) {
    char fname[62];
    sprintf(fname, "flow_%d.dump", nb_pic);
    dump(file = fname);
    sprintf(fname, "caps_%d.dump", nb_pic);
    dump_capsules(fname);
  }
#endif

/**
We generate pictures featuring the velocity norm and the capsule, as well as the
z-vorticity field.
*/
event pictures (OUTPUT_FREQ) {
  scalar omega_x[];
  scalar omega_y[];
  scalar omega_z[];
  scalar un[];
  foreach() un[] = norm(u);

  char name[64];
  view(fov = 18.9*.5, bg = {1,1,1},
    tx = -1./8, ty = 1./8,
    width = 2400, height = 2400);
  clear();
  squares("un", n = {0,0,1}, map = cool_warm, min = 0., max = 3*U_AVG, linear = true);
  #if NCAPS > 0
    draw_lags(lw = .5, edges = true, facets = true);
  #endif
  #if INITIAL_FLOW_FIELD
    sprintf(name, "iff_un_%d.png", nb_pic);
  #else
    sprintf(name, "un_%d.png", nb_pic);
  #endif
  save(name);

  view(fov = 18.9*.2, bg = {1,1,1},
    tx = -.06, ty = .06,
    width = 2400, height = 2400);
  clear();
  squares("un", n = {0,0,1}, map = cool_warm, min = 0., max = 3*U_AVG, linear = true);
  #if NCAPS > 0
    draw_lags(lw = .5, edges = true, facets = true);
  #endif
  #if INITIAL_FLOW_FIELD
    sprintf(name, "iff_un_zoom_%d.png", nb_pic);
  #else
    sprintf(name, "un_zoom_%d.png", nb_pic);
  #endif
  save(name);

  nb_pic += 1.;
}

/**
When $t = t_{end}$, the simulation is stopped. If we were generating the initial
flow field, now is the time to save it to a dump file.
*/
event end (t = T_END) {
  fclose(fperf);
  #if INITIAL_FLOW_FIELD
    fclose(iff_stats);
    char dump_name[64];
    sprintf(dump_name, "initial_flow_field_re%s_lvl%s.dump", STRING(REw),
      STRING(MAXLEVEL));
    dump(file = dump_name);
  #endif
  return 0.;
}

/**
## Results

The data is analized in [Julia](https://en.wikipedia.org/wiki/Julia_(programming_language)) using Pluto notebooks. The graphs are plotted in [post_processing_corner.jl](https://damienhuet.github.io/images/basilisk_sandbox_files/corner_post_processing/post_processing_corner.jl.html), which depends on [caps_toolbox.jl](https://damienhuet.github.io/images/basilisk_sandbox_files/corner_post_processing/caps_toolbox.jl.html).

## References
~~~bib
@Article{zhu2015motion,
  author    = {Zhu, Lailai and Brandt, Luca},
  journal   = {Journal of Fluid Mechanics},
  title     = {The motion of a deforming capsule through a corner},
  year      = {2015},
  pages     = {374--397},
  volume    = {770},
  file      = {:the-motion-of-a-deforming-capsule-through-a-corner.pdf:PDF},
  groups    = {Biological flows},
  publisher = {Cambridge University Press},
}
~~~
*/
