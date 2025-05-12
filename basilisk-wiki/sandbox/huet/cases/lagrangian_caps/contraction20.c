/**
# Flow of a capsule through a 20:1 channel contraction at $Re = 2000$

<p align="center">
<iframe width="500" height="500" src="https://www.youtube.com/embed/2wirYpP1kaI" title="3D flow of a capsule through a contraction at Re=2000" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" allowfullscreen></iframe>
<iframe width="500" height="500" src="https://www.youtube.com/embed/nx775DXj__M" title="3D flow of a capsule through a 20:1 contraction (full domain)" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" allowfullscreen></iframe>
</p>

## Definition of the non-dimensional numbers and solver parameters
*/

/** Geometry defitions */
#define RADIUS 1.
#define HEIGHT (5*RADIUS)
#define CURVATURE_RADIUS (RADIUS)
#define DEPTH (2*HEIGHT)
#define L0 (40*HEIGHT)

/** Fluid properties */
#define U_AVG 1.
#define RHO 1.
#define REc 2000.
#define MU (RHO*U_AVG*2*HEIGHT/REc)
#define TAU_ADVECTIVE (40*HEIGHT/U_AVG)
#define STOKES false

/** Capsule properties */
#define NCAPS 1
#define CAPILLARY (0.015)
#define E_S (MU*U_AVG/CAPILLARY)
#define ND_EB 0.005
#define E_B (ND_EB*E_S*sq(RADIUS))
#define INITIAL_POSITION_X (-2*HEIGHT)
#define INITIAL_POSITION_Y (2*HEIGHT)
#define INITIAL_POSITION_Z (0.)

/** Solver properties */
#define T_END (TAU_ADVECTIVE + 100.)
#define DT_MAX (1.e-2)
#define U_TOL (0.05*U_AVG)
#define LAGLEVEL 4
#define MINLEVEL 5
#define MAXLEVEL 11
#define MY_TOLERANCE (1.e-3*U_AVG)
#define JACOBI 1
#define OUTPUT_FREQ ( t += .5 )
#define DUMP_OUTPUT_FREQ ( t += 5. )
#define PARAVIEW_CAPSULE 1
#define PARAVIEW_FLOW_FIELD 0
#define RESTART_CASE 0
#define RESTORE_FRAME -1

/** Some additional definitions are required for the capsule-free initial simulation */
#ifndef INITIAL_FLOW_FIELD
  #define INITIAL_FLOW_FIELD 0
#endif
#if INITIAL_FLOW_FIELD
  #undef T_END
  #define T_END (TAU_ADVECTIVE)
  #undef DT_MAX
  #define DT_MAX (T_END/100.)
  #undef OUTPUT_FREQ
  #define OUTPUT_FREQ (i+=10)
  #undef NCAPS
  #define NCAPS 0
#endif
#define STR(s) #s
#define STRING(s) STR(s)

#include "grid/octree.h"
//#include "lagrangian_caps/adapt2.h"
#include "embed.h"
#include "navier-stokes/centered.h"
#include "lagrangian_caps/capsule-ft.h"
#include "lagrangian_caps/neo-hookean-ft.h"
#include "lagrangian_caps/bending-ft.h"
#include "lagrangian_caps/common-shapes-ft.h"
#include "lagrangian_caps/view-ft.h"


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
FILE* foutput = NULL;
int nb_pic;


int main(int argc, char* argv[]) {
    origin(-.5*L0 + 1.e-11, -.5*L0 + 1.e-11, -.5*L0 + 1.e-11);
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

/** Boundary conditions */
u.n[left] = dirichlet(U_AVG/sq(20));
u.t[left] = dirichlet(0.);
u.r[left] = dirichlet(0.);
u.n[right] = neumann(0.);
u.t[right] = dirichlet(0.);
u.r[right] = dirichlet(0.);
u.n[top] = dirichlet(0.);
u.n[bottom] = dirichlet(0.);
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
p[left] = neumann(0.);
pf[left] = neumann(0.);
p[right] = dirichlet(0.);
pf[right] = dirichlet(0.);


#define large_channel (-x)
#define small_channel (intersection(HEIGHT - fabs(z), HEIGHT - fabs(y)))
#define xc1 (CURVATURE_RADIUS)
#define yc1 (HEIGHT + CURVATURE_RADIUS)
#define zc1 (HEIGHT + CURVATURE_RADIUS)
#define horizontal_curvatures (intersection(intersection(sq(x - xc1) + sq(fabs(y) - yc1) - sq(CURVATURE_RADIUS), xc1 - x), yc1 - fabs(y)))
#define vertical_curvatures (intersection(intersection(sq(x - xc1) + sq(fabs(z) - zc1) - sq(CURVATURE_RADIUS), xc1 - x), zc1 - fabs(z)))
#define curvatures (intersection(horizontal_curvatures, vertical_curvatures))
#define rounded_edge ((fabs(y) > HEIGHT || fabs(z) > HEIGHT) && \
  fabs(y) < HEIGHT + CURVATURE_RADIUS && fabs(z) <  HEIGHT + CURVATURE_RADIUS && \
  x > 0 && x < CURVATURE_RADIUS)

scalar mcs[];
int last_introduced_caps;
event init (i = 0) {
    fperf = fopen("perf.csv", "w");
    fprintf(fperf, "total_nb_cells, nb_iter_viscous, resb_viscous, resa_viscous, nb_relax_viscous, nb_iter_pressure, resb_pressure, resa_pressure, nb_relax_pressure\n");

    for (scalar s in {u})
    s.third = true;

    #if INITIAL_FLOW_FIELD
        iff_stats = fopen("iff_stats.csv", "w");
        fprintf(iff_stats, "time iteration number_of_cells avg_velocity_change max_velocity_change max_x_velocity outlet_flow_rate avg_output_vel\n");

        astats ss;
        int ic = 0;
        do {
        fprintf(stderr, "Initial mesh refinement: iteration %d\n", ic);
        ic++;
        solid (cs, fs, union(union(small_channel, large_channel), curvatures));
        foreach() mcs[] = cs[];

        tag_ibm_stencils();
        // foreach() {
        //   if (rounded_edge) stencils[] = noise();
        // }

        ss = adapt_wavelet ((scalar*){mcs, stencils}, (double[]){1.e-30, 1.e-30}, MAXLEVEL, MINLEVEL);
        } while ((ss.nf || ss.nc) && ic < 100);
    #else
        char dump_name[64];
        #if RESTART_CASE
            sprintf(dump_name, "flow_%s.dump", STRING(RESTORE_FRAME));
            restore(dump_name);
            sprintf(dump_name, "caps_%s.dump", STRING(RESTORE_FRAME));
            restore_capsules(dump_name);
            last_introduced_caps = NCAPS - 1;
	    foutput = fopen("output.csv", "a");
        #else
	    foutput = fopen("output.csv", "w");
            sprintf(dump_name, "initial_flow_field_re%s_lvl%s.dump", 
                STRING(REw), STRING(MAXLEVEL));
            restore(dump_name);
            activate_spherical_capsule(&CAPS(0), level = LAGLEVEL,
                radius = RADIUS);
            last_introduced_caps = 0;
            for(int k=0; k<CAPS(0).nln; k++) {
                CAPS(0).nodes[k].pos.x += INITIAL_POSITION_X;
                CAPS(0).nodes[k].pos.y += INITIAL_POSITION_Y;
                CAPS(0).nodes[k].pos.z += INITIAL_POSITION_Z;
            }
            generate_lag_stencils(no_warning = true);
            comp_centroid(&CAPS(0));
        #endif
            /**
        With the membrane properly created to an initially strain-free spherical
        shape, we now proceed to refining the Eulerian mesh around the membrane.
            */
        astats ss;
        int ic = 0;
        do {
            fprintf(stderr, "Initial mesh refinement: iteration %d\n", ic);
            ic++;
            solid (cs, fs, union(union(small_channel, large_channel), 
                curvatures));
            foreach() mcs[] = cs[];
            tag_ibm_stencils();
            // foreach() if (rounded_edge) stencils[] = noise();
            ss = adapt_wavelet ((scalar*){mcs, stencils}, (double[]){1.e-30, 1.e-30}, MAXLEVEL, MINLEVEL);
            generate_lag_stencils(no_warning = true);
        } while ((ss.nf || ss.nc) && ic < 100);


    /**
    We also erase the content of the files where the membrane position and 
    triangle areas will be stored, and store the initial configuration.
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
  tag_ibm_stencils();
  foreach() if (rounded_edge) stencils[] = noise();
  adapt_wavelet ((scalar*) {mcs, stencils, u}, (double[]) {1.e-30, 1.e-30,
      U_TOL,  U_TOL, U_TOL}, MAXLEVEL, MINLEVEL);
  #if NCAPS>0
    generate_lag_stencils();
  #endif
}

vector prev_u[];
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
          if (fabs(x - L0/2.) < Delta/2.)
            flow_rate += u.x[]*sq(Delta)*fs.x[];
        }
      }
      if (nbc > 0) avg_deltau /= nbc;
      double avg_output_vel = flow_rate/sq(2*HEIGHT);
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
  scalar un[];
  foreach() un[] = norm(u);

  char name[64];
  view(fov = 18.9, bg = {1,1,1},
    // tx = -1./8, ty = 1./8,
    width = 2400, height = 2400);
  clear();
  cells(n = {0., 0., 1.});
  squares("un", n = {0,0,1}, map = cool_warm, min = 0., max = 2*U_AVG, linear = true);
  #if NCAPS > 0
    draw_lags(lw = .5, edges = true, facets = true);
  #endif
  #if INITIAL_FLOW_FIELD
    sprintf(name, "iff_un_%d.png", nb_pic);
  #else
    sprintf(name, "un_%d.png", nb_pic);
  #endif
  save(name);

  view(fov = 18.9*.25, bg = {1,1,1},
    tx = -.04,
    width = 2400, height = 2400);
  clear();
  cells(n = {0., 0., 1.});
  squares("un", n = {0,0,1}, map = cool_warm, min = 0., max = 2*U_AVG, linear = true);
  #if NCAPS > 0
    draw_lags(lw = .5, edges = true, facets = true);
  #endif
  #if INITIAL_FLOW_FIELD
    sprintf(name, "iff_un_zoom_%d.png", nb_pic);
  #else
    sprintf(name, "un_zoom_%d.png", nb_pic);
  #endif
  save(name);

  view(fov = 18.9*.25, bg = {1,1,1},
    tx = -.04,
    width = 2400, height = 2400);
  clear();
  squares("un", n = {0,0,1}, map = cool_warm, min = 0., max = 2*U_AVG, linear = true);
  #if NCAPS > 0
    draw_lags(lw = .2, edges = true, facets = true);
  #endif
  #if INITIAL_FLOW_FIELD
    sprintf(name, "iff_un_zoom_nogrid_%d.png", nb_pic);
  #else
    sprintf(name, "un_zoom_nogrid_%d.png", nb_pic);
  #endif
  save(name);

  nb_pic += 1.;

  fprintf(foutput, "%g, %d, %g, %g, %g\n", t, i, CAPS(0).centroid.x,
            CAPS(0).centroid.y, CAPS(0).centroid.z);
  fflush(foutput);
}

event end (t = T_END) {
    fclose(fperf);
    #if INITIAL_FLOW_FIELD
        fclose(iff_stats);
        char dump_name[64];
        sprintf(dump_name, "initial_flow_field_re%s_lvl%s.dump", 
            STRING(REw), STRING(MAXLEVEL));
        dump(file = dump_name);
    #else
        fclose(foutput);
    #endif
    return 0.;
}
