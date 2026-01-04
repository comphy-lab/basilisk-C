//Geometric parameters of the jet
#define D (0.0079)
#define RADIUS (D/2.)
#define LENGTH (2*RADIUS/3.3)
#define AREAI (pi*sq(RADIUS))
#define Qi 1e-4 //0.6 4.1666667e-5 1.5 2.5e-5 2.5 4.1666667e-5 3.5 5.8333333e-5 4.5 7.5e-5 6 1e-4
#define VELOCITY (Qi/AREAI)
#define VOFTHR 1e-6
#define REMOVE_DROPLETS 1
#define ARCELOR 1
#define RESTART_MASS_TRANSFER 0
#define SUBGRID 1

//used in load balancing scheme. implies that there is N times more work in oil-water interfacial cells than in other cells
#define INTERFACE_WEIGHTED 30


// Physical parameters 20°C Re=116


// Physical parameters of fluids 20°C
#define RHOG 1.225*10
#define RHOL 998
#define RHOO 920.0
#define MUG 1.85e-5
#define MUL 1.002e-3
#define MUO 0.079
#define NUL (MUL/RHOL)
#define NUO (MUO/RHOO)

#define SIGMAGL 72e-3
#define SIGMALO 25.5e-3
#define SIGMAGO 31.7e-3

#define ScL 1.
#define ScO 10.

#define ScL10 10.
#define ScO10 100.

#define ScL40 40.
#define ScO40 400.

#define ScL1500 1.48e3
#define ScO1500 1.3e7

coord interface_normal (Point point, scalar c);
#define interface_normal(point, c) interface_normal (point, c)


scalar f10[], f11[];


//scalar triple_point_tag[];

#define rho(f1,f2,f3) ((clamp(f3,0.,1.) + clamp(f2,0.,1.) + clamp(f1,0.,1.)) > 0.0 ? ((clamp(f3,0.,1.)*rho3 + clamp(f2,0.,1.)*rho2 + clamp(f1,0.,1.)*rho1)/(clamp(f3,0.,1.) + clamp(f2,0.,1.) + clamp(f1,0.,1.))) : rho2)

//weighted harmonic mean
#define mu(f1,f2,f3) ((clamp(f3,0.,1.) + clamp(f2,0.,1.) + clamp(f1,0.,1.)) > 0.0 ? (((clamp(f3,0.,1.) + clamp(f2,0.,1.) + clamp(f1,0.,1.))/(clamp(f3,0.,1.)*rho3/mu3 + clamp(f2,0.,1.)*rho2/mu2 + clamp(f1,0.,1.)*rho1/mu1))*rho(f1,f2,f3)) : mu2)


#include "grid/my_octree.h"
#include "navier-stokes/mycentered.h"
#include "navier-stokes/perfs.h"

#include "three_phase_NN_v2.h"
#include "navier-stokes/conserving3f.h"
#include "tension_csf_3f_3D.h"
//#include "tension_3f.h"
#include "species_transfer_two_field_ARCELOR_NN.h"

#include "tag.h"
#include "view.h"

#include "VOFI/include/vofi.h"
#pragma autolink -L$HOME/basilisk/src/VOFI/lib -lvofi

#define WIDTH 0.27
#define Ho 0.2
#define Ha 0.20667
#define MAXTIME 20.


static double xc = 0.0, zc = 0.0, Rc = RADIUS;

static double cylinder (creal p[dimension])
{
  return sq(p[0] - xc) + sq(p[2] - zc) - sq(Rc);
}

static void vofi (scalar c, int levelmax)
{
  double fh = Get_fh (cylinder, NULL, 1./(1 << levelmax), dimension, 0);
  foreach() {
    creal pp[3] = {x - Delta/2., y - Delta/2., z - Delta/2.};
    c[] = Get_cc (cylinder, pp, Delta, fh, dimension);
  }
}


scalar zeros[];


scalar Ts_l1500[];
scalar delta_b_Tsl1500[];
scalar c_boundary_Ts_l1500[];


scalar Ts_o1500[];
scalar delta_b_Tso1500[];
scalar c_boundary_Ts_o1500[];


vector h1[], h2[], h3[];

vector h10[], h11[];

/**
  The default maximum level of refinement is 8 and the error threshold
  on velocity is 0.01. */

int maxlevel = 9;
int minlevel = 1;
double uemax = 0.01;
double cmax = 0.01;
scalar f1_cumt[], f3_cumt[], f1_cumt2[], f3_cumt2[];
vector cumt_u[], cumt2_u[];


double inflow_summ = 0.0;
double inflow_before = 0.0;
scalar inflow_measure[]; 

//Boundary conditions
f1[bottom]   = f10[];
f2[bottom]   = 1.-f10[];
f3[bottom]   = 0.;
inflow_measure[bottom] = f10[];

h1.n[bottom] = h10.y[];
h1.t[bottom] = h10.z[];
h1.r[bottom] = h10.x[];
h2.n[bottom] = h11.y[];
h2.t[bottom] = h11.z[];
h2.r[bottom] = h11.x[];

u.n[bottom]  = dirichlet(VELOCITY*f10[]);
u.t[bottom]  = dirichlet(0.0);
u.r[bottom]  = dirichlet(0.0);
uf.n[bottom]  = VELOCITY*f10[];

p[bottom]   = neumann(0.);


u.t[left] = dirichlet(0);
u.r[left] = dirichlet(0);
u.n[left] = dirichlet(0);
uf.n[left]  = 0.0;

u.t[front] = dirichlet(0);
u.r[front] = dirichlet(0);
u.n[front] = dirichlet(0);
uf.n[front]  = 0.0;

u.t[back] = dirichlet(0);
u.r[back] = dirichlet(0);
u.n[back] = dirichlet(0);
uf.n[back]  = 0.0;

u.t[right] = dirichlet(0);
u.r[right] = dirichlet(0);
u.n[right] = dirichlet(0);
uf.n[right]  = 0.0;

u.n[top] = u.n[] > 0. ? neumann(0) : 0.0;
p[top]    = dirichlet(0.);
inflow_measure[top] = 0.0;

f1[top] = 1.0;
f2[top] = 0.0;
f3[top] = 0.0;


scalar f1_old[], f2_old[], f3_old[];

int main (int argc, char * argv[])
{

  sprintf (NN_filepath, "/ccc/work/cont003/gen14629/maarekja/delim_file.csv");
  read_NN_matrix();

  size (WIDTH);
  origin (-L0/2, 0, -L0/2);
  init_grid (1 << 7);
  rho1 = RHOG;
  rho2 = RHOL;
  rho3 = RHOO;
  mu1 = MUG;
  mu2 = MUL;
  mu3 = MUO;

  f1.refine = f1.prolongation = fraction_refine;
  f2.refine = f2.prolongation = fraction_refine;
  f3.refine = f3.prolongation = fraction_refine;
  f10.refine = f10.prolongation = fraction_refine;
  f11.refine = f11.prolongation = fraction_refine;

  f1.injection = 1.0;
  f2.injection = 0.0;
  f3.injection = 0.0;

  f1.height = h1;
  f2.height = h2;
  f3.height = h3;
  f10.height = h10;
  f11.height = h11;

  f1.sigma = (SIGMAGO+ SIGMAGL- SIGMALO)/2.;
  f2.sigma = (-SIGMAGO+ SIGMAGL+ SIGMALO)/2.;
  f3.sigma = (SIGMAGO- SIGMAGL+ SIGMALO)/2.;

  inflow_measure.gradient = zero;
  f1.tracers = {inflow_measure};

  f2.bl_tracers_NN = {Ts_l1500};


  f3.bl_tracers = {Ts_o1500};

  #if SUBGRID
	fprintf(ferr, "SUBGRID 1\n");
  #else
	fprintf(ferr, "SUBGRID 0\n");
  #endif


  run();
}

event defaults (i = 0) {


  Ts_l1500.gradient = minmod2;
  Ts_l1500.c_boundary = c_boundary_Ts_l1500;
  Ts_l1500.delta_b = delta_b_Tsl1500;
  Ts_l1500.c_inf = 1.0;
  Ts_l1500.interface1 = f3;
  Ts_l1500.phase_diffusivity = (NUL/ScL1500);

  Ts_o1500.gradient = minmod2;
  Ts_o1500.c_boundary = c_boundary_Ts_o1500;
  Ts_o1500.delta_b = delta_b_Tso1500;
  Ts_o1500.c_inf = 0.0;
  Ts_o1500.interface1 = f2;
  Ts_o1500.phase_diffusivity = (NUO/ScO1500);

}



event init (t = 0) {
  if (!restore (file = "restart")) {
    scalar f1bis[];
    refine (y < 1.2*LENGTH && sq(x) + sq(z) < 2.*sq(RADIUS) && level < maxlevel);


  vofi (f10, maxlevel);
  restriction ({f10}); // for boundary conditions on levels

  foreach()
    f11[] = 1-f10[];
  restriction({f11});

  heights (f10, h10);
  heights (f11, h11);

  boundary((scalar *){h10, h11});

  boundary({f10});
  restriction({f10});

  boundary({f11});
  restriction({f11});

  scalar ftemp1[];
  scalar ftemp2[];

  fraction (ftemp1, (Ha-y));
  fraction (ftemp2, (Ho-y));
  foreach(){
    f3[] = ftemp1[] - ftemp2[];
  }

  boundary({f3});
  restriction({f3});

  fraction (f1bis, (y-Ha));

  foreach(){
    f1[] = f1bis[];
    f2[] = ftemp2[];
    if (y < LENGTH){
      f1[] += f10[];
      f2[] -= f10[];
    }
    u.y[] = ((y < LENGTH) ? f1[]*VELOCITY: 0.0);
    }

      foreach(){
    delta_b_Tsl1500[] = -1;

    delta_b_Tso1500[] = -1;
    zeros[] = 0.0;
  }

  foreach(){

    Ts_l1500[] = clamp(f2[],0.,1.);
    Ts_o1500[] = clamp(f3[],0.,1.);
  }

  }
  else{

  fprintf(ferr, "restart\n");
    vofi (f10, maxlevel);
    restriction ({f10}); // for boundary conditions on levels

    foreach()
      f11[] = 1-f10[];
    restriction({f11});

    heights (f10, h10);
    heights (f11, h11);

    boundary((scalar *){h10, h11});

    boundary({f10});
    restriction({f10});

    boundary({f11});
    restriction({f11});
  }

  if (f1.height.x.i){
    heights (f1, f1.height);
  }

  if (f2.height.x.i){
    heights (f2, f2.height);
  }

  if (f3.height.x.i){
    heights (f3, f3.height);
  }

  boundary ((scalar *){f1,f2,f3});


  foreach(){
    f1_old[] = 0.0;
    f2_old[] = 0.0;
    f3_old[] = 0.0;
    inflow_measure[] = 0.0;
  }

  event("properties");


    foreach(){
    delta_b_Tsl1500[] = -1;
    delta_b_Tso1500[] = -1;
    zeros[] = 0.0;
  }

  #if RESTART_MASS_TRANSFER

  fprintf(ferr, "restart mass transfer\n");

  foreach(){


    Ts_l1500[] = clamp(f2[],0.,1.);
    Ts_o1500[] = 0.0;
  }

  #endif

}

static inline void myrestrict(Point point, scalar s){
  s[] = 0;
}

static inline void myprolongation(Point point, scalar s){
  foreach_child()
    s[] = 0;
}

scalar interface[];

event tracer_diffusion(i++){

  f1.refine = f1.prolongation = fraction_refine;
  f2.refine = f2.prolongation = fraction_refine;
  f3.refine = f3.prolongation = fraction_refine;

  f1_old.refine = f1_old.prolongation = fraction_refine;
  f2_old.refine = f2_old.prolongation = fraction_refine;
  f3_old.refine = f3_old.prolongation = fraction_refine;

  interface.prolongation = myprolongation;
  interface.restriction = myrestrict;
  interface.refine = refine_injection;


  //save old volume fractions to for local mass conservation calculation
  foreach(){
  f1_old[] = f1[];
  f2_old[] = f2[];
  f3_old[] = f3[];}

  //remove small artifacts
  foreach(){
    f1[] = (f1[] > VOFTHR ? f1[] : 0.0);
    f2[] = (f2[] > VOFTHR ? f2[] : 0.0);
    f3[] = (f3[] > VOFTHR ? f3[] : 0.0);
  }

  boundary ((scalar *){f1,f2,f3});

  if (f1.height.x.i){
    heights (f1, f1.height);
  }

  if (f2.height.x.i){
    heights (f2, f2.height);
  }

  if (f3.height.x.i){
    heights (f3, f3.height);
  }

  boundary ((scalar *){f1,f2,f3});


  scalar ctemp1500_sl[];
  scalar ctemp1500_so[];
  double interfacial_cells = 0.0;
  double water_clipped = 0.0;
  double oil_clipped = 0.0;
  foreach(reduction(+:oil_clipped) reduction(+:water_clipped) reduction(+:interfacial_cells)){

      ctemp1500_sl[] = (f2[] > VOFTHR ? Ts_l1500[]/f2[] : 0.0);
      ctemp1500_so[] = (f3[] > VOFTHR ? Ts_o1500[]/f3[] : 0.0);
      if ((f2[] > 0. && f2[] < 1.) && (f3[] > 0. && f3[] < 1. )){
        interfacial_cells += 1.0;
	water_clipped += (ctemp1500_sl[] < 0.0) ? 1.0 : 0.0;
	oil_clipped += (ctemp1500_so[] < 0.0) ? 1.0 : 0.0;
	ctemp1500_sl[] = max(ctemp1500_sl[], 0.0);
	ctemp1500_so[] = max(ctemp1500_so[], 0.0);
      }
  }
  fprintf(ferr, "diffusion clipped %g %g\n", water_clipped/interfacial_cells, oil_clipped/interfacial_cells);


#if REMOVE_DROPLETS
  //remove small droplets
  //fprintf(ferr, 'remove_droplets');
  if (i%50 == 0 && i > 0){

//taken from remove droplets function found in tag.h. Instead of just removing the droplet we
//replace the volume with the phase which is a maximum in the nearest neighbors



    scalar d[];
    double threshold = 1e-6;
    int min_diameter = 2;
    int minsize = pow (min_diameter, dimension);

    //f1
    foreach()
      d[] = f1[] > threshold;
    int n1 = tag (d), size1[n1];
    for (int i = 0; i < n1; i++)
      size1[i] = 0;
    foreach (serial)
      if (d[] > 0)
        size1[((int) d[]) - 1]++;
    #if _MPI
    MPI_Allreduce (MPI_IN_PLACE, size1, n1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    #endif
    foreach()
      if (d[] > 0 && size1[((int) d[]) - 1] <= minsize){
        double f1_sum = 0;
        double f2_sum = 0;
        double subtract = 0.0;
        foreach_neighbor(2){
          f1_sum += f2[];
          f2_sum += f3[];
        }
      if (f1_sum > f2_sum){
        subtract = f1[];
        f1[] -= subtract;
        f2[] += subtract;
      }
      else{
        subtract = f1[];
        f1[] -= subtract;
        f3[] += subtract;
      }
    }

    boundary ((scalar *){f1,f2,f3});

    if (f1.height.x.i){
      heights (f1, f1.height);
    }

    if (f2.height.x.i){
      heights (f2, f2.height);
    }

    if (f3.height.x.i){
      heights (f3, f3.height);
    }

    boundary ((scalar *){f1,f2,f3});

   //f3
    foreach()
      d[] = f3[] > threshold;
    int n3 = tag (d), size3[n3];
    for (int i = 0; i < n3; i++)
      size3[i] = 0;
    foreach (serial)
      if (d[] > 0)
        size3[((int) d[]) - 1]++;
    #if _MPI
    MPI_Allreduce (MPI_IN_PLACE, size3, n3, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    #endif
    foreach()
      if (d[] > 0 && size3[((int) d[]) - 1] <= minsize){
        double f1_sum = 0;
        double f2_sum = 0;
        double subtract = 0.0;
        foreach_neighbor(2){
          f1_sum += f1[];
          f2_sum += f2[];
        }
      if (f1_sum > f2_sum){
        subtract = f3[];
        f3[] -= subtract;
        f1[] += subtract;
      }
      else{
        subtract = f3[];
        f3[] -= subtract;
        f2[] += subtract;
      }
    }

    boundary ((scalar *){f1,f2,f3});

    if (f1.height.x.i){
      heights (f1, f1.height);
    }

    if (f2.height.x.i){
      heights (f2, f2.height);
    }

    if (f3.height.x.i){
      heights (f3, f3.height);
    }

    boundary ((scalar *){f1,f2,f3});

    //f2
    foreach()
      d[] = f2[] > threshold;
    int n2 = tag (d), size2[n2];
    for (int i = 0; i < n2; i++)
      size2[i] = 0;
    foreach (serial)
      if (d[] > 0)
        size2[((int) d[]) - 1]++;
    #if _MPI
    MPI_Allreduce (MPI_IN_PLACE, size2, n2, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    #endif
    foreach()
      if (d[] > 0 && size2[((int) d[]) - 1] <= minsize){
        double f1_sum = 0;
        double f2_sum = 0;
        double subtract = 0.0;
        foreach_neighbor(2){
          f1_sum += f1[];
          f2_sum += f3[];
        }
      if (f1_sum > f2_sum){
        subtract = f2[];
        f2[] -= subtract;
        f1[] += subtract;
      }
      else{
        subtract = f2[];
        f2[] -= subtract;
        f3[] += subtract;
      }
    }

  boundary ((scalar *){f1,f2,f3});

  if (f1.height.x.i){
    heights (f1, f1.height);
  }

  if (f2.height.x.i){
    heights (f2, f2.height);
  }

  if (f3.height.x.i){
    heights (f3, f3.height);
  }

  boundary ((scalar *){f1,f2,f3});
  }
#endif

foreach() {
  double summ = (f1[]+f2[]+f3[]);
  if (summ > 0.) {
    f1[] = clamp(f1[]/summ,0.,1.);
    f2[] = clamp(f2[]/summ,0.,1.);
    f3[] = clamp(f3[]/summ,0.,1.);
  }
  else{
    f1[] = 1.0;
    f2[] = 0.0;
    f3[] = 0.0;}
}

boundary ((scalar *){f1,f2,f3});

if (f1.height.x.i){
  heights (f1, f1.height);
}

if (f2.height.x.i){
  heights (f2, f2.height);
}

if (f3.height.x.i){
  heights (f3, f3.height);
}

boundary ((scalar *){f1,f2,f3});
restriction ((scalar *){f1,f2,f3});


foreach(){

  Ts_l1500[] = ctemp1500_sl[]*f2[];
  Ts_o1500[] = ctemp1500_so[]*f3[];
}

foreach(){
    interface[] = 0.0;
    double flagg = 0.0;
    if ( (f2[] > 0. && f2[] < 1.) && (f3[] > 0. && f3[] < 1.) )
	flagg = 1.0;
    bool flag2 = (level  == grid->maxdepth);
    foreach_neighbor(2){
      if ( (f2[] > 0. && f2[] < 1.) && (f3[] > 0. && f3[] < 1. ) && flag2)
        flagg = 1.0;
    }
    flag2 = (level  == grid->maxdepth-1);
    foreach_neighbor(1){
      if ( (f2[] > 0. && f2[] < 1.) && (f3[] > 0. && f3[] < 1. ) && flag2)
        flagg = 1.0;
    }
    if (flagg)
      interface[] = 1.0;
}


event("properties");


if (i%50 == 0){

static FILE * fp1 = fopen ("fluids.dat", "a");
static FILE * fp2 = fopen ("local_mass_conservation.dat", "a");


  double f1_change = 0;
  double f2_change = 0;
  double f3_change = 0;

  foreach (reduction(+:f1_change) reduction(+:f2_change) reduction(+:f3_change)){
    f1_change += fabs(f1[] - f1_old[])*dv();
    f2_change += fabs(f2[] - f2_old[])*dv();
    f3_change += fabs(f3[] - f3_old[])*dv();
  }

  stats p = statsf (f1);
  stats q = statsf (f2);
  stats r = statsf (f3);
  fprintf (fp2, "%d %12e %12e %12e %12e\n", i, dt, f1_change, f2_change, f3_change);
  fprintf (fp1, "%12e %12e %12e %12e %12e %12e %12e %12e %12e %12e\n", t, p.min, p.max, p.sum, q.min, q.max, q.sum, r.min, r.max, r.sum);
  fflush (fp1);
  fflush (fp2);

  static FILE * fp3 = fopen ("oil_water.dat", "a");
  double oil_water_count = 0.0;
  double oil_water_changed = 0.0;
  double oil_water_change_avg = 0.0;
  double oil_water_change_avg2 = 0.0;
  double oil_water_area = 0.0;

  foreach(reduction(+:oil_water_count) reduction(+:oil_water_changed) reduction(+:oil_water_change_avg) reduction(+:oil_water_change_avg2) reduction(+:oil_water_area)){
    if (f2[] > 0.0 && f2[]< 1.0 && f3[] > 0.0 && f3[] < 1.0){
      oil_water_count += 1.0;
      oil_water_changed += (fabs(f2[]- f2_old[]) > 0.0 ? 1 : 0.0);
      oil_water_change_avg += ((f2[]- f2_old[]) - (f3[]- f3_old[]))/2;
      oil_water_change_avg2 += (fabs(f2[] - f2_old[]) + fabs(f3[]- f3_old[]))/2;

      coord n_temp2, n_temp3, p2, p3;
      n_temp2 = interface_normal(point, f2);
      n_temp3 = interface_normal(point, f3);

      //oil_water_area += sqrt(plane_area_center(n_temp2, plane_alpha(f2[], n_temp2), &p2)*plane_area_center(n_temp3, plane_alpha(f3[], n_temp3), &p3))*sq(Delta);
  oil_water_area += 0.5*(plane_area_center(n_temp2, plane_alpha(f2[], n_temp2), &p2) + plane_area_center(n_temp3, plane_alpha(f3[], n_temp3), &p3))*sq(Delta);

    }
  }

  fprintf (fp3, "%g %g %g %g %g\n", oil_water_area, oil_water_count, oil_water_changed, oil_water_change_avg/oil_water_count, oil_water_change_avg2/oil_water_count);
  fflush (fp3);
}

/*f should be the oil phase and f2 should be the water phase*/

//henry_coef should be 350

  double r_sum_output = species_transfer_two_field_henry_sgs(Ts_o1500, Ts_l1500, f3, f2, dt, 350.0, 1e-3*sqrt(Ts_l1500.phase_diffusivity));


  static FILE * fp4 = fopen ("tracers_v2.dat", "a");
  
  fprintf(fp4, "%12e %12e %12e\n", t, dt, r_sum_output);

  fflush (fp4);

  //MPI_Barrier(MPI_COMM_WORLD);
  //fprintf(ferr, "end sgs diffusion\n");
  //MPI_Barrier(MPI_COMM_WORLD);
}

event viscous_term (i++,last){
  TOLERANCE = 1e-6;
}

event acceleration (i++){
  face vector av = a;
  foreach_face(y)
    av.y[] -= 9.81;
  TOLERANCE = Qi*1e-1; //1e-6 for Qi = 1e-5
}



/*
//should not be necessary but saw issues when restart simulations with boundary conditions on pressure field
event projection(i++){
  boundary({p});
}*/




event cumulative (i++) {
  f1_cumt.gradient = minmod;
  f1_cumt2.gradient = minmod;
  f3_cumt.gradient = minmod;
  f3_cumt2.gradient = minmod;
  cumt_u.x.gradient = minmod;
  cumt2_u.x.gradient = minmod;
  cumt_u.y.gradient = minmod;
  cumt2_u.y.gradient = minmod;
  cumt_u.z.gradient = minmod;
  cumt2_u.z.gradient = minmod;
  if (t >=4.){
    foreach(){
      f1_cumt[] += dt*f1[];
      f1_cumt2[] += dt*sq(f1[]);
      f3_cumt[] += dt*f3[];
      f3_cumt2[] += dt*sq(f3[]);
      foreach_dimension(){
       cumt_u.x[] += dt*(uf.x[] + uf.x[1])/2;
       cumt2_u.x[] += dt*sq((uf.x[] + uf.x[1])/2);
      }
    }
  }
  else {
    foreach(){
      f1_cumt[] = 0.;
      f1_cumt2[] = 0.;
      f3_cumt[] = 0.;
      f3_cumt2[] = 0.;
      foreach_dimension(){
        cumt_u.x[] = 0.;
        cumt2_u.x[] = 0.;
      }
    }
  }
}


event statistic (t += 0.005){

  double v_gas =0.0, v_liq =0.0, v_oil =0.0;
  foreach (reduction(+:v_gas) reduction(+:v_liq) reduction(+:v_oil)) {
    v_gas += f1[]*dv();
    v_liq += f2[]*dv();
    v_oil += f3[]*dv();}

  double ket = 0., kef1 = 0.,kef2 = 0., kef3 = 0., voltot = 0.;
  //  double vd=0., er1= 0., er2= 0., er3= 0.;
  foreach(reduction(+:ket) reduction(+:kef1) reduction(+:kef2) reduction(+:kef3)
      reduction(+:voltot)) {
    double  vo = f3[]*dv();
    double  vw = f2[]*dv();
    double  va = f1[]*dv();
    double  vt = (f1[]+f2[]+f3[])*dv();
    // double w2=0.;
    foreach_dimension() {
      //Total kinetic energy
      ket += rho(f1[],f2[],f3[])*sq(u.x[])*vt;
      //Phase specific kinetic energy
      kef1 += rho1*sq(u.x[])*va;
      kef2 += rho2*sq(u.x[])*vw;
      kef3 += rho3*sq(u.x[])*vo;
      voltot += vt;
      // Enstrophy
      //   er1 += dv()*f1[]*sq(u.x[0,1] - u.x[0,-1] - u.y[1,0] + u.y[-1,0])/sq(2.*Delta);
      // viscous dissipation
      //      vd += dv()*(sq(u.x[1] - u.x[-1]) + sq(u.x[0,1] - u.x[0,-1]) + sq(u.x[0,0,1] - u.x[0,0,-1]))/sq(2.*Delta);
    }
    //     w2 /= sq(2.*Delta);
    //     er1 += dv()*f1[]*w2;
  }
  ket /= 2.;
  kef1 /= 2.;
  kef2 /= 2.;
  kef3 /= 2.;

  double u_avg_x = 0.0, u_avg_y = 0.0, u_avg_z = 0.0;
  foreach(reduction(+:u_avg_x) reduction(+:u_avg_y) reduction(+:u_avg_z)){
    u_avg_x += (uf.x[0] + uf.x[1,0,0])/2*(dv()/pow(L0, dimension));
    u_avg_y += (uf.y[0] + uf.y[0,1,0])/2*(dv()/pow(L0, dimension));
    u_avg_z += (uf.z[0] + uf.z[0,0,1])/2*(dv()/pow(L0, dimension));
  }

  static  FILE * fp5 = fopen ("output.dat", "a");
  fprintf (fp5, "%g %12e %12e %12e %12e %12e %12e %12e %12e %12e %12e %12e %12e %12e %12e %12e\n",
      t, ket/voltot, kef1/v_gas, kef2/v_liq, kef3/v_oil, ket, kef1, kef2, kef3,
      voltot, interface_area(f3), interface_area(f1), interface_area(f2), u_avg_x, u_avg_y, u_avg_z);
  fflush (fp5);


  static FILE * fp3 = fopen ("tracers.dat", "a");



  double Tsl_1500_sum = 0.0;
  double Tso_1500_sum = 0.0;


  foreach(reduction(+:Tsl_1500_sum), reduction(+:Tso_1500_sum)){


    Tsl_1500_sum += f2[] > 0.0 ? (1.-Ts_l1500[]/f2[])*f2[]*dv() : 0.0;
    Tso_1500_sum += f3[] > 0.0 ? (Ts_o1500[]/f3[])*f3[]*dv() : 0.0;
    }

  fprintf (fp3, "%g %12e %12e %12e %12e %12e\n", t, v_gas, v_liq, v_oil, Tsl_1500_sum, Tso_1500_sum);
  fflush (fp3);
}


event inflow(i++){
  double inflow_temp1 = 0.0;
  double inflow_temp2 = 0.0;
  foreach(reduction(+:inflow_temp1)){
    if (y > 0.05 || f1[] <= 0.0){
      inflow_temp1 += inflow_measure[]*dv();
      inflow_measure[] -= inflow_measure[];
    }
  }

  foreach(reduction(+:inflow_temp2)){
      inflow_temp2 += inflow_measure[]*dv();
    }

  inflow_summ += inflow_temp2 - inflow_before + inflow_temp1;
  inflow_before = inflow_temp2;

  static FILE * fp5 = fopen ("output_inflow.dat", "a");
  fprintf (fp5, "%12e %12e\n", t, inflow_summ);
  fflush (fp5);
}

event outputfiles (t += 0.1) {
  char name[80];
  sprintf (name, "snapshot-%g", t);
  dump(name, list = {f1, f2, f3, u.x, u.y, u.z, f1_cumt, f3_cumt, f1_cumt2, f3_cumt2, cumt_u.x, cumt_u.y, cumt_u.z, cumt2_u.x, cumt2_u.y, cumt2_u.z, Ts_l1500, Ts_o1500, triple_point_tag});
}

event adapt (i++) {
  //adapt_wavelet ({oil_water, f1,f2,f3, u.x, u.y, u.z, Tl, To, Ts_l, Ts_o}, (double[]){0.1, cmax,cmax,1e-6,uemax,uemax,uemax, 1e-3,1e-3,1e-3,1e-3}, maxlevel);
  adapt_wavelet ({interface, triple_point_tag, f1,f2,f3, u.x, u.y, u.z, Ts_l1500, Ts_o1500}, (double[]){0.1, 0.1, 1e-6, 1e-6,1e-6,uemax,uemax,uemax, 1e-3,1e-2}, maxlevel);

  vofi (f10, maxlevel);
  restriction ({f10}); // for boundary conditions on levels

  foreach()
    f11[] = 1-f10[];
  restriction({f11});

  heights (f10, h10);
  heights (f11, h11);


  scalar ctemp1500_sl[];
  scalar ctemp1500_so[];
  foreach(){

      ctemp1500_sl[] = (f2[] > VOFTHR ? Ts_l1500[]/f2[] : 0.0);
      ctemp1500_so[] = (f3[] > VOFTHR ? Ts_o1500[]/f3[] : 0.0);
  }

  foreach(){
    double summ = (f1[]+f2[]+f3[]);
    if (summ > 0.) {
      f1[] = clamp(f1[]/summ,0.,1.);
      f2[] = clamp(f2[]/summ,0.,1.);
      f3[] = clamp(f3[]/summ,0.,1.);
    }
    else{
      f1[] = 1.0;
      f2[] = 0.0;
      f3[] = 0.0;}
  }

  boundary ((scalar *){f1,f2,f3});

  if (f1.height.x.i){
    heights (f1, f1.height);
  }

  if (f2.height.x.i){
    heights (f2, f2.height);
  }

  if (f3.height.x.i){
    heights (f3, f3.height);
  }

  boundary ((scalar *){f1,f2,f3});
  restriction ((scalar *){f1,f2,f3});


  foreach(){
    Ts_l1500[] = ctemp1500_sl[]*f2[];
    Ts_o1500[] = ctemp1500_so[]*f3[];
  }
}



event movie (t =0; t += 0.01) {
  clear();
  view (fov = 19.5755, quat = {0.0784591,0,0,0.996917}, tx = 0.0223444, ty = -0.502485, bg = {0.3,0.4,0.6}, width = 1600, height = 1200);
  box();
  draw_vof("f1",fc= {0,1,0});
  draw_vof("f3",fc= {1,0,0});
  save ("f1f3.mp4");

  clear();
  box();
  draw_vof ("f1");
  save ("f1.mp4");

  clear();
  box();
  draw_vof("f3",fc= {1,0,0});
  save ("f3.mp4");
}

event movie2 (t =0; t += 0.01) {
  clear();
  view (fov = 19.5755, quat = {0.0784591,0,0,0.996917}, tx = 0.0223444, ty = -0.502485, bg = {0.3,0.4,0.6}, width = 1600, height = 1200);

  clear();
  box();
  cells(alpha=1e-6);
  save ("cells.mp4");

      clear();
    box();
    squares("cumt_u.y",alpha=1e-6, linear = true);
    save ("cum_uy_midplane.mp4");

  //  clear();
  //  box();
  //  squares("f1_cum",alpha=1e-6,linear = true);
  //  // squares("f1_cum",n= {1,0,0},alpha=1e-6,linear = true);
  //  save ("f1_cum.mp4");
  //
  clear();
  view (fov = 32.7291, quat = {0.707107,0,0,-0.707107}, tx = -0.00214276, ty = -0.0118471, bg = {0.3,0.4,0.6}, width = 1600, height = 1200);
  box();
  draw_vof("f3",fc= {1,0,0});
  save ("f3top.mp4");

  clear();
  box();
  draw_vof("f1",fc= {0,1,0});
  draw_vof("f3",fc= {1,0,0});
  save ("f1f3top.mp4");

  //  view (fov = 26.032, quat = {-0.0724867,0.381504,0.030025,0.921031}, tx = -0.0671838, ty = -0.423359, bg = {0.3,0.4,0.6}, width = 1600, height = 1200);
  //  box();
  //  squares("f3_cum",alpha=1e-6,linear = true);
  //  squares("f3_cum",n= {1,0,0},alpha=1e-6,linear = true);
  //  save ("f3_cum_crossplane.mp4");
  //
}

static double _maxruntime = 86400.0;

event runtime (i += 3) {
  mpi_all_reduce (perf.t, MPI_DOUBLE, MPI_MAX);
  if (perf.t >= _maxruntime - 1200.0) { // we allow 5 minutes for termination
    dump(file = "restart", list = {f1, f2, f3, u.x, u.y, u.z, f1_cumt, f3_cumt, f1_cumt2, f3_cumt2, cumt_u.x, cumt_u.y, cumt_u.z, cumt2_u.x, cumt2_u.y, cumt2_u.z, Ts_l1500, Ts_o1500, triple_point_tag});
    return 1; // exit
  }
}


event end (t = MAXTIME) {
  printf ("i = %d t = %g\n", i, t);
}

/*
#if TRACE > 1
event profiling (i += 10) {
  static FILE * fp = fopen ("profiling", "w");
  trace_print (fp, 1);
}
#endif
*/