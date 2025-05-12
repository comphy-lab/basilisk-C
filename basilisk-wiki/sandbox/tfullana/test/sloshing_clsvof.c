// #define rho(f) (clamp(f,0,1)*(rho1 - rho2) + rho2)
// #define mu(f)  (1./(clamp(f,0,1)*(1./mu1 - 1./mu2) + 1./mu2))
#include "grid/multigrid3D.h"
// #include "grid/octree.h"
#include "embed.h"
#include "navier-stokes/centered.h"
#include "mixed-embed.h"
// #include "two-phase-clsvof.h"
#include "two-phase.h"
#include "redistance.h"
// #include "navier-stokes/conserving.h"
#include "tension.h"
#include "contact-embed2.h"
#include "navier-stokes/perfs.h"
#include "vof-tracer-particles.h"
#include "view.h"
#include "scatter2.h"

/* --------------------------- */
/* Dimensional properties (SI) */
/* --------------------------- */ 

/* Oil or Water in fluid 1 */ 
double nu = 50e-6;

double rho_1 = 899.1786;
double sig = 21e-3; 
#define mu_1 (rho_1*nu)

/* Air in fluid 2 */
double rho_2 = 1.2047;           
double mu_2 = 1.8205e-5;         

double grav = 9.8067;                                         
double radius = 51.2e-3;                             
double initial_height = 111e-3;                      
#define total_height (1.15*initial_height)
double z_measurements = 100e-3;

double w_1 = 180.; //rpm

double ratio_f = 0.67;
double eps = 0.057;

#define period (1./(ratio_f*w_1/60.))
#define omega (2*pi*ratio_f*w_1/60.)
#define amplitude (eps*radius)

#define u_ref (omega*radius)                        /* Reference velocity */     

int n_periods = 110;
#define tend (n_periods*period + 0.25*period)                                 /* Final time */ 
#define tout 0.1                                  /* Extract data files every */
#define tdump period                                 /* Extract dump files every */
#define tmovie 0.1                                /* Extract movie every */

/* ----------------------------- */
/* Dimensionless properties (SI) */ 
/* ----------------------------- */ 

#define MU_unit (rho_1*u_ref*radius)              /* Scaling viscosity */
#define SIGMA_unit (rho_1*radius*u_ref*u_ref)     /* Scaling surface tension */
#define T_unit (radius/u_ref)                     /* Characteristic time scale */
#define G_unit (radius/(T_unit*T_unit))           /* Scaling gravity */
#define RHO_1 1.                                  /* Dimensionless density 1 */
#define RHO_2 (rho_2/rho_1)                       /* Dimensionless density 2 */
#define MU_1 (mu_1/MU_unit)                       /* Dimensionless viscosity 1 */
#define MU_2 (mu_2/MU_unit)                       /* Dimensionless viscosity 2 */
#define SIG (sig/SIGMA_unit)                      /* Dimensionless surface tension */
#define GRAV (grav/G_unit)                        /* Dimensionless gravity */
#define RADIUS 1.                                 /* Dimensionless diameter */
#define U 1.                                      /* Dimensionless velocity */
#define AMP (eps)                                 /* Dimensionless stage amplitude */
#define PERIOD (period/T_unit)                      /* Dimensionless stage frequency */
#define OMEGA (omega*T_unit)                      /* Dimensionless stage frequency */
#define TEND (tend/T_unit)                        /* Dimensionless final time */ 
#define TOUT (tout/T_unit)                        /* Dimensionless extraction time data*/
#define TDUMP (tdump/T_unit)                      /* Dimensionless extraction time dump*/
#define TMOVIE (tmovie/T_unit)                    /* Dimensionless extraction time movie*/
#define L0 (total_height/radius)                  /* Dimensionless domain size */
#define INITIAL_HEIGHT (initial_height/radius)    /* Dimensionless initial position */
#define Z_MEASUREMENTS (z_measurements/radius)

/* --------------------- */
/* Tank movement */ 
/* --------------------- */

#define TANK_VEL_x(t) AMP*OMEGA*sin(OMEGA*t)     /* Dimensionless stage velocity x direction */
#define TANK_VEL_y(t) -AMP*OMEGA*cos(OMEGA*t)     /* Dimensionless stage velocity y direction */

#define TANK_ACC_x(t) AMP*OMEGA*OMEGA*cos(OMEGA*t)     /* Dimensionless stage acceleration x direction */
#define TANK_ACC_y(t) AMP*OMEGA*OMEGA*sin(OMEGA*t)     /* Dimensionless stage acceleration y direction */

/* --------------------- */
/* Adaptivity parameters */ 
/* --------------------- */ 

int LEVEL = 7;                                    /* Level of refinement */
double angle_eq = 90.;
/* ------------------- */
/* Boundary conditions */ 
/* ------------------- */

u.n[back] = dirichlet(0);
u.n[front] = dirichlet(0);
u.n[top] = dirichlet(0);
u.n[bottom] = dirichlet(0);
u.n[left] = dirichlet(0);
u.n[right] = dirichlet(0);

u.t[back] = dirichlet(0);
u.t[front] = dirichlet(0);
u.t[top] = dirichlet(0);
u.t[bottom] = dirichlet(0);
u.t[left] = dirichlet(0);
u.t[right] = dirichlet(0);

u.r[back] = dirichlet(0);
u.r[front] = dirichlet(0);
u.r[top] = dirichlet(0);
u.r[bottom] = dirichlet(0);
u.r[left] = dirichlet(0);
u.r[right] = dirichlet(0);

double l_cl = 100. [0];
double l_d = 1e-4 [0];
double _pf = 5. [0]; // Prefactor for slip region width
double _num_slip = 1. [0]; // Grid dependent slip everywhere
int _vert = 1; // Flag for vertical slip
int _tanh = 1; // Flag for tanh localization slip function
#define delta_slip (1./sqrt(RHO_1*GRAV*RADIUS*RADIUS/SIG))

double exp_slip(double z, double l_cl, double l_d, double _pf) {
    return (l_cl > l_d ? l_cl*exp(-(z/(_pf*delta_slip))*log(l_d/l_cl)) : 0.);
}

double tanh_slip(double z, double l_cl, double _pf){
  return l_cl*(1. - pow(tanh(z/(_pf*delta_slip)), 2));
}

const scalar c[] = angle_eq*pi/180.;
scalar slip_length[];
const scalar u_wall[] = 0.;

u.n[embed] = (_vert == 1 ? dirichlet(0.) : dirichlet(mx.x[]));
u.t[embed] = (_vert == 1 ? dirichlet(0.) : dirichlet(mx.y[]));
u.r[embed] = dirichlet(mx.z[]);

// u.n[embed] = dirichlet(0.);
// u.t[embed] = dirichlet(0.);
// u.r[embed] = dirichlet(0.);

vector mean_u[];
vector wave_u[];
scalar ff[];

Particles line_trackers;
Particles hori_trackers;
Particles vert_trackers;
int nb_particles = 100;
#define lstep (1./(double)nb_particles)
#define hstep (2./(double)nb_particles)
#define nb_vert_particles (int)nb_particles*((int)INITIAL_HEIGHT)
#define vstep ((double)INITIAL_HEIGHT/(double)nb_vert_particles)

// Cantor pairing function
unsigned int cantor_pairing(unsigned int x, unsigned int y) {
    return (x + y) * (x + y + 1) / 2 + y;
}

// Inverse Cantor pairing function
void inverse_cantor_pairing(unsigned int val, unsigned int *x, unsigned int *y) {
    unsigned int w = (unsigned int)(floor((sqrt(8 * val + 1) - 1) / 2));
    *y = val - (w * (w + 1)) / 2;
    *x = w - *y;
}

int count_eulerian = 0;
int count_lagrangian = 0;
int count_waveflow = 0;

scalar d[];

int main(int argc, char *argv[]) {
    if (argc > 1)
      nu = atof(argv[1])*1e-6;
    if (argc > 2)
      eps = atof(argv[2]);
    if (argc > 3)
      ratio_f = atof(argv[3]);
    if (argc > 4)
      angle_eq = atof(argv[4]);
    if (argc > 5)
      LEVEL = atoi(argv[5]);
    if (argc > 6)
      nb_particles = atoi(argv[6]);
    if (argc > 7)
      l_cl = atof(argv[7]);
    if (argc > 8)
      _pf = atof(argv[8]);
    if (argc > 9)
      _num_slip = atof(argv[9]);
    if (argc > 10)
      _vert = atoi(argv[10]);
    if (argc > 11)
      _tanh = atoi(argv[11]);

    size (L0);
    origin(-L0/2, -L0/2, 0.);
    init_grid (1 << LEVEL);

    rho1 = RHO_1, mu1 = MU_1;
    rho2 = RHO_2, mu2 = MU_2;

    f.sigma = SIG;

    DT = HUGE [0];
    // TOLERANCE = 1e-3 [*];

    contact_angle = c;
    lambda = slip_length;
    U0 = u_wall;
    run();

}

event init (t = 0) {
    if (!restore(file = "restart"))
    {
    char name[360];
    sprintf (name, "properties_nu%g_eps%g_omega%g_angle%g_LEVEL%d_l%g_delta%g_num%g_vert%d_tanh%d", nu*1e6, eps, ratio_f, angle_eq, LEVEL, l_cl, _pf, _num_slip, _vert, _tanh);
    FILE * fp = fopen(name, "w");

    fprintf(fp, "Re = %g, Re_nodim = %g\n", rho_1*omega*radius*radius/mu_1, RHO_1*OMEGA*RADIUS*RADIUS/MU_1);
    fprintf(fp, "Bo = %g, Bo_nodim = %g\n", rho_1*grav*radius*radius/sig, RHO_1*GRAV*RADIUS*RADIUS/SIG);
    fprintf(fp, "Fr = %g, Fr_nodim = %g\n", amplitude*omega*omega/grav, AMP*OMEGA*OMEGA/GRAV);
    fprintf(fp, "AMP OMEGA OMEGA %g\n", AMP*OMEGA*OMEGA);
    fprintf(fp, "GRAV %g\n", GRAV);
    fprintf(fp, "OMEGA %g\n", OMEGA);
    fprintf(fp, "T_unit %g\n", T_unit);
    fprintf(fp, "u_ref %g\n", u_ref);
    fprintf(fp, "PERIOD %g\n", PERIOD);
    fprintf(fp, "TEND %g\n", TEND);
    fprintf(fp, "Grid (m) %g\n", total_height/pow(2,LEVEL)/npe());
    fprintf(fp, "Z_MEASUREMENTS %g\n", Z_MEASUREMENTS);
    
    fclose(fp);

    solid (cs, fs, (sq(RADIUS) - sq(y) - sq(x)));
    fractions_cleanup (cs, fs);

    fraction (f, -(z - INITIAL_HEIGHT));
    foreach()
      d[] = -(z - INITIAL_HEIGHT);

    line_trackers = new_vof_tracer_particles (0, 1);
    hori_trackers = new_vof_tracer_particles (0, 1);
    vert_trackers = new_vof_tracer_particles (0, 1);
    }else{
    solid (cs, fs, (sq(RADIUS) - sq(y) - sq(x)));
    fractions_cleanup (cs, fs);
    }
}

event acceleration (i++){
  face vector av = a;
  foreach_face(x){
    av.x[] += TANK_ACC_x(t);
  }
  foreach_face(y){
    av.y[] += TANK_ACC_y(t);
  }
  foreach_face(z){
    av.z[] -= GRAV;
  }
  foreach(){
    if (cs[] == 1.){
      ff[] = f[];
    }
  }
}

#if TREE
event adapt(i++)
{
  adapt_wavelet({cs, f, u}, (double[]){0., 0., sq(eps), sq(eps), sq(eps)}, LEVEL);
}
#endif

event output_trace (i++, t <= TEND){
  char tracefile[360];
  sprintf (tracefile, "trace_nu%g_eps%g_omega%g_angle%g_LEVEL%d_l%g_delta%g_num%g_vert%d_tanh%d", 
      nu*1e6, eps, ratio_f, angle_eq, LEVEL, l_cl, _pf, _num_slip, _vert, _tanh);
  static FILE *fp1 = fopen(tracefile, "w");

  vector h[];
  heights(f, h);

  double max_height = 0., u0_mag = 0., u0_x = 0., u0_y = 0.;
  foreach (reduction(max:max_height)){
    if ((f[] < (1. - 1e-6)) & (f[] > (1e-6)) & (cs[] == 1.0) & (h.z[] != nodata)){
      double tmp = z + Delta * height(h.z[]);
      if (tmp > max_height) {
        max_height = tmp;
      }
    }
  }

  foreach_point (0.0, 0.0, Z_MEASUREMENTS, reduction(max:u0_mag) reduction(max:u0_x) reduction(max:u0_y)){
    u0_mag = sqrt(sq(u.x[]) + sq(u.y[]));
    u0_x = fabs(u.x[]);
    u0_y = fabs(u.y[]);
  }

  fprintf(fp1, "%-12g %-12g %-12g %-12g %-12g %-12g\n",
          t, t/PERIOD, max_height - INITIAL_HEIGHT, u0_mag, u0_x, u0_y);
  fflush(fp1);


  double weight = 0.1;
  foreach()
    if (f[] > 1e-6 && f[] < 1. - 1e-6) {
      coord n = interface_normal (point, f);
      normalize (&n);
      double alpha = plane_alpha (f[], n);
      d[] = (1. - weight)*d[] + weight*Delta*alpha;
    }
  
  redistance (d, imax = 3);

  foreach()
    slip_length[] = _num_slip*Delta + Delta*(_tanh == 0 ? exp_slip(-fabs(d[]), l_cl, l_d, _pf) : tanh_slip(-fabs(d[]), l_cl, _pf));
}

event Eulerian_meanflow (t = pi/2., t += PERIOD/1000.){
  foreach(){
    foreach_dimension(){
      mean_u.x[] += u.x[];
    }
  }
  count_eulerian += 1;
  if (count_eulerian == 1000){ // End of the period
    char tmp_Eulerian_line[360], Eulerian_line[360], tmp_Eulerian_hori[360], Eulerian_hori[360], tmp_Eulerian_vert[360], Eulerian_vert[360];
    sprintf (tmp_Eulerian_line, "tmp_Eulerian_line_t%g_pid%d", t, pid());
    sprintf (tmp_Eulerian_hori, "tmp_Eulerian_hori_t%g_pid%d", t, pid());
    sprintf (tmp_Eulerian_vert, "tmp_Eulerian_vert_t%g_pid%d", t, pid());
    FILE *fp1 = fopen(tmp_Eulerian_line, "w");
    FILE *fp2 = fopen(tmp_Eulerian_hori, "w");
    FILE *fp3 = fopen(tmp_Eulerian_vert, "w");
    sprintf (Eulerian_line, "Eulerian_line_nu%g_eps%g_omega%g_angle%g_LEVEL%d_l%g_delta%g_num%g_vert%d_tanh%d_t%g", nu*1e6, eps, ratio_f, angle_eq, LEVEL, l_cl, _pf, _num_slip, _vert, _tanh, (t - pi/2. + 2.*pi/1000.)/(2.*pi));
    sprintf (Eulerian_hori, "Eulerian_hori_nu%g_eps%g_omega%g_angle%g_LEVEL%d_l%g_delta%g_num%g_vert%d_tanh%d_t%g", nu*1e6, eps, ratio_f, angle_eq, LEVEL, l_cl, _pf, _num_slip, _vert, _tanh, (t - pi/2. + 2.*pi/1000.)/(2.*pi));
    sprintf (Eulerian_vert, "Eulerian_vert_nu%g_eps%g_omega%g_angle%g_LEVEL%d_l%g_delta%g_num%g_vert%d_tanh%d_t%g", nu*1e6, eps, ratio_f, angle_eq, LEVEL, l_cl, _pf, _num_slip, _vert, _tanh, (t - pi/2. + 2.*pi/1000.)/(2.*pi));
    foreach(){
      if ((z > (Z_MEASUREMENTS - Delta/2.)) & (z < (Z_MEASUREMENTS + Delta/2.)) & (y < Delta) & (y > 0) & (fabs(x) <= 1.) & (x > 0.)){
        double xvel = ((mean_u.x[]+mean_u.x[0,-1,0])/2.)/1000.;
        double yvel = ((mean_u.y[]+mean_u.y[0,-1,0])/2.)/1000.;
        double zvel = ((mean_u.z[]+mean_u.z[0,-1,0])/2.)/1000.;
        fprintf(fp1, "%-12g %-12g %-12g %-12g %-12g %-12g\n", x, y, z, xvel, yvel, zvel);
      }
    }
    foreach(){
      if ((z > (Z_MEASUREMENTS - Delta/2.)) & (z < (Z_MEASUREMENTS + Delta/2.)) & (cs[] > 0.)){
        fprintf(fp2, "%-12g %-12g %-12g %-12g %-12g %-12g\n", x, y, z, mean_u.x[]/1000., mean_u.y[]/1000., mean_u.z[]/1000.);
      }
    }
    foreach(){
      if ((y < Delta) & (y > 0) & (cs[] > 0.)){
        double xvel = ((mean_u.x[]+mean_u.x[0,-1,0])/2.)/1000.;
        double yvel = ((mean_u.y[]+mean_u.y[0,-1,0])/2.)/1000.;
        double zvel = ((mean_u.z[]+mean_u.z[0,-1,0])/2.)/1000.;
        fprintf(fp3, "%-12g %-12g %-12g %-12g %-12g %-12g\n", x, y, z, xvel, yvel, zvel);
      }
    }
    fclose(fp1);
    fclose(fp2);
    fclose(fp3);
    MPI_Barrier(MPI_COMM_WORLD);
    if (pid() == 0) {
      char command1[1000];
      sprintf(command1, "cat tmp_Eulerian_line_t%g_pid* > %s", t, Eulerian_line);
      system(command1);
  
      char command2[1000];
      sprintf(command2, "cat tmp_Eulerian_hori_t%g_pid* > %s", t, Eulerian_hori);
      system(command2);

      char command3[1000];
      sprintf(command3, "cat tmp_Eulerian_vert_t%g_pid* > %s", t, Eulerian_vert);
      system(command3);

      char command4[1000];
      sprintf(command4, "rm tmp_Eulerian*");
      system(command4);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    count_eulerian = 0;
    foreach(){
      foreach_dimension(){
        mean_u.x[] = 0.;
      }
    }
  }
}

event Lagrangian_meanflow (t = pi/2., t += PERIOD/1000.){
  if (count_lagrangian == 0){ // Put Lagrangian particles at the begining of the period
    for (int x_tag = 0; x_tag <= nb_particles - 1; x_tag += 1) {
      particle p = {.x = 0. + lstep/2. + lstep*x_tag, .y = 0., .z = Z_MEASUREMENTS, .tag = x_tag};
      set_a_particle_attributes (&p);
      add_particle (p, line_trackers);
    }
    for (int x_tag = 0 ; x_tag <= nb_particles - 1; x_tag += 1) {
      for (int y_tag = 0 ; y_tag <= nb_particles - 1; y_tag += 1) {
        unsigned int cantor_tag = cantor_pairing(x_tag, y_tag);
        particle p = {.x = -1. + hstep/2. + hstep*x_tag, .y = -1. + hstep/2. + hstep*y_tag, .z = Z_MEASUREMENTS, .tag = cantor_tag};
        set_a_particle_attributes (&p);
        add_particle (p, hori_trackers);
      }
    }
    remove_particles (hori_trackers, sqrt(sq(x) + sq(y)) > 1.); // Remove particles outside of the cylinder 
    for (int x_tag = 0 ; x_tag <= nb_particles - 1; x_tag += 1) {
      for (int z_tag = 0 ; z_tag <= nb_vert_particles - 1; z_tag += 1) {
        unsigned int cantor_tag = cantor_pairing(x_tag, z_tag);
        particle p = {.x = -1. + hstep/2. + hstep*x_tag, .y = 0., .z = 0. + vstep/2. + vstep*z_tag, .tag = cantor_tag};
        set_a_particle_attributes (&p);
        add_particle (p, vert_trackers);
      }
    }
    remove_particles (vert_trackers, sqrt(sq(x) + sq(y)) > 1.); // Remove particles outside of the cylinder
  }
  count_lagrangian += 1;
  if (count_lagrangian == 1001){ // End of the period
    // fprintf(stdout, "%d %g %g\n", count_lagrangian, (t - pi/2.), (t - pi/2.)/(2.*pi));
    char tmp_Lagrangian_line[360], Lagrangian_line[360], tmp_Lagrangian_hori[360], Lagrangian_hori[360], tmp_Lagrangian_vert[360], Lagrangian_vert[360];
    sprintf (tmp_Lagrangian_line, "tmp_Lagrangian_line_t%g_pid%d", t, pid());
    sprintf (tmp_Lagrangian_hori, "tmp_Lagrangian_hori_t%g_pid%d", t, pid());
    sprintf (tmp_Lagrangian_vert, "tmp_Lagrangian_vert_t%g_pid%d", t, pid());
    FILE *fp1 = fopen(tmp_Lagrangian_line, "w");
    FILE *fp2 = fopen(tmp_Lagrangian_hori, "w");
    FILE *fp3 = fopen(tmp_Lagrangian_vert, "w");
    sprintf (Lagrangian_line, "Lagrangian_line_nu%g_eps%g_omega%g_angle%g_LEVEL%d_l%g_delta%g_num%g_vert%d_tanh%d_t%g", nu*1e6, eps, ratio_f, angle_eq, LEVEL, l_cl, _pf, _num_slip, _vert, _tanh, (t - pi/2.)/(2.*pi));
    sprintf (Lagrangian_hori, "Lagrangian_hori_nu%g_eps%g_omega%g_angle%g_LEVEL%d_l%g_delta%g_num%g_vert%d_tanh%d_t%g", nu*1e6, eps, ratio_f, angle_eq, LEVEL, l_cl, _pf, _num_slip, _vert, _tanh, (t - pi/2.)/(2.*pi));
    sprintf (Lagrangian_vert, "Lagrangian_vert_nu%g_eps%g_omega%g_angle%g_LEVEL%d_l%g_delta%g_num%g_vert%d_tanh%d_t%g", nu*1e6, eps, ratio_f, angle_eq, LEVEL, l_cl, _pf, _num_slip, _vert, _tanh, (t - pi/2.)/(2.*pi));
    foreach_particle_in(line_trackers){
      double xp = 0. + lstep/2. + lstep*p().tag;
      double yp = 0.;
      double zp = Z_MEASUREMENTS;
      fprintf(fp1, "%-12g %-12g %-12g %-12g %-12g %-12g %-12d\n", xp, yp, zp, (x - xp)/(2.*pi), (y - yp)/(2.*pi), (z - zp)/(2.*pi), p().tag);
    }
    foreach_particle_in(hori_trackers){
      unsigned int x_tag, y_tag;
      inverse_cantor_pairing(p().tag, &x_tag, &y_tag);
      double xp = -1. + hstep/2. + hstep*x_tag;
      double yp = -1. + hstep/2. + hstep*y_tag;
      double zp = Z_MEASUREMENTS;
      fprintf(fp2, "%-12g %-12g %-12g %-12g %-12g %-12g %-12d\n", xp, yp, zp, (x - xp)/(2.*pi), (y - yp)/(2.*pi), (z - zp)/(2.*pi), p().tag);
    }
    foreach_particle_in(vert_trackers){
      unsigned int x_tag, z_tag;
      inverse_cantor_pairing(p().tag, &x_tag, &z_tag);
      double xp = -1. + hstep/2. + hstep*x_tag;
      double yp = 0.;
      double zp = 0. + vstep/2. + vstep*z_tag;
      fprintf(fp3, "%-12g %-12g %-12g %-12g %-12g %-12g %-12d\n", xp, yp, zp, (x - xp)/(2.*pi), (y - yp)/(2.*pi), (z - zp)/(2.*pi), p().tag);
    }
    fclose(fp1);
    fclose(fp2);
    fclose(fp3);
    MPI_Barrier(MPI_COMM_WORLD);
    if (pid() == 0) {
      char command1[1000];
      sprintf(command1, "cat tmp_Lagrangian_line_t%g_pid* > %s", t, Lagrangian_line);
      system(command1);
  
      char command2[1000];
      sprintf(command2, "cat tmp_Lagrangian_hori_t%g_pid* > %s", t, Lagrangian_hori);
      system(command2);

      char command3[1000];
      sprintf(command3, "cat tmp_Lagrangian_vert_t%g_pid* > %s", t, Lagrangian_vert);
      system(command3);

      char command4[1000];
      sprintf(command4, "rm tmp_Lagrangian*");
      system(command4);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    count_lagrangian = 0;
    remove_particles (line_trackers, z < 1e15);
    remove_particles (hori_trackers, z < 1e15);
    remove_particles (vert_trackers, z < 1e15);
  }
  
  if (count_lagrangian == 0){ // Put Lagrangian particles at the begining of the period
    for (int x_tag = 0; x_tag <= nb_particles - 1; x_tag += 1) {
      particle p = {.x = 0. + lstep/2. + lstep*x_tag, .y = 0., .z = Z_MEASUREMENTS, .tag = x_tag};
      set_a_particle_attributes (&p);
      add_particle (p, line_trackers);
    }
    for (int x_tag = 0 ; x_tag <= nb_particles - 1; x_tag += 1) {
      for (int y_tag = 0 ; y_tag <= nb_particles - 1; y_tag += 1) {
        unsigned int cantor_tag = cantor_pairing(x_tag, y_tag);
        particle p = {.x = -1. + hstep/2. + hstep*x_tag, .y = -1. + hstep/2. + hstep*y_tag, .z = Z_MEASUREMENTS, .tag = cantor_tag};
        set_a_particle_attributes (&p);
        add_particle (p, hori_trackers);
      }
    }
    remove_particles (hori_trackers, sqrt(sq(x) + sq(y)) > 1.); // Remove particles outside of the cylinder 
    for (int x_tag = 0 ; x_tag <= nb_particles - 1; x_tag += 1) {
      for (int z_tag = 0 ; z_tag <= nb_vert_particles - 1; z_tag += 1) {
        unsigned int cantor_tag = cantor_pairing(x_tag, z_tag);
        particle p = {.x = -1. + hstep/2. + hstep*x_tag, .y = 0., .z = 0. + vstep/2. + vstep*z_tag, .tag = cantor_tag};
        set_a_particle_attributes (&p);
        add_particle (p, vert_trackers);
      }
    }
    remove_particles (vert_trackers, sqrt(sq(x) + sq(y)) > 1.); // Remove particles outside of the cylinder
    count_lagrangian += 1;
  }
}

event wave_flow(t = pi/2., t += PERIOD){
  coord us = {0., 0., 0.};
  char tmp_Wavefile[360], Wavefile[360];
  sprintf (tmp_Wavefile, "tmp_Waveflow_t%g_pid%d", t, pid());
  FILE *fp1 = fopen(tmp_Wavefile, "w");
  sprintf (Wavefile, "Waveflow_nu%g_eps%g_omega%g_angle%g_LEVEL%d_l%g_delta%g_num%g_vert%d_tanh%d_t%g", nu*1e6, eps, ratio_f, angle_eq, LEVEL, l_cl, _pf, _num_slip, _vert, _tanh, (t - pi/2.)/(2.*pi));
  foreach(){
    if ((z > (Z_MEASUREMENTS - Delta/2.)) & (z < (Z_MEASUREMENTS + Delta/2.)) & (cs[] > 0.)){
        fprintf(fp1, "%-12g %-12g %-12g %-12g %-12g %-12g\n", x, y, z, u.x[] + us.x, u.y[] + us.y, u.z[] + us.z);
    }
  }
  fclose(fp1);
  MPI_Barrier(MPI_COMM_WORLD);
  if (pid() == 0) {
    char command1[1000];
    sprintf(command1, "cat tmp_Waveflow_t%g_pid* > %s", t, Wavefile);
    system(command1);

    char command2[200];
    sprintf(command2, "rm tmp_Waveflow*");
    system(command2);
  }
  MPI_Barrier(MPI_COMM_WORLD);
}

event snapshot(t = 2.*pi, t += 2.*pi, t <= TEND)
{
    char name[400];
    sprintf (name, "snap_nu%g_eps%g_omega%g_angle%g_LEVEL%d_l%g_delta%g_num%g_vert%d_tanh%d_t%g", 
      nu*1e6, eps, ratio_f, angle_eq, LEVEL, l_cl, _pf, _num_slip, _vert, _tanh, t/(2.*pi));
    // scalar pid[];
    // foreach()
    //   pid[] = fmod(pid()*(npe() + 37), npe());
    dump(name);
}