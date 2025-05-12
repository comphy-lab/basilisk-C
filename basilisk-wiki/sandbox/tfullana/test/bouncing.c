/* Flags */
#define SAVEDUMP 1
#define SAVEFACETS 1
#define SAVEMOVIE 1
#define HARMONICS 1

/* Header files */
#define mu(f)  (1./(clamp(f,0,1)*(1./mu1 - 1./mu2) + 1./mu2))
#include "axi.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "navier-stokes/conserving.h"
#include "tension.h"
#include "view.h"

/* --------------------------- 
Dimensional Oil/Air properties (SI)
--------------------------- */

double rho_1 = 0.95e3 [1, -3, 0]; // Oil density
double mu_1 = 20e-3 [1, -1, -1]; // Oil viscosity
double rho_2 = 1.2047 [1, -3, 0]; // Air density
double mu_2 = 1.8205e-5 [1, -1, -1]; // Air viscosity
double sig = 21e-3 [1, 0, -2]; // Oil to Air surface tension
double grav = 9.8067 [0, 1, -2]; // Gravity
double initial_width = 1.52e-3 [0, 1, 0], initial_height = 1.52e-3 [0, 1, 0];
double radius = 0.76e-3 [0, 1, 0];
#define INITIAL_WIDTH initial_width/radius
#define INITIAL_HEIGHT initial_height/radius

/* --------------------------- 
Time and velocity scalings
--------------------------- */

#define T_unit sqrt((rho_1*cube(radius))/sig)
#define U_unit sqrt(sig/(rho_1*radius))

/* --------------------------- 
Dimensionless numbers
--------------------------- */

#define Bo (rho_1*grav*sq(radius)/sig)
#define Bo2 (rho_2*grav*sq(radius)/sig)
#define Oh (mu_1/sqrt(sig*rho_1*radius))
#define Oh2 (mu_2/sqrt(sig*rho_1*radius))
#define We (rho_1*radius*sq(U_unit)/sig)

#define Bo_nd (1.*Bo*sq(1.)/1.)
#define Bo2_nd ((rho_2/rho_1)*Bo*sq(1.)/1.)
#define Oh_nd (Oh/sqrt(1.*1.*1.))
#define Oh2_nd ((mu_2/mu_1)*Oh/sqrt(1.*1.*1.))
#define We_nd (1.*1.*sq(1.)/1.)


/* --------------------------- 
Stage properties
--------------------------- */

double freq = 40.0 [0, 0, -1], amp = 0.112e-3 [0, 1, 0];
#define OMEGA (freq*(2.*pi*T_unit))
#define GAMMA (amp*sq(2.*pi*freq)/grav)
#define period (1./freq)

#define stage_pos(t) (GAMMA/sq(OMEGA))*grav*sq(T_unit)*cos(OMEGA*(t) + phase*pi)                               
#define stage_vel(t) -(GAMMA/OMEGA)*grav*T_unit*sin(OMEGA*(t) + phase*pi)                
#define stage_acc(t) -GAMMA*grav*cos(OMEGA*(t) + phase*pi) 

#define STAGE_POS(t) (GAMMA/sq(OMEGA))*Bo*cos(OMEGA*(t) + phase*pi)                              
#define STAGE_VEL(t) -(GAMMA/OMEGA)*Bo*sin(OMEGA*(t) + phase*pi)                
#define STAGE_ACC(t) -GAMMA*Bo*cos(OMEGA*(t) + phase*pi) 


/* --------------------------- 
Final and output times
--------------------------- */

double tend = 0.3 [0, 0, 1];
double tout = 1./3000. [0, 0, 1], tdump = 1./3000. [0, 0, 1], tmovie = 1./3000. [0, 0, 1]; 

#define TEND (tend/T_unit)     
#define TOUT (tout/T_unit)     
#define TDUMP (tdump/T_unit)   
#define TMOVIE (tmovie/T_unit) 


/* --------------------------- 
Domain size and initial conditions
--------------------------- */

double L1 = 8. [0];      
double initial_pos = 2.337443255e-3 [0, 1, 0]; 
double initial_vel = -0.044039997 [0, 1, -1]; 
#define INITIAL_POS ((initial_pos/radius) + STAGE_POS(0))     
#define INITIAL_VEL ((initial_vel/U_unit) + STAGE_VEL(0)) 
double phase = 0.000 [0];

/* --------------------------- 
Energy budget
--------------------------- */

double viscousD1 = 0. [0], viscousD2 = 0. [0], stage_energy = 0. [0], reaction_force = 0. [0], E0 = 0. [0], Etot = 0. [0];


/* --------------------------- 
Adaptivity parameters
--------------------------- */

double u_err = 1e-4 [0, 1, -1], f_err = 1e-3; 
int LEVEL = 8;                     


/* --------------------------- 
Boundary conditions
--------------------------- */

uf.n[bottom] = 0.;

p[right] = dirichlet(0);
pf[right] = dirichlet(0);

u.t[left] = dirichlet(0);
u.n[left] = dirichlet(0);

f[left] = 0.;

int main(int argc, char *argv[])
{
  if (argc > 1)
    radius = atof(argv[1])*1e-3;
  if (argc > 2)
    amp = atof(argv[2])*1e-3;
  if (argc > 3)
    phase = atof(argv[3]);
  if (argc > 4)
    LEVEL = atoi(argv[4]);

  size(L1);

  origin(0., 0.);

  init_grid(1 << LEVEL);

  DT = HUGE [0];

  rho1 = 1. [0], mu1 = Oh;
  rho2 = (rho_2/rho_1), mu2 = Oh2;

  f.sigma = 1. [0];
  TOLERANCE = 1e-4 [*];
  run();
}

event init(t = 0)
{
  char name[360];
  sprintf(name, "properties_pos%.2f_radius%.2f_OMEGA%.2f_GAMMA%.2f_freq%.2f_amp%.2f_phi%.3f_tend%.2f_LEVEL%d",
          initial_pos*1e3, radius*1e3, OMEGA, GAMMA, freq, amp*1e6, phase, tend, LEVEL);
  FILE *fp = fopen(name, "w");

  fprintf(fp, "%-12s %-12s %-12s %-12s %-12s\n", "Bo", "Oh", "We", "Bo_air", "Oh_air");
  fprintf(fp, "%-12.4f %-12.4f %-12.4f %-12.4f %-12.4f\n\n", Bo, Oh, We, Bo2, Oh2);
  fprintf(fp, "%-12s %-12s %-12s %-12s %-12s\n", "Bo_nd", "Oh_nd", "We_nd", "Bo_air_nd", "Oh_air_nd");
  fprintf(fp, "%-12.4f %-12.4f %-12.4f %-12.4f %-12.4f\n\n", Bo_nd, Oh_nd, We_nd, Bo2_nd, Oh2_nd);
  fprintf(fp, "%-12s %-12s\n", "T_unit (s)", "U_unit (m/s)");
  fprintf(fp, "%-12.2e %-12.2e\n\n", T_unit, U_unit);
  fprintf(fp, "%-12s %-12s\n", "Domain (m)", "Domain_nd");
  fprintf(fp, "%-12.2e %-12.2f\n\n", radius*L1, L1);
  fprintf(fp, "%-12s %-12s\n", "Grid (m)", "Grid_nd");
  fprintf(fp, "%-12.2e %-12.2e\n\n", L1*radius/pow(2, LEVEL), L1/pow(2, LEVEL));

  fclose(fp);

  fraction(f, -sq((x - INITIAL_POS)/(INITIAL_HEIGHT/2.)) - sq(y/(INITIAL_WIDTH/2.)) + sq(1.));

  foreach ()
    u.x[] = f[]*INITIAL_VEL;

  system("mkdir -p energies");
  system("mkdir -p positions");
  system("mkdir -p velocities");
  system("mkdir -p volumes");
  system("mkdir -p facets");
  system("mkdir -p dumps");
  system("mkdir -p harmonics");
}

event adapt(i++)
{
  adapt_wavelet({f, u}, (double[]){f_err, u_err/U_unit, u_err/U_unit, u_err/U_unit}, LEVEL);
}

event acceleration(i++)
{
  face vector av = a;
  foreach_face(x)
      av.x[] -= (Bo - STAGE_ACC(t));

  // Stops if the interface exits the domain
	foreach_boundary(right){
		if (f[] != 0)
			return 1.;
	}
}

event compute_time_integrals(i++)
{ 
  double vd1 = 0., vd2 = 0., instantaneous_stage_energy = 0.;
  reaction_force = 0.;

  foreach (reduction(+:vd1) reduction(+:vd2))
  {
    double dv1 = f[]*dv();
    double dv2 = (1. - f[])*dv();
    double vd = 0.;
    double dudr = sq(u.y[0,1] - u.y[0,-1])/sq(2.*Delta);
    double dudy = sq(u.y[1,0] - u.y[-1,0])/sq(2.*Delta);
    double dvdr = sq(u.x[0,1] - u.x[0,-1])/sq(2.*Delta);
    double dvdy = sq(u.x[1,0] - u.x[-1,0])/sq(2.*Delta);
    double ur = sq(u.y[]/y);
    vd = 2.*(dudr + dvdy + ur) + (dvdr + dudy);
    vd1 += mu1*dv1*vd*dt;
    vd2 += mu2*dv2*vd*dt;
  }

  foreach_boundary(left, reduction(+:reaction_force) reduction(+:instantaneous_stage_energy))
  {
    if (f[] > 1e-6 && f[] < 1. - 1e-6)
    {
      coord n = interface_normal(point, f), p_;
      double alpha = plane_alpha(f[], n);
      double area_ = 2.*pi*y*pow(Delta, dimension - 1)*plane_area_center(n, alpha, &p_);
      reaction_force += -p[]*area_*n.x;
      instantaneous_stage_energy += p[]*area_*n.x*(STAGE_VEL(t))*dt; 
    }
  }

  stage_energy += instantaneous_stage_energy;
  viscousD1 += vd1;
  viscousD2 += vd2;

}

event logfile(t += TOUT, t <= TEND)
{
  char datafile1[360], datafile2[360], datafile3[360], datafile4[360], datafile5[360], datafile6[360], datafile7[360];
  double xc = 0., yc = 0., uc = 0., vc = 0., wd = 0., hd = 0.;
  double vol1 = 0., ke1 = 0., gpe1 = 0.;
  double vol2 = 0., ke2 = 0., gpe2 = 0.;
  double x2 = 0., x1 = 0., y1 = 0., xmin = 100.;
  double area = 0.;

  coord us = {STAGE_VEL(t), 0., 0.};

  sprintf(datafile1, "./positions/POS_pos%.2f_radius%.2f_OMEGA%.2f_GAMMA%.2f_freq%.2f_amp%.2f_phi%.3f_tend%.2f_LEVEL%d",
          initial_pos*1e3, radius*1e3, OMEGA, GAMMA, freq, amp*1e6, phase, tend, LEVEL);
  sprintf(datafile2, "./positions/pos_pos%.2f_radius%.2f_OMEGA%.2f_GAMMA%.2f_freq%.2f_amp%.2f_phi%.3f_tend%.2f_LEVEL%d",
          initial_pos*1e3, radius*1e3, OMEGA, GAMMA, freq, amp*1e6, phase, tend, LEVEL);
  sprintf(datafile3, "./velocities/VEL_pos%.2f_radius%.2f_OMEGA%.2f_GAMMA%.2f_freq%.2f_amp%.2f_phi%.3f_tend%.2f_LEVEL%d",
          initial_pos*1e3, radius*1e3, OMEGA, GAMMA, freq, amp*1e6, phase, tend, LEVEL);
  sprintf(datafile4, "./velocities/vel_pos%.2f_radius%.2f_OMEGA%.2f_GAMMA%.2f_freq%.2f_amp%.2f_phi%.3f_tend%.2f_LEVEL%d",
          initial_pos*1e3, radius*1e3, OMEGA, GAMMA, freq, amp*1e6, phase, tend, LEVEL);
  sprintf(datafile5, "./energies/ENERGY_pos%.2f_radius%.2f_OMEGA%.2f_GAMMA%.2f_freq%.2f_amp%.2f_phi%.3f_tend%.2f_LEVEL%d",
          initial_pos*1e3, radius*1e3, OMEGA, GAMMA, freq, amp*1e6, phase, tend, LEVEL);
  sprintf(datafile6, "./energies/energy_pos%.2f_radius%.2f_OMEGA%.2f_GAMMA%.2f_freq%.2f_amp%.2f_phi%.3f_tend%.2f_LEVEL%d",
          initial_pos*1e3, radius*1e3, OMEGA, GAMMA, freq, amp*1e6, phase, tend, LEVEL);
  sprintf(datafile7, "./volumes/volume_pos%.2f_radius%.2f_OMEGA%.2f_GAMMA%.2f_freq%.2f_amp%.2f_phi%.3f_tend%.2f_LEVEL%d",
          initial_pos*1e3, radius*1e3, OMEGA, GAMMA, freq, amp*1e6, phase, tend, LEVEL);

  static FILE *fp1 = fopen(datafile1, "w");
  static FILE *fp2 = fopen(datafile2, "w");
  static FILE *fp3 = fopen(datafile3, "w");
  static FILE *fp4 = fopen(datafile4, "w");
  static FILE *fp5 = fopen(datafile5, "w");
  static FILE *fp6 = fopen(datafile6, "w");
  static FILE *fp7 = fopen(datafile7, "w");

  if (t == 0)
  {
    fprintf(fp1, "%-12s\n", "#data file non_dimensional unit");
    fprintf(fp1, "%-12s %-12s %-12s %-12s %-12s %-12s %-12s\n\n",
            "#1 time", "#2 x_stage", "#3 x_center", "#4 y_center", "#5 width", "#6 height", "#7 air layer height");

    fprintf(fp2, "%-12s\n", "#data file dimensional unit (SI)");
    fprintf(fp2, "%-12s %-12s %-12s %-12s %-12s %-12s %-12s\n\n",
            "#1 time", "#2 x_stage", "#3 x_center", "#4 y_center", "#5 width", "#6 height", "#7 air layer height");

    fprintf(fp3, "%-12s\n", "#data file non_dimensional unit");
    fprintf(fp3, "%-12s %-12s %-12s %-12s %-12s %-12s\n\n",
            "#1 time", "#2 ux_stage", "#3 gx_stage", "#4 ux_center", "#5 uy_center", "#6 Webber");

    fprintf(fp4, "%-12s\n", "#data file dimensional unit (SI)");
    fprintf(fp4, "%-12s %-12s %-12s %-12s %-12s %-12s\n\n",
            "#1 time", "#2 ux_stage", "#3 gx_stage", "#4 ux_center", "#5 uy_center", "#6 Webber");

    fprintf(fp5, "%-12s\n", "#data file non_dimensional unit");
    fprintf(fp5, "%-12s %-12s %-12s %-12s %-12s %-12s %-12s %-12s %-12s %-12s %-12s\n\n",
            "#1 time", "#2 kinetic_1", "#3 kinetic_2", "#4 potential_1", "#5 potential_2", "#6 viscousD_1", "#7 viscousD_2", "#8 Interfacial Energy", "#9 stage energy", "#10 reaction force", "#11 Etot");

    fprintf(fp6, "%-12s\n", "#data file dimensional unit (SI)");
    fprintf(fp6, "%-12s %-12s %-12s %-12s %-12s %-12s %-12s %-12s %-12s %-12s %-12s\n\n",
            "#1 time", "#2 kinetic_1", "#3 kinetic_2", "#4 potential_1", "#5 potential_2", "#6 viscousD_1", "#7 viscousD_2", "#8 Interfacial Energy", "#9 stage energy", "#10 reaction force", "#11 Etot");

    fprintf(fp7, "%-12s\n", "#Volume conservation");
    fprintf(fp7, "%-12s %-12s %-12s %-12s %-12s %-12s %-12s\n\n",
            "#1 time", "#2 area", "#3 vol1", "#4 vol2", "#5 time_dim", "#6 vol1_dim", "#7 vol2_dim");
  }

  foreach (reduction(+:xc) reduction(+:yc) reduction(+:uc) reduction(+:vc)
    reduction(+:vol1) reduction(+:ke1) reduction(+:gpe1)
    reduction(+:vol2) reduction(+:ke2) reduction(+:gpe2)
    reduction(+:area))
  {
    double dv1 = f[]*dv();
    double dv2 = (1. - f[])*dv();
    double norm2 = 0.;
    foreach_dimension()
    {
      norm2 += sq(u.x[] - us.x);
    }
    uc += (u.x[] - us.x)*dv1;
    vc += u.y[]*dv1;
    xc += x*dv1;
    yc += y*dv1;
    vol1 += dv1;
    ke1 += rho1*norm2*dv1;
    gpe1 += rho1*Bo*(x - STAGE_POS(t))*dv1;
    vol2 += dv2;
    ke2 += rho2*norm2*dv2;
    gpe2 += rho2*Bo*(x - STAGE_POS(t))*dv2;
    if (f[] > 1e-6 && f[] < 1. - 1e-6)
      {
        coord n = interface_normal(point, f), p_;
        double alpha = plane_alpha(f[], n);
        area += 2.*pi*y*pow(Delta, dimension - 1)*plane_area_center(n, alpha, &p_);
      }
  }

  ke1 /= 2.;
  ke2 /= 2.;

  if (t == 0)
    E0 = ke1 + ke2 + gpe1 + 1.*4.*pi*sq(1.);

  vector h[];
  heights(f, h);

  foreach (reduction(max:y1) reduction(min:xmin))
  {
    if ((h.y[] != nodata) & (f[] < (1. - 1e-6)) & (f[] > (1e-6)) & (x < (xc/vol1 + Delta)) & (x > (xc/vol1 - Delta)))
    {
      double pos = y + Delta*height(h.y[]);
      if (pos > yc/vol1)
        y1 = pos;
    }
    if ((h.x[] != nodata) & (f[] < (1. - 1e-6)) & (f[] > (1e-6)))
    {
      double pos = x + Delta*height(h.x[]);
      if (pos < xmin)
        xmin = pos;
    }
  }
  foreach_boundary(bottom, reduction(max:x1) reduction(max:x2))
  {
    if ((h.x[] != nodata) & (f[] < (1. - 1e-6)) & (f[] > (1e-6)))
    {
      double pos = x + Delta*height(h.x[]);
      if (pos > xc/vol1)
        x1 = pos;
      if (pos < xc/vol1)
        x2 = pos;
    }
  }

  wd = 2.*fabs(y1);
  hd = fabs(x1 - x2);
  
  Etot = ke1 + gpe1 + 1.*area + viscousD1 + viscousD2 - stage_energy;

  fprintf(fp1, "%-12g %-12g %-12g %-12g %-12g %-12g %-12g\n",
          t, STAGE_POS(t), ((xc/vol1)), yc/vol1, wd, hd, xmin);
  fprintf(fp2, "%-12g %-12g %-12g %-12g %-12g %-12g %-12g\n",
          t*T_unit, stage_pos(t), radius*((xc/vol1)), radius*yc/vol1, radius*wd, radius*hd, radius*xmin);

  fprintf(fp3, "%-12g %-12g %-12g %-12g %-12g %-12g\n",
          t, STAGE_VEL(t), STAGE_ACC(t), uc/vol1, vc/vol1, rho1*sq(uc/vol1)*1./1.);
  fprintf(fp4, "%-12g %-12g %-12g %-12g %-12g %-12g\n",
          t*T_unit, stage_vel(t), stage_acc(t), U_unit*uc/vol1, U_unit*vc/vol1, rho_1*sq(U_unit*uc/vol1)*radius/sig);

  fprintf(fp5, "%-12g %-12g %-12g %-12g %-12g %-12g %-12g %-12g %-12g %-12g %-12g\n",
          t, ke1/E0, ke2/E0, gpe1/E0, gpe2/E0, viscousD1/E0, viscousD2/E0, 1.*area/E0, stage_energy/E0, reaction_force, Etot/E0);
  fprintf(fp6, "%-12g %-12g %-12g %-12g %-12g %-12g %-12g %-12g %-12g %-12g %-12g\n",
          t*T_unit, ke1, ke2, gpe1, gpe2, viscousD1, viscousD2, 1.*area, stage_energy, reaction_force, Etot);

  fprintf(fp7, "%-12g %-12g %-12g %-12g %-12g %-12g %-12g\n",
          t, area, 2.*pi*vol1, 2.*pi*vol2, t*T_unit, 2.*pi*vol1*cube(radius), 2.*pi*vol2*cube(radius));

  fflush(fp1);
  fflush(fp2);
  fflush(fp3);
  fflush(fp4);
  fflush(fp5);
  fflush(fp6);
  fflush(fp7);
}


#if SAVEMOVIE
event moviefile(t += TMOVIE)
{
  char vofmovie[360];

  sprintf(vofmovie, "vof_pos%.2f_radius%.2f_OMEGA%.2f_GAMMA%.2f_freq%.2f_amp%.2f_phi%.3f_tend%.2f_LEVEL%d.mp4",
          initial_pos*1e3, radius*1e3, OMEGA, GAMMA, freq, amp*1e6, phase, tend, LEVEL);

  view (quat = {0.000, 0.000, -0.707, 0.707},
      fov = 30, near = 0.01, far = 1000,
      tx = -0.006, ty = -0.451, tz = -1.743,
      width = 1024, height = 1024);
  clear();
  squares("u.x");
  draw_vof(c = "f", fc = {0, 0, 0}, lc = {0, 0, 0}, lw = 2);
  mirror ({0,1}) {
      squares ("f");
      cells();
  }
  save(vofmovie);
}
#endif

#if HARMONICS
// Spherical harmonics
#define l_max 10   // Number of coefficients to compute
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_sf_legendre.h>

void read_facet_file(const char* filename, double** x, double** y, int* num_points) {
    FILE *file = fopen(filename, "r");
    if (!file) {
        perror("Failed to open file");
        exit(EXIT_FAILURE);
    }

    // Count number of lines/points
    int count = 0;
    double temp_x, temp_y;
    while (fscanf(file, "%lf %lf", &temp_x, &temp_y) == 2) {
        count++;
    }
    rewind(file);

    // Allocate memory for x and y arrays
    *x = (double*)malloc(count*sizeof(double));
    *y = (double*)malloc(count*sizeof(double));
    *num_points = count;

    // Read the file data into arrays
    int index = 0;
    while (fscanf(file, "%lf %lf", &(*x)[index], &(*y)[index]) == 2) {
        index++;
    }

    fclose(file);
}

void compute_centroid(double* x, double* y, int num_points, double* cx, double* cy) {
    *cx = 0;
    *cy = 0;
    for (int i = 0; i < num_points; i++) {
        *cx += x[i];
        *cy += y[i];
    }
    *cx /= num_points;
    *cy /= num_points;
}

event dumpfile(t += TDUMP)
{
  if (pid() == 0) {
  char dumpname[360], facetname[360], name[360], datafile1[360], datafile2[360];

  sprintf(dumpname, "./dumps/dump_pos%.2f_radius%.2f_OMEGA%.2f_GAMMA%.2f_freq%.2f_amp%.2f_phi%.3f_tend%.2f_LEVEL%d_t%.0f",
          initial_pos*1e3, radius*1e3, OMEGA, GAMMA, freq, amp*1e6, phase, tend, LEVEL, t/TDUMP);
  sprintf(facetname, "./facets/points_pos%.2f_radius%.2f_OMEGA%.2f_GAMMA%.2f_freq%.2f_amp%.2f_phi%.3f_tend%.2f_LEVEL%d_t%.4f",
          initial_pos*1e3, radius*1e3, OMEGA, GAMMA, freq, amp*1e6, phase, tend, LEVEL, t*T_unit);
  sprintf(name, "interfaces_t%.4f_pid%d", t*T_unit, pid());
  sprintf(datafile1, "./harmonics/HAR_pos%.2f_radius%.2f_OMEGA%.2f_GAMMA%.2f_freq%.2f_amp%.2f_phi%.3f_tend%.2f_LEVEL%d",
          initial_pos*1e3, radius*1e3, OMEGA, GAMMA, freq, amp*1e6, phase, tend, LEVEL);
  sprintf(datafile2, "./harmonics/har_pos%.2f_radius%.2f_OMEGA%.2f_GAMMA%.2f_freq%.2f_amp%.2f_phi%.3f_tend%.2f_LEVEL%d",
          initial_pos*1e3, radius*1e3, OMEGA, GAMMA, freq, amp*1e6, phase, tend, LEVEL);
#if SAVEDUMP
  dump(dumpname);
#endif
  FILE *fp1 = fopen(name, "w");
  static FILE *fp2 = fopen(datafile1, "w");
  static FILE *fp3 = fopen(datafile2, "w");

  if (t == 0){
    fprintf(fp2, "%-12s\n", "#data file non_dimensional unit");
    fprintf(fp2, "%-12s %-12s %-12s %-12s %-12s %-12s %-12s %-12s %-12s %-12s %-12s %-12s\n\n",
            "#1 time", "#2 l = 0", "#2 l = 1", "#2 l = 2", "#2 l = 3", "#2 l = 4", "#2 l = 5", "#2 l = 6", "#2 l = 7", "#2 l = 8", "#2 l = 9", "#2 l = 10");

    fprintf(fp3, "%-12s\n", "#data file dimensional unit (SI)");
    fprintf(fp3, "%-12s %-12s %-12s %-12s %-12s %-12s %-12s %-12s %-12s %-12s %-12s %-12s\n\n",
            "#1 time", "#2 l = 0", "#2 l = 1", "#2 l = 2", "#2 l = 3", "#2 l = 4", "#2 l = 5", "#2 l = 6", "#2 l = 7", "#2 l = 8", "#2 l = 9", "#2 l = 10");
  }

  foreach(){
    if (f[] > 1e-6 && f[] < 1. - 1e-6)
      {
        coord n = interface_normal(point, f), p;
        double alpha = plane_alpha(f[], n);
        plane_area_center(n, alpha, &p);
        fprintf (fp1, "%g %g\n", x + Delta*p.x, y + Delta*p.y);
        fprintf (fp1, "%g %g\n", x + Delta*p.x, -(y + Delta*p.y)); // Plane symmetry in y
      }
  }

  fclose(fp1);

  if (pid() == 0) {
      char command[1000];
      sprintf(command, "cat interfaces_t%.4f_pid* > %s", t*T_unit, facetname);
      system(command);

      char command2[200];
      sprintf(command2, "rm interfaces*");
      system(command2);
  }

  double *xp, *yp;
  int num_points;

  // Read the facet file
  read_facet_file(facetname, &xp, &yp, &num_points);

  double cx, cy;
  compute_centroid(xp, yp, num_points, &cx, &cy);

  gsl_matrix *A = gsl_matrix_alloc(num_points, l_max);
  gsl_vector *b = gsl_vector_alloc(num_points);
  gsl_vector *x = gsl_vector_alloc(l_max);

  // Fill the matrix A and vector b
  for (int i = 0; i < num_points; ++i) {
      double r = sqrt(sq(xp[i] - cx) + sq(yp[i] - cy)) - 1.;
      gsl_vector_set(b, i, r);
      double theta = atan2((xp[i] - cx), (yp[i] - cy)) + pi/2.;  // Angle from y-axis
      for (int j = 0; j < l_max; ++j) {
          double P_lj = gsl_sf_legendre_Pl(j, cos(theta));
          gsl_matrix_set(A, i, j, P_lj);
      }
  }

  free(xp);
  free(yp);
  // Solve the linear system A*x = b
  gsl_matrix *AT = gsl_matrix_alloc(l_max, num_points);
  gsl_matrix_transpose_memcpy(AT, A);
  gsl_matrix *ATA = gsl_matrix_alloc(l_max, l_max);
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, AT, A, 0.0, ATA);

  gsl_vector *ATb = gsl_vector_alloc(l_max);
  gsl_blas_dgemv(CblasNoTrans, 1.0, AT, b, 0.0, ATb);

  gsl_permutation *p = gsl_permutation_alloc(l_max);
  int s;
  gsl_linalg_LU_decomp(ATA, p, &s);
  gsl_linalg_LU_solve(ATA, p, ATb, x);

  double COEF[l_max], coef[l_max];
  // // Print the coefficients
  for (int i = 0; i < l_max; ++i) {
      COEF[i] = gsl_vector_get(x, i);
      coef[i] = radius*gsl_vector_get(x, i);
  }

  fprintf (fp2, "%g %g %g %g %g %g %g %g %g %g %g\n", t, COEF[0], COEF[1], COEF[2], COEF[3], COEF[4], COEF[5], COEF[6], COEF[7], COEF[8], COEF[9]);
  fprintf (fp3, "%g %g %g %g %g %g %g %g %g %g %g\n", t*T_unit, coef[0], coef[1], coef[2], coef[3], coef[4], coef[5], coef[6], coef[7], coef[8], coef[9]);

  fflush(fp2);
  fflush(fp3);

  // Free the allocated memory
  gsl_matrix_free(A);
  gsl_matrix_free(AT);
  gsl_matrix_free(ATA);
  gsl_vector_free(b);
  gsl_vector_free(ATb);
  gsl_vector_free(x);
  gsl_permutation_free(p);
#if !SAVEFACETS
  if (pid() == 0) {
      char command3[200];
      sprintf(command3, "rm ./facets/points*");
      system(command3);
  }
#endif
  }
}
#endif