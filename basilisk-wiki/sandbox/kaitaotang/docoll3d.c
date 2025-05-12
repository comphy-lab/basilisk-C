#include "grid/octree.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "navier-stokes/conserving.h"
#include "tension.h"
#include "tag.h"
#include "view.h"

/**
   Include profiling information. */

#include "navier-stokes/perfs.h"
#include "profiling.h"

#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_fft_real.h>
#include <gsl/gsl_fft_halfcomplex.h>


#include <sys/stat.h>
#include <sys/types.h>

#define LEVEL 11
#define MINLEVEL 7
#define MAXTIME 0.03
//#define R0 0.05 // Initial bubble radius
#define R0 0.05
#define WE 280. //Weber number
//#define OH 0.14 //Ohnesorge number
#define OH 0.01 // Try this
#define OFFSET 0.95
#define ETA 0.06
#define RATIO 1./833.
//#define MURATIO 17.4e-6/8.9e-4 //dynamic viscosity ratio, water to air
#define MURATIO 1./55. //viscosity ratio water to air, not chosen according to Pierson but to previous bubble-in-turbulence studies; expect it will not have much effect, especially for values below 0.1
#define ADDNOISE 1 //boolean noise 
#define CONSTANT 0.01 // noise scaling

double raw_series_l[4096]; double raw_series_r[4096]; 
double raw_series_l_old[4096]; double raw_series_r_old[4096]; 
int series_len;

u.n[top] = neumann(0.);
u.t[top] = neumann(0.);
u.r[top] = neumann(0.);

FILE * fk = NULL;
FILE * fe = NULL;


//----------------------MAIN-----------------------//
//Begin main function
int main() {
  L0 = 2.;
  periodic(front);
  origin (-L0/2., 0., 0.);
  rho1 = 1.0;
  rho2 = RATIO;
  /**
     For calculating viscosities, interpret Reynolds number as defined at depth of unity. */
  f.sigma = 1.; //0.01 is very stable
  mu1 = OH*sqrt(rho1*f.sigma*2*R0); //double check this definition
  mu2 = mu1*MURATIO;
  /**
     Use depth length scale for Bond number as well. */
  //TOLERANCE = HUGE;
  //NITERMIN = 20;
  
  mkdir("./interface", 0777);
  mkdir("./frag_stats", 0777);
  
  char name[80];
  sprintf (name, "3Denergy-%d.dat", LEVEL);
  fe = fopen (name, "w");

  sprintf (name, "3Dkinetics-%d.dat", LEVEL);
  fk = fopen (name, "w");
  
  series_len = (1 << LEVEL); 
  
  init_grid(1 << (MINLEVEL+1));
  run();
  
  fclose (fe); fclose (fk);
}

//---------------------INITIALIZATION------------------------//
void set_noise_signal() {
  gsl_fft_real_wavetable * real;
  gsl_fft_halfcomplex_wavetable * hc;
  gsl_fft_real_workspace * work;

  if (pid() == 0) {
    srand(time(NULL));
    for (int i=0; i<series_len; i++) {
      raw_series_l[i] = (double)rand() / (double)RAND_MAX - 0.5;
      raw_series_r[i] = (double)rand() / (double)RAND_MAX - 0.5;
      raw_series_l_old[i] = raw_series_l[i]; raw_series_r_old[i] = raw_series_r[i]; 
    }
    
    work = gsl_fft_real_workspace_alloc(series_len);
    real = gsl_fft_real_wavetable_alloc(series_len);

    gsl_fft_real_transform (raw_series_l, 1, series_len, real, work);
    gsl_fft_real_transform (raw_series_r, 1, series_len, real, work);
    gsl_fft_real_wavetable_free (real);

    for (int i = 121; i<series_len; i++) {raw_series_l[i] = 0; raw_series_r[i] = 0;}
    for (int i = 1; i<120; i++) if (i % 2 == 1) {raw_series_l[i] = 0; raw_series_r[i] = 0;}

    hc = gsl_fft_halfcomplex_wavetable_alloc (series_len);

    gsl_fft_halfcomplex_inverse (raw_series_l, 1, series_len, hc, work);
    gsl_fft_halfcomplex_inverse (raw_series_r, 1, series_len, hc, work);
    gsl_fft_halfcomplex_wavetable_free (hc);
    gsl_fft_real_workspace_free (work);
  
  
    double series_ave_l = 0., series_var_l = 0., series_ave_r = 0., series_var_r = 0.;
    for (int i=0; i<series_len; i++) {
      series_ave_l += raw_series_l[i]/(double)(series_len); series_var_l += sq(raw_series_l[i])/(double)(series_len);
      series_ave_r += raw_series_r[i]/(double)(series_len); series_var_r += sq(raw_series_r[i])/(double)(series_len);
    }
  
    series_var_l -= sq(series_ave_l); series_var_r -= sq(series_ave_r);
    
    for (int i=0; i<series_len; i++) {
      raw_series_l[i] = (raw_series_l[i] - series_ave_l) / sqrt(series_var_l);
      raw_series_r[i] = (raw_series_r[i] - series_ave_r) / sqrt(series_var_r);
    }
  }
  
  MPI_Bcast (&(raw_series_l[0]), series_len, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast (&(raw_series_r[0]), series_len, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  
  char resultname[100];
  sprintf (resultname, "perturbation_series.dat");
  
  if (pid() == 0) {
    FILE * fseries = fopen(resultname, "w");
    
    for (int i=0; i<series_len; i++) {
      fprintf(fseries, "%d %g %g %g %g\n", i, raw_series_l_old[i], raw_series_l[i], raw_series_r_old[i], raw_series_r[i]);
      fflush(fseries);
    }
    
    fclose(fseries);
  }
}



double dropletValue(double x, double y, double z, double amp){
//  double sgn = -1.;
//  double sum = 0.; // initializing summation

  // wavenumber for loop
  // for (int k_i = 0; k_i < 25; k_i++) {
    //double a_i = forcing_const/k_i;
    //sum = sum + noise_sig_coeff[k_i] * sin(2*pi/L0 * (k_i+1)*z);
  //}
  //double amp = sum; // assigning the amplitude to the summation
  double stp = L0/(double)(series_len);
  int index = (int)((z - stp/2. - Z0)/stp);
  if (index < 0) index = 0; if (index > series_len) index = series_len;
  
  if (x < 0.) return -sq(x + R0*OFFSET) - sq(y) + sq(R0 + amp*R0*raw_series_l[index]);
    else return -sq(x - R0*OFFSET) - sq(y) + sq(R0 + amp*R0*raw_series_r[index]);
}

//
/*Write initialization event*/
event init (i = 0)
{
  if (!restore("restart")){
    set_noise_signal();
//    int jdx = 0;
    double dropvel = sqrt(WE*f.sigma/(2.*R0)/rho1)/2.;
//    double femax = 2e-4;
//    double uemax = 2e-2;
//    do {
//      jdx += 1;
//      fprintf(ferr, "Initializing, %d\n", jdx );
      //fraction (f, -sq(x - R0*OFFSET) - sq(y) + sq(R0));
      refine (fabs(x) < 2.2*R0 && fabs(y) < 1.1*R0 && level < (LEVEL-1));
      fraction (f, dropletValue(x, y, z, ETA));
      boundary ((scalar *){f,u});
      foreach()
        u.x[] = -f[]*dropvel*sign(x);
//    }
//    while (adapt_wavelet ({f,u}, (double[]){femax,uemax, uemax, uemax}, LEVEL, MINLEVEL).nf);
  }
}

//-------------------ADAPTIVITY---------------------//
/*Adapt once on error in volume fraction, velocity field, and beach fraction*/
event adapt(i++) {
  double femax = 1e-2;
  double uemax = 5e-2;
  adapt_wavelet ({f,u.x,u.y,u.z}, (double[]){femax,uemax,uemax,uemax}, LEVEL, 5);
}

//------------------------IMAGING--------------------------
/*
event intoutput(t += 0.0001) {
  char resultname[100];
  char tmpname[100];
  sprintf( tmpname, "interface/_res-%d.dat", pid() );
  if (pid() == 0)
    sprintf (resultname, "interface/results%4.4f.dat", t );
  
  FILE * fp = fopen(tmpname, "w");
  scalar xpos[];
  scalar ypos[];
  scalar zpos[];
  
  position (f, xpos, {1, 0, 0});
  position (f, ypos, {0, 1, 0});
  position (f, zpos, {0, 0, 1});
  
  foreach()
    if (xpos[] != nodata)
      fprintf (fp, "%g %g %g\n", xpos[], ypos[], zpos[]);
  fclose (fp);
  
  char copycat[200];
  if (pid() == 0){
    sprintf (copycat, "cat interface/_res-*.dat > %s", resultname);
    system (copycat);
    system ("rm interface/_res-*.dat");
  }   
}
*/



int dissipation_rate (vector u, double* rates)
{
  double rateWater = 0.0;
  double rateAir = 0.0;
  foreach (reduction (+:rateWater) reduction (+:rateAir)) {
    double dudx = (u.x[1] - u.x[-1])/(2.0*Delta);
    double dudy = (u.x[0,1] - u.x[0,-1])/(2.*Delta);
    double dudz = (u.x[0,0,1] - u.x[0,0,-1])/(2.*Delta);
    double dvdx = (u.y[1]   - u.y[-1]  )/(2.*Delta);
    double dvdy = (u.y[0,1] - u.y[0,-1])/(2.0*Delta);
    double dvdz = (u.y[0,0,1] - u.y[0,0,-1])/(2.*Delta);
    double dwdx = (u.z[1]     - u.z[-1]    )/(2.*Delta);
    double dwdy = (u.z[0,1]   - u.z[0,-1]  )/(2.*Delta);
    double dwdz = (u.z[0,0,1] - u.z[0,0,-1])/(2.0*Delta);
    double SDeformxx = dudx;
    double SDeformxy = 0.5*(dudy + dvdx);
    double SDeformxz = 0.5*(dudz + dwdx);
    double SDeformyx = SDeformxy;
    double SDeformyy = dvdy;
    double SDeformyz = 0.5*(dvdz + dwdy);
    double SDeformzx = SDeformxz;
    double SDeformzy = SDeformyz;
    double SDeformzz = dwdz; 
    double sqterm = 2.0*dv()*(sq(SDeformxx) + sq(SDeformxy) + sq(SDeformxz) +
			      sq(SDeformyx) + sq(SDeformyy) + sq(SDeformyz) +
			      sq(SDeformzx) + sq(SDeformzy) + sq(SDeformzz));
    
    rateWater += mu1*f[]*sqterm; //water
    rateAir   += mu2*(1. - f[])*sqterm; //air
  }

  rates[0] = rateWater;
  rates[1] = rateAir;
  return 0;
}



event logfile (t = 1e-6; t <= MAXTIME; t += 2e-5) {
  double ke_d = 0., ke_a = 0., area = 0.;

  foreach(reduction(+:ke_d) reduction(+:ke_a)) {
    // Kinetic energy of the droplet and ambient airflow
    ke_d += 0.5*dv()*(sq(u.x[]) + sq(u.y[]) + sq(u.z[])) * rho1*f[];
    ke_a += 0.5*dv()*(sq(u.x[]) + sq(u.y[]) + sq(u.z[])) * rho2*(1-f[]);
  }

  // Droplet surface area
  area = interface_area(f);



  // Viscous dissipation rate
  double rates[2];
  dissipation_rate(u, rates);
  double diss_w = rates[0], diss_a = rates[1];
  
  
  
  double delta = L0 / (double)(1 << LEVEL);
  double h_contact_ave = 0.0; int n_contact = 0;
  double h_contact_max = 0.0;
  
  foreach(reduction(+:h_contact_ave) reduction(+:n_contact) reduction(max:h_contact_max))
    if (f[] > 1e-6 && f[] < 1.-1e-6 && -0.5*delta <= x && x <= 0.5*delta) {
      h_contact_ave += y; h_contact_max = fmax(h_contact_max, y);
      n_contact++;
    }
  h_contact_ave /= n_contact;

  fprintf (fk, "%.8f %.8f %.8f %.8f\n", t, h_contact_ave, h_contact_max, area);
  fflush (fk);

  fprintf (fe, "%.8f %.8f %.8f %.8f %.8f\n", t, ke_d, ke_a, diss_w, diss_a);
  fflush (fe);
}



event count_droplets (t = 0; t <= MAXTIME; t += 2e-5) {
  scalar m_cd[];
  foreach() m_cd[] = f[] > 0.5;
  int n = tag (m_cd);

  double v[n], ux_drop[n], uy_drop[n], uz_drop[n], area[n], ke[n];
  coord b[n];

  for (int j = 0; j < n; j++)
    v[j] = b[j].x = b[j].y = b[j].z = ux_drop[j] = uy_drop[j] = uz_drop[j] = area[j] = ke[j] = 0.;

  foreach_leaf()
    if (m_cd[] > 0) {
      int j = m_cd[] - 1;
      v[j] += dv()*f[];

      coord p = {x,y,z};
      foreach_dimension()
	    b[j].x += dv()*f[]*p.x;

      ux_drop[j] += dv()*f[]*u.x[];
	  uy_drop[j] += dv()*f[]*u.y[];
	  uz_drop[j] += dv()*f[]*u.z[];
	  
      ke[j] += 0.5 * rho1*dv()*f[] * (sq(u.x[]) + sq(u.y[]) + sq(u.z[]));
	  
      // interfacial area
      if (f[] > 1e-4 && f[] < 1. - 1e-4) {
        coord m = mycs (point, f);
        double alpha = plane_alpha (f[], m);
        coord p;
        area[j] += sq(Delta) * plane_area_center(m, alpha, &p);
      }
    }

#if _MPI
  MPI_Allreduce (MPI_IN_PLACE, v, n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce (MPI_IN_PLACE, area, n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce (MPI_IN_PLACE, ke, n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce (MPI_IN_PLACE, b, 3*n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce (MPI_IN_PLACE, ux_drop, n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce (MPI_IN_PLACE, uy_drop, n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce (MPI_IN_PLACE, uz_drop, n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif

  if (n > 1 && pid() == 0) {
    char resultname[100];
    sprintf(resultname, "./frag_stats/frag_stat_%4.4f.dat", t);
    FILE * fdrop = fopen(resultname, "w");
  
    for (int j=0; j<n; j++) {
      fprintf (fdrop, "%g %d %g %g %g %g %g %g %g %g %g\n", t, j, v[j], b[j].x/v[j], b[j].y/v[j], b[j].z/v[j], ux_drop[j]/v[j], uy_drop[j]/v[j], uz_drop[j]/v[j], ke[j], f.sigma*area[j]);
    }
    
    fclose(fdrop);
  }
}



#define POPEN(name, mode) fopen (name ".ppm", mode)

event movie (t = 1e-6; t += 2e-5) {
#if TREE
  scalar l[];
  foreach()
    l[] = level;
#endif
#if dimension == 2
  static FILE * fp = POPEN ("f", "w");
  output_ppm (f, fp, n=1024);
  static FILE * fp4 = POPEN ("l", "w");
  output_ppm (l, fp4, min=1, max=LEVEL);
#else //dimension == 3
//  view (width = 800, height = 600, tx=-0.5, ty=-0.5, theta=-pi/4., phi=-pi/6., fov=10);
  view (camera = "iso");
  view (tx = 0.4, ty=0.3, fov=20);
  clear();
  box();
  draw_vof("f");
  //squares("l", n = {1,0,0}, alpha=0.);
  static FILE *fp = POPEN ("movie", "w");
  save (fp = fp);
  //view (theta=pi/4., phi=pi/6., tx = 0.5, ty=0.5);
  clear();
  view (tx = 0.0001, ty=0.0001);
  view (camera = "right");
  view (tx = 0.5, ty = -0.5);
  draw_vof("f");
  static FILE *fpr = POPEN ("right", "w");
  save (fp = fpr);
#endif //dimension
}

event snapshot (t += 0.002) { // dump files
  char dname[100];
  sprintf(dname, "dump%g", t);
  dump(dname);
}

event snapshot_restart (t += 1e-4) { // dump files
  char dname[100];
  sprintf(dname, "restart");
  dump(dname);
}


event end (t=MAXTIME) { 
  printf ("i = %d t=%g\n", i, t);
  dump ("end");
}
