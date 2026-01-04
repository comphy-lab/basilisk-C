#include "grid/octree.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "navier-stokes/conserving.h"
#include "tension.h"
#include "view.h"
#include "tag.h"
#include "lambda2.h"
#include "navier-stokes/perfs.h"

#include "turb-inlet-imposed.h"
//#include "adapt_wavelet_limited.h"
//#include "frac-dist.h"

#include "maxruntime.h"

#define SQUARES
#include "signature_rescaled_new_3D.h"

#include <sys/stat.h>
#include <sys/types.h>



#define POPEN(name, mode) fopen (name ".ppm", mode)

face vector av[];
face vector muv[];

int MINLEV = 5;
int MAXLEV_TURB = 8; // Max refinement level if not use limited refinement
int MAXLEV_DROP = 11;
int SIGN_LEV = 11;


double WE = 15.0; // Default Weber number
double OH = 0.005; // Default Reynolds number

#define INSERTTIME 50.2 
#define MAXTIME 600.

#define RHOR 833.
#define MUR 20.
#define WIDTH 20.
#define Xc 5
#define R0 1.
#define U0 1.
#define L_turb 1.
#define u_rms 0.25



FILE * fe = NULL; FILE * fk = NULL; FILE * fturb = NULL;
FILE * fslice_x = NULL; FILE * fslice_y = NULL; FILE * fslice_z = NULL; 

/**
   We need to store the variable forcing field av[]. */

scalar phii[], sign[], M[];
scalar l2[], omegay[], omegax[], pp[];

const int max_change = 2000; 
double T_MD = 0.01;
double epsilon = 0.02;
bool large = true;



u.n[right] = neumann(0.);
u.t[right] = neumann(0.);
u.r[right] = neumann(0.);

p[right] = dirichlet(0.);
p[left]  = neumann(0.);



/**
   Some paramters are passed in from the command line. */
int main(int argc, char *argv[]) {

  size(WIDTH);
  origin (-Xc, -L0/2, -L0/2.);
  init_grid (1 << (MINLEV));
  

  // According to https://basilisk.fr/Basilisk%20C#boundary-conditions
  // for top, u.n = u.y, u.t = u.z, u.r = u.x
  // periodic (top); periodic (back);
  
  rho1 = 1.;
  rho2 = 1./RHOR;
  
  f.sigma = rho2*sq(U0)*(2.*R0)/WE;
  
  mu1 = OH*sqrt(rho1*f.sigma*2.*R0);
  mu2 = mu1/MUR;
  
  a = av; mu = muv;
  
  double Re = MUR*sqrt(WE) / sqrt(RHOR)/OH;
  
  TOLERANCE = 1e-6;
  
  fprintf (stderr, "INSERTTIME = %g\n", INSERTTIME);
  fprintf (stderr, "Re = %g, Delta = %g\n", Re, WIDTH/(1 << MAXLEV_TURB));
  
  char name[80];

  sprintf (name, "drop-energy-%d.dat", MAXLEV_DROP);
  fe = fopen (name, "w");

  sprintf (name, "drop-kinetics-%d.dat", MAXLEV_DROP);
  fk = fopen (name, "w");
  
  sprintf (name, "turb-kinetics-%d.dat", MAXLEV_TURB);
  fturb = fopen (name, "w");
  
  mkdir("./chirco_debugs", 0777);
  mkdir("./injection_debugs", 0777);
  mkdir("./interface_profiles", 0777);
  mkdir("./frag_stats", 0777);
  mkdir("./field", 0777);
  mkdir("./profs", 0777);
  
  sprintf (name, "./field/fluc-ux-xpos-%d.dat", MAXLEV_TURB); fslice_x = fopen (name, "w");
  sprintf (name, "./field/fluc-uy-xpos-%d.dat", MAXLEV_TURB); fslice_y = fopen (name, "w");
  sprintf (name, "./field/fluc-uz-xpos-%d.dat", MAXLEV_TURB); fslice_z = fopen (name, "w");
  
  run();
}



event set_turb(i++) {
  set_inlet(MAXLEV_TURB, L_turb, U0, u_rms);
}



/**
   Random noise gets killed by adaptive mesh refinement. Therefore, instead of initializing with a mean log profile plus random noise, we initialize with only the mean, and we wait for the instability to naturally develop (because of the shear profile). */

event init (i = 0) {
  if (!restore("restart")) {
    do {      
      foreach() {
        f[] = 0.0;
        u.x[] = U0; u.y[] = 0.; u.z[] = 0.;
      }
      
      set_inlet(MAXLEV_TURB, L_turb, U0, u_rms);
#if TREE
    // No need for adaptation when starting 
    } while (0);
#else
    while (0);
#endif
  }
}



// Insert the drop
event insert_drop (t = INSERTTIME) {
  refine (sq(x) + sq(y) + sq(z) - sq(1.05*R0) < 0 && sq(x) + sq(y) + sq(z) - sq(0.95*R0) > 0 && level < MAXLEV_DROP);
  fraction (f, -pow(sq(x) + sq(y) + sq(z), 0.25) + sqrt(R0));
  
  foreach() {
    foreach_dimension() {u.x[] = u.x[] * (1.-f[]);}
  }
  boundary ((scalar *){u});
}



event count_droplets (t = INSERTTIME + 50; t <= MAXTIME; t += 0.02) {
  scalar m_cd[];
  foreach() m_cd[] = f[] > 0.5;
  int n = tag (m_cd);

  /**
  Once each cell is tagged with a unique droplet index, we can easily
  compute the volume *v* and position *b* of each droplet. Note that
  we use *foreach (serial)* to avoid doing a parallel traversal when
  using OpenMP. This is because we don't have reduction operations for
  the *v* and *b* arrays (yet). */

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

 /**
 When using MPI we need to perform a global reduction to get the
 volumes and positions of droplets which span multiple processes. */

#if _MPI
  MPI_Allreduce (MPI_IN_PLACE, v, n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce (MPI_IN_PLACE, area, n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce (MPI_IN_PLACE, ke, n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce (MPI_IN_PLACE, b, 3*n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce (MPI_IN_PLACE, ux_drop, n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce (MPI_IN_PLACE, uy_drop, n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce (MPI_IN_PLACE, uz_drop, n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif

  /**
  Finally we output the volume and position of each droplet to
  standard output. */

  if (n > 1 && pid() == 0) {
    char resultname[100];
    sprintf(resultname, "./frag_stats/frag_stat_%4.2f.dat", t);
    FILE * fdrop = fopen(resultname, "w");
  
    for (int j=0; j<n; j++) {
      fprintf (fdrop, "%g %d %g %g %g %g %g %g %g %g %g\n", t, j, v[j], b[j].x/v[j], b[j].y/v[j], b[j].z/v[j], ux_drop[j]/v[j], uy_drop[j]/v[j], uz_drop[j]/v[j], ke[j], f.sigma*area[j]);
    }
    
    fclose(fdrop);
  }
}



event neck_detect(t = INSERTTIME + 80; t <= MAXTIME; t += T_MD){
  foreach(){
    phii[] = 2*f[] - 1;
    sign[] = 7;
  }
  
  int l_sign = SIGN_LEV;
  
  for (int ilev = depth() - 1; ilev >= l_sign; ilev--)  
    foreach_level(ilev){
      if(is_refined(cell))
      restriction_average(point, phii);
    }
  
  compute_signature_neigh_level (f, phii, sign, l_sign);
  
  if (pid()==0)  
    printf("time %g level used for moments %d and depth is %d \n", t, l_sign, depth()); 
  
  for (int ilev = l_sign; ilev < depth(); ilev++)  
    foreach_level(ilev){
      sign.prolongation = phii.prolongation = refine_injection;
      if(is_refined(cell)){
        sign.prolongation (point, sign);
        phii.prolongation (point, phii);
      }
    }

  change_topology (f, sign, M, l_sign, max_change, large, epsilon, WIDTH/(1 << l_sign), T_MD);
}



/** 
   Output video and field. We first do simple movies of the volume fraction, level of refinement fields. In 3D, these are in a $z=0$ cross-section. */
#define POPEN(name, mode) fopen (name ".ppm", mode)

event movie (t = INSERTTIME + 0.001; t += 0.05)
{
  view (quat = {0.000, 0.707, 0.000, 0.707}, fov = 60, near = 0.01, far = 1000, tx = 0.001, ty = -0.002, tz = -0.642, width = 720, height = 720);
  clear();
  draw_vof ("f");
  box();

  static FILE * fp = POPEN ("movie", "w");
  save (fp = fp);
}



event dumpforrestart (t = 0.; t <= MAXTIME; t += 0.1) {
  char dname[100];
  //pair.nodump = true;
  
  sprintf (dname, "restart");
  dump (dname);
}


event dumpstep (t = 0.; t <= MAXTIME; t += 1) {
  char dname[100];
  
  lambda2 (u, l2);
  foreach() pp[] = p[];
  
  sprintf (dname, "dump-%g", t);
  dump (dname);
}



void dissipation_rate (vector u, double* rates) {
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

  rates[0] = rateWater; rates[1] = rateAir;
}



event logfile_drop (t = INSERTTIME; t <= MAXTIME; t += 0.05) {
  double xb = 0., yb = 0., zb = 0., sb = 0., vbx = 0., vby = 0., vbz = 0., ke_d = 0., ke_a = 0., a = 0., diss_w = 0, diss_a = 0.;

  foreach(reduction(+:xb) reduction(+:yb) reduction(+:zb) 
      reduction(+:vbx) reduction(+:vby) reduction(+:vbz) 
      reduction(+:sb) reduction(+:ke_d) reduction(+:ke_a) reduction(+:a)) {
    double dv_l = f[]*dv();
    
    xb += x*dv_l; yb += y*dv_l; zb += z*dv_l;
    vbx += u.x[]*dv_l; vby += u.y[]*dv_l; vbz += u.z[]*dv_l;
    sb += dv_l;

    // Kinetic energy of the droplet and ambient airflow
    ke_d += 0.5*dv()*(sq(u.x[]) + sq(u.y[]) + sq(u.z[])) * rho1*f[];
    ke_a += 0.5*dv()*(sq(u.x[]) + sq(u.y[]) + sq(u.z[])) * rho2*(1-f[]);
  }

  // Droplet surface area
  a = interface_area(f);



  // Viscous dissipation rate
  double rates[2];
  dissipation_rate(u, rates);
  diss_w = rates[0], diss_a = rates[1];
  
  // Droplet CM's location and velocity
  xb = xb/sb; yb = yb/sb; zb = zb/sb; vbx = vbx/sb; vby = vby/sb; vbz = vbz/sb;

 double h1_max = -100.0, h1_min = 100.0, h2_max = -100.0, h2_min = 100.0, t_max = -100.0, t_min = 100.0;
  foreach(reduction(max:h1_max) reduction(min:h1_min) reduction(max:h2_max) reduction(min:h2_min) reduction(max:t_max) reduction(min:t_min))
    if (f[] > 1e-6 && f[] < 1.-1e-6) {
      h1_max = fmax(h1_max, y); h1_min = fmin(h1_min, y);
      h2_max = fmax(h2_max, z); h2_min = fmin(h2_min, z); 
      t_max = fmax(t_max, x); t_min = fmin(t_min, x);
    }

  fprintf (stderr,
	   "%.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f\n", 
	   t, sb, a, xb, yb, zb,
	   vbx, vby, vbz, 0.5*(h1_max-h1_min+h2_max-h2_min), t_max-t_min);
  fflush (stderr);

  fprintf (fk,
	   "%.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f\n", 
	   t, sb, a, xb, yb, zb, vbx, vby, vbz, 0.5*(h1_max-h1_min+h2_max-h2_min), t_max-t_min);
  fflush (fk);

  double ke_d_t = 0.5*rho1*sb*(sq(vbx)+sq(vby)+sq(vbz));
  fprintf (fe, 
          "%.8f %.8f %.8f %.8f %.8f %.8f %.8f\n", 
          t, ke_d_t, ke_d-ke_d_t, ke_a, diss_w, diss_a, (a-(4*M_PI*sq(R0)))*f.sigma);
  fflush (fe);
}



event logfile_turb (t = 0; t <= MAXTIME; t += 0.1) {
  coord ubar; double vol_air = 0.;
  
  double ubarx = 0., ubary = 0., ubarz = 0.; 
  foreach(reduction(+:vol_air)) vol_air += dv() * (1.-f[]);
  
  foreach(reduction(+:ubarx) reduction(+:ubary) reduction(+:ubarz)) {
    ubarx += u.x[] * dv() * (1.-f[]) / vol_air;
    ubary += u.y[] * dv() * (1.-f[]) / vol_air;
    ubarz += u.z[] * dv() * (1.-f[]) / vol_air;
  }
  
  ubar.x = ubarx; ubar.y = ubary; ubar.z = ubarz;
  
  double ke = 0., vd = 0.;
  foreach(reduction(+:ke) reduction(+:vd)) {
    //foreach_dimension() {
      // mean fluctuating kinetic energy
      ke += dv()*sq(u.x[] - ubar.x) * (1. - f[]);
      ke += dv()*sq(u.y[] - ubar.y) * (1. - f[]);
      ke += dv()*sq(u.z[] - ubar.z) * (1. - f[]);
      
      // viscous dissipation
      vd += dv()*(sq(u.x[1] - u.x[-1]) + sq(u.x[0,1] - u.x[0,-1]) + sq(u.x[0,0,1] - u.x[0,0,-1]))/sq(2.*Delta) * (1. - f[]);
      vd += dv()*(sq(u.y[1] - u.y[-1]) + sq(u.y[0,1] - u.y[0,-1]) + sq(u.y[0,0,1] - u.y[0,0,-1]))/sq(2.*Delta) * (1. - f[]);
      vd += dv()*(sq(u.z[1] - u.z[-1]) + sq(u.z[0,1] - u.z[0,-1]) + sq(u.z[0,0,1] - u.z[0,0,-1]))/sq(2.*Delta) * (1. - f[]);
    //}
  }
  
  ke *= 0.5 * rho2 / vol_air;
  vd *= mu2 / rho2 / vol_air;

  if (t == 0) {
    fprintf (stderr, "t vol_turb dissipation_norm energy Reynolds ubar\n");
    fprintf (fturb, "t vol_turb dissipation energy Reynolds ubar\n");
    fflush(fturb);
  }

  fprintf (stderr, "%g %g %g %g %g %g\n", t, vol_air, vd, ke, 2./3. * ke / (mu2/rho2) * sqrt(15.*mu2/rho2/vd), ubar.x);
  fprintf (fturb, "%g %g %g %g %g %g\n", t, vol_air, vd, ke, 2./3. * ke / (mu2/rho2) * sqrt(15.*mu2/rho2/vd), ubar.x);
  fflush(fturb);
}



event outputInterface(t = INSERTTIME; t <= MAXTIME; t += 1) {
  char resultname[100];
  sprintf(resultname, "./interface_profiles/results%4.2f_%d.dat", t, pid());
  FILE * fp = fopen(resultname, "w");

  scalar xpos[], ypos[], zpos[];
  position (f, xpos, {1, 0, 0});
  position (f, ypos, {0, 1, 0});
  position (f, zpos, {0, 0, 1});
  foreach()
      if (xpos[] != nodata)
        fprintf (fp, "%g %g %g\n", xpos[], ypos[], zpos[]);
  
  //output_facets (f, fp);
  fclose (fp);
  
  /*
  char command[80];
  sprintf(command, "LC_ALL=C  cat results* > interface_%g.dat", t);
  system(command);
  sprintf(command, "LC_ALL=C  rm results*");
  system(command);
  */
}



event u_rms_slice (t = 0.; t <= INSERTTIME; t += 0.5) {
  int Nslice = 50;
  double xslice = -Xc + L0/2./Nslice;
  int nn = 1 << MAXLEV_TURB;
  double r_ave[3]; double r_sq_ave[3]; double r_var[3];
  
  double ** ux_field = matrix_new (nn+1, nn+1, sizeof(double));
  double ** uy_field = matrix_new (nn+1, nn+1, sizeof(double));
  double ** uz_field = matrix_new (nn+1, nn+1, sizeof(double));
  
  for (int i=0; i<Nslice-1; i++) {
    xslice += L0/Nslice;
    
    sliceYZ (ux_field, u.x, xslice, MAXLEV_TURB);
    sliceYZ (uy_field, u.y, xslice, MAXLEV_TURB);
    sliceYZ (uz_field, u.z, xslice, MAXLEV_TURB);
    
    for (int dir_iter = 0; dir_iter < 3; dir_iter++) {
      r_ave[dir_iter] = 0.; r_sq_ave[dir_iter] = 0.;
    }

    for (int j=0; j<nn; j++)
      for (int k=0; k<nn; k++) {
        r_ave[0] += ux_field[j][k]; r_sq_ave[0] += sq(ux_field[j][k]);
        r_ave[1] += uy_field[j][k]; r_sq_ave[1] += sq(uy_field[j][k]);
        r_ave[2] += uz_field[j][k]; r_sq_ave[2] += sq(uz_field[j][k]);
      }
    
    for (int dir_iter = 0; dir_iter < 3; dir_iter++) {
      r_ave[dir_iter] /= sq(nn); r_sq_ave[dir_iter] /= sq(nn); r_var[dir_iter] = r_sq_ave[dir_iter] - sq(r_ave[dir_iter]);  
    }
 
    if (pid() == 0) {
      fprintf (fslice_x, "%g ", sqrt(r_var[0])); fprintf (fslice_y, "%g ", sqrt(r_var[1])); fprintf (fslice_z, "%g ", sqrt(r_var[2]));
      fflush(fslice_x); fflush(fslice_y); fflush(fslice_z);
    }
  }     
  
  matrix_free (ux_field); matrix_free (uy_field); matrix_free (uz_field);
  fprintf (fslice_x, "\n"); fprintf (fslice_y, "\n"); fprintf (fslice_z, "\n"); 
}



/** 
   Adaptive function. uemax is tuned to be 0.3Ustar in most cases. We need a more strict criteria for water speed once the waves starts moving. */ 

#if TREE
event adapt (i++) {
  vector u_a[], u_w[]; 
  scalar levfield[];
  
  vector gf[];
  gradients ({f}, {gf});

  foreach() {
    foreach_dimension() u_a.x[] = u.x[] * (1.-f[]);
    foreach_dimension() u_w.x[] = u.x[] * f[];
  }
  
  foreach() {
    if (sq(gf.x[]) + sq(gf.y[]) + sq(gf.z[]) > 0.01) levfield[] = MAXLEV_DROP;
     else if (f[] > 0.99) levfield[] = MINLEV;
       else levfield[] = MAXLEV_TURB;
  }

  double uemax = 0.01;
  double femax = 0.01;
  
  if (t <= INSERTTIME) adapt_wavelet ({u_a.x,u_a.y,u_a.z}, (double[]){uemax,uemax,uemax}, MAXLEV_TURB, MINLEV);
  //if (t > INSERTTIME) adapt_wavelet_limited ({gf.x, gf.y, gf.z, u_a.x, u_a.y, u_a.z}, (double[]){femax,femax,femax,uemax,uemax,uemax}, levfield, MINLEV);
  if (t > INSERTTIME) adapt_wavelet ({gf.x, gf.y, gf.z}, (double[]){femax,femax,femax}, MAXLEV_DROP, MAXLEV_TURB);
}

#endif