//#include "grid/multigrid3D.h"
#include "grid/octree.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "navier-stokes/conserving.h"
#include "tension.h"
#include "reduced.h"  // reduced gravity
#include "view.h"
#include "tag.h"
#include "lambda2.h"
#include "navier-stokes/perfs.h"

//#include "adapt_wavelet_limited.h"
//#include "sandbox/frac-dist.h" // Extra headerfiles used in profiling function
//#include "sandbox/profile6.h"  // From Antoon

#define POPEN(name, mode) fopen (name ".ppm", mode)

double RELEASETIME = 200; 
int MAXLEVEL = dimension == 2 ? 10 : 5; // Max refinement level if not use limited refinement

double RE_tau = 720.; // Air side friction Reynolds number  
double BO = 200; // Default Bond number
double RE = 10000.; // Wave Reynolds number that is dependent on c
double k_ = 1;
double g_ = 1;
double h_ = 1;
double ak = 0.1;
double c_;
double UstarRATIO = 0.5;
double Ustar;

double uemax = 0.01;
double femax = 0.0001;
double uwemax = 0.001;
double uemaxRATIO = 0.01;

#define RATIO 1.225/1000. // Density ratio, air to water
#define MURATIO 18.31e-6/10.0e-4 // Dynamic viscosity ratio, air to water
// kinematic viscosity air = 16*water
double alter_MU = 1.; // An extra factor if deviating from MURATIO. Default is 1, which means MURATIO is used.
                      // We use alter_MU > 1 for an artificially smaller water viscosity to prevent fast damping. 

vector u_water[]; // Velocity field only used with uwemax

/**
   We need to store the variable forcing field av[]. */
double amp_force = 0.1; // Amplitude of the forcing
face vector av[];

/**
   Some paramters are passed in from the command line. */
int main(int argc, char *argv[]) {
  if (argc > 1)
    RE_tau  = atof (argv[1]);
  if (argc > 2)
    BO = atof(argv[2]);
  if (argc > 3)
    MAXLEVEL = atoi(argv[3]);
  if (argc > 4)
    g_ = atof(argv[4]);
  if (argc > 5)
    ak = atof(argv[5]);
  if (argc > 6)
    RELEASETIME = atof(argv[6]);
  if (argc > 7)
    uemaxRATIO = atof(argv[7]);
  if (argc > 8)
    alter_MU = atof(argv[8]);

  L0 = 2*pi;
  h_ = 1; // Water depth
  k_ = 4; // Four waves per box
  origin (-L0/2., 0, -L0/2.);
  // According to http://basilisk.fr/Basilisk%20C#boundary-conditions
  // for top, u.n = u.y, u.t = u.z, u.r = u.x
  u.r[top] = neumann(0);
  u.r[bottom] = neumann(0);
  u.n[top] = dirichlet(0); // This is supposed to be neumann 
  u.n[bottom] = dirichlet(0);
  u.t[top] = neumann(0);
  u.t[bottom] = neumann(0);
  // TO-DO: Test if setting to neumann change 
  periodic (right);
  periodic (front);
  rho1 = 1.;
  rho2 = RATIO;
  f.sigma = g_*rho1/(BO*sq(k_));
  G.y = -g_;
  c_ = sqrt(g_/k_+f.sigma/rho1*k_);
  // Ustar = c_*UstarRATIO; // Obsolete
  Ustar = 0.25; // Pick a fixed value
  mu2 = Ustar*rho2*(L0-h_)/RE_tau;
  mu1 = mu2/(MURATIO)/alter_MU;
  RE = rho1*c_*(2*pi/k_)/mu1; // RE now becomes a dependent Non-dim number on c 
  fprintf (stderr, "g = %g, c = %g, Ustar = %g, MURATIO = %g, mu_w = %g, rho_w = %g, mu_a = %g, rho_a = %g, sigma = %g, Bo = %g, RE = %g, Re_tau = %g\n", g_, c_, Ustar, MURATIO*alter_MU, mu1, rho1, mu2, rho2, f.sigma, BO, RE, RE_tau);
  // Give the address of av to a so that acceleration can be changed
  a = av;
  init_grid (1 << 7);
  // Refine according to 
  uemax = uemaxRATIO*Ustar;
  uwemax = 0.001*c_;
  fprintf (stderr, "RELEASETIME = %g, uemax = %g \n", RELEASETIME, uemax);
  run();
}

/** 
   Specify the interface shape. Since we are considering ak=0.1/0.2/0.3, use Stokes instead of linear profile. */

double WaveProfile_linear (double x, double z) {
  double a_ = ak/k_;
  double eta1 = a_*cos(k_*x);
  return eta1 + h_;
}

double WaveProfile (double x, double y)
{
  double a_ = ak/k_;
  double eta1 = a_*cos(k_*x);
  double alpa = 1./tanh(k_*h_);
  double eta2 = 1./4.*alpa*(3.*sq(alpa) - 1.)*sq(a_)*k_*cos(2.*k_*x);
  double eta3 = -3./8.*(cube(alpa)*alpa - 
			3.*sq(alpa) + 3.)*cube(a_)*sq(k_)*cos(k_*x) + 
    3./64.*(8.*cube(alpa)*cube(alpa) + 
	    (sq(alpa) - 1.)*(sq(alpa) - 1.))*cube(a_)*sq(k_)*cos(3.*k_*x);
  return eta1 + eta2 + eta3 + h_;
}

/**
   Random noise gets killed by adaptive mesh refinement. Therefore, instead of initializing with a mean log profile plus random noise, we initialize with only the mean, and we wait for the instability to naturally develop (because of the shear profile). */

event init (i = 0) {
  if (!restore("restart")){
    double rand = 0;
    double ytau = (mu2/rho2)/Ustar;
    do {
      fraction (f, WaveProfile(x,z)-y);
      foreach() {
	// Initialize with a fairly accurate profile
	if ((y-WaveProfile(x,z)) > 0.05)
	  u.x[] = (1-f[])*(log((y-WaveProfile(x,z))/ytau)*Ustar/0.41);
	else
	  u.x[] = 0.;
	u.y[] = 0.;
	u.z[] = 0.;
      }
      boundary ((scalar *){u});
    }
#if TREE
    // No need for adaptation when starting 
    while (0);
#else
    while (0);
#endif
  }
}

/**
  Set the wave velocity 0, until the RELEASETIME. */
event set_wave(i=0; i++; t<RELEASETIME) {
  fraction (f, WaveProfile(x,z)-y);
  foreach(){
    u.x[] = (1.0 - f[])*u.x[];
    u.y[] = (1.0 - f[])*u.y[];
    u.z[] = (1.0 - f[])*u.z[];
  }
  boundary ((scalar *){u});
}

/**
   Release the wave at RELEASETIME. We don't do any adaptation at this step. 
   And we use linear wave instead of stokes. (UPDATE: use stokes for 0.1/0.2/0.3) */

double u_x_linear (double x, double y) {
  return ak*c_*cos(x*k_)*exp(y*k_);
}
double u_y_linear (double x, double y) {
  return ak*c_*sin(x*k_)*exp(y*k_);
}

// TO-DO: maybe include g correction by Bond number. High Bond should be fine.
double u_x (double x, double y)
{
  double alpa = 1./tanh(k_*h_);
  double a_ = ak/k_;
  double sgma = sqrt(g_*k_*tanh(k_*h_)*
		     (1. + k_*k_*a_*a_*(9./8.*(sq(alpa) - 1.)*
					(sq(alpa) - 1.) + sq(alpa))));
  double A_ = a_*g_/sgma;
  return A_*cosh(k_*(y + h_))/cosh(k_*h_)*k_*cos(k_*x) +
    3.*ak*A_/(8.*alpa)*(sq(alpa) - 1.)*(sq(alpa) - 1.)*
    cosh(2.0*k_*(y + h_))*2.*k_*cos(2.0*k_*x)/cosh(2.0*k_*h_) +
    1./64.*(sq(alpa) - 1.)*(sq(alpa) + 3.)*
    (9.*sq(alpa) - 13.)*
    cosh(3.*k_*(y + h_))/cosh(3.*k_*h_)*ak*ak*A_*3.*k_*cos(3.*k_*x);
}

double u_y (double x, double y)
{
  double alpa = 1./tanh(k_*h_);
  double a_ = ak/k_;
  double sgma = sqrt(g_*k_*tanh(k_*h_)*
		     (1. + k_*k_*a_*a_*(9./8.*(sq(alpa) - 1.)*
					(sq(alpa) - 1.) + sq(alpa))));
  double A_ = a_*g_/sgma;
  return A_*k_*sinh(k_*(y + h_))/cosh(k_*h_)*sin(k_*x) +
    3.*ak*A_/(8.*alpa)*(sq(alpa) - 1.)*(sq(alpa) - 1.)*
    2.*k_*sinh(2.0*k_*(y + h_))*sin(2.0*k_*x)/cosh(2.0*k_*h_) +
    1./64.*(sq(alpa) - 1.)*(sq(alpa) + 3.)*
    (9.*sq(alpa) - 13.)*
    3.*k_*sinh(3.*k_*(y + h_))/cosh(3.*k_*h_)*ak*ak*A_*sin(3.*k_*x);
}

event start(t = RELEASETIME) {
  // A slightly changed version of the test/stokes.h wave function as y = 0 at the bottom now so y+h -> y
  fraction (f, WaveProfile(x,z)-y);
  foreach () {
    u.x[] += u_x(x, y-h_)*f[];
    u.y[] += u_y(x, y-h_)*f[];
  }
  boundary ((scalar *){u});
}

/**
   Forcing term equivalent to the pressure gradient in x direction. */
event acceleration (i++) {
  double ampl = sq(Ustar)/(L0-h_);
  foreach_face(x)
    av.x[] += ampl*(1.-f[]);
}

/** 
   Output video and field. We first do simple movies of the volume fraction, level of refinement fields. In 3D, these are in a $z=0$ cross-section. */
#  define POPEN(name, mode) fopen (name ".ppm", mode)
event movies (t += 0.1) {
  scalar omega[];
  vorticity (u, omega);
  clear();
  view (fov = 40, camera = "iso", ty = -0.25,
  width = 600, height = 600, bg = {1,1,1}, samples = 4);
  squares ("u.x", linear = true, n = {0,0,1}, alpha = -3.1415);
  squares ("omega", linear = true, n = {1,0,0}, alpha = -3.1415);
  cells (n = {1,0,0}, alpha = -3.1415);
  draw_vof ("f", color = "u.x");
  char s[80];
  sprintf (s, "t = %0.1f", t);
  draw_string (s, size = 30);
  {
    static FILE * fp = POPEN ("3D", "a");
    save (fp = fp);
  }
  scalar l2[];
  lambda2 (u, l2);
  view (fov = 40, camera = "iso", ty = -0.25,
  width = 600, height = 600, bg = {1,1,1}, samples = 4);
  squares ("u.y", linear = true, n = {0,0,1}, alpha = -3.1415);
  squares ("u.x", linear = true, n = {1,0,0}, alpha = -3.1415);
  cells (n = {1,0,0}, alpha = -3.1415);
  draw_vof ("f", color = "u.x");
  isosurface ("l2", -1);
  draw_string (s, size = 30);
  {
    static FILE * fp = POPEN ("vortex", "a");
    save (fp = fp);
  }
}

/**
   Generate averaged profile in y direction on the fly (not used in post-processing, just for quick diagnosis). Uses Antoon's */
//event profile_output (t += 0.5) {
//  char file[99];
//  sprintf (file, "prof_%g", t);
//  scalar uxuy[],uxux[],uyuy[],uzuz[];
//  foreach () {
//    uxuy[] = u.x[]*u.y[];
//    uxux[] = u.x[]*u.x[];
//    uyuy[] = u.y[]*u.y[];
//    uzuz[] = u.z[]*u.z[];
//  }
//  vertex scalar phi[];
//  foreach_vertex ()
//    phi[] = y;
//  // default phi is y
//  profiles ({u.x, u.y, u.z, uxuy, uxux, uyuy, uzuz}, phi, rf = 0.5, fname = file, min = 0.8, max = 2.*pi);
//}

/** 
   Outputting slices on the fly. Fields are interpolated onto a 2d (x-y) structure grid of maxlevel, at z=zp. */
void sliceXY(char * fname, scalar s, double zp, int maxlevel){
  FILE *fpver =fopen (fname,"w"); 
  int nn = (1<<maxlevel);
  double ** field = matrix_new (nn, nn, sizeof(double));
  double stp = L0/(double)nn;
  for (int i = 0; i < nn; i++)
    {
      double xp = stp*i + X0 + stp/2.;
      for (int j = 0; j < nn; j++) 
	{
	  double yp = stp*j + Y0 + stp/2.;
	  Point point = locate (xp, yp,zp);
	  field[i][j] = point.level >= 0 ? s[] : nodata;
	}
    }
  if (pid() == 0){ // master
#if _MPI
    MPI_Reduce (MPI_IN_PLACE, field[0], sq(nn), MPI_DOUBLE, MPI_MIN, 0,
		MPI_COMM_WORLD);
#endif
    for (int i = 0; i < nn; i++) {
      for (int j = 0; j < nn; j++) {
	fprintf (fpver, "%g\t", field[i][j]);
      }
      fputc ('\n', fpver);
    }
    fflush (fpver);
  }
#if _MPI
  else // slave
    MPI_Reduce (field[0], NULL, nn*nn, MPI_DOUBLE, MPI_MIN, 0,
		MPI_COMM_WORLD);
#endif
  matrix_free (field);
}


/**
   Output eta on the fly (an obsolete function). */
void output_twophase (double snapshot_time) {
  scalar pos[];
  coord G = {0.,1.,0.}, Z = {0.,0.,0.};
  position (f, pos, G, Z);
  char etaname[100];
  sprintf (etaname, "./eta/eta_t%g_%d", snapshot_time, pid());
  FILE * feta = fopen (etaname, "w");
  fprintf(feta, "x,z,pos,p,p_p1,p_p2\n");
  // printing out quantities: p_p1 for p at plus 1, p_m1 for p at minus 1 etc.
  foreach(){
    if (interfacial (point, f)){
      fprintf (feta, "%g,%g,%g,%g,%g,%g\n", 
	       x, z, pos[], p[], p[0,1], p[0,2]);
      /* tau.x[], tau.y[], u.x[], u.y[], n.x/norm_2, n.y/norm_2); */
    }
  }
  fclose (feta);
} 

/**
   Output eta on the fly. This function first find the interfacial grid points by interfacial (point, f), then try to locate the (leaf?) cell above it, offset by stp. It is not guaranteed to find it so a `if (point.level > 0)` filtering is needed. The outputs are x, z, eta, slope, pressure, dudy, dvdx, dudx, dvdy. */
void output_twophase_locate (double snapshot_time) {
  scalar pos[];
  coord G = {0.,1.,0.}, Z = {0.,0.,0.};
  position (f, pos, G, Z);
  char etaname[100];
  sprintf (etaname, "./eta/eta_loc_t%g_%d", snapshot_time, pid());
  FILE * feta = fopen (etaname, "w");
  fprintf(feta, "x,z,pos,epsilon,p,dudy,dvdx,dudx,dvdy\n");
  // printing out quantities: p_p1 for p at plus 1, p_m1 for p at minus 1 etc.
  double stp = L0/256.;
  foreach(){
    if (interfacial (point, f)){
      if (point.level == MAXLEVEL) {
	coord n = mycs (point, f);
	double norm_2 = sqrt(sq(n.x) + sq(n.y));
	//fprintf (stderr, "Interface cell x = %g, y = %g, Delta = %g, ux = %g \n", x, y, Delta, u.x[]);
	//fflush (stderr);
	double eta = pos[];
	double yp = y + stp;
	point = locate (x, yp, z);
	// Getting the local normal vector
	// n is norm 1 and has to be normalized to norm 2
	// double dudx = (u.x[1] - u.x[-1])/(2.*Delta);
	// fprintf (stderr, "Inline 1! \n");
	if (point.level > 0) {
	  POINT_VARIABLES;
	  /* fprintf (stderr, "Inline 2! \n"); */
	  /* fprintf (stderr, "Above cell x = %g, y = %g \n", x, y); */
	  /* fprintf (stderr, "Above cell Delta = %g \n", Delta); */
	  /* fprintf (stderr, "Above cell ux = %g \n", u.x[]); */
	  double dudy1 = (u.x[0,1] - u.x[0,-1])/(2.*Delta);
	  //double dudy2 = (u.x[0,2] - u.x[0,0])/(2.*Delta);
	  //double dudz = (u.x[0,0,1] - u.x[0,0,-1])/(2.*Delta);
	  double dvdx1 = (u.y[0,1] - u.y[0,-1])/(2.*Delta);
	  //double dvdx2 = (u.y[1,1] - u.y[1,-1])/(2.*Delta);
	  double dudx1 = (u.x[0,1] - u.x[0,-1])/(2.*Delta);
	  //double dudx2 = (u.x[1,1] - u.x[1,-1])/(2.*Delta);
	  double dvdy1 = (u.y[0,1] - u.y[0,-1])/(2.*Delta);
	  //double dvdy2 = (u.y[0,2] - u.y[0,0])/(2.*Delta);
	  //tau.x[] = 2*mu2*(SDeform.x.x[]*n.x + SDeform.y.x[]*n.y)/norm_2; 
	  //tau.y[] = 2*mu2*(SDeform.x.y[]*n.x + SDeform.y.y[]*n.y)/norm_2;  
	  //double tau1 = 2*mu2*(dudy1 + dvdx); 
	  //double tau2 = 2*mu2*(dudy2 + dvdx); 
	  fprintf (feta, "%g,%g,%g,%g,%g,%g,%g,%g,%g\n",
		   x, z, eta, -n.x/n.y, p[0,0], dudy1, dvdx1, dudx1, dvdy1);
	/* tau.x[], tau.y[], u.x[], u.y[], n.x/norm_2, n.y/norm_2); */
	}
      }
    }
  }
  fflush (feta);
  fclose (feta);
}

event eta_output (t += 0.1) {
  if (t > RELEASETIME) // Event does not take variable time condition 
    output_twophase_locate (t);
}

/**
   We can wither take slices in the spanwise direction or the vertical direction. We choose to take 256 slices of size 512x512 in the spanwise direction.
   
/* event output_slice (t += 0.05) */
/* { */
/*   char filename[100]; */
/*   double zslice = 0.; */
/*   sprintf (filename, "./field/ux_t%g_center", t); */
/*   sliceXY (filename,u.x,zslice,MAXLEVEL-1); */
/*   sprintf (filename, "./field/uy_t%g_center", t); */
/*   sliceXY (filename,u.y,zslice,MAXLEVEL-1); */
/*   sprintf (filename, "./field/uz_t%g_center", t); */
/*   sliceXY (filename,u.z,zslice,MAXLEVEL-1); */
/*   sprintf (filename, "./field/f_t%g_center", t); */
/*   sliceXY (filename,f,zslice,MAXLEVEL-1); */
/* } */

scalar pair[];
event turbulence_stat (t += 0.5) {
  char filename[100];
  int Nslice = 256;
  double L0 = 2*pi;
  double zslice = -L0/2+L0/2./Nslice;
  for (int i=0; i<Nslice; i++) {
    zslice += L0/Nslice;
    sprintf (filename, "./field/ux_t%g_slice%d", t, i);
    sliceXY (filename,u.x,zslice,9);
    sprintf (filename, "./field/uy_t%g_slice%d", t, i);
    sliceXY (filename,u.y,zslice,9);
    sprintf (filename, "./field/uz_t%g_slice%d", t, i);
    sliceXY (filename,u.z,zslice,9);
    sprintf (filename, "./field/f_t%g_slice%d", t, i);
    sliceXY (filename,f,zslice,9);
  }
}

event p_stat (t += 0.2) {
  char filename[100];
  int Nslice = 256;
  double L0 = 2*pi;
  double zslice = -L0/2+L0/2./Nslice;
  foreach () {
    pair[] = p[]*(1-f[]);
  }
  for (int i=0; i<Nslice; i++) {
    zslice += L0/Nslice;
    sprintf (filename, "./field/pair_run_t%g_slice%d", t, i);
    sliceXY (filename,pair,zslice,9);
  }  
}
/**
   Dump at the end, and dump regularly for restarting the run. */

event end (t = 1000.) {
  fprintf (fout, "i = %d t = %g\n", i, t);
  dump ("end");
}

event dumpforrestart (t += 1) { // There seems to be problem with tiger with dumping too frequently. Otherwise could do t+=0.1
  char dname[100];
  u_water.x.nodump = true;
  u_water.y.nodump = true;
  u_water.z.nodump = true;
  // p.nodump = false;
  sprintf (dname, "restart");
  dump (dname);
}

event dumpstep (t += 0.5) {
  char dname[100];
  u_water.x.nodump = true;
  u_water.y.nodump = true;
  u_water.z.nodump = true;
  pair.nodump = true;
  // p.nodump = false;
  sprintf (dname, "dump%g", t);
  dump (dname);
}

/** 
   Adaptive function. uemax is tuned to be 0.3Ustar in most cases. We need a more strict criteria for water speed once the waves starts moving. */ 
#if TREE
event adapt (i++) {
  if (i == 5)
    fprintf(stderr, "uemaxRATIO = %g\n", uemaxRATIO);
  if (t < RELEASETIME)
    adapt_wavelet ({f,u}, (double[]){femax,uemax,uemax,uemax}, MAXLEVEL);
  if (t >= RELEASETIME) {
    foreach () {
      foreach_dimension ()
	u_water.x[] = u.x[]*f[];
    }
    adapt_wavelet ({f,u,u_water}, (double[]){femax,uemax,uemax,uemax,uwemax,uwemax,uwemax}, MAXLEVEL);
  }
}
#endif
