/**
 This code simulates granular discharge from a silo under oscillatory  conditions.
 The system includes gravity and centrifugal acceleration, with optional Euler and Coriolis forces to account for rotational effects.
 */

#include "grid/multigrid.h"
#include "navier-stokes/centered.h"
#include "vof.h"

// ==== Acceleration term toggles ====
#define USE_CENTRIPETAL  1   // Centrifugal acceleration
#define USE_EULER        0   // Euler acceleration
#define USE_CORIOLIS     0   // Coriolis acceleration
#define USE_GRAVITY      1   // Gravity

#define LY 4.      // height of Domain
#define LX 1.      // width of Domain

#define RHOF 1e-4
#define mug  1e-5


double  H0,R0,D,W,tmax,Q,Wmin,muwall,WDOMAIN,LEVEL;
double  thetamax,T;

char s[80];
FILE * fpf, *fwq;
scalar f[];
scalar * interfaces = {f};
face vector alphav[];
face vector muv[];
face vector av[];
scalar rhov[];
scalar mu_I[];
scalar eta[];

p[top] = dirichlet(0);
u.n[top] = neumann(0);
u.t[bottom] =  fabs(x-LX/2.)<= W/2 ? neumann(0):  dirichlet(0);
u.n[bottom] =  fabs(x-LX/2.)<= W/2 ? neumann(0):  dirichlet(0);
p[bottom]   =  fabs(x-LX/2.)<= W/2 ? dirichlet(0): neumann(0);
u.n[right] = dirichlet(0);
u.n[left] = dirichlet(0);
u.t[right] = dirichlet(0);
u.t[left] = dirichlet(0);

int main(int argc, char ** argv) {
    
    LEVEL    = 6;
    T        = 6;                // rocking period
    thetamax = 45 * (M_PI/180.0);// rocking amplitude
    W        = 0.1875;           // orifice width
    tmax     = 20;
    dt       = 0.008;
    
    dimensions (nx = 1, ny = 4);
    L0 = LX;

    WDOMAIN = 0.25;
    
    N = ldexp(L0,LEVEL);
    muwall = 0.1;
    TOLERANCE = 1e-3;
    H0 = LY - 0.1;
    R0 = LX;
    D = 1./90.;        // Grain size
    fwq = fopen("outWQW", "w");
    fclose(fwq);

    DT = dt;
    a = av;
    alpha = alphav;
    mu = muv;
    rho = rhov;
    
    Q = 0;
    fpf = fopen ("interface.txt", "w");
    run();
    fclose (fpf);
    fprintf (stdout,"\n");
    fwq = fopen ("outWQ", "a");
    fprintf(fwq," %lf %lf \n", W, Q);
    fclose (fwq);
}
/**
   initial heap, a rectangle
*/
event init (t = 0) {
  scalar phi[];
  foreach_vertex()
    phi[] = min(H0 - y, R0 - x);
  fractions (phi, f);

  foreach()
    p[] = (fabs(x-LX/2.)<= W/2 && fabs(y)<= .1) ?  0 : max(H0 - y,0) ;
  boundary({p});
}

// Stay still for the first 2 seconds
#define QUIET_T 2.0
static inline void kinematics(double t, double *theta, double *theta_dot, double *theta_ddot)
{
  if (T <= 0.) {
    *theta     = thetamax;
    *theta_dot = 0.;
    *theta_ddot= 0.;
    return;
  }
  if (t < QUIET_T) {
    *theta     = 0.;
    *theta_dot = 0.;
    *theta_ddot= 0.;
    return;
  }
  const double omega = 2.*M_PI/T;
  const double tau   = t - QUIET_T;
  const double s     = sin(omega*tau);
  const double c     = cos(omega*tau);
  *theta      = thetamax * s;
  *theta_dot  = thetamax * omega * c;
  *theta_ddot = -thetamax * omega * omega * s;
}

static inline void accel_components_cell(
  double x, double y, double ux, double uy,
  double x0, double y0, double theta, double theta_dot, double theta_ddot, double gphys,
  double *a_cen_x, double *a_cen_y, double *a_eul_x, double *a_eul_y,
  double *a_cor_x, double *a_cor_y, double *g_x, double *g_y)
{
  const double dx = x - x0, dy = y - y0;
  const double r  = sqrt(dx*dx + dy*dy);

  if (USE_GRAVITY) {
    *g_x = -sin(theta);
    *g_y = -cos(theta);
  } else {
    *g_x = 0., *g_y = 0.;
  }

  if (USE_CENTRIPETAL && r > 0.) {
    const double cen_mag = (theta_dot*theta_dot) * r / gphys;
    *a_cen_x = cen_mag * (dx / r);
    *a_cen_y = cen_mag * (dy / r);
  } else *a_cen_x = 0., *a_cen_y = 0.;

  if (USE_EULER) {
    *a_eul_x = ( theta_ddot * dy)/gphys;
    *a_eul_y = (-theta_ddot * dx)/gphys;
  } else *a_eul_x = 0., *a_eul_y = 0.;

  if (USE_CORIOLIS) {
    *a_cor_x = ( 2.0 * theta_dot * uy)/gphys;
    *a_cor_y = (-2.0 * theta_dot * ux)/gphys;
  } else *a_cor_x = 0., *a_cor_y = 0.;
}

event acceleration (i++) {
  const double gphys = 9.81;
  const double x0 = LX/2.;
  const double y0 = 0.75*LY;

  double theta, theta_dot, theta_ddot;
  kinematics(t, &theta, &theta_dot, &theta_ddot);

  scalar axc[], ayc[];
  foreach () {
    double a_cen_x, a_cen_y, a_eul_x, a_eul_y, a_cor_x, a_cor_y, g_x, g_y;
    accel_components_cell(x, y, u.x[], u.y[], x0, y0,
                          theta, theta_dot, theta_ddot, gphys,
                          &a_cen_x, &a_cen_y, &a_eul_x, &a_eul_y,
                          &a_cor_x, &a_cor_y, &g_x, &g_y);

    axc[] = a_cen_x + a_eul_x + a_cor_x + g_x;
    ayc[] = a_cen_y + a_eul_y + a_cor_y + g_y;
  }
  boundary({axc, ayc});

  foreach_face(x) av.x[] = 0.5*(axc[] + axc[-1,0]);
  foreach_face(y) av.y[] = 0.5*(ayc[] + ayc[0,-1]);


  boundary((scalar *){av});
  return 0;
}

event inertial_dump (t = QUIET_T; t += 0.25; t <= tmax) {
  static int made = 0;
  if (!made) { system("mkdir -p inertial"); made = 1; }

  const double gphys = 9.81;
  const double x0 = LX/2.;
  const double y0 = 0.75*LY;

  double theta, theta_dot, theta_ddot;
  kinematics(t, &theta, &theta_dot, &theta_ddot);

  char fname[128];
  snprintf(fname, sizeof(fname), "inertial/inertial-%.2f.txt", t);
  FILE *fp = fopen(fname, "w");
  if (!fp) return 0;

  fprintf(fp,
    "# t=%.6g\n"
    "# x y "
    "a_cen_x a_cen_y a_cen "
    "a_eul_x a_eul_y a_eul "
    "a_cor_x a_cor_y a_cor "
    "g_x g_y g_abs "
    "sum_x sum_y\n", t);

  foreach() {
    double a_cen_x, a_cen_y, a_eul_x, a_eul_y, a_cor_x, a_cor_y, g_x, g_y;
    accel_components_cell(x, y, u.x[], u.y[], x0, y0,
                          theta, theta_dot, theta_ddot, gphys,
                          &a_cen_x, &a_cen_y, &a_eul_x, &a_eul_y,
                          &a_cor_x, &a_cor_y, &g_x, &g_y);

    double a_cen = sqrt(a_cen_x*a_cen_x + a_cen_y*a_cen_y);
    double a_eul = sqrt(a_eul_x*a_eul_x + a_eul_y*a_eul_y);
    double a_cor = sqrt(a_cor_x*a_cor_x + a_cor_y*a_cor_y);
    double g_abs = sqrt(g_x*g_x + g_y*g_y);
    double sum_x = a_cen_x + a_eul_x + a_cor_x + g_x;
    double sum_y = a_cen_y + a_eul_y + a_cor_y + g_y;

    fprintf(fp,
      "%.9g %.9g "
      "%.9g %.9g %.9g "
      "%.9g %.9g %.9g "
      "%.9g %.9g %.9g "
      "%.9g %.9g %.9g "
      "%.9g %.9g\n",
      x, y,
      a_cen_x, a_cen_y, a_cen,
      a_eul_x, a_eul_y, a_eul,
      a_cor_x, a_cor_y, a_cor,
      g_x, g_y, g_abs,
      sum_x, sum_y
    );
  }
  fclose(fp);
  return 0;                  
}

#define rho(f) ((f) + RHOF*(1. - (f)))
event properties (i++) {
    trash({alphav});
    foreach() {
        eta[] = mug;
        if (p[] > 0.) {
            double D2 = 0.;
            foreach_dimension() {
                double dxx = u.x[1,0] - u.x[-1,0];
                double dxy = (u.x[0,1] - u.x[0,-1] + u.y[1,0] - u.y[-1,0])/2.;
                D2 += sq(dxx) + sq(dxy);
            }
            if (D2 > 0.) {
                D2 = sqrt(D2)/(2.*Delta);  // this is D2
                double sD2 = sqrt(2.)*D2; // this sD2 is (sqrt(2) D2)
                double In = sD2*D/sqrt(p[]);
                double muI = 0.4 + 0.28*In/(0.4 + In);
                mu_I[] = muI;
                double etamin = sqrt(D*D*D);
                eta[] = max((muI*p[])/sD2, etamin);
                eta[] = min(eta[], 100);
            }
        }
    }
    boundary({eta,mu_I});
    scalar fa[];
    foreach()
        fa[] = (4.*f[] +
                2.*(f[-1,0] + f[1,0] + f[0,-1] + f[0,1]) +
                f[1,1] + f[-1,1] + f[1,-1] + f[-1,-1])/16.;
    boundary({fa});
    foreach_face() {
        double fm = (fa[] + fa[-1,0])/2.;
        muv.x[] = (fm*(eta[] + eta[-1,0])/2. + (1. - fm)*mug);
        // mu.x[] = 1./(2.*fm/(eta[] + eta[-1,0]) + (1. - fm)/mug);
        alphav.x[] = 1./rho(fm);
    }
    foreach()
     rhov[] = rho(fa[]);
     boundary ({muv,alphav,rhov});
}

void mg_print (mgstats mg)
{
    if (mg.i > 0 && mg.resa > 0.)
        fprintf (stderr, "#   %d %g %g %g\n", mg.i, mg.resb, mg.resa,
                 exp (log (mg.resb/mg.resa)/mg.i));
}

event logfile (i++) {
    stats s = statsf(f);
    fprintf(stderr, "%g %d %g %g %g %g\n", t, i, dt, s.sum, s.min, s.max - 1.);
    mg_print(mgp);
    mg_print(mgpf);
    mg_print(mgu);
    fflush(stderr);
}

/**
   wall friction
   $$\frac{du}{dt} = \frac{-2 \mu_w p u}{W ||u||}$$
   equation with discretization:
   $$\frac{u^{n+1}-u^n}{\delta t} = \frac{-2 \mu_w p u^{n+1}}{W ||u^n||} $$
*/


/**
event friction (i++) {
    foreach() {
            double m = 2. * muwall * dt * p[] / WDOMAIN;
            double U = norm(u);
            foreach_dimension()
                u.x[] = U > 0 ? max(U - m, 0) * u.x[] / U : 0;
    }
}
 */

event interface (t = 2. ; t+=1. ; t<= tmax) {
#if dimension == 2
  output_facets (f, fpf);
#endif
  char s[80];
  sprintf (s, "field-%g.txt", t);
  FILE * fp = fopen (s, "w");
  output_field ({f,p,u,uf,pf,mu_I,eta}, fp, linear = true);
  fclose (fp);
}

event debit (t += 0.05; t <= tmax ) {
    static double Vold,V=1,Qinst=0;
    Vold=V;
    V=0;
    foreach()
        V += f[] * Delta * Delta;
    Qinst = -(V-Vold)/0.05;
    if(t>=.1) fprintf (stdout,"%lf %lf %lf \n",t,V/L0/H0,Qinst);
    fflush (stdout);
}

#if 1

event movie (t += 0.05) {
    
    scalar l[];
    
    foreach()
        l[] = f[] * (1 + sqrt(sq(u.x[]) + sq(u.y[])));
    boundary({l});
    static FILE * fp1 = popen("ppm2mpeg > velo.mpg", "w");
    output_ppm(l, fp1, min = 0, max = 2., linear = true,n = 256, box = {{0,0},{LX,LY}});


    static FILE * fp2 = popen("ppm2mpeg > p.mpg", "w");
    foreach()
        l[] = p[];
    output_ppm(l, fp2, min = 0, linear = true, n = 256, box = {{0,0},{LX,LY}});
}

event pictures (t == 4) {
    output_ppm(f, file = "f.png", min = 0, max = 2, spread = 2, n = 2048, linear = true,
               box = {{0,0},{2,2}});
}
#endif

/**
 
 #to run
  qcc  -g -O2 -o granular_sandglass_muw granular_sandglass_muw.c -lm
  ./granular_sandglass_muw > out

   Plot of pressure at time 5
 
   ~~~gnuplot pressure
 
   reset
   set terminal pngcairo size 700,1400
   set output 'pressure_t5.png'
   set pm3d; set palette rgbformulae 22,13,-31; unset surface
   set ticslevel 0.
   unset border
   unset xtics
   unset ytics
   unset ztics
   unset colorbox

   set xrange [0:1]
   set yrange [0:4]
   set size ratio -1
   set view 0,0
   splot 'field-5.txt' using 1:2:4

   ~~~


   Plot of velocity at time 5
 
   ~~~gnuplot velocity
 
   reset
   set terminal pngcairo size 700,1400
   set output 'velocity_t5.png'
   set pm3d; set palette rgbformulae 22,13,-31; unset surface
   set ticslevel 0.
   unset border
   unset xtics
   unset ytics
   unset ztics
   unset colorbox
   set xrange [0:1]
   set yrange [0:4]
   set size ratio -1
   set view 0,0
   splot 'field-5.txt' using 1:2:(sqrt($5*$5 + $6*$6))

   ~~~
 
   Authors: M.Li, Y. Zhou & P.-Y. Lagrée
*/
