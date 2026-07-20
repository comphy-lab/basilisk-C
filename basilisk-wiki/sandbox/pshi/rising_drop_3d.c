#include "navier-stokes/centered.h"

#define mu(f) (1./(clamp(f,0,1.)*1./mu1 + clamp(1-f,0,1)*1./mu2)) 
#define rho(f) (clamp(f,0.,1.)*rho1  + clamp(1-f,0,1)*rho2) 

#include "two-phase.h"
#include "tension.h"
#include "tag.h"
#include "output_htg.h"

#define RHOR 1.156
#define MUR 1.618

#define eq_radius (1.5) // mm 
#define MAXTIME 1000.

#define WIDTH 960.0
#define Z0 14.5
#define D0 2.0
#define S1 WIDTH/2.0

//scalar l2[], omegay[];

//#define half_domain

int LEVELMAX = 16;

double Ga, Bo;

int main () {
  
  double N = 128;
  size (WIDTH);
  #ifdef half_domain
  origin (0, 0, 0);
  #else
	  origin (0, 0, -L0/2);
  #endif
  init_grid (N);
  rho2 = 1.0;
  rho1 = rho2 / RHOR;

  Ga = 14.9895 * pow(eq_radius, 1.5);
  Bo = 0.037745 * pow(eq_radius, 2.0);

  mu2 = (1.0 - rho1) / Ga;
  mu1 = mu2 / MUR;

  alpha = alphav;
  rho = rhov;
  
  f.sigma = (rho2-rho1)/Bo;
  TOLERANCE = 1e-5;
  run();
}

event init (t = 0) {
	if (!restore (file = "restart")) {
    refine (sq(x-S1) + sq(y - Z0) + sq(z) - sq(4.0*D0) < 0 && level < (14));
    fraction(f, -(sq(x-S1) + sq(y - Z0) + sq(z) - sq(D0/2.)));
	}
}

event acceleration (i++) {
  face vector av = a;
  foreach_face(y){
    av.y[] -= 1.;
  }
	
}

//event stability (i++) {CFL = 0.05; }

event adapt (i++) {
	
	double femax = 1e-3;
    double uemax = 5.0e-3;
	scalar g1f[],g2f[];
	foreach(){
		g1f[] = sqrt(sq(f[1,0,0] - f[-1,0,0])+sq(f[0,1,0] - f[0,-1,0])+sq(f[0,0,1] - f[0,0,-1]));
	}
	boundary ({g1f});
	foreach(){
        g2f[] = sqrt(sq(g1f[1,0,0] - g1f[-1,0,0])+sq(g1f[0,1,0] - g1f[0,-1,0])+sq(g1f[0,0,1] - g1f[0,0,-1]));
	}
	boundary ({g2f});
    adapt_wavelet ({g2f,u},(double []){femax,uemax,uemax,uemax},LEVELMAX,6);
}
/**
event remove_drops ( i+=20)
{
  scalar m[];
  foreach()
    m[] = f[] > 0.;
  int n = tag (m);
  double v[n];
  int remove_flag[n];
  for (int j = 0; j < n; j++) {
    v[j] = 0.;
    remove_flag[j] = 0;
  }
  foreach_leaf()
    if (m[] > 0) {
      int j = m[] - 1;
      v[j] += dv()*f[];
    }

#if _MPI
  MPI_Allreduce (MPI_IN_PLACE, v, n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif

  for (int j = 0; j < n; j++)
    if ( v[j] < 0.01 )
       remove_flag[j] = 1;

  foreach()
    if (m[] > 0) {
      int j = m[] - 1;
      if ( remove_flag[j] == 1 )
         f[]=0.;
    }
	  boundary ({f});
}
**/

event information (t = 0; t <= MAXTIME; t=t+0.05) {
//event information (i++) {

  double xb = 0., yb = 0., zb = 0., sb = 0.;
  double vbx = 0., vby = 0., vbz = 0.;
  double rmax = -HUGE, rmin = HUGE;
  double area = 0.; 

  double E_i = 0., E_e = 0.; 
  foreach(reduction(+:xb) reduction(+:yb) reduction(+:zb) 
	  reduction(+:vbx) reduction(+:vby) reduction(+:vbz) 
	  reduction(+:sb)) {
      double dv = (f[])*dv();
      xb += x*dv;
      yb += y*dv;
      vbx += u.x[]*dv;
      vby += u.y[]*dv;
      sb += dv;
	  #ifdef half_domain
	  zb += 0.*dv;
	  vbz += 0.*dv;
	  
		#else
      zb += z*dv;
      vbz += u.z[]*dv;
		#endif
	
  }

  foreach(reduction(+:E_i) reduction(+:E_e)) {
      double dv_i = (f[])*dv();
	  double dv_e = (1.0 - f[])*dv()*(Delta > 4. ? 0. : 1.0);
	  double r_e = sqrt((x-xb/sb) * (x-xb/sb) + (z-zb/sb) * (z-zb/sb));
	  double e_theta_x = -(z-zb/sb) / r_e;
      double e_theta_z = (x-xb/sb) / r_e;
	  double u_theta = (u.x[] * e_theta_x + u.z[] * e_theta_z);
	  E_i += rho(f[])*dv_i*u_theta*u_theta;
	  E_e += rho(f[])*dv_e*u_theta*u_theta;
  }
E_i = E_i/sb;
E_e = E_e/sb;

  foreach(reduction(max:rmax) reduction(min:rmin) reduction(+:area)) {
    if (f[] > 1.0e-4 && f[] < 1) {
		double fmax = -HUGE, fmin = HUGE;
			for (int i = -1; i <= 1; i++)
				for (int j = -1; j <= 1; j++) 
					for (int k = -1; k <= 1; k++) {
					if (i == 0 && j == 0 && k == 0) continue;
					fmax = max(fmax, f[i,j,k]);
					fmin = min(fmin, f[i,j,k]);
				}

    if (fmax > 0.9 && fmin < 0.1) {
      coord p;
      coord n = mycs (point, f); // n is a vector normal to the interface 
      double alpha = plane_alpha (f[], n);// and alpha is the intercept (http://basilisk.fr/src/geometry.h).
      double s = plane_area_center (n, alpha, &p);//s denotes the surface fractions i.e. the fractions of the faces of the cell which are inside the interface.(http://basilisk.fr/src/fractions.h#mycs)
      double rad  = sqrt(sq(x + Delta*p.x-xb/sb) + sq(y + Delta*p.y-yb/sb) + sq(z + Delta*p.z-zb/sb)); 
      if (rad > rmax)
		 rmax = rad;
      if (rad < rmin)
	     rmin = rad;
      area += pow(Delta,2)*s;
	}
	}
  }
  

  // 1-time 2-volume 3-x 4-y 5-z 6-vx 7-vy 8-vz 9-dmax/2.0 10-dmin/2.0 11-xmin 

  fprintf (stderr,
	   //"%.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %03d %.8f\n", 
	   "%.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.3d %.8f\n", 
	   // 1    2    3    4    5    6    7    8    9    10   11   12   13   14   15   16   17   18    19   20    21   22   23    24   25    26   27); 
	   //t, sb, xb/sb, yb/sb, zb/sb, vbx/sb, vby/sb, vbz/sb, rmax, rmin, xmin, Ga, Bo, LEVELMAX, area, incAngle
	   t, sb, xb/sb, yb/sb, zb/sb, //1-5
	   vbx/sb, vby/sb, vbz/sb, rmax, rmin, //6-10
	   E_i, E_e, Ga, Bo, LEVELMAX, area// 11-15
	   );//25-30

  fflush (stderr);

}

event snapshot1 (t = 0; t <= MAXTIME; t=t+1.0)
{
  char name[80];
  //sprintf (name, "dump-%03d", (int) t);
  sprintf (name, "dump-%g", t);
  dump (file = name);
}

event snapshot2 (t = 0; t <= MAXTIME; t=t+1.)

{
  char path[]="htg"; // no slash at the end!!
  char prefix[80];

  scalar omegax[], omegay[], omegaz[];
    foreach(){
    omegax[] = (u.z[0,1,0] - u.z[0,-1,0] - u.y[0,0,1] + u.y[0,0,-1])/(2.*Delta);
    omegay[] = (u.x[0,0,1] - u.x[0,0,-1] - u.z[1,0,0] + u.z[-1,0,0])/(2.*Delta);
	omegaz[] = (u.y[1,0,0] - u.y[-1,0,0] - u.x[0,1,0] + u.x[0,-1,0])/(2.*Delta);
	}
	boundary ({omegax, omegay,omegaz});

  sprintf(prefix, "data-%g", t);
  output_htg((scalar *) {f,p,omegax,omegay,omegaz,u.x,u.y,u.z},(vector *){uf}, path, prefix, i, t);
}