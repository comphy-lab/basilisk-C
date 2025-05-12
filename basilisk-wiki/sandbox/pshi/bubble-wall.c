#include "navier-stokes/centered.h"

#define mu(f) (1./(clamp(f,0,1.)*1./mu1 + clamp(1-f,0,1)*1./mu2)) 
#define rho(f) (clamp(f,0.,1.)*rho1  + clamp(1-f,0,1)*rho2) 

#include "two-phase.h"
#include "tension.h"
#include "tag.h"
#include "output_htg.h"

#define RHOR 1000.
#define MUR 100.

# define Ga 30.
# define Bo 0.25
   
# define MAXTIME 200.

#define WIDTH 480.0
#define Z0 14.5
#define D0 2.0
#define S1  2.  //Initial distance

//scalar l2[], omegay[];

#define half_domain

int LEVELMAX = 15;
int LEVELMAX_a = 15;
int LEVELMAX_b = 16;
int LEVELMAX_c = 17;

f[left] = 0.;
u.t[left]  = dirichlet(0);
u.n[left]  = dirichlet(0);
u.r[left]  = dirichlet(0);

int main () {
  
  double N = 128;
  size (WIDTH);
  #ifdef half_domain
  origin (0, 0, 0);
  #else
	  origin (0, 0, -L0/2);
  #endif
  init_grid (N);

  rho2 =1.;
  mu2 = 1./Ga;
  rho1 = 1./RHOR;
  mu1 = 1./(MUR*Ga);
  alpha = alphav;
  rho = rhov;
  
  f.sigma =  1./Bo;
  TOLERANCE = 1e-4;
  run();
}

event init (t = 0) {
	if (!restore (file = "restart")) {
    refine (sq(x-S1) + sq(y - Z0) + sq(z) - sq(2.0*D0) < 0 && level < (12));
    fraction(f, -(sq(x-S1) + sq(y - Z0) + sq(z) - sq(D0/2.)));
	}
}

event acceleration (i++) {
  face vector av = a;
  foreach_face(y)
    av.y[] -= 1.;
}

event adapt (i++) {
	
	double femax = 1e-3;
    double uemax = 5.0e-3;
	double xmin_ad = HUGE;
	foreach(reduction(min:xmin_ad)) {
       if (f[] > 1.0e-4 && f[] < 1 && x < xmin_ad) {
		   xmin_ad = x;
        }
	}
	if (xmin_ad > 0.15)
		LEVELMAX = LEVELMAX_a;
	else{
		if (xmin_ad > 0.075)
			LEVELMAX = LEVELMAX_b;
		else
			LEVELMAX = LEVELMAX_c;
	}
	adapt_wavelet ({f,u}, (double[]){femax,uemax,uemax,uemax}, LEVELMAX, 7);

}

event information (t = 0; t <= MAXTIME; t=t+0.05) {
  double xb = 0., yb = 0., zb = 0., sb = 0.;
  double vbx = 0., vby = 0., vbz = 0.;
  double rmax = -HUGE, rmin = HUGE;
  double xmin = HUGE;
  double area = 0.; 

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

  foreach(reduction(max:rmax) reduction(min:rmin) reduction(min:xmin) reduction(+:area)) {
    if (f[] > 1.0e-4 && f[] < 1) {
      coord p;
      coord n = mycs (point, f); // n is a vector normal to the interface 
      double alpha = plane_alpha (f[], n);// and alpha is the intercept (http://basilisk.fr/src/geometry.h).
      double s = plane_area_center (n, alpha, &p);//s denotes the surface fractions i.e. the fractions of the faces of the cell which are inside the interface.(http://basilisk.fr/src/fractions.h#mycs)
      double rad  = sqrt(sq(x + Delta*p.x-xb/sb) + sq(y + Delta*p.y-yb/sb) + sq(z + Delta*p.z)); 
      if (rad > rmax)
		 rmax = rad;
      if (rad < rmin)
	     rmin = rad;
	 double radx  = x + Delta*p.x; 
      if (radx < xmin)
	     xmin = radx;
      area += pow(Delta,2)*s;
    }
 }  
  

  // 1-time 2-volume 3-x 4-y 5-z 6-vx 7-vy 8-vz 9-dmax/2.0 10-dmin/2.0 11-xmin 
  fprintf (stderr,
	   //"%.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %03d %.8f\n", 
	   "%.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.3d %.8f\n", 
	   // 1    2    3    4    5    6    7    8    9    10   11   12   13   14   15   16   17   18    19   20    21   22   23    24   25    26   27); 
	   //t, sb, xb/sb, yb/sb, zb/sb, vbx/sb, vby/sb, vbz/sb, rmax, rmin, xmin, Ga, Bo, LEVELMAX, area, incAngle
	   t, sb, xb/sb, yb/sb, zb/sb, //1-5
	   vbx/sb, vby/sb, vbz/sb, rmax, rmin, //6-10
	   xmin, Ga, Bo, LEVELMAX, area// 11-15
	   );//25-30

  fflush (stderr);

}

event snapshot1 (t = 0; t <= MAXTIME; t=t+0.05)
{
  char name[80];
  //sprintf (name, "dump-%03d", (int) t);
  sprintf (name, "dump-%g", t);
  dump (file = name);
}

event snapshot2 (t = 0; t <= MAXTIME; t=t+0.25)
{
  char path[]="htg"; // no slash at the end!!
  char prefix[80];
  scalar omegax[], omegay[], omegaz[];
  
    omegax[] = (u.z[0,1,0] - u.z[0,-1,0] - u.y[0,0,1] + u.y[0,0,-1])/(2.*Delta);
    omegay[] = (u.x[0,0,1] - u.x[0,0,-1] - u.z[1,0,0] + u.z[-1,0,0])/(2.*Delta);
    omegaz[] = (u.y[1,0,0] - u.y[-1,0,0] - u.x[0,1,0] + u.x[0,-1,0])/(2.*Delta);
	}
	boundary ({omegax, omegay,omegaz});

  sprintf(prefix, "data-%g", t);
  output_htg((scalar *) {f,p,omegax,omegay,omegaz,u.x,u.y,u.z},(vector *){uf}, path, prefix, i, t);
}

