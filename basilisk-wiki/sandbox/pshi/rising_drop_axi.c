#include "axi.h"
#include "navier-stokes/centered.h"

#define mu(f) (1./(clamp(f,0,1.)*1./mu1 + clamp(1-f,0,1)*1./mu2)) // mu = mu1 if f = 1
#define rho(f) (clamp(f,0.,1.)*rho1  + clamp(1-f,0,1)*rho2) 

#include "two-phase.h"
#include "tension.h"
#include "tag.h"
#include "output_htg.h"

// liquid-liquid system in dimentional form
#define RHOR_drop 862.3 // kg m^-3
#define RHOR_fluid 997.02 // kg m^-3
#define RHOR (RHOR_drop/RHOR_fluid)

#define MU_drop (0.55e-3) // Pa s 
#define MU_fluid (0.89e-3) // Pa s 
#define MUR (MU_drop/MU_fluid)

#define SurTen (35.0e-3) // N/m 
#define a_grav 9.8067 // m s^-2 
#define eq_radius (1.2e-3) // m 
#define v_grav (sqrt(fabs(1.0-RHOR) * eq_radius * a_grav)) // m s^-1 

# define Ga (RHOR_fluid * eq_radius * v_grav / MU_fluid)
# define Bo (fabs(RHOR_fluid - RHOR_drop) * a_grav * eq_radius * eq_radius / SurTen)

//#define RHOR 1.156
//#define MUR 1.618

//# define Ga 19.7
//# define Bo 0.0543
# define MAXTIME 500.

#define WIDTH 480.
#define D0 2.0
#define S1 15.

//scalar l2[], omegay[];
scalar omega_i[];

//#define half_domain

int LEVELMAX = 15;
//int LEVELMAX_a = 14;
//int LEVELMAX_b = 15;
//int LEVELMAX_c = 16;

f[right] = 0.;
//u.t[left]  = dirichlet(0);
//u.n[left]  = dirichlet(0);
//u.r[left]  = dirichlet(0);
u.t[right]  = dirichlet(0);
u.n[right]  = dirichlet(0);
uf.n[bottom] = 0.;
uf.n[top] = 0.;

int main () {
  
  double N = 128;
  size (WIDTH);
  //#ifdef half_domain
  origin (0, 0);
  //#else
//	  origin (0, 0, -L0/2);
//  #endif
  init_grid (N);
  rho2 = 1.; // liquid density
  rho1 = rho2*RHOR;
  //mu2 = (rho2-rho1)/Ga; // liquid dynamic viscosity
  mu2 = sqrt(fabs(1.0-RHOR))*rho2/Ga; // liquid dynamic viscosity
  mu1 = mu2*MUR;
  alpha = alphav;
  rho = rhov;
  f.sigma = fabs(rho2-rho1)/Bo;
  TOLERANCE = 1e-6;
  run();
}

event init (t = 0) {
	if (!restore (file = "restart")) {
	refine (sq(x - S1  ) + sq(y)  - sq(1.5*D0) < 0 && level < (LEVELMAX));
    fraction (f, - (sq(x - S1 ) + sq(y) - sq(1)));
	}
}

event acceleration (i++) {
  face vector av = a;
  foreach_face(x)
    av.x[] -= 1.;
}

event adapt (i++) {
	
	double femax = 1e-3;
    double uemax = 5.0e-3;
	/**
	double xmax_ad = -HUGE;
	foreach(reduction(max:xmax_ad)) {
       if (f[] > 1.0e-4 && f[] < 1 && x > xmax_ad) {
		   xmax_ad = x;
        }
	}
	if (xmax_ad < 239.85)
		LEVELMAX = LEVELMAX_a;
	else{
		if (xmax_ad < 239.925)
			LEVELMAX = LEVELMAX_b;
		else
			LEVELMAX = LEVELMAX_c;
	}
	**/
	//adapt_wavelet ({f,u}, (double[]){femax,uemax,uemax,uemax}, LEVELMAX, 6);
	adapt_wavelet ((scalar *){f,u.x,u.y},(double []){femax,femax,femax},LEVELMAX,LEVELMAX-8);
	//adapt_wavelet ((scalar *){g2f,grad_u,grad_v,grad_w},(double []){femax,uemax,uemax,uemax},LEVELMAX,6);
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
  foreach(){
   if(f[]<1.0e-5)
     f[]=0.0;
  }
 //remove_droplets (f,true);
 // DT = 2e-5;
}

event properties (i++){
	vorticity (u, omega_i);  
	scalar g1f[],  g2f[];
	foreach(){
	g1f[] = sqrt(sq(f[1,0,0] - f[-1,0,0])+sq(f[0,1,0] - f[0,-1,0])+sq(f[0,0,1] - f[0,0,-1]));
	}
	boundary ({g1f});
	foreach(){
	g2f[] = sqrt(sq(g1f[1,0,0] - g1f[-1,0,0])+sq(g1f[0,1,0] - g1f[0,-1,0])+sq(g1f[0,0,1] - g1f[0,0,-1]));
	}
	boundary ({g2f});
	foreach(){
		omega_i[] = g2f[]>0?0:omega_i[];
	}
}

event information2 (t = 0; t <= MAXTIME; t=t+0.05) {
//event information (i++) {
  double xb = 0., yb = 0., sb = 0.;
  double vbx = 0., vby = 0.;
  double re_i = 0., re_o = 0., cd = 0.;
  double rmax_x = -HUGE, rmin_x = HUGE, rmax_y = -HUGE, omega_i_max = -HUGE, omega_o_max = -HUGE;
  double area = 0., ek = 0.; //sum_yp is the positive Omega_y, sum_yn is the negative Omega_y, sum_o is the total Omega, area is the interface area

  foreach(reduction(+:xb) reduction(+:yb)
	  reduction(+:vbx) reduction(+:vby)
	  reduction(+:sb)) {
    double dv = (f[])*2.0*pi*y*sq(Delta);
    xb += x*dv;
    yb += y*dv;
    vbx += u.x[]*dv;
    vby += u.y[]*dv;
    sb += dv;
  }
/**
  foreach(reduction(+:ek)) {
    double dv = (f[])*dv()*0.001+(1-f[])*dv();
    ek += 0.5*(u.x[]*u.x[] + u.y[]*u.y[])*dv;
    //ek + = dv();
  }
**/
 foreach(reduction(max:rmax_x) reduction(min:rmin_x) reduction(max:rmax_y) reduction(+:area)) {
    if (f[] > 1e-4 && f[] < 1) {
      coord p;
      coord n = mycs (point, f);
      double alpha = plane_alpha (f[], n);
      double s = plane_area_center (n, alpha, &p);
     //double rad  = sqrt(sq(x + Delta*p.x-xb/sb) + sq(y + Delta*p.y-yb/sb) + sq(z + Delta*p.z-zb/sb)); 
        double rad_x = x + Delta*p.x-xb/sb;
        double rad_y = y + Delta*p.y;  

      if (rad_x > rmax_x)
	rmax_x = rad_x;
      if (rad_x < rmin_x)
	rmin_x = rad_x;
      if (rad_y > rmax_y)
	rmax_y = rad_y;
	area += s*Delta*2.0*pi*y;
    }

 }
 foreach(reduction(max:omega_i_max)) {
    if (f[] > 0.99) {
	  double vor_i = omega_i[];
      if (vor_i > omega_i_max)
	omega_i_max = vor_i;
    }
 }
 foreach(reduction(max:omega_o_max)) {
    if (f[] < 1e-4) {
	  double vor_o = omega_i[];
      if (vor_o > omega_o_max)
	omega_o_max = vor_o;
    }
 }
 re_o = 2.0 * Ga * (vbx/sb) / sqrt(fabs(1.0-RHOR));
 re_i = re_o * RHOR / MUR;
 cd = 32.0 / 3.0 * Ga * Ga / re_o / re_o;
 omega_i_max = omega_i_max / (vbx/sb);
 omega_o_max = omega_o_max / (vbx/sb);

// 1-time 2-volume(*M_PI) 3area(*2M_PI) 4-x 5-y 6-vx 7-vy 8-dmax/2.0 9-dmin/2.0 
  fprintf (stderr,
	   "%.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.3d\n", 
	   t, sb, area, xb/sb, vbx/sb, (2*rmax_y)/(rmax_x - rmin_x), re_o, re_i, cd, omega_i_max, omega_o_max, Ga, Bo, LEVELMAX);
  fflush (stderr);

}

event snapshot1 (t = 0; t <= MAXTIME; t=t+1)
{
  char name[80];
  sprintf (name, "dump-%g", t);
  dump (file = name);
}

event snapshot2 (t = 0; t <= MAXTIME; t=t+1)

{
  char path[]="htg"; // no slash at the end!!
  char prefix[80];
  scalar omegaz[], max_f[];
  
  vorticity (u, omegaz);  
  /**
    foreach(){
		max_f[] = f[];

    for (int i = -1; i <= 1; i++) {
		for (int j = -1; j <= 1; j++) {
				max_f[] = max(f[i,j],max_f[]);
		}
    }
    omegaz[] = omegaz[]*(max_f[]>0?0:1);
	}
	boundary ({omegaz});
  **/
  sprintf(prefix, "data-%g", t);
  output_htg((scalar *) {f,p,omegaz,u.x,u.y},(vector *){uf}, path, prefix, i, t);
}

