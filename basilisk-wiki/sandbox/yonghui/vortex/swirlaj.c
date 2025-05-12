/**
# Rayleigh-Plateau instability for an axisymmetric air jet with rotation

We are interested in the offect of rotation on RPI, especially to a air column.

We use a linear code to solve its normal mode, and then implement the initial perturbation (velocity and interface) into this code as initial condition.
*/
#include <gsl/gsl_fit.h>  //needed gsl lib
#pragma autolink -lgsl -lgslcblas
//#include "grid/multigrid.h"
#include "axi.h"
#include "navier-stokes/centered.h"
#include "navier-stokes/swirl.h"
#include "two-phase.h"
#include "tension.h"
#include "view.h"
/**
Here we define the Ohnesorge number as 
$Oh = \frac{\mu_l }{\sqrt { \rho_l \gamma R_0}}$, the ratio of $\frac{R_0}{\lambda}$ as $LL$, the ratio of density $\rho_r$ and viscosity $\mu_r$ for air/liquid. the $q$ is the maximum vorticity*/


#define LL 7.85 	// 	wave length L / R_0
#define Oh 0.1		//	1 / Re
#define mur 1.e-3 	// 	ratio of viscosity 	gas / liquid
#define rhor 1.e-3 	//	ratio of density 	gas / liquid
#define tend 2.5
#define qq 0.5		//	maximum vorticity
#define READ 1		//	add perturbation
#define LIQUID 1	//	1:solid rotation 	0: vortex cone
#define ROBIN 1		//	1:Robin condition

/**
##main event 

We set a periodic condition at left and right, the size of box is set as $\lambda$ with the initial mean radius $R_0 = 1$, the interface perturbation is of 1\% grid size.
*/
FILE *fp, *fp2, *fp3;
double fa1,fa2,fa3,fap,bbx,fc1,fc2,fc3,fcp,ddd;
double prmax=-1.;
double hh[999999],tt[999999],epsilon;
scalar velow;
int kkk;
int main() {
  fp=fopen("w","w");
  fp2=fopen("wmax","w");
  fp3=fopen("wtest","w");
  size(LL);
  periodic(right);
  rho2 = 1., mu2 = Oh;
  rho1 = rhor, mu1 = Oh*mur;
  f.sigma=1.;  
  N = 1024;
  kkk = 0;	
  epsilon = 0.1*LL/N;
  fprintf(stderr,"epsilon = %g L0 = %g\n" ,epsilon,L0);
  fprintf (fp,"\n  \n");
  fprintf (fp2,"\n  \n");	
  run(); 
}

/**
## init event 

We can refine the mesh by proposing that the error is in the form $e = u''(y) \Delta x^2$, with $u (y)$ the theoretical profile in the form of Bessel I for the inner liquid and Bessel K for the external liquid.

We then define the base state as follows: $U_x = 0$, $V_r = 0$ and $W_ \theta$ and $\delta$ the vortex sheet thickness .*/
event init (t = 0) {
  fraction (f, (1.+ epsilon*cos(2.*M_PI*x/LL)) -y ); 
  view(fov=10.,tx=-0.5,ty=-0.2,width=1280,height=640);
  draw_vof("f", lw = 1.5); 
  cells();
  save ("fim.png");
  #if LIQUID
  //initial velocity for liquid/gas
  double dd2 = 0.1;
  foreach(){
    velow[] = y < 1. ? (qq*y) :		\
      (qq/y*(1.+sq(dd2)*(1.- exp(-sq((y-1.)/dd2)))  \
             + dd2*sqrt(M_PI)*erf((y-1.)/dd2)));	
    w[] = velow[];
  }
  //robin parmeters
  double Wrmax, Wrmaxp;	
  Wrmax = (qq/L0*(1.+sq(dd2)*(1.- exp(-sq((L0-1.)/dd2)))  \
                  + dd2*sqrt(M_PI)*erf((L0-1.)/dd2)));	
  Wrmaxp = qq*(-(1.+sq(dd2)) - dd2*sqrt(M_PI)*erf((L0-1.)/dd2) \
               + (sq(dd2)+2.*sq(L0))*exp(-sq((L0-1.)/dd2)))/sq(L0);
  #else
  //initial velocity for gas/liquid
  double dd1 = 0.3;
  double dd2 = 3.0;
  double rr1,rr2;
  double ww = qq/2.*(dd1*sqrt(M_PI) - sq(dd1)); // W_1 at interface
  double ot = qq*mur + 2.*(1.-mur)*ww ; //omega_2
  fprintf(stderr,"dd1 %g, dd2 %g, ww %g, omega2 %g\n" ,dd1,dd2,ww, ot);
  foreach(){
    rr1 = (y-1.)/dd1;
    rr2 = (y-1.)/dd2; 
    velow[] = y < 1. ? \
      ( qq/(2.*y)*( dd1*sqrt(M_PI)*(erf(rr1)+1.)	\
                   -sq(dd1)*exp(-sq(rr1)))	) :		\
      ((ww+ot*( dd2*sqrt(M_PI)/2.*erf(rr2)	\
               +sq(dd2)/2.*(1.-exp(-sq(rr2))) ))/y);	
    w[] = velow[];
  }	
  //robin parmeters
  double Wrmax,Wrmaxp;
  Wrmax = ((ww+ot*( dd2*sqrt(M_PI)/2.*erf((L0-1.)/dd2)	\
                   +sq(dd2)/2.*(1.-exp(-sq((L0-1.)/dd2))) ))/L0);	
  Wrmaxp = ((-ww+ot*(-dd2*sqrt(M_PI)/2.*erf((L0-1.)/dd2)	\
                     - sq(dd2)/2. + (sq(dd2)/2.+sq(L0)) \
                     *exp(-sq((L0-1.)/dd2)) ))/sq(L0));
  #endif
  fprintf(stderr,"W(L0) %g, W'(L0) %g\n" ,Wrmax, Wrmaxp);
 
  fa1 = (0.7 + 0.5/L0);
  fa2 = (0.7 + 0.5/L0);
  fa3 = (0.7 + 1.5/L0);
  fap = (0.7 + 0.5/L0);
  aax = fa1;
  bbx = 1.;
  fc1 = 0.;
  fc2 = 0.;
  fc3 = Wrmax*fa3 + Wrmaxp*bbx;
  fcp = prmax*fap + rho2*sq(Wrmax)/L0;
  ddd = L0/N;
  fprintf(stderr,"WRmax%g,Pm %g", 2.*fc3*ddd,prmax);
  //test of azimuthal velocity
  //  the azimuthal velocity plot
  view(fov=20.,tx=-0.5,ty=-0.5,width=640,height=640);
  clear();
  squares("w");
  draw_vof("f", lw = 1.5, max =0.6, min = -1.2); 
  save ("winit.png");
}

/**
## Boundary conditions 
Robin or free outflow condition
*/

#if ROBIN
u.t[top] = robin(fa1,bbx,fc1);
u.n[top] = robin(fa2,bbx,fc2);
w[top] = robin(fa3,bbx,fc3);
p[top] = robin(fap,bbx,fcp);
#else
u.n[top] = neumann(0.);
u.t[top] = neumann(0.);
w[top] = neumann(0.);
p[top] = dirichlet(prmax);
#endif
/**
We implement the initial disterbance correspoing to the initial interface disterbance.
The files contains the data correspond to a N=256 regular mesh.

This trick can help us avoid the transition regime at start as much as possible.
*/
#if READ
event readfile (t = 0) {
  int Ncell = 256;
  double initvelo[Ncell][4];
  for(int i=0; i < Ncell; i++){
    for(int j=0; j < 4; j++){
      initvelo[i][j] = 0.;
    }
  }

  FILE *file ;
  file = fopen("velobas.txt", "r");
  for(int i=0; i < Ncell; i++){
    for(int j=0; j < 4; j++){
      fscanf(file, "%lf",&initvelo[i][j]);
    }
  }
  fclose(file);

  foreach(){
    int pos = floor(y / L0 * Ncell);
    u.x[] = -1.*epsilon*initvelo[pos][1]*sin(2.*M_PI*x/LL);
    u.y[] = -1.*epsilon*initvelo[pos][2]*cos(2.*M_PI*x/LL);
    w[] += -1.*epsilon*initvelo[pos][3]*cos(2.*M_PI*x/LL);
  }	
}
#endif




/**
movies: we output some movies to see the motion.
*/
event moviex(t += tend/600.){
  scalar omega[];
  view(fov=20.,tx=-0.5,ty=-0.5,width=640,height=640);
  clear();
  cells();
  squares("w");
  draw_vof("f", lw = 1.5); 
  save ("wtheta.mp4");

  view(fov=20.,tx=-0.5,ty=-0.5,width=640,height=640);  clear();
  cells();
  squares("omega");
  draw_vof("f", lw = 1.5); 
  save ("omega.mp4");
}

/**
Here we save the total kinematic energy for liquid 1 and the position of maximum radius,for the output data, it will be rescaled to dimensionless values. We can find out that the perturbation in time 0.3 tend to 0.4 tend is almost in the linear zone, we calculate the growthrate in this region*/
event slope (i += 3 ; t <= tend) {
  scalar pos2[];
  position(f,pos2,{0,1});
  double hh2 = (statsf(pos2).max - statsf(pos2).min)/2.;
  double scale = fabs((statsf(pos2).max + statsf(pos2).min)-2.)/hh2;
  double tt2 =t;
  fprintf(stderr," %g %g %g %g %g %g %d\n", tt2, hh2, scale,
          (statsf(pos2).max-1.)*N, epsilon, L0, N);
  double endeta;
  //We start at t=0.1 to avoid the transition regime. 
  if (t > 0.1){
    hh[kkk] = log(hh2);	
    tt[kkk] = tt2;
    fprintf(fp2,"%d %d %g %g %g %g\n", N, kkk, t, hh2, tt2,scale);
    kkk += 1;
    endeta = (statsf(pos2).max-1.)*N;
  }
  //The endtime is set to 1.5, the scale < 1e-2 indicate the linear regime
  if( (t > 0.15) && (scale < 1e-2)){
    double c0, c1, cov00, cov01, cov11, sumsq;
    gsl_fit_linear (tt, 1, hh, 1, kkk-2, 
                    &c0, &c1, &cov00, &cov01, &cov11, 
                    &sumsq);
    fprintf (stderr,"#best fit: ln(h) = %g + %g t\n", c0, c1);
    fprintf (stderr,"#covariance matrix:\n");
    fprintf (stderr,"#[ %g, %g\n  #%g, %g]\n", 
             cov00, cov01, cov01, cov11);
    fprintf (stderr,"#sumsq = %g\n", sumsq);
    fprintf (fp,"%g %g %d %g %g\n",  epsilon, L0, N, c1, endeta );
    exit(0);
  }
}


