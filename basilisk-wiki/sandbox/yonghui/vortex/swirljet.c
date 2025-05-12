/**
# Rayleigh-Plateau instability for an axisymmetric liquid jet with rotation

We are interested in the rotating effect on a capillary jet.
We use a linear code to solve the problem, and then implement the initial perturbation (velocity and interface) of the normal mode into this DNS code.

A rough test of Robin conditino can be used here.
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
$Oh = \frac{\mu_l }{\sqrt { \rho_l \gamma R_0}}$, the ratio of $\frac{\lambda}{R_0}$ as $LL$, the ratio of density $rhor$ and viscosity $mur$ for water/air. the $qq$ is the maximum angular velocity*/

#define LL 9. // wave length L /R_0
#define Oh 0.001
#define mur 0.02 // ratio of mu2 / mu1
#define rhor 1.e-3 // ratio of rho2 / rho
#define tend 1.
#define qq 9.8
#define robin 0

/**
##main event 

We set a periodic condition at left and right, the size of box is set as $\lambda$ with the initial mean radius $R_0 = 1$, the interface perturbation is of $1\%$ grid size.
*/
FILE *fp, *fp2, *fp3;
double hh[999999],tt[999999],epsilon;
scalar velow;
int kkk;
int main() {
  fp=fopen("all","w");
  fp2=fopen("linear","w");
  fp3=fopen("growth","w");
  size(LL);
  periodic(right);
  rho2 = rhor, mu2 = Oh*mur;
  rho1 = 1., mu1 = Oh;
  f.sigma=1.;  
  N = 256;
  kkk=0;	
  epsilon = 0.01*LL/N;
  fprintf(stderr,"now mu1 = %g,mu2 = %g\n" ,mu1,mu2);
  fprintf (fp,"\n  \n");
  fprintf (fp2,"\n  \n");	
  run(); 
}

/**
## boundary conditions 
we can apply a Robin conditions or a normal one.
*/
double fa1,fa2,fb,fc,ddd;
#if robin
u.t[top] = (2.*fb-fa1)*u.t[]/(2.*fb+fa1*ddd) ;
u.n[top] = (2.*fb-fa2)*u.n[]/(2.*fb+fa2*ddd) ;
#else
u.n[top] = dirichlet(0.);
#endif
p[top]   = neumann(0.);
w[top] = neumann(0.);

/**
## init event 
We then define the initial velocity as follows: $u_x = 0 + u_x'$, $v_r = 0 +v_r'$ and $w_\theta = W_\theta + w_\theta'$ as a vortex layer near the interface, with $dd$ the initial vortex size .

The values of initial perturbations as a function of radius is saved in "velobas.txt".

*/
event init (t = 0) {
  fraction (f, (1.+ epsilon*cos(2.*M_PI*x/LL)) -y ); 
  view(fov=10.,tx=-0.5,ty=-0.2,width=1280,height=640);
  draw_vof("f", lw = 1.5); 
  cells();
  save ("fim.png");

  //initial velocity
  double dd = 0.1;
  double ww = qq*dd*sqrt(M_PI)/2.;
  fprintf(stderr,"dd %g,ww%g\n" ,dd,ww);
  foreach(){
    velow[] = y < 1. ? (qq*y) :		\
      (qq/y*(1.+sq(dd)*(1.- exp(-sq((y-1.)/dd)))  \
             + dd*sqrt(M_PI)*erf((y-1.)/dd)));	
    w[] = velow[];
  }	
  //robin parmeters
  fa1 = (0.7 + 0.5/9.);
  fa2 = (0.7 + 0.5/9.);
  fb = 1.;
  fc = 0.;
  ddd = L0/N;
}
/**
We read the perturbation. The files contains the data correspond to a N=256 regular mesh.
*/
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

/**
movies: we output some movies to see the effect.
*/
event moviex(t += tend/600.){
  scalar omega[];
  vorticity (u, omega);
  cells();
  squares("p");
  draw_vof("f", lw = 1.5); 
  save ("pressure.mp4");

  clear();
  cells();
  squares("w");
  draw_vof("f", lw = 1.5); 
  save ("wtheta.mp4");

  clear();
  cells();
  squares("omega");
  draw_vof("f", lw = 1.5); 
  save ("omega.mp4");
}


/**
Here we save the total kinematic energy for liquid 1 and the position of maximum radius,for the output data, it will be rescaled to dimensionless values. we calculate the growthrate in the linear regime.*/

event slope (i += 3 ; t <= tend) {
  scalar pos2[];
  position(f,pos2,{0,1});
  double scqle = (statsf(pos2).max + statsf(pos2).min)-2.;
  double hh2 =statsf(pos2).max-1.;
  double tt2 =t;
  //we start calculate the slope
  //fprintf(stderr,"%d %g %g\n", i,t,scqle);
  fprintf(stderr," %g %g %g %g %g %g %d\n", tt2, hh2, scqle,
          (statsf(pos2).max-1.)*N, epsilon, L0, N);
  double endeta;
  if ((t > 0.2) && (fabs(scqle) < 1e-3) ){
    hh[kkk] = log(hh2);	
    tt[kkk] = tt2;
    fprintf(fp2,"%d %d %g %g %g %g\n", N, kkk, t, hh2, tt2,scqle);
    kkk += 1;
    endeta = (statsf(pos2).max-1.)*N;
  }

  if( ((t > 0.2) && (fabs(scqle) > 1e-3)) || (t > 0.8)){
    double c0, c1, cov00, cov01, cov11, sumsq;
    gsl_fit_linear (tt, 1, hh, 1, kkk-2, 
                    &c0, &c1, &cov00, &cov01, &cov11, 
                    &sumsq);
    fprintf (stderr,"#best fit: ln(h) = %g + %g t\n", c0, c1);
    fprintf (stderr,"#covariance matrix:\n");
    fprintf (stderr,"#[ %g, %g\n  #%g, %g]\n", 
             cov00, cov01, cov01, cov11);
    fprintf (stderr,"#sumsq = %g\n", sumsq);
    fprintf (fp3,"%g %g %d %g %g\n",  epsilon, L0, N, c1, endeta);
    exit(0);
  }
}
