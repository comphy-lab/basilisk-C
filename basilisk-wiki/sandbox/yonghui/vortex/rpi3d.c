/**
# Rayleigh-Plateau instability for different initial amplitude

The dimensionless momuntun equation can be write in the form
$$ \bar \rho  \frac{D \bar u}{D \bar t }
= - \nabla \bar P + Oh \sqrt{ \frac{R_0}{\lambda} } \bar \mu \Delta \bar u 
+ \bar \gamma  \kappa \delta_s  n 
$$
with 
$Oh = \frac{\mu_0 U_0}{L_0^2}$ 
and the scaling law
$u_{sim}= \sqrt{\frac{\lambda}{R_0}} u_d$,
$l_{sim} = \frac{R_0}{\lambda} l_d$,
$t_{sim} = \sqrt{\frac{R_0^3}{\lambda^3}} t_d$,
$p_{sim} = \frac{\lambda}{R_0} p_d $.
axi.h will be used here as in the axisymmetric case
*/
//#include <gsl/gsl_fit.h>  //needed gsl lib
//#pragma autolink -lgsl -lgslcblas
//#include "grid/multigrid.h"
#include "grid/octree.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "tension.h"
#include "view.h"
/**
Here we define the amplitude of perturbation as $\epsilon$, $\frac{R_0}{\lambda}$ as LL, possible to change ratio of density and viscosity. a fixe mesh at $2^{LEV+1}$ will be used.*/

#define LL 9. // wave length/R_0
#define Oh 0.00001
#define epsilon 0.01
#define mur 1. // ratio of mu2 / mu1
#define rhor 1 // ratio of rho2 / rho  
#define tend 0.6
FILE *fp, *fp2, *fp3;
/**
##main event 

We set a periodic condition at left and right to simulate a axisymmetric cylinder with infinity length, here we varies the values of $\epsilon$ to get the results for different initial amplitude.
*/
double hh[99999],tt[99999];
int kkk,jjj;
int main() {
  fp=fopen("energy","w");
  fp2=fopen("r_data","w");
  fp3=fopen("r_slope","w");
  fprintf (fp3,"#alpha slope N oh\n");
  N = 64; 	
  origin (0.,-L0/2.,-L0/2.);
  periodic(right);
  rho1 = 1., mu1 = Oh*sqrt(1. / LL);
  rho2 = rhor, mu2 = mu1*mur;
  f.sigma=1.;  
  kkk=0;	
  jjj =0;
  fprintf(stderr,"now N= %d,test kkk %d jjj%d\n" ,N,kkk, jjj);
  run(); 
}

/**
## boundary conditions 
simply no penetrating wall at top */
u.n[top]  = dirichlet(0.);
u.n[bottom]  = dirichlet(0.);
u.n[front]  = dirichlet(0.);
u.n[back]  = dirichlet(0.);

p[top]   = neumann(0.);
p[bottom]   = neumann(0.);
p[front]   = neumann(0.);
p[back]   = neumann(0.);


/**
## init event 

We use the fraction to seperate the 2 liquid. for f(x)>0, it defines as liquid 1 with f = 1. some useless picture can be output for fun.*/
event init (t = 0) {
  refine( (sq(y)+sq(z) > sq(1./(LL + 1.)) ) && (sq(y)+sq(z) < sq(1./(LL - 1.)) )  && (level < 9)  );
  fraction (f,  sq((1.+ epsilon*cos(2.*M_PI*x))/LL ) - sq(y) -sq(z) ); 
  //general fraction field view
  output_ppm (f, file = "finit.png", linear = true, n=512, min = 0, max = 1);
  //contains details of meshes
  view(fov=10.,tx=-0.5,ty=-0.2,width=1280,height=640);
  draw_vof("f", lw = 1.5); 
  cells();
  save ("fim.png");

}

/**
## output 
The movies of surface variation will be output every $t / 400$ for fun.*/
event movies (t += tend/400; t <= tend) {
  char legend[2000]; //time
  sprintf(legend, "t = %0.2g, \n LL = %g", t, LL);
  //general view
  view(fov=10.,tx=-0.5,ty=-0.2,width=1280,height=640);
  squares("p");
  draw_vof("f", lw = 1.5); 
  draw_string(legend, 1, size = 30., lw = 2.);
  save ("pressure.mp4");
  //movie 2 for detatch details
  view(fov=3.,tx=-0.6,ty=-0.05,width=1280,height=640);
  box();
  squares("p",min=0,max=1);
  draw_vof("f", lw = 1.5); 
  draw_string(legend, 1, size = 30., lw = 2.);
  save ("detail.mp4");
} 

/**
Here we save the total kinematic energy for liquid 1 and the position of maximum radius,for the output data, it will be rescaled to dimensionless values. We can find out that the perturbation in time 0.3 tend to 0.4 tend is almost in the linear zone, we calculate the growthrate in this region*/

event datas (t += tend / 100.; t <= tend) {
  double etotal =0.;
  foreach()
    etotal += f[]*sq(Delta)*(sq(u.x[])+sq(u.y[]))/2.;

  scalar pos[];
  position(f,pos,{0,1});
  double hhh =statsf(pos).max;
  fprintf(stderr,"%d %g %g %g\n", i, t*pow(LL, 1.5), etotal*LL, hhh*LL-1.);
  fprintf(fp,"%d %g %g %g\n", N, t*pow(LL, 1.5), etotal*LL, hhh*LL-1.);
}


event slope (i += 1; t <= tend) {
  //we start calculate the slope
  if ((t*pow(LL, 1.5) > 10.) && (t*pow(LL, 1.5) < 14.)){
    scalar pos2[];
    position(f,pos2,{0,1});
    double hh2 =statsf(pos2).max*LL-1.;
    double tt2 =t*pow(LL, 1.5);

    hh[kkk] = log(hh2);	
    tt[kkk] = tt2;

    fprintf(fp2,"%d %d %g %g %g %g %g\n", N,kkk, t, hh2, tt2,hh[kkk],tt[kkk]);
    kkk += 1;	
  }

  view(fov=10.,tx=-0.5,ty=-0.2,width=1280,height=640);
  draw_vof("f", lw = 1.5); 
  cells();
  save ("fnow.png");

  if((t*pow(LL, 1.5) > 14.)  && (jjj ==0)){
    double c0, c1, cov00, cov01, cov11, sumsq;
    gsl_fit_linear (tt, 1, hh, 1, kkk-2, 
                    &c0, &c1, &cov00, &cov01, &cov11, 
                    &sumsq);
    fprintf (stdout,"#best fit: ln(h) = %g + %g t\n", c0, c1);
    fprintf (stdout,"#covariance matrix:\n");
    fprintf (stdout,"#[ %g, %g\n  #%g, %g]\n", 
             cov00, cov01, cov01, cov11);
    fprintf (stdout,"#sumsq = %g\n", sumsq);
    fprintf (fp3,"%g %g %d %g\n", 2.*M_PI/LL, c1, N, Oh);	
    jjj += 1;
  }

}

int k=0;
event snapshot (t += tend /20.)
{
  k += 1;
  char name[40];
  sprintf (name, "dump-%d", k);
  dump (file = name);
}

/**
event adapt (i++) {
  double tolerance = 1e-4;
  adapt_wavelet ({f}, &tolerance, 9);
  event ("coefficients");
}
*/
