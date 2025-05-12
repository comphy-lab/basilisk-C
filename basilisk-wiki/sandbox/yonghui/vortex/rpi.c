/**
# Rayleigh-Plateau instability 

we want to study the growthrate of a axisymmetric Rayleigh-Plateau instability in the linear regime.
The dimensionless NS equation can be written in the form
$$ \bar \rho  \frac{D \bar u}{D \bar t }
= - \nabla \bar P + Oh \sqrt{ \frac{R_0}{\lambda} } \bar \mu \Delta \bar u 
+ \bar \gamma  \kappa \delta_s  n 
$$
with he Ohnesorge number
$Oh = \frac{\mu_1}{\sqrt {\rho_1 \gamma R_0}}$ as the reverse $Re$ number, 
here we proposed the caracteristic length $L_0 = \lambda$.

The values $f_{sim}$ of this simulation can be scaled to $f_d$ as
$u_{sim}= \sqrt{\frac{\lambda}{R_0}} u_d$,
$l_{sim} = \frac{R_0}{\lambda} l_d$,
$t_{sim} = \sqrt{\frac{R_0^3}{\lambda^3}} t_d$
if we propose $L_0 = R_0$
*/
#include <gsl/gsl_fit.h>  //need gsl lib
#pragma autolink -lgsl -lgslcblas
//#include "grid/multigrid.h"
#include "axi.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "tension.h"
#include "view.h"
// AMR switch
#define STATIC 1
#define ADAPT 0
#define LVMAX 9
#define LL 9. // $\lambda$ / $R_0$
#define Oh 1e-3 
#define mur 1. // ratio of mu2 / mu1
#define rhor 1. // ratio of rho2 / rho  
#define tend 0.55 // pinchoff moment ~ tend 1.5
FILE *fp, *fp2, *fp3;
/**
##main event 

The simulation is periodic at left and right.

We set the initial perturbation amplitude as $\epsilon = 0.01 \frac{\Delta x}{R_0}$ as a optimal choice of error study.
*/
double hh[99999],tt[99999],epsilon;
int kkk,jjj;
int main() {
  fp=fopen("fposition","w");
  fp2=fopen("lineardata","w");
  fp3=fopen("slope","w");
  fprintf (fp3,"#alpha slope N oh\n");
  
  #if STATIC
  N = 32; 
  #else
  N=pow(2,LVMAX);
  #endif  
  periodic(right);
  epsilon = 0.01*LL/(pow(2,LVMAX));
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
classic condition of no penetrating wall at top */
u.n[top]  = dirichlet(0.);
p[top]   = neumann(0.);


/**
## init event 

The two phases are seperated by the initial interface $r = r_0 (1+\epsilon cos(ks))$.

With the fraction field $f$, the inner liquid noted liquid 1 is presented by $f = 1$, and the external fluid noted liquid 2 is presented by $f = 0$.

We assume the error is related to the Taylor expension of the radial velocity at order $\alpha$, as presented by equation:
$$
 L = L_0 - int\left[\frac{1}{\alpha}\,log_2\left(\frac{{\frac{\partial^{\alpha}u}{\partial r^{\alpha}}}(r=r_0)}{\frac{\partial^{\alpha}u}{\partial r^{\alpha}}(r)}\right)\right]
$$
The theorytic profil of the radial velocity perturbation of RPI in the linear regime is under the form of Bessel function, such that $u_{r1}' = C_1 I_1(kr)$ and $u_{r2}' = C_2 K_1(kr)$. 

The grid redistribution with $\alpha = 2$ is implemented with the switch "STATIC".
*/
event init (t = 0) {
  #if STATIC
  double position[] = {0.7,1.13,1.78,2.66,3.70,4.85,6.05};
  for(int GG= 1 ; GG <= 6 ; GG++){
    refine(  (y >0.03 ) && (y < position[GG]/(2.*M_PI) ) 
           && (level <(LVMAX+ 1 -GG))  ); }
  refine ((y < 0.03 ) && (level <(LVMAX-1)));
  #endif  

  fraction (f,  1./LL*(1.+ epsilon*cos(2.*M_PI*x)) -y ); 
  //general fraction field view
  view(fov=10.,tx=-0.5,ty=-0.2,width=1280,height=640);
  draw_vof("f", lw = 1.5); 
  cells();
  save ("fim.png");
  fprintf(stderr,"now cell = %ld\n" , grid->tn);
}

/**
## output 
The movies of interface and pressure variation.*/
event movies (t += tend/100; t <= tend) {
  char legend[2000]; //time
  sprintf(legend, "t = %0.2g, \n LL = %g", t, LL);
  //general view
  view(fov=10.,tx=-0.5,ty=-0.2,width=1280,height=640);
  squares("p");
  draw_vof("f", lw = 1.5); 
  draw_string(legend, 1, size = 30., lw = 2.);
  save ("pressure.mp4");
  //movie 2 for detatch details
  //view(fov=3.,tx=-0.6,ty=-0.05,width=1280,height=640);
  //box();
  //squares("p",min=0,max=1);
  //draw_vof("f", lw = 1.5); 
  //draw_string(legend, 1, size = 30., lw = 2.);
  //save ("detail.mp4");
}

/**
Here we save the position of maximum radius.
A control parameter "scale" is used to verify if the simulation is still in the linear regime, that $\frac{\eta_1 -\eta_2}{r_0} = \frac{r_{max} + r_{min}}{r_0} -2 < 0.02$, with $\eta$ the radius perturbation at max/min radius. The simulation will stoped if the linear condition is not satisfied.

If the small perturbation theory is satisfied, the perturbation growth with time $\eta(r,t) = \epsilon(r) e^{st}$, with $s$ the growthrate.

We use a LLS as a linear regression method to calculate the growthrate
$$
ln \, (\frac{\lambda}{r_0} \eta) = C_{st} + \bar{s} ( t (\frac{\lambda}{r_0})^{\frac{3}{2}})
$$
*/
event slope (i += 1; t <= tend) {
  scalar pos2[];
  position(f,pos2,{0,1});
  double scale = (statsf(pos2).max + statsf(pos2).min)*LL-2.;
  double hh2 =statsf(pos2).max*LL-1.;
  double tt2 =t*pow(LL, 1.5);
  fprintf(fp," %g %g %g %g %g %g %d\n", tt2, hh2, scale,
          (statsf(pos2).max-1./LL)*N, epsilon, L0, N);
  double endeta;
  if ((t*pow(LL, 1.5) > 8.) && (fabs(scqle) < 1e-2) ){
    hh[kkk] = log(hh2);	
    tt[kkk] = tt2;
    fprintf(fp2,"%d %d %g %g %g %g\n", N, kkk, t, hh2, tt2,scale);
    kkk += 1;
    endeta = (statsf(pos2).max-1./LL)*N;
  }

  if( ((t*pow(LL, 1.5) > 8.) && (fabs(scale) > 1e-2)) || (t*pow(LL, 1.5) > 14.)){
    double c0, c1, cov00, cov01, cov11, sumsq;
    gsl_fit_linear (tt, 1, hh, 1, kkk-2, 
                    &c0, &c1, &cov00, &cov01, &cov11, 
                    &sumsq);
    fprintf (stderr,"#best fit: ln(h) = %g + %g t\n", c0, c1);
    fprintf (stderr,"#covariance matrix:\n");
    fprintf (stderr,"#[ %g, %g\n  #%g, %g]\n", 
             cov00, cov01, cov01, cov11);
    fprintf (stderr,"#sumsq = %g\n", sumsq);
    fprintf (fp3,"%g %g %d %g %g\n",  epsilon, L0, N, c1, endeta );	
    return 1;
  }

}
/**
Adapt wavelet function.

As for now, the non-static AMR seems worse than the static AMR in the linear regime.
*/
#if ADAPT
event adaptf(i++){
  doubel TOL = 1e-4;
  if (i>200){
    fprintf(fp2, "%d %g	%g %ld\n", i, t, TOL, grid->tn);
    adapt_wavelet ({f,u}, (double[]){TOL,TOL,TOL}, LVMAX,4);}
}
#endif
