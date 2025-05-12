#include "axi.h"
#include "../src/tpccontact.h"
#include "../src/compressible-tension.h"
#include "view.h"

int MINLEVEL = 7;
int MAXLEVEL = 11;
double CFLac = 0.4;

double rhoL = 1., rhoG = 0.001;
double pg0;
double p0 = 1./10.;
double pinf = 1.;
double tend = 1.5;

double Rbub = 1.;
double lambda = 100.;
double Reynolds = 10000.*0.1;
double Web = 100./72.*1000.*0.1;

scalar keliq[];

uf.n[right] = neumann(0.);
p[right]    = dirichlet(pinf);
q.n[right]  = neumann(0.);
f[right]    = dirichlet(1);
fE1[right] = neumann(0.);

uf.n[top] = neumann(0.);
p[top]    = dirichlet(pinf);
q.n[top]  = neumann(0.);
f[top]    = dirichlet(1.);
fE1[top] = neumann(0.);

uf.n[bottom] = 0.;
uf.t[bottom] = dirichlet(0.);

void careful_refinement(){
  for(int ii = MINLEVEL ; ii <= MAXLEVEL; ii++){
    refine(level < ii /*&& sqrt(sq(x) + sq(y)) > (Rbub - 4.*sqrt(2.)*lambda/(1<<(ii-1)))*/ && sqrt(sq(x) + sq(y)) < (2.5*Rbub + 4.*sqrt(2.)*lambda/(1<<(ii-1))));
  }
}

event stability(i++)
{
  double cspeed;
  foreach()
    {
      double fc = clamp (f[],0.,1.);
      double invgammaavg = fc/(gamma1 - 1.) + (1. - fc)/(gamma2 - 1.);
      double PIGAMMAavg = (fc*PI1*gamma1/(gamma1 - 1.) +
			   (1. - fc)*PI2*gamma2/(gamma2 - 1.));
      
      double cspeedsq = (p[]*(invgammaavg + 1.) + PIGAMMAavg)/invgammaavg/(frho1[]+frho2[]);
      if (cspeedsq > 0.)
	cspeed = sqrt(cspeedsq);
      else
	cspeed = sqrt(gamma1*(pinf + PI1));
      double dtmaxac = CFLac*Delta/cspeed;
      dtmax = min(dtmax,dtmaxac);
    }
}

int main() {

  tend *= 0.915;
  pg0 = p0 + 2./Web;
  
  f.gradient = zero;
  f.sigma= 1./Web;
  
  mu1 = 1./Reynolds;
  mu2 = mu1*0.01;
    
  gamma2 = 1.4;
  gamma1 = 7.14;
  PI1 = 1./sq(0.007)/7.14 - pinf;
    
  L0 = lambda;  
  
  X0 = 0.;
  init_grid(1<<MINLEVEL);

  run();
  
  free_grid();
}

event init (i = 0) {
  if (!restore (file = "restart"))
    {
      //outside c=1, inside 0
      careful_refinement();

      vertex scalar phi[];
      foreach_vertex() {
	phi[] = sq(x) + sq(y) - Rbub*Rbub;
      }
      boundary ({phi});
      fractions (phi, f);
      boundary({f});
      
      foreach() {
	frho1[]  = f[]*rhoL;
	frho2[]  = (1. - f[])*rhoG;
	
	double pL = pinf*(1. - Rbub/sqrt(sq(x) + sq(y))) + (pg0 - 2*f.sigma/Rbub/Rbub)*Rbub/sqrt(sq(x) + sq(y));
	double pg = pg0;
	
	fE1[]   = f[]*(pL/(gamma1 - 1.) + PI1*gamma1/(gamma1 - 1.));
	fE2[]   = (1.-f[])*pg/(gamma2 - 1.);
	
	double invgammaavg = f[]/(gamma1 - 1.) + (1. - f[])/(gamma2 - 1.);
	double PIGAMMAavg = (f[]*PI1*gamma1/(gamma1 - 1.) +
			     (1. - f[])*PI2*gamma2/(gamma2 - 1.));
	
	p[] = (fE1[] + fE2[] - PIGAMMAavg)/invgammaavg;
	q.x[] = 0.;
	q.y[] = 0.;
      }
      boundary ((scalar *){q,frho1,frho2,p,fE1,fE2});
      dump("init");
    }
}

// scalar unorm[],fscalar[];
// event adapt (i++) {  
//   foreach(){
//     unorm[] = norm(q)/(frho1[]+frho2[]);
//     fscalar[] = f[];
//   }
//   boundary ((scalar *){fscalar,unorm});
 
//   adapt_wavelet((scalar *){fscalar,unorm},(double[]){0.001,0.002},maxlevel = MAXLEVEL);
// }

event logfile (i++) {

  scalar pgas[];
  double volume = 0.;
  double ekmax = 1.e-20;
  
  foreach (reduction(+:volume) reduction(max:ekmax)){

   pgas[] = p[]*(1. - f[]);
   
   double Ek = 0.;
   foreach_dimension()
     Ek += sq(q.x[]);
   
   keliq[] = (Ek/(frho1[] + frho2[]))*f[];
   
   ekmax = max(ekmax,keliq[]);
   volume += dv()*(1. - f[]);
 }
 boundary({keliq,pgas});

 if(i == 0)
   fprintf(ferr,"#i \t t \t volume \t statsf(keliq).sum \t statsf(pgas).sum/volume \n");
   
 fprintf(ferr,"%d %10.9f %10.9f %10.9f %10.9f\n",i,t,volume,statsf(keliq).sum,statsf(pgas).sum/volume);

}

event movie (t += tend/90.) {
  char s[80];
  sprintf (s, "t = %.2f T_RP", t/0.915);
  clear();
  view(fov = 1., ty = -0.02, quat = {0,0,-cos(pi/4.),cos(pi/4.)}, width = 1980, height = 1980);
  draw_vof("f",lw=4);
  squares("p", min = 1./5., max = 1., linear=true, map = cool_warm);
  mirror({0,1}) {
    draw_vof("f",lw=4);
    squares("keliq", min = 0, max = 0.5, linear=true, map = cool_warm);
    draw_string (s, pos = 2, size = 100, lc = {255,255,255}, lw = 4);
    draw_string ("P/P0", pos = 3, size = 100, lc = {0,0,0}, lw = 4);
    draw_string ("KE/rho UC**2", size = 100, lc = {0,0,0}, lw = 4);
  }
  save("bubble.mp4");
}

/**

<center>
<video width="426" height="660" controls>
<source src="collapse/bubble.mp4" type="video/mp4">
</video>
</center>

 */

event end (t = tend){}

/**
Comparison with the Rayleigh-Plesset solution
~~~gnuplot Bubble radius as a function of time
set xlabel 't/T_{RP}'
set ylabel 'R/R_0'
plot "log" every 100 u 2:($3*3.)**(1./3.) w p pt 6 t "Basilisk",\
"../RP.dat" w l lw 2 t "RP"
~~~
*/