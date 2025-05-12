#include "axi.h"
#include "compressible/two-phase.h"
#include "compressible/Mie-Gruneisen.h"
#include "contact.h"
#include "tension.h"
#include "compressible/tension.h"
#include "view.h"

int minlevel = 8;
int maxlevel = 12;

double rhoL = 1., rhoG = 0.001;
double pg0;
double p0 = 1.;
double pinf = 10.;
double tend = 1.;

double Rbub = 1.;
double lambda = 100.;
double Reynolds = 10000.*0.1;
double Web = 100./72.*1000.*0.1;

vector h[];
scalar keliq[];
double theta0;

p[right]    = dirichlet(pinf);
q.n[right]  = neumann(0.);

p[top]    = dirichlet(pinf);
q.n[top]  = neumann(0.);

p[left]    = neumann(0);
q.n[left]  = dirichlet(0.);
q.t[left]  = dirichlet(0.);
h.t[left] = contact_angle(theta0);

uf.n[bottom] = 0.;
uf.t[bottom] = dirichlet(0.);

void careful_refinement(){
  double xc = cos(theta0);
  refine ( level <= (10 - (sqrt(x*x + y*y + z*z)/(3.) > 1)*30*pow(sqrt(x*x + y*y + z*z)/(3.) - 1.,1) ));
  refine(level < maxlevel && sqrt(sq(x - xc) + sq(y)+ z*z) > (1. - 10.*sqrt(2.)*100./(1<<12)) && sqrt(sq(x - xc) + sq(y)+ z*z) < (1. + 10.*sqrt(2.)*100./(1<<12)));
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

FILE * fp, *fpr;

double findval(int lvl,int I,int J,long int nocells,int datal[]
	       ,int datai[],int dataj[],double datap[]){

  double pval;
  bool found =0;
  for(long int ii=0; ii<nocells; ++ii){
    if(lvl == datal[ii])
      if(I == datai[ii])
	if(J == dataj[ii])
	  {
	    pval = datap[ii];
	    found =1;
	    break;
	  }
  }
  if(found == 0)
    printf("MISTAKE");
  if(found == 0)
    fprintf(ferr,"%d %d %d\n",lvl,I,J);
  assert(found);
  return pval;
}

int main() {

  pg0 = p0 + 2./Web;
  tend = 1.5*0.915*sqrt (1./fabs(pinf - pg0));
  
  f.gradient = zero;
  f.height = h;
  f.sigma= 1./Web;
  
  mu1 = 1./Reynolds;
  mu2 = mu1*0.01;

  CFLac = 10;
    
  gamma2 = 1.4;
  gamma1 = 7.14;
  PI1 = 1./sq(0.007)/7.14 - pinf;
    
  L0 = lambda;  
  
  X0 = 0.;
  init_grid(1<<minlevel);

#if CASE2
  theta0=12.*pi/18.;
#else  
  theta0=7.5*pi/18.;
#endif
  run();

}

event init (i = 0) {
  careful_refinement();
  
  double xc = cos(theta0);
  vertex scalar phi[];
  foreach_vertex() {
    phi[] = sq(x-xc) + sq(y) - Rbub*Rbub;
  }
  fractions (phi, f);
  
  char namew[50];
  sprintf(namew,"../embedpotential/data-%2.2f.dat",theta0);
  fpr = fopen(namew,"r");
  
  long int nocells;
  fscanf(fpr,"%ld",&nocells);
  
  fprintf(ferr,"#%s %ld %ld\n",namew,nocells,grid->tn);
  
  int *datal = (int*)malloc(nocells * sizeof(int));
  int *datai = (int*)malloc(nocells * sizeof(int));
  int *dataj = (int*)malloc(nocells * sizeof(int));
  double *datap = (double *)malloc(nocells * sizeof(double));
  
  for(long int ii=0; ii<nocells; ii++){
    fscanf(fpr,"%d %d %d %lf",&datal[ii],&datai[ii],&dataj[ii],&datap[ii]);
  }
  
  foreach() {
    frho1[]  = f[]*rhoL;
    frho2[]  = (1. - f[])*rhoG;
    
    double pL = findval(point.level,point.i,point.j,nocells,datal,datai,dataj,datap);
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
}

scalar unorm[],fscalar[];
event adapt (i++) {  
  foreach(){
    unorm[] = norm(q)/(frho1[]+frho2[]);
    fscalar[] = f[];
  }
  boundary ((scalar *){fscalar,unorm});
 
  adapt_wavelet((scalar *){fscalar,unorm},(double[]){0.001,0.002},maxlevel = maxlevel);
}

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
  squares("p", min = 1., max = 10., linear=true, map = cool_warm);
  mirror({0,1}) {
    draw_vof("f",lw=4);
    squares("keliq", min = 0, max = 0.5, linear=true, map = cool_warm);
    draw_string (s, pos = 2, size = 100, lc = {255,255,255}, lw = 4);
    draw_string ("P/P0", size = 100, lc = {0,0,0}, lw = 4);
    draw_string ("KE/rho UC**2", pos = 3, size = 100, lc = {0,0,0}, lw = 4);
  }
  save("bubble.mp4");
}

/**
Bubble expansion and collapse near solid boundary

<video width="426" height="660" controls>
<source src="collapsewall/bubble.mp4" type="video/mp4">
</video>

<video width="426" height="660" controls>
<source src="collapsewall2/bubble.mp4" type="video/mp4">
</video>

 */

event facetwriting(t += tend/8.){
  output_facets(f,stdout);
}
  
event end (t=tend){}

/**
~~~gnuplot Facets at different times for 75 degrees contact angle
unset border

set style line 1 linecolor rgb '#ce2d4f' dt 1 linewidth 4 ps 0.2 pt 7

set key top outside horizontal font "sans-serif,10"
set size ratio -1

unset xtics
unset ytics
set arrow from -2,0 to 2,0 nohead lw 6 lc rgb "#000000"
set lmargin 1
set bmargin 2.5
set tmargin 3
set rmargin 2

p "out" u 2:1 w l ls 1 t "{/Symbol a} = 75^o",'' u (-$2):1 w l ls 1 notitle
     
~~~

~~~gnuplot Facets at different times for 120 degrees contact angle
unset border

set style line 2 linecolor rgb '#5aae61' dt 1 linewidth 4 ps 0.2 pt 7

set key top outside horizontal font "sans-serif,10"
set size ratio -1

unset xtics
unset ytics
set arrow from -2,0 to 2,0 nohead lw 6 lc rgb "#000000"
set lmargin 1
set bmargin 2.5
set tmargin 3
set rmargin 2

p "../collapsewall2/out" u 2:1 w l ls 2 t "{/Symbol a} = 120^o",'' u (-$2):1 w l ls 2 notitle
     
~~~

~~~gnuplot Kinetic energy

reset
set style line 1 linecolor rgb '#ce2d4f' dt 1 linewidth 4 ps 0.2 pt 7
set style line 2 linecolor rgb '#5aae61' dt 1 linewidth 4 ps 0.2 pt 7

set key top right horiz font "sans-serif,16"

set lmargin 9
set bmargin 4
set tmargin 3
set rmargin 2

set xlabel "t/t_c" font "sans-serif,16"
set ylabel "KE/p_{inf} V_0" font "sans-serif,16"

p "log" u ($2/0.915*sqrt(9)):(0.5*$4/10./0.459790698) w l ls 1 t "{/Symbol a} = 75^o",\
"../collapsewall2/log" u ($2/0.915*sqrt(9)):(0.5*$4/10./0.104141846) w l ls 2 t "{/Symbol a} = 120^o"
~~~
 */
