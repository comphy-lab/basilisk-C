/**
# Rising bubble

An axisymmetric bubble is released in a square box and raises
under the influence of buoyancy.

We solve the incompressible, variable-density, Navier--Stokes
  equations with interfaces and variable surface tension. */

#include "axi.h"

#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "marangoni.h"
#include "surfactant-transport.h"
#include "navier-stokes/perfs.h"
#include <sys/stat.h>


#define RHOR 1000.
#define MUR 100.
double Pe = 100.[*] ;

int LEVEL = 8;
double Ga = 1.;
double Ma = 1.;
double Bo = 0.545;
double MAXTIME = 20.;
double R0 = 0.5;    //radius of the bubble
const double WIDTH = 20.;
const double Zi = -9.;
double cutoff = 0.1; //minimum surface tension 


int main(int argc, char * argv[]) {

  LEVEL = 9; 
  Ga = 100.;
  Bo = 0.5;
  Ma = 1.;
  /**
    The domain will span $[0:2]\times[0:0.5]$ and will be resolved with
    $256\times 64$ grid points. */

  size (WIDTH);
  origin (-L0/2, 0.);

  periodic(right);

  N = 1 << 6;
  
  rho1 = 1.;
  rho2 = rho1/RHOR;

  mu1 = 1./Ga;
  mu2 = mu1/(MUR);
  /**
    We reduce the tolerance on the Poisson and viscous solvers to
    improve the accuracy. */
  
    TOLERANCE = 1e-4 [*];

  D_s = 1./Pe; //Surface Diffusivity of surfactant D_s

  surfactant = 1;
  addmarangoni = 1; //when using without surfactant or Ma = 0


  run();
}

double gamma0 = 0.1, gamma_inf = 1.; // Tritox X fdhila et. al.


event init (i = 0) {
  
  /*  The bubble is centered on (Zi,0) and has a radius of R0 */
    if (!restore (file = "restart")) {
      refine (sq(x - Zi) + sq(y)   - sq(0.75) < 0 && level < LEVEL);
      fraction (f, sq(x - Zi) + sq(y)  - sq(R0));
     

    

    event("properties2"); //call this so that pfield can be populated and then c1 can be initialized	
  foreach(){

    double deltas = (pfield[]*(1. - pfield[]))/EPSILON;


    c1[] = gamma0*deltas;
    if (deltas > 1.e-2){
      double gamma = c1[]*4.*EPSILON;
      sigmaf[] =   1./Bo*(1.   - tanh(Ma* Bo/Ga * gamma/gamma0));
 
    }
    else
      sigmaf[] = 0.;
      
  }
  boundary({sigmaf});
 }
}
//Gravitation

event acceleration (i++) {
  
  face vector av = a;
  
  foreach_face(x)
    av.x[] += (-1. + alphav.x[]*rho1/max(y,1e-20));
  
  boundary ((scalar *){av});
}

event stability (i++){

  foreach(){
       
      double xi = (pfield[]*(1. - pfield[]))/EPSILON;

    if (xi > 1.e-1 && surfactant == 1){
	double gamma = c1[]*4.*EPSILON; 
      sigmaf[] =   1./Bo*(1.   - tanh(Ma* Bo/Ga * gamma/gamma0));
    }
    else
    sigmaf[] = 0.;
      
  }
  boundary({sigmaf});
}


event adapt (i++) {

double femax = 0.001;
 double uxemax = 0.001;
 double uyemax = 0.001;

 double c1emax = 0.001;
 double pfemax = 0.001;


  adapt_wavelet ({f,u,pfield,c1}, (double[]){femax,uxemax,uyemax,pfemax,c1emax}, LEVEL); 

}

int count = 0;
double tprev = 0., yprev = 0., xprev = 0., zprev = 0.;
event output_files(t = 0.; t +=0.01; t <= 0.01 ){

  
    foreach(){
      gamma2[] = c1[]*4.*EPSILON;
    }
    char filename[100];


    sprintf(filename, "D2dump-fdhila--Ga%g-Bo%g-Ma%g-%d-%g",Ga,Bo,Ma,(1 << LEVEL),t);
    dump(filename);

    double xb = 0., yb = 0., zb = 0., vb = 0.;
    foreach (reduction(+:xb) reduction(+:yb) reduction(+:zb) reduction(+:vb)) {
      double dvb = (1. - f[])*dv();
      vb += dvb;          // volume of the bubble
      xb += x*dvb;	// location of the bubble
      yb += y*dvb;	// location of the bubble
      zb += z*dvb;	// location of the bubble

     
    }
    
    sprintf(filename, "D2surf-fdhilla-velstats-Ga%g-Bo%g-Ma%g-%d.dat",Ga,Bo,Ma,(1 << LEVEL));

    static FILE * fileptr2 = fopen(filename,"w");

    double dxb_dt = 0., dyb_dt = 0.,dzb_dt = 0., tm = 0.;
    if (i > 0){
        dxb_dt = (xb/vb - xprev)/(t - tprev);
       dyb_dt = (yb/vb - yprev)/(t - tprev);
       dzb_dt = (zb/vb - zprev)/(t - tprev);
       tm = (t + tprev)/2.;
    }
    
    fprintf (fileptr2,"%.12e %.12e %.12e %.12e %.12e %.12e %.12e %.12e\n",t,xb/vb,yb/vb,zb/vb,tm,dxb_dt,dyb_dt,dzb_dt);

    fflush(fileptr2);

    xprev = xb/vb;
    yprev = yb/vb;
    zprev = zb/vb;
    tprev = t;


    count++;
  }
