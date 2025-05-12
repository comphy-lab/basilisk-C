/**
# Plateau-Rayleigh instability for topological changes
*/


#include "axi.h" 
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "marangoni.h"
int LEVEL = 8;

#include "surfactant-transport.h"

#define RHOR 1.
#define MUR 0.01
double MUD = 1.;
#define Db 0.1
double Oh = 0.2;
#define sigmas 100.
#define Pe 10. //defination of Pe = gammadot*(4.*db)^2/D_s
#define gamma0 0.1
double beta = 0.8;

double MAXTIME = 110.;
double gammadot = 1.; //shin et al. shear rate 1/s
 double a0 = 0.02; 
    double kk = 2.*M_PI;
double ndtime = 1.;


int main(int argc, char * argv[]) {

  LEVEL = atoi(argv[1]);
  Oh = atof(argv[2]);
  beta = atof(argv[3]);


  N = 1 << 6;
  
    reinit_skip_steps = 10;

  MUD = Oh*pow(rho2*sigmas*Db, 0.5); //Shear rate

  rho1 = 0.01;      //outside jet
  rho2 = 1.;  //inside jet

  mu2 = MUD; //inside jet
  mu1 = MUD*MUR; //outside jet
  /**
    We reduce the tolerance on the Poisson and viscous solvers to
    improve the accuracy. */
  
    TOLERANCE = 1e-4;
    ndtime = kk/pow(sigmas*pow(kk,3)/rho2,0.5);

    D_s = pow(sigmas/kk/rho2,0.5) / Pe; //Surface Diffusivity of surfactant D_s
    // printf("gd = %g %g",gammadot,D_s);
    surfactant = 1;
    addmarangoni = 1; //when using without surfactant or Ma = 0

    

    run();
}


event init (i = 0) {
  // exit(0);

printf("ndtime = %g",ndtime);
  /**
    The droplet is centered on (0,0) and has a radius of Db */
    if (!restore (file = "restart")) {
      
     
      refine (y - (1.5*(Db/2. - a0*cos(kk*(x - L0/2.)))) < 0 && level < LEVEL);

      fraction (f, (y)  - (Db/2. - a0*cos(kk*(x-L0/2.))));
      //epsilon = (L0/(1 << grid->maxdepth))*0.75;



      event("properties2"); //call this so that pfield can be populated and then c1 can be initialized	

      foreach(){


	double deltas = (pfield[]*(1. - pfield[]))/EPSILON;


	c1[] = gamma0*(2 + cos(kk*(x - L0/2.)))*deltas;
    
	if (deltas > 1.e-1 && surfactant == 1 ){

	  double gamma = c1[]*4.*EPSILON;
	  sigmaf[] =   sigmas*(1.   - (beta * gamma/gamma0));

	}
	else
	  sigmaf[] = sigmas;
      

	boundary({sigmaf});
      }
    }
}


event stability (i++){
   printf("\n%d %g",i,t); fflush(stdout);

  foreach(){
       
    double deltas = (pfield[]*(1. - pfield[]))/EPSILON;

    if (deltas > 1.e-1 && surfactant == 1){
      double gamma = c1[]*4.*EPSILON;
	  sigmaf[] =   sigmas*(1.   - (beta * gamma/gamma0));
    }
    else
      sigmaf[] = sigmas;
      
  }
  boundary({sigmaf});
}


event adapt (i++) {

  double femax = 0.001;
  double uxemax = 0.01;
  double uyemax = 0.01;
//   double uzemax = 0.01;
  double c1emax = 0.001;
  double pfemax = 0.001;
 
 

  adapt_wavelet ({f,u,pfield,c1}, (double[]){femax,uxemax,uyemax,pfemax,c1emax}, LEVEL); 
   printf("\n%d %g",i,t); fflush(stdout);
}
#include "tag.h"
int count = 0;
double tprev = 0., yprev = 0., xprev = 0., zprev = 0.;
event vtk(t = 0.; t +=0.0001; t <= 10.*ndtime ){
  //event vtk( i++; t <= 56. ){
  
    foreach(){
      gamma2[] = c1[]*4.*EPSILON;
    }

      scalar m1[];
     foreach()
          m1[] = f[] < 0.99;

        int n1 = tag(m1);
        double v1[n1]; //volume of bubbles

        
    
    count++;
  }

