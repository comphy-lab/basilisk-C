/**
# Rising bubble

An axisymmetric bubble is released in a square box and raises
under the influence of buoyancy.

We solve the incompressible, variable-density, Navier--Stokes
  equations with interfaces and variable surface tension.  We can used standard or "reduced"
  gravity. */


#include "grid/octree.h" 
#include "navier-stokes/centered.h"
#include "two-phase.h"
// #include "two-phase-clsvof.h"
//# include "reduced.h" //keep reduced above marangoni

#include "marangoni.h"
// vector iforce[];
// #include "integral.h"

#include "surfactant-transport.h"
#include "navier-stokes/perfs.h"
//#include <dirent.h>

#include <sys/stat.h>
//#include <unistd.h>
//#include "vtknew_cell.h"

/**
The boundary conditions are slip lateral walls (the default) and
no-slip on the right and left walls. */

// #if MOMENTUM
// q.t[right] = dirichlet(0);
// q.t[left]  = dirichlet(0);
// #else
// u.t[right] = dirichlet(0);
// u.t[left]  = dirichlet(0);

// #endif

/**
We make sure there is no flow through the top and bottom boundary,
  otherwise the compatibility condition for the Poisson equation can be
  violated. */
/*

c1[left] = dirichlet(0.);
c1[right] = dirichlet(0.);
c1[top] = dirichlet(0.);
c1[bottom] = dirichlet(0.);
c1[front] = dirichlet(0.);
c1[back] = dirichlet(0.);

pfield[left] = dirichlet(0.);
pfield[right] = dirichlet(0.);
pfield[top] = dirichlet(0.);
pfield[bottom] = dirichlet(0.);
pfield[front] = dirichlet(0.);
pfield[back] = dirichlet(0.);

f[left] = dirichlet(1.);
f[right] = dirichlet(1.);
f[top] = dirichlet(1.);
f[bottom] = dirichlet(1.);
f[front] = dirichlet(1.);
f[back] = dirichlet(1.);
*/

//Bottom boundary is symmetric
// c1[bottom] = dirichlet(0.);

#define RHOR 1000.
#define MUR 100.
int LEVEL = 8;
double Ga = 1.;
double Ma = 1.;
double Bo = 0.545;
double MAXTIME = 110.;
double R0 = 0.5;    //radius of the bubble
const double WIDTH = 20.;
const double Zi = -9.;
double cutoff = 0.1; //minimum surface tension 

// scalar sigmaf[];

int main(int argc, char * argv[]) {

  LEVEL = atoi(argv[1]);
  Ga = atof(argv[2]);
  Bo = atof(argv[3]);
  Ma = atof(argv[4]);
  /**
    The domain will span $[0:2]\times[0:0.5]$ and will be resolved with
    $256\times 64$ grid points. */

  size (WIDTH);
  origin (-L0/2, -L0/2, -L0/2);
  // init_grid (64);
foreach_dimension()
  periodic(right);

  N = 1 << 6;
  
  rho1 = 1.;
  rho2 = rho1/RHOR;
  // mu1 = pow(Mo/pow(Bo,3.) , 0.25);
   mu1 = 1./Ga;
  mu2 = mu1/(MUR);
  /**
    We reduce the tolerance on the Poisson and viscous solvers to
    improve the accuracy. */
  
    TOLERANCE = 1e-4;
	//NITERMAX = 20;
  // #if REDUCED
 // G.x = -1.;
 // Z.x = 1.;
  // #endif
  D_s = 0.1; //Surface Diffusivity of surfactant D_s
  surfactant = 1;
  addmarangoni = 1; //when using without surfactant or Ma = 0
  //  d.sigmaf = sigmaf;
   struct stat st = {0};
    char name[80];
    sprintf (name, "dat");
    char newname[80];
    sprintf (name, "dat");

    if (stat(name, &st) == -1)
      {
	mkdir(name, 0700);
      }

    /**
      Preserve old stats restoring the simulation with a timestamp suffixed. */		
      sprintf (name, "stats.dat");
	
    if (stat(name, &st) == 0)
      {
	sprintf (newname, "stats-%ld.dat",time(0));
	rename(name, newname);
      }
    sprintf (name, "perfs");
	
    if (stat(name, &st) == 0)
      {
	sprintf (newname, "perfs-%ld",time(0));
	rename(name, newname);
      }

    //Uncomment this line if using clsvof
  // tracers = list_concat(tracers, {c1,d2,pfield});

  run();
}

double gamma0 = 0.01, gamma_inf = 1.; // Tritox X fdhila et. al.


event init (i = 0) {
    // exit(0);
  
/*  struct dirent **namelist;
  int n = scandir("./dat/",&namelist,0,alphasort);
  char last_file_name[200];
  sprintf(last_file_name,namelist[n-1]->d_name);
  printf("\n %d %s",n,last_file_name);
  exit(0);*/
  /**
    The bubble is centered on (Zi,0) and has a radius of R0 */
    if (!restore (file = "restart")) {
      refine (sq(x - Zi) + sq(y) + sq(z)  - sq(0.75) < 0 && level < LEVEL);
      fraction (f, sq(x - Zi) + sq(y) + sq(z)  - sq(R0));
      //epsilon = (L0/(1 << grid->maxdepth))*0.75;

    

    event("properties2"); //call this so that pfield can be populated and then c1 can be initialized	
  foreach(){
     //d2[] = R0 - sqrt (sq(x - Zi) + sq(y) + sq(z));
     //pfield[] = 0.5*(1. - tanh((d2[])/2./epsilon));
     //pfield[] = clamp(pfield[], 0., 1.);	
    double deltas = (pfield[]*(1. - pfield[]))/EPSILON;

    // double gamma0 = gamma_inf*0.01;
    c1[] = gamma0*deltas;
      if (deltas > 1.e-2){
	double gamma = c1[]*4.*EPSILON;
	// sigmaf[] =   1./Bo * (1. + Ma_gd*pow(Mo*Bo,0.25) * log(max(1. - gamma/gamma_inf, 1.e-6)));
  	// sigmaf[] =   1./Bo   + Ma_gd*pow(Mo*Bo,0.25) * gamma/gamma0;
    sigmaf[] =   1./Bo*(1.   - tanh(Ma* Bo/Ga * gamma/gamma0));
//    if (sigmaf[] < cutoff/Bo ){
//      sigmaf[] = cutoff/Bo;
//    }
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
    av.x[] += (-1. + alphav.x[]*rho1);
  
  boundary ((scalar *){av});
}

event stability (i++){

  foreach(){
       
      double deltas = (pfield[]*(1. - pfield[]))/EPSILON;

      if (deltas > 1.e-1 && surfactant == 1){
	double gamma = c1[]*4.*EPSILON;
  sigmaf[] =   1./Bo*(1.   - tanh(Ma* Bo/Ga * gamma/gamma0));
//	if (sigmaf[] < cutoff/Bo ){
//	  sigmaf[] = cutoff/Bo;
//	}
      }
      else
        sigmaf[] = 0.;
      
  }
  boundary({sigmaf});
}


event adapt (i++) {

 double femax = 0.001;
 double uxemax = 0.01;
 double uyemax = 0.01;
 double uzemax = 0.01;
 double c1emax = 0.001;
 double pfemax = 0.001;
 
 

  adapt_wavelet ({f,u,pfield,c1}, (double[]){femax,uxemax,uyemax,uzemax,pfemax,c1emax}, LEVEL); 
//  printf("\n%d %g",i,t); fflush(stdout);
}

int count = 0;
double tprev = 0., yprev = 0., xprev = 0., zprev = 0.;
event vtk(t = 0.; t +=0.1; t <= MAXTIME ){
//event vtk( i++; t <= 56. ){
  
    foreach(){
      gamma2[] = c1[]*4.*EPSILON;
    }
    char filename[100];
    // sprintf(filename, "vtk/D3fdhila--Bo%g-Ma%g-Ga%g-%d-%d.vtk",Bo,Ma,Ga,(1 << LEVEL),count);
    // FILE * fileptr = fopen(filename,"w");
    
    // scalar * list = {f,pfield,c1,u.x,u.y,sigmaf,gamma2,d2,iforce.x,iforce.y};
    // output_vtk( list , fileptr);
    // fclose(fileptr);
    // get dump files too

    sprintf(filename, "dat/D3dump-fdhila--Ga%g-Bo%g-Ma%g-%d-%g",Ga,Bo,Ma,(1 << LEVEL),t);
    dump(filename);

    double xb = 0., yb = 0., zb = 0., vb = 0.;
    foreach (reduction(+:xb) reduction(+:yb) reduction(+:zb) reduction(+:vb)) {
      double dvb = (1. - f[])*dv();
      vb += dvb;          // volume of the bubble
      xb += x*dvb;	// location of the bubble
      yb += y*dvb;	// location of the bubble
      zb += z*dvb;	// location of the bubble

     
    }
    
    //sprintf(filename, "D3surf-fdhilla-velstats-Ga%g-Bo%g-Ma%g-%d.dat",Ga,Bo,Ma,(1 << LEVEL));

    static FILE * fileptr2 = fopen("stats.dat","w");

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

    //      printf("\n%d %g",i,t); fflush(stdout);

    count++;
  }

