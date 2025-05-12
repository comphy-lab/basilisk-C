/**
# 2D Kelvin-Helmholtz 

This is a shear flow instability simulation created by flows with different velocities and densities : 
the Kelvin-Helmholtz instability in 2D. The goal is to reproduce some of the tests that were made in the article : " Arrufat, 
T., Crialesi-Esposito, M., Fuster, D., Ling, Y., Malan, L., Pal, S., Scardovelli, R., Tryggvason, G. and Zaleski, S., 2021. 
A mass-momentum conserving, Volume-of-Fluid method for incompressible flow on staggered grids. Computers & Fluids, 215, p.104785.",
to test the Basilisk solving methods
*/

#include "navier-stokes/centered.h"
#include <complex.h>
#include "vof.h"
#include "tag.h"
#include "two-phase.h"

/**
Uncomment for conserving method, which implements the momentum-conserving VOF advection of the velocity components.
*/

/*#include "navier-stokes/conserving.h" */

/**
We define air and water densities in this simulation.
*/

#define rho_1 (1e1)
#define rho_2 (1e0)
#define L0 (2.)

double Pi = 3.141592653589793238;
double U_init = 0.5;
/*double nwnr = 2.;
double lambda = 2./ nwnr; */
double kwnr = 6.28318530718 * (L0/2.); // 2*pi/lambda 
double Ah = 1e-4;
double r = rho_2/rho_1;
  
int maxlevel;
double vis, Re; 
FILE * fpn;
char fnamev[99];
char fnameg[99];

/*char name[99];*/

scalar omega[], lev[], k[];
scalar psi[];

/*
Function reuturning half the size of a cell according to the mesh chosen.
*/
double halfcell ()
  {
    double hshift = 0.;
    foreach ()
      hshift = Delta/2.;
    return hshift; 
  }


/**2
The domain is periodic in $x$ and resolved using 128$^2$ points. */

int main(){
  periodic(left); 

  //Psi boundaries
/*  L0 = 1.;*/
  psi[bottom] = dirichlet(0.5*(L0/2)); // -U_init * 0.5 or 1  L0/2
  psi[top] = dirichlet(0.5*(L0/2)); // U_init * 0.5 or 1
  X0 = -0.5*L0;
  Y0 = -0.5*L0; // flow position 
  size (L0); // domain size 

  rho1 = rho_1 ;
  rho2 = rho_2 ;
  maxlevel =  4; 

/*  display_control (maxlevel, 4, 12);  //jview maxlevel control*/

  N = 64*L0;

  run();

}

double shape (double x, double y)
  {
    double omi = 2. * kwnr * U_init * sqrt(r) / (1.+r);  //  imaginary part
    double omr = (r-1.) * kwnr * U_init/(1. + r); // real part
    double complex smallomega = omr + omi * I;
    return -y + halfcell() + Ah * (cos(-kwnr*x) + I * sin(-kwnr*x))*(cos(-smallomega*t) + I * sin(-smallomega*t)); 
  }

/**
We don't involve viscosity effects.
*/

event init(t = 0.){ // time initialization
  Re = 0.;
  vis = 0.;
  const face vector muc[] = {0., 0.};
  mu = muc ;

/*  display_control (r, 1e-15, 1.);  //jview ratio control*/



//a = 0

// frequency definition
  double omi = 2. * kwnr * U_init * sqrt(r) / (1.+r);  //  imaginary part
  double omr = (r-1.) * kwnr * U_init/(1. + r); // real part
  double complex smallomega = omr + omi * I;

// Variables

  double complex A1 = ( U_init - (smallomega/kwnr) )*Ah ; 
  double complex B2 = -( U_init + (smallomega/kwnr) )*Ah;

/**
# Psi : stream function definition
*/

  foreach() 

  {

  // modexp function
    double ymax = L0/2;
    double ymin = -L0/2;
    double modexp1 = (exp(-kwnr*y) - exp(-kwnr*ymax))/(1. - exp(-kwnr*(ymax-ymin))) ;
    double modexp2 = (exp(kwnr*y) - exp(-kwnr*ymax))/(1. - exp(-kwnr*(ymax-ymin))) ;

  // Psi function 
/**
Only the real part of the complex number is taken into account in the exponential with the <complex.h> package,
therefore we are using the cosinus and sinus form.
*/

    double complex psiup = U_init * y + creal (A1 * (cos(-kwnr*x) + I * sin(-kwnr*x))) * modexp1; // 0.0001 * (cos(4Pi)+Isin(4Pi))=0
    double complex psidown = -U_init * y + creal(B2 * (cos(-kwnr*x) + I * sin(-kwnr*x))) * modexp2;

/* 
We leave the separation at 0 because a cell can't have two different values, not like the interface.
*/
    psi[] = (y > 0. ?  psiup : psidown);

}

  boundary ({psi});


/**
# The stream function Psi allows to define velocity fields.
*/

  foreach() {

    u.x[] =  (psi[0, 0] - psi[0, -1])/(Delta);
    u.y[] = -(psi[0, 0] - psi[-1, 0])/(Delta);

  }

  boundary({u.x});
  boundary({u.y});

  fraction (f,shape(x,y));

  sprintf (fnamev, "KH%g.mp4", Re);
  sprintf (fnameg, "KHlev%g.mp4", Re);
  char name1[99];
  sprintf (name1, "nrofvortices%g", Re);
  fpn = fopen (name1, "w");

}



event measurements (t += 0.001)
{
    stats func = statsf(u.y);
    double vmax  = func.max;
    double thgrowth = Ah *  exp((kwnr * U_init )* t);     
    fprintf (stderr,"%g %g %g\n", t, vmax,thgrowth);

}


/**
We compute the vorticity.
*/

event output (t = 0.; t += 0.001){

  int n = 0;
  double m = 0;
  foreach(reduction(+:n) reduction(max:m)){
    n++;
    lev[] = level;
    omega[] = (u.x[0,1]-u.x[0,-1] - (u.y[1,0]-u.y[-1,0]))/(2*Delta);
    if (fabs(omega[]) > m)
      m = fabs(omega[]);
  }
  boundary({omega});
  output_ppm (omega, file = fnamev, n = 512, min = -10, max = 10, linear = true);
  output_ppm (lev, file = fnameg, n = 512, min = 2, max = 11);
  foreach()
    k[] = (omega[] > m/3.); //1 or 0
  int nrv = tag(k);
  fprintf (fpn,"%g\t%d\t%g\t%d\t%d\t%d\n",t, nrv, m,
           n,((1 << (maxlevel*dimension)))/(n), i);
}


/**
#  Kinetic energy
*/

/*static double energy()
{
  double ke = 0.;
    foreach(reduction(+:ke))
      ke += (sq(u.x[]) + sq(u.y[]))/2.*sq(Delta);
  return ke;
}

event kinetic_energy (i += 25)
{
  fprintf (stderr, "%g %g\n", t, energy());
}*/

event end (t = 2.5)
  fclose(fpn);


   