/**
# Rain in saint-venant : Torrential test case with rain and Darcy

In this example, we reproduce the MacDonald benchmark$^{[1]}$ with
rain. We test the friction and the rain in a torrential flow. The
leading equations can be written as :

$$ 
\begin{array}{cc} \partial_t h + \partial_x q = R(x,t)
\\ \partial_t q + \partial_x \left(\frac{q^2}{h}+ \frac{1}{2} g h^2
\right) = -gh (S_0 + S_f) \end{array} 
$$ 

$S_0$ is the slope, $S_f$ is the friction term in its kinematic
form and R(x,t) the rain intensity.*/
/**   
## Declarations

We call the saint venant solver and we add the Darcy friction term.*/
#include "grid/cartesian1D.h"
#include "saint-venant.h"
#include "saintvenant/darcy.h"
#include "saintvenant/rain.h"

int LEVEL;

scalar e[];
double e1 = 0., e2 = 0., emax = 0.;
int ne = 0;
double pause, tmax, q0 = 2.5, z0, h0, tb = 300, intens = 1e-3 ;

// Q is not cst
double q(double x){
  return q0+x*intens;
}

// Analytical solution for h and dh/dx
double hex(double x){
  return pow(4/G,1/3.)*(1 - 0.2*exp(-36*pow(-0.5 + x/1000.,2)));
}
double dhex(double x){
  return pow(4/G,1/3.)*0.0144*exp(-36*sq(-0.5 + x/1000.))*(-0.5 + x/1000.);
}

// Darcy friction term in kinematic formulation
double sf(double x){
  return  -f/(8*G)*q(x)*q(x)/pow(hex(x),3);
}

// Z and dz/dx
double dzex(double x) {
  return (q(x)*q(x)/(G*pow(hex(x),3))-1)*dhex(x)-2*q(x)*intens/(G*sq(hex(x)))+sf(x);
}
double zex(double x, double z) {
  double dx = L0/N;  return z + dx/3.*(dzex(x-dx)+dzex(x-0.5*dx)+dzex(x));
}

/**
## Parameters

Definition of parameters and calling of the saint venant subroutine run().*/
int main()
{
  f = 0.065;
  pause=0.2;
  L0 = 1000.;
  X0 = 0;
  G = 9.81;
  tmax = 400;
  for( LEVEL = 6; LEVEL <= 10; LEVEL++){
    e1 = e2 = emax = 0.;
    ne = 0;
    N = 1 << LEVEL;
    run();
    fprintf (stderr, "%d %g %g %g\n", N, e1/ne, sqrt(e2)/ne, emax);
  }
}

/** 
## Boundary condition*/

h[left] = dirichlet(max(hex(0),0));
eta[left] =  dirichlet(max(hex(0)+zb[],zb[]));
u.n[left] = dirichlet(max(q0/hex(0),0));

h[right] = dirichlet(max(hex(1000),0));
eta[right] =  dirichlet(max(hex(1000)+zb[],zb[]));
u.n[right] = dirichlet(max(q(L0)/hex(1000),0));

/**
## Initial conditions
*/
event init (i = 0)
{
  DT = 1e-2;
  z0=0;
  foreach(){
    rain[] = intens;
    zb[] = zex(x,z0);
    z0=zb[];
    u.x[] = 0;
    h[]=0;
  }
  boundary(all);
}

/**
## Computing error
*/
event error (i++; t<=tmax) {
  foreach()
      e[] = (h[] - hex(x));
  norm no = normf (e);
  if(t > tb){
    e1 += no.avg;
    e2 += no.rms*no.rms;
    ne++;
    if (no.max > emax)
      emax = no.max;
  }
   if( N == 256){
    static FILE * fp1 = fopen("ErrorN512.dat","w");
    fprintf(fp1,"%g \t %g \t %g \t %g \n",t,no.avg,no.rms,no.max);
  }
}
/**
## Gnuplot output

We can use gnuplot to produce an animation of the water surface.*/
/*
event plot ( t <= tmax; t += 15 ) {
  if( N == 64 || N == 1024 ){
    printf("set title 'Manning friction fluvial ----- t= %.3g , N = %i '\n"
	   "p[%g:%g][-60:1]  '-' u 1:($2+$4) t'free surface' w p pt 1,"
	   "'' u 1:4 t'topo' w l lt 4,"
	   "'' u 1:5 t'Analytical' w l lt -1 \n",t,N,X0,X0+L0);
    
    foreach()
      printf (" %g %g %g %g %g %g\n", x, h[], u.x[], zb[], hex(x)+zb[], t);
    
    printf ("e\n"
	    "pause %.4g \n\n",pause);
  }
}
*/
/**
Print the water profile along the channel at final time.
*/

event printprofile ( t = tmax-1){
    char name[100];
    FILE * fp;
    sprintf(name,"profil-%i.dat",N);
    fp=fopen(name,"w");
    foreach() 
      fprintf(fp,"%g  \t %g  \t %g \t %g \t %g \t %g  \t %g \n"
	      ,x,h[],zb[],hex(x),u.x[],u.x[]*h[], h[] > dry ? u.x[]/(sqrt(G*h[])) : 0);
    fclose(fp);
}

/**
## Results

![Error convergence](/rain_darcy_tor_error.png)
![Error time](/rain_darcy_tor_errt.png)
![Water depth profiles](/rain_darcy_tor_height.png)


*/
/**
## References 
[1] I. MacDonald, M. Baines, N. Nichols, and P. G. Samuels, “Analytic Benchmark Solutions for Open-Channel Flows,” no. November, pp. 1041–1045, 1997.
*/
