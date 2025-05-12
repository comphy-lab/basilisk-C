/**
# Rain in saint-venant : Fluvial Test case with rain and Darcy friction

In this example, we reproduce the MacDonald benchmark$^{[1]}$ and we add rain. 
We test the friction ans the rain in a fluvial flow. The leading equations can
be written as : $$ \begin{array}{cc} \partial_t h + \partial_x q = R(x,t)
\\ \partial_t q + \partial_x \left(\frac{q^2}{h}+ \frac{1}{2} g h^2
\right) = -gh (S_0 + S_f) \end{array} $$ $S_0$ is the slope, $S_f$
is the friction term in its kinematic form and R(x,t) the rain intensity.*/

/**
## Declarations
We call the saint venant solver in 1D and we add the Darcy friction term and the rain.*/

#include "grid/cartesian1D.h"
#include "saint-venant.h"
#include "saintvenant/darcy.h"
#include "saintvenant/rain.h"

int LEVEL;

scalar e[];
double e1 = 0., e2 = 0., emax = 0.;
int ne = 0;
double pause, tmax, q0 = 1, z0, zf, h0, tb = 2000, intens = 1e-3;

// Q is not cst
double q(double x){
  return q0+x*intens;
}
// Analytical solution for h and dh/dx
double hex(double x) {
  return pow(4/G,1/3.)*(1 + 0.5*exp(-16*pow(-0.5 + x/1000.,2)));
}

double dhex(double x) {
  return -0.016*pow(4/G,1/3.)*exp(-16*pow(-0.5 + x/1000.,2))*(-0.5 + x/1000.);
}

// Darcy friction term in kinematic formulation
double sf(double x){
  return  -f/(8*G)*q(x)*q(x)/pow(hex(x),3);
}

// Z and dz/dx
double dzex(double x) {
  return (q(x)*q(x)/(G*pow(hex(x),3))-1)*dhex(x)-2*q(x)*intens/(G*sq(hex(x)))+sf(x);
}
double zex(double x, double z){
  double dx = L0/N;
  return z + dx/4.*(dzex(x-dx)+2*dzex(x-0.5*dx)+dzex(x));
}

/**
## Parameters

Definition of parameters and calling of the saint venant subroutine run().*/

int main()
{
  f = 0.093;
  pause=0.2;
  L0 = 1000.;
  X0 = 0;
  G = 9.81;
  tmax = 2500;
  for( LEVEL = 4; LEVEL <= 8; LEVEL++){  
    e1 = e2 = emax = 0.;
    ne = 0;
    N = 1 << LEVEL;
    run();
    fprintf (stderr, "%d %g %g %g\n", N, e1/ne, sqrt(e2)/ne, emax/ne);
  }
}
/** 
## Boundary condition

We fix h and q (u) at both boundaries (fluvial case).*/
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
  // Because the slope is initially dry, we fix an artificial time-step. 
  DT = 1e-2;
  z0=0;
  foreach(){
    rain[] = intens;
    zb[] = zex(x,z0);
    z0=zb[];
    u.x[] = 0;
    h[]=0;
    zf=z0;
  }
  boundary(all);
}

/**
## Computing error

Noticing $tb$ the time when the flow is already stationary, we define
the following relative error norms [2]:

$$
|h|_1 =  \frac{\int_{tb}^T \int_0^L  |h() - hex()| dx dt}{(T-tb)*L} 
$$
$$
|h|_2 =  \frac{\sqrt{\int_{tb}^T \int_0^L  (h() - hex()^2 dx dt}}{(T-tb)*L}
$$
$$
|h|_{max} = \frac{\int_{tb}^T max(h() - hex())}{T-tb} 
$$

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
   if( N == 128){
    static FILE * fp1 = fopen("ErrorN128.dat","w");
    fprintf(fp1,"%g \t %g \t %g \t %g \n",t,no.avg,no.rms,no.max);
  }
}

/**
## Gnuplot output

We can use gnuplot to produce an animation of the water surface.*/
/*
event plot ( t <= tmax; t += 10 ) {
    if( N == 32 || N == 128 ){
      printf("set title 'Manning friction fluvial ----- t= %.3g , N = %i '\n"
	     "p[%g:%g][-5:1]  '-' u 1:($2+$4) t'free surface' w p pt 1,"
	     "'' u 1:4 t'topo' w l lt 4,"
	     "'' u 1:5 t'Analytical' w l lt -1 \n",t,N,X0,X0+L0);
      foreach()
	printf (" %g %g %g %g %g %g\n", x, h[], u.x[], zb[],hex(x)+zb[], t);
      printf ("e\n"
	      "pause %.4g \n\n",pause);
    }
}
*/
/**
Print the water profile along the channel at final time.*/

event printprofile ( t = tmax-1 ){
    char name[100];
    FILE * fp;
    sprintf(name,"profil-%i.dat",N);
    fp=fopen(name,"w");
    foreach() {
      fprintf(fp,"%g  \t %g  \t %g \t %g \t %g   \n"
	      ,x,h[],zb[],hex(x),u.x[]);
    }
    fclose(fp);
}

/**
## Results

![Error convergence](/rain_fluv_darcy_error.png)
![Error time](/rain_fluv_darcy_errt.png)
![Water depth profiles](/rain_darcy_fluv_height.png)

 
*/
/**
## References 
[1] I. MacDonald, M. Baines, N. Nichols, and P. G. Samuels, “Analytic Benchmark Solutions for Open-Channel Flows,” no. November, pp. 1041–1045, 1997.

[2] S. Popinet, “Quadtree-adaptive tsunami modelling,” Ocean Dyn., vol. 61, no. January, pp. 1261–1285, 2011.
*/
