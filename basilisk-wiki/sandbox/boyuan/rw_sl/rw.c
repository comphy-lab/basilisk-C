/**
# Roll-wave Simulation Using the Single-layer Non-hydrostatic Model

A trial to reproduce Brock's experimental results (1967). Code modified from [this test case](/src/test/gaussian.c) */

#include "grid/multigrid1D.h"
#include "layered/hydro.h"
#define phi q
#include "layered/nh.h"
#include "layered/remap.h"
#include "layered/perfs.h"

/**
The primary parameters defining the normal flow and some other parameters */
// problem-sepcific parameters
double So = 0.05011;
double normalDepth = 0.00798;
double normalVelocity = 1.0377;
double gravityCoeff = 9.81;
double disMag = 0.01;
double disPeriod = 0.933;
double simTime = 50.0;
double Qin;
double froudeNum=3.71;
/**
Balance between bed-friction and gravity:
$$
\frac{c_fU^2}{2H}=g\sin\theta
$$
$c_f$ evaluated in Init based on SWEs*/
double cf; // = gravityCoeff * So * 2. * normalDepth / (sq(normalVelocity));

/** **main params**

full-length domain is sufficiently long for roll-wave generation */
int main()
{
  L0 = 40.;
  G = 9.81;
  N = 4096; // less damping with 1024
  /**
No viscosity for now (single-layer model)*/
  
  nu = 0.0;
  nl = 1;
  
  // max_slope = 1.
 
  NITERMIN = 2;
  
  // use a smaller CFL number (based on hydrostatic models)
  CFL_H = 0.375;
  
  /** <span style="color:red"> **need tweak breaking parameter b?** </span> */

  run();
}

/**
## Initialisation and boundary conditions
Note that we use a single-layer model for now.

The inflow is a simple time-dependent Dirichlet BC. Outflow is Neumann zero consistent with supersonic outlet. */

event init (i = 0)
{
  // B.C.'s and I.C.
   Qin = normalDepth*normalVelocity;
   h[left] = dirichlet(normalDepth + disMag * normalDepth * sin(2. * pi * t / disPeriod));
   // For simplicity, BC for $u_n$ is set to constant
   u.n[left] = dirichlet(normalVelocity);
   // u.n[left] = dirichlet(froudeNum*(sqrt(gravityCoeff*(normalDepth + disMag * normalDepth * sin(2. * pi * t / disPeriod)))));
   
   //outlet B.C.
   u.n[right] = neumann(0.);
   h[right] = neumann(0.);
   
   cf = gravityCoeff * So * 2. * normalDepth / (sq(normalVelocity));
   
   // I.C.
  foreach() {
    zb[] = -So * x;
    h[] = normalDepth;
    u.x[] = normalVelocity;
  }

  /**
  In the non-hydrostatic case, a boundary condition is required for
  the non-hydrostatic pressure $\phi$ of each layer. */
  
  phi[right] = dirichlet(0.);
}

/**
## Quadratic Friction (source term integration)
A simple backward Euler integration for the friction term.*/

event friction(i++)
{
    // double uMed;
     
     foreach ()
     {
          double a = h[] < dry ? HUGE : 1. + (cf / (2.)) * dt * norm(u) / h[];
          // double a = 1. + (cf / (2.)) * dt * u.x[] / h[];
          foreach_dimension(){
//                   u.x[] /= a;
//               u.x[] = (u.x[] + G * So * dt) / (1. + (cf / (2.)) * dt * norm(u) / (h[]));
              // u.x[] = (u.x[] + G * So * dt) / a;
              u.x[] /= a;
          
//           uMed = u.x[] + dt * chezyBedFriction(u.x[], h[], cf, G, So);
//           uMed = (3. / 4.) * u.x[] + (1. / 4.) * uMed + (1. / 4.) * dt * chezyBedFriction(uMed, h[], cf, G, So);
//           u.x[] = (1. / 3.) * u.x[] + (2. / 3.) * uMed + (2. / 3.) * dt * chezyBedFriction(uMed, h[], cf, G, So);
              
        }
     }
     boundary((scalar *){u.x}); // note that the input should be a list (at least for 1d)
}

/**
## Viscosity and Miscellaneous
We can optionally add horizontal viscosity. <span style="color:red"> (not considered for now) </span> */
#if 0
event viscous_term (i++)
{
  // add horizontal viscosity (small influence)
#if HYDRO
  scalar * list = {u.x};
#else
  scalar * list = {u.x,w};
#endif
  scalar d2u[];
  foreach_layer()
    for (scalar s in list) {
      foreach()
	d2u[] = (u.x[1] + u.x[-1] - 2.*u.x[])/sq(Delta);
      foreach()
	u.x[] += dt*nu*d2u[];
    }
  boundary (list);
}
#endif

/**
*plot_layer* meaningless for single-layer model */

// #if 0
// #include "plot_layers.h"
// #endif

/**
## Outputs

* Water depth,  velocity, NL pressure profiles output by out-\*.txt.
* Gnuplot depth profiles in t%.0f.png graphs.*/

event output(t = 0; t <= simTime; t += 1)
{
     char name[80];
     sprintf(name, "out-%.0f", t);
     FILE *fp = fopen(name, "w");
     double H = 0.;
     foreach ()
     {
         foreach_layer()
             H += h[];
         fprintf(fp, "%g %g %g %g %g\n", x, H, u.x[], w[], phi[]);
         H = 0.;
          
     }
     // fprintf(fp, "\n");
     fclose(fp);
}

void plot_profile(double t, FILE *fp)
{
     fprintf(fp,
             "set term pngcairo enhanced size 800,600 font \",10\"\n"
             "set output 't%.0f.png'\n"
             "set title 't = %.2f'\n"
             "set xlabel 'x(m)'\n"
             "set ylabel 'h(m)'\n"
             "plot [0:][0:]'-' u 1:2 w l lw 2\n",
             t, t);
     foreach ()
     {
          fprintf(fp, "%g %g\n", x, h[]);
     }
     fprintf(fp, "e\n\n");
     fflush(fp);
}

event gnuplot(t = 0; t <= simTime; t += 1)
{
     static FILE *fp = popen("gnuplot 2> /dev/null", "w");
     plot_profile(t, fp);
}

/** 
## results 
<span style="color:red"> **TODO** </span> 

![](~/Desktop/res1.png)*/