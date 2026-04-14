/**
# Axisymetric buoyant plume

We want to study the dynamics of a distribued buoyant plume, near the injection site.

![Tracer](buoyant_plume_axi/buoyant_plume_axi/c.mp4)

Using the Boussinesq approximation, eqations for the problem are 

$$
\tilde\nabla\cdot\mathbf{\tilde u} = 0,
$$
$$
\partial_t\mathbf{\tilde u}+\tilde \nabla\cdot(\mathbf{\tilde u}\otimes\mathbf{\tilde u}) = 
- \tilde \nabla  \tilde p + \dfrac{1}{\mathrm{Re}} \tilde \Delta \mathbf{\tilde u} + \dfrac{1}{\mathrm{Fr}^2} \tilde c,
$$ 

$$
\partial_t \tilde c + \mathbf{\tilde u} \cdot \tilde \nabla \tilde c = \dfrac{1}{\mathrm{Re} \mathrm{Sc}} \tilde \nabla c.
$$


We will need the Navier-Stokes solver with advection-diffusion of tracer in cylindrical coordinates.
*/

#include "axi.h"
#include "navier-stokes/centered.h"
#include "tracer.h"
#include "diffusion.h"

/**
The dimensionless parameters are defined with respect to the injection:

$$ \mathrm{Re} = \dfrac{2 Q}{\pi r_0 \nu} $$
$$ \mathrm{Fr} = \dfrac{Q}{\pi r_0^2 \sqrt{g' 2 r_0}}$$
$$ \mathrm{Sc} = \dfrac{\nu}{D} $$
*/ 
// Dimensionless parameters
double Re = 1.;
double Sc = 1000.;
double Fr = 0.1;


// Domain definition
int maxlevel = 7;
double r0 = 1. [1];
double L_box = 5. [1];


// Scalar definition
scalar c[];
scalar * tracers = {c};

/** 
For the boundary conditions, we impose a Poiseuil flow with constant concentration at the input.
*/
c[left]   = dirichlet(fabs(y)<r0);
c[right]  = neumann(0);
u.n[left] = dirichlet((fabs(y)<r0)*(1-y*y/sq(r0))*1 [1,-1]);
u.t[left] = dirichlet(0);
u.n[right] = neumann(0.);

p[left]    = neumann(0.);
pf[left]   = neumann(0.);

p[right]   = dirichlet(0.);
pf[right]  = dirichlet(0.);


// Vector fields
face vector av[];
face vector muc[];
face vector D[];


//
// Main Function
//
int main(){
  size(L_box);
  init_grid(128);
  a = av;
  mu = muc;
  run();
}

event adapt(i++){
  adapt_wavelet ((scalar *){c,u}, (double[]){2.e-3,2.e-3,2.e-3},maxlevel,4);
}

event init(i=0){
  foreach()
    c[] = 0;
}


/**
Implementation of the viscosity, gravitational effect and diffusion.
We have,
$$ \mathrm{Re} = \dfrac{2 Q}{\pi r_0 \nu} 
\Rightarrow \mu = \dfrac{2 Q \rho}{\pi r_0 \mathrm{Re}} = \dfrac{2}{\pi \mathrm{Re}}$$

$$ \mathrm{Fr} = \dfrac{Q}{\pi r_0^2 \sqrt{g' 2 r_0}}
\Rightarrow g' = \dfrac{ Q^2 }{2 \pi^2 r_0^5 \mathrm{Fr}^2} = \dfrac{1}{2 \pi^2 \mathrm{Fr}^2}$$



*/
event acceleration (i++) {
  foreach_face(x)
    av.x[] +=  (c[] + c[-1])/2. * 1/ (2*sq(pi*Fr));
}

event properties (i++) {
  foreach_face(){
    muc.x[] = fm.x[]*2/pi/Re;}
}

// Diffusion du tracer
event tracer_diffusion (i++) {
    foreach_face()
        D.x[] = 1./(Re*Sc);
    diffusion(c, dt, D);
}

/**
Outputs of the simulation : 
*/
event output(t+=0.1){
  char name[80];
  //sprintf(name, "c-Re%g-Fr%g-maxlevel%d.mp4",Re,Fr,maxlevel);
  sprintf(name, "c.mp4");
   //output_ppm (c);
  output_ppm(c,file = name);
  printf("t= %f\n",t);
}

event end(t=2){
  printf("end");
}


/*
event data_dump (t += 2) {
  char name[80];
  sprintf(name, "c-Re%g-Fr%g-maxlevel%d-t%g.csv",Re,Fr,maxlevel, t);

  FILE * fp = fopen(name, "w");

  foreach() {
    fprintf(fp, "%g,%g,%g\n", x, y, c[]);
  }

  fclose(fp);
 }*/