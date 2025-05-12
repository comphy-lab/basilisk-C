/**
# Rayleight Benard instability in 2D

My first code is inspired by the work of Andrés Castillo. I simulate a fluid heated from below and cooled from above to approach the convection of the Earth's mantle.

## General set-up
From the Basilisk source code, we include multigrid.h and perfs.h. We use multigrid.h because we have a non-cubic domains and MPI parallelism. 
*/
#include "grid/multigrid.h"
/**
For run the code in parallel with make using 4 process :
CC='mpicc -D_MPI=4' make rbi2d.tst
([Parallel informations](http://basilisk.dalembert.upmc.fr/src/Tips)).

## Model equations

We see that we obtain a system of equation dependent on $\sqrt{Ra}$ because we have adimensioned the velocity as follows: $[U] = \frac{\kappa}{H}\sqrt{\frac{g\beta\Delta T H^3}{\kappa\nu}} =  \frac{\kappa}{H}\sqrt{Ra}$

We use the same dimensionless
Boussinesq equations of Andrés :
$$
\nabla \cdot \vec{u} = 0
$$
$$
\partial_t \vec{u} + \nabla \cdot \left(\vec{u} \otimes \vec{u}
\right) = -\nabla p + \nabla \cdot \left( PrRa^{0.5} \nabla \vec{u} \right) +
Pr\theta\vec{e}_y
$$
$$
\partial_t \theta + \nabla \cdot \left( \vec{u} \theta \right) = \nabla \cdot
\left( Ra^{-0.5} \nabla \theta \right)
$$
Under the Boussinesq approximation our system
is fully defined in terms of two dimensionless
parameters: the *Rayleigh* and *Prandtl* numbers,
in addition to the geometry and boundary conditions.

To use the Boussinesq approximationq we use predefined Basiliks functions: Navier-Stokes centered equation and the advection-diffusion equation.
*/

	#include "convection_boussinesq.h"

/**

For our model we calculate the number of nusselt at the top and bottom of our domain.
Morevoer we use profil5.h to have acces of averaged profiles scalarfields in the y-direction. Then perfs.h is use for monitor the performance statistics.

*/

        #include "global_nusselt.h"
        #include "profil5.h"
        #include "navier-stokes/perfs.h"

/**

## Dimensionless parameters

MINLEVEL and MAXLEVEL are variables used for adaptive meshing but here we use a non-cubic domain (multigrid.h). Non-cubic and adaptive simulations are not possible yet. That is why they are equal. the function npe() gives the number of processor used and dimension(ny()) imposes the number of processes along y axis ([more details here](http://basilisk.dalembert.upmc.fr/src/Tips)). For this example we have a $2^{6} grid$.

DT is the maximum time step to help run the code.
Then TOLERANCE define the minimum resisdus to be reached ([Multigrid Poisson–Helmholtz solvers](http://basilisk.dalembert.upmc.fr/src/poisson.h#multigrid-solver)).

Here Ra=1e5 and Pr=1.

*/
#define MINLEVEL 8
#define MAXLEVEL 8

double EndTime= 300.;
int main() {
  size (npe());
  origin (-0.5, -0.5);
  dimensions (ny = 1);
  DT = 0.1;
  TOLERANCE = 1e-6;
  Ra = 1e5; Pr = 1.; N = 1 << MINLEVEL ; run();
}
/**

## Boundary conditions

We impose our temperature at the top and bottom of the domain and then for the side walls we have a neumann condition.
In addition to no-slip walls. 

*/
T[top] = dirichlet(-0.5);
T[left] = neumann(0.);
T[right] = neumann(0.);
T[bottom] = dirichlet(0.5);

u.n[top] = dirichlet(0.);
u.n[bottom] = dirichlet(0.);
u.n[right] = dirichlet(0.);
u.n[left] = dirichlet(0.);
u.t[top] = dirichlet(0.);
u.t[bottom] = dirichlet(0.);
u.t[right] = dirichlet(0.);
u.t[left] = dirichlet(0.);
/**

## Initial conditions

Initial conditions correspond to a linear temperature profile and no motion.

*/
event init (t=0) {
  foreach(){
    T[] = - y;
    foreach_dimension()
      u.x[] = 0.;
  }
  boundary ({T,u});
}
/**

## Outputs

We write in the log file all statistical quantities.
statsf() function returns the minimum, maximum, volume sum, standard deviation and volume for a field ([utils.h](http://basilisk.dalembert.upmc.fr/src/utils.h#simple-field-statistics)).
We write in the data file all our physical quantities allowing us to characterize our physical problem.

*/
event logfile (t += 1.0; t <= EndTime) {
  scalar div[];
  foreach() {
    div[] = 0.;
    foreach_dimension()
      div[] += u.y[1] - u.y[];
    div[] /= Delta;
  }
  stats s0 = statsf (div);
  fprintf (ferr, "%f %.9g %.9g %d %d\n",
	   t, s0.sum/s0.volume, s0.max, mgT.i, mgT.nrelax);
  output_ppm (T, file="temperature.mp4", n = 1024, box = {{-0.5,-0.5},{-0.5 + L0, 0.5}});
  double nu_vol=0. , nu_t=0. ,nu_b=0. ;
  nu_vol = nusselt_vol(T,u);
  nu_t = nusselt_top(T);
  nu_b = nusselt_bot(T);
  static FILE * fout = fopen ("data", "w");
  if (t<1.){
    fprintf (fout, "[1]Ra [2]Pr [3]N [4]nutop [5]nubot [6]nuvol [7]umin [8]umax [9]vmin [10]vmax [11]t\n");
  }
  stats velx = statsf (u.x), vely = statsf (u.y);
  fprintf (fout, "%.9g %.9g %d %.9g %.9g %.9g %.9g %.9g %.9g %.9g %f \n",
	   Ra,Pr,N,nu_t,nu_b,nu_vol, velx.min, velx.max, vely.min, vely.max, t);

}
/**

Furthermore, in the last time step of simulation averaged profile of temperature are calculated, we do this by evaluating the solution on an regular grid within the domain. To see if we have enought points about the boundary layer.

*/
event tempfile (t=EndTime) {
  profile({T},-0.5,0.5,t,1);
}
/**

## Results

~~~gnuplot
set output "Nu=f(t).png"
set xlabel 'times t'
set ylabel 'Nusselts'
titre=system("awk 'NR==2 {print $1}' data")
set title sprintf("Ra=%s",titre)
plot 'data' u 11:4 w l title 'Nu-top',\
      '' u 11:5 w l title 'Nu-bot'
~~~

![Temperature field.](rbi2d/temperature.mp4)

## References

*/