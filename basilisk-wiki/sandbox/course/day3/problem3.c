/**
# Advection of a sharp function

We wish to compare the performance of various ways to compute the flux 
when advecting a sharp interface.

For simplicity, we will solve for the 1D problem
*/

#include "grid/cartesian1D.h"
#include "utils.h"

/** Assuming a compressible substance, in 1D the velocity has to be uniform.
In this case we will assume that is also constant*/
double uadv = 1.;

/** We define different ways to compute the flux knowing that u>0*/

void flux_centered (scalar f, scalar df, double dt)
{
  foreach()
    df[] = (f[1] + f[])/2.*uadv*dt;                                                           
}

void flux_upwind (scalar f, scalar df, double dt)
{
  foreach()
    df[] = f[]*uadv*dt;
}

void flux_upwind_second_order (scalar f, scalar df, double dt)
{
  foreach()
    df[] = (3*f[] - f[-1])/2.*uadv*dt;
}

void flux_VOF (scalar f, scalar df, double dt)
{
  foreach() {
    if (dt <= (1 - f[])*Delta/uadv) {
      df[] = 0.;
    } else {
      df[] = (f[] - 1.)*Delta + uadv*dt;
    }
  }
}

/** We are ready to write the main core of the program */

int main() {

  /** We initialize the main variables */

  int n = 100;        //grid points
  double CFL = 0.1;   //CFL number
  double t = 0., tend = 0.2; //end of the simulation     
  double x0 = 0.3;  //  initial interface possition

  double dt = CFL/n/uadv;
  init_grid (n);

  scalar f0[], f[], fexact[], df[], e[];
  foreach() {
    f0[] = (x < x0);
    fexact[] = (x < (x0 + uadv*tend));
    f[] = f0[];
  }

  /** we start advecting the interface using centered flux calculation */
  for (int i = 1; t <= tend; i++) {

    // we obtain the flux
    flux_centered (f, df, dt);

    // numerical solution at t+dt
    foreach()
      f[] -= (df[]-df[-1])/Delta;

    t += dt;
  }

  /** we write the solution at the end of the advection */
  FILE * fp    = fopen ("centered.dat", "w");
  foreach () 
    fprintf(fp, "%g %g %g \n", x, f[], fexact[]);

  fprintf(fp, " \n");
  fclose(fp);


  /** We try now the 1st order Upwind scheme

    We reinitialize the variables and then we proceed with the advection */
  t = 0.;
  foreach()
    f[] = f0[];

  for (int i = 1; t <= tend; i++) {

    flux_upwind (f, df, dt);

    // numerical solution at t+dt
    foreach()
      f[] -= (df[]-df[-1])/Delta;
    t += dt;
  }

  //we write the solution at the end of the advection
  fp    = fopen ("upwind.dat", "w");
  foreach ()  
    fprintf(fp, "%g %g %g \n", x, f[], fexact[]);
  fprintf(fp, " \n");
  fclose(fp);

  /** We try now the 2st order Upwind scheme

    We reinitialize the variables and then we proceed with the advection */

  t = 0.;
  foreach()
    f[] = f0[];
  for (int i = 1; t <= tend; i++) {

    // we obtain the flux
    flux_upwind_second_order (f, df, dt);

    // numerical solution at t+dt
    foreach()
      f[] -= (df[]-df[-1])/Delta;
    
    t += dt;
  }

  //we write the solution at the end
  fp    = fopen ("upwind_2nd.dat", "w");
  foreach ()  
    fprintf(fp, "%g %g %g \n", x, f[], fexact[]);
  fprintf(fp, " \n");
  fclose(fp);

  /** We try now the Volume of Fluid Method (VOF)

    We reinitialize the variables and then we proceed with the advection */

  t = 0.;
  foreach()
    f[] = f0[];
  for (int i = 1; t <= tend; i++) {

    // we obtain the flux
    flux_VOF (f, df, dt);

    // numerical solution at t+dt
    foreach()
      f[] -= (df[]-df[-1])/Delta;

    t += dt;
  }

  //we write the solution at the end
  fp    = fopen ("VOF.dat", "w");
  foreach ()  
    fprintf(fp, "%g %g %g \n", x, f[], fexact[]);
  fprintf(fp, " \n");
  fclose(fp);
}

/** After the program has finished we plot the function and the error

~~~gnuplot Function after advection
set output 'function.png'
set xlabel 'x'
set ylabel 'f(x)'
set key right
set grid
p "centered.dat" u 1:2 t 'centered' w lp, "upwind.dat" u 1:2 t 'upwind 1st' w lp, "upwind_2nd.dat" u 1:2 t 'upwind 2nd' w l, "VOF.dat" u 1:2 t 'VOF' w lp, "centered.dat" u 1:3 not w l lc 0 lw 2
~~~ 

~~~gnuplot Spatial distribution of the error
set output 'error.png'
set xlabel 'x'
set ylabel 'f-f_{exact}'
p "centered.dat" u 1:(abs($2-$3)) t 'centered' w l, "upwind.dat" u 1:(abs($2-$3)) t 'upwind 1st' w l, "upwind_2nd.dat" u 1:(abs($2-$3)) t 'upwind 2nd' w l, "VOF.dat" u 1:(abs($2-$3)) t 'VOF' w lp
~~~ 

The VOF method solves the numerical issues related to the advection of a sharp interface naturally. No numerical diffusion 
and accurate interface representation.
*/  
