/**
The two-phase Navier-Stokes solver is included.
*/

#include "axi.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "tension.h"
#define L_domain 50.0
#define L0 40.0
/** The geometrical shape of the filament is defined. This function is positive outside the filament and negative inside the filament.*/
double geometry(double x, double y)
{
    double Line_up = 0.5 - y;
    double Line_right = L0-0.5-x;
    double retrangle = min(Line_right,Line_up);
    double circle = sq(0.5) - sq(x-(L0-0.5)) - sq(y);
    double shape_final = max(retrangle,circle);
    return -shape_final;
}

int main()
{
    size(L_domain);
    origin(0.,0.);
    init_grid(1<<10);
/** density ratio and viscosity ratio is fixed.

Basilisk solves Navier-Stokes equations with surface tension terms which read:

$$\rho(\partial_{t}\bm{u}+\bm{u}\cdot\nabla\bm{u})=-\nabla p+\nabla\cdot(2\mu\bm{D})+\sigma\kappa\delta_{s}\bm{n}$$

Diemnsionless equation is as below (details of non-dimensionalization can be found in the published paper):

$$\rho^{*}\frac{D\bm{u}}{Dt}= -\nabla p + Oh\mu^{*} \Delta \bm{u} + \delta^{*}_{s} \kappa^{*}\bm{n}$$
*/
    rho1 = 1;
    rho2 = 1./828.4;
  /** By varying the Ohnesorge number $\mu_{1}$, the cases in the paper can be reproduced. */
    mu1 = 0.1;
    mu2 = mu1/55.4;
    f.sigma = 1.;
    run();
}

event init(i=0)
{
  /** The area around the interface is refined at the initial step.*/
    refine(x<L0-0.5 && y<0.6 && y>0.4 && level < 12);
    refine(x>L0-0.5 && sq(y)+sq(x-(L0-0.5))>sq(0.4) && sq(y)+sq(x-(L0-0.5))<sq(0.6) && level < 12);
    fraction(f,geometry(x,y));
}

event adapt(i++)
{
    adapt_wavelet({f,u},(double []){1e-6,1e-3,1e-3},12,8);
}
char file_name[100];
/** Simulation output and we do post-processing after the simulation*/
event outputfile(t+=0.01)
{
    p.nodump = false;
    sprintf(file_name,"dumpfile_%.2f",t);
    dump(file=file_name);
}

event end(t=20.)
{
    
}
