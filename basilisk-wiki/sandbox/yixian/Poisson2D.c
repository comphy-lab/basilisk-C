/**

# The poisson equation
 The poisson equation is a diffusion Fick equation of a specie $c$
 $$
 \nabla^2 c = s
 $$
 we solve it in 2D.
 we compare the solution with the analytical solution of equation
 $$
 \nabla^2 c = 1
 $$
 with boundary condition dirichlet 1 at every boundary
 the solution is 
 $$
ce(x,y) = \sum\limits_{l=1}^{\infty}\sum\limits_{m=1}^{\infty}U_{lm}sin((2l-1)\pi x)sin((2m-1)\pi y)
 $$
 with
 $$
 U_{lm}=\frac{-16}{((2l-1)^2+(2m-1)^2)\pi ^4 (2l-1)(2m-1)}
 $$
 This system can be solved with the poisson solver. */
#include "grid/multigrid.h"
#include "poisson.h"
#include "utils.h"
#define pi 3.1415926
/** Concentration source  and flux, and some obvious variables
 */
scalar c[],s[],flux[];
double dx,dy;
/**
 the exact solution without mask
*/
double cexact(double x, double y)
{ double ce;
    ce =0.;
    //ce = (bi*y+1)/(bi*L0+1);
    for (double l = 1 ; l <= 200; l += 1){
        for (double m= 1 ; m<= 200; m += 1){
            ce += -16*sin((2*l-1)*pi*x)*sin((2*m-1)*pi*y)/(((2*l-1)*(2*l-1)+(2*m-1)*(2*m-1))*pi*pi*pi*pi*(2*l-1)*(2*m-1));
        }
    }
    return ce;
}

/** Boundary conditions: dirichlet 1
 */
c[top]    =  dirichlet(0);
c[bottom] =  dirichlet(0);
c[right]  =  dirichlet(0);
c[left]   =  dirichlet(0);
/**

## Parameters

The size of the domain `L0`.
 $0<x<L0$ $0<y<L0$
 */
int main() {
    L0 = 1.;
    X0 = 0.;
    Y0 = 0.;
    N =  32;
    dx = L0/N;
    dy = L0/N;
    init_grid (N);
    /**
    
 ##Solve Laplace
     Initialisation of the solution of laplace equation with `poisson`
     $$
     \nabla^2 c = s
     $$
     with a  source term $s=1$
     */
    foreach()
    c[] = 0.;
    foreach()
    s[]= 1.;
    boundary ({c, s});
    poisson (c,s);
    /**
     Computation of the flux
     */
    foreach()
    flux[] =  ( c[0,0] - c[0,-1] )/Delta;
    boundary ({flux});
    /**
     Save the results
     
*/
    /**
    analytical solution
     */
    FILE *  fpx = fopen("Fce.txt", "w");
    for (double x = 0 ; x <=1; x += dx){
        for (double y = 0; y <= 1; y += dy){
            fprintf (fpx, "%g %g %g \n",
                     x, y,  cexact (x,y));}
        fprintf (fpx,"\n");}
    fclose(fpx);
    /**
    numerical solution
     */
    FILE *  fpc = fopen("Fcs.txt", "w");
    for (double x = 0 ; x <1-dx; x += dx){
        for (double y = 0; y < 1-dy; y += dy){
            fprintf (fpc, "%g %g %g \n",
                     x, y,   interpolate (c, x, y));}
        fprintf (fpc,"\n");}
    fclose(fpc);
    /**
     Error
     */
    FILE *  fpe = fopen("Fcerror.txt", "w");
    static double abscemax=0., cemax=0., abscsmax=0., csmax=0.;
    for (double x = 0 ; x <=1; x += dx){
        for (double y = 0; y <= 1; y += dy){
            cemax = fabs(cexact (x,y));
            if(cemax > abscemax) abscemax = cemax;
            csmax = fabs( interpolate (c, x, y));
            if(csmax > abscsmax) abscsmax = csmax;
            }}
    fprintf (fpe, " %g \n", fabs((abscemax-abscsmax)/abscemax));
    fclose(fpe);
    fprintf(stdout,"Err = ");
    fprintf(stdout," %g\n",  fabs((abscemax-abscsmax)/abscemax));
    fprintf(stdout," end\n");
    }



/**

##Run
 Then compile and run:
 
~~~bash
 qcc -O2 -Wall -o Poisson2D Poisson2D.c -lm; ./Poisson2D
~~~
 
##Results
 
 
~~~gnuplot  3D plot of analycal solution and numerical solution
set hidden3d
sp[][][:] 'Fcs.txt' , 'Fce.txt' 
~~~
 

##Bibliography

* http://www.ufrmeca.univ-lyon1.fr/~buffat/COURS/COURSDF_HTML/node32.html
 
 Version 1: december 2015
 */