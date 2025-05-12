/**
 # we solve a case of 2D silo with orifice at the bottom by using the poisson equation
 $$
 \nabla^2 c = 0
 $$
 with a flux at the top of the silo
 This system can be solved with the poisson solver. */
#include "grid/multigrid.h"
#include "poisson.h"
#include "utils.h"
/** Concentration source, and some obvious variables
 */
scalar c[],s[];
double dx,dy,W;


/** Boundary conditions: 
 a flux at the top of silo
 an orifice with size W at the bottom of the silo
 */
c[top]    =  neumann(1);
c[bottom] =  fabs(x)<= W/2. ? dirichlet(0): neumann(0);
c[right]  =  neumann(0);
c[left]   =  neumann(0);
/**
 ## Parameters
 The size of the domain `L0`.
 $0<x<L0$ $0<y<L0$
 */
int main() {
    L0 = 1.;
    X0 = -L0/2.;
    Y0 = -L0/2.;
    N =  64;
    dx = L0/N;
    dy = L0/N;
    W = 0.25;
    init_grid (N);
    /**
     ##Solve Laplace
     Initialisation of the solution of laplace equation with `poisson`
     $$
     \nabla^2 c = s
     $$
     with a  source term $s=0$
     */
    foreach()
    c[] = s[] = 0.;
    boundary ({c, s});
    poisson (c,s);

    /**
     Save the results
     */
    
    /**
    numerical solution
     */
    FILE *  fpc = fopen("Fcsilo.txt", "w");
    for (double x =-L0/2 ; x <= L0/2; x += dx){
        for (double y = -L0/2 ; y <= L0/2; y += dy){
            fprintf (fpc, "%g %g %g \n",
                     x, y,   interpolate (c, x, y));}
        fprintf (fpc,"\n");}
    fclose(fpc);
    fprintf(stdout," end\n");
    }



/**
 ##Run
 Then compile and run:
 
 ~~~bash
 qcc -O2 -Wall -o silo silo.c -lm; ./silo
 ~~~
 ##Results
 
 
 ~~~gnuplot  3D plot of analycal solution and numerical solution
 set hidden3d
sp 'Fcsilo.txt'
 ~~~
 
 ~~~gnuplot  with colors
 reset
 set pm3d map
 set size square
 set palette rgbformulae 22,13,-31
 sp 'Fcsilo.txt'
 ~~~

 ##Bibliography
 
 Version 1: 4 december 2015/ 14:58
 */