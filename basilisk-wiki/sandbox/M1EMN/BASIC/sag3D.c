/**
# The SAG equation
 The SAG equation is a diffusion Fick equation of a specie $c$:

 $$
   \nabla^2 c = 0
 $$
 with boundary condition of cristal growth on the substrate
 $\frac{\partial c}{\partial y} = bi\; c$ on $y=0$ ($bi$ is a kind of Biot number)
 and no growth $\frac {\partial c}{\partial y}=0$ on $y=0$ on the mask
 and far from the wall, there is a constant arrival of species $c(x,top)=1$,
 right and left are periodic conditions (or infinite domain),
 we solve it in 3D.
 
 This system can be solved with the poisson solver. */
//#include "grid/multigrid.h"
#include "grid/multigrid3D.h"
#include "poisson.h"
#include "utils.h"

/** Concentration source  and flux, and some obvious variables
 */
scalar c[],s[],flux[];
double bi,dx;
/** 
 the exact solution without mask
 */
double cexact(double y)
{ double ce;
    ce = (bi*y+1)/(bi*L0+1);
    return ce;
}
/**
 the  mask function
 */
double masque(double x, double z)
{
    return ((fabs(x)<1 && z<0)? 1:0);
}
/** Boundary conditions: a given concentration at the top,
 no reaction on the mask (no flux: $\partial c(t)/\partial y=0$), a complete reaction on the cristal
 $\partial c(t)/\partial y = bi\;  c(t)$
 this mixed condition is written
 `(c[0,0]-c[bottom])/Delta = bi (c[0,0] + c[bottom])/2`
 */
c[top]    =  dirichlet(1);
c[bottom] = ( masque(x,z))?   neumann(0) :c[]*(2.-bi*Delta)/(2.+bi*Delta)  ;
c[right]  =  neumann(0);
c[left]   =  neumann(0);
c[back]   =  neumann(0);
c[front]  =  neumann(0);


flux[bottom] = neumann(0);

/**
## Parameters
 The size of the domain `L0`. it is a cube of reacting surface on the plane
 $-L_0/2<x<L0/2$ $-L_0/2<z<L0/2$, the gaz is $0<y<L0$
*/
int main() {
    L0 = 10.;
    X0 = -L0/2;
    Y0 = 0;
    Z0 = -L0/2;
    N =  256/2;
    dx = L0/N;
    init_grid (N);
/**
   parameter of flux
*/
    bi = 1.2;
/**
##Solve Laplace
  Initialisation of the solution of laplace equation with `poisson`
 $$
 \nabla^2 c = s
 $$
 with a zero source term $s=0$ 
*/
    foreach()
      c[] = s[] = 0.;
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
    double c00=1./(1+bi*L0);
    FILE *  fpx = fopen("Fcxz.txt", "w");
    for (double x = -L0/2 ; x < L0/2; x += dx){
        for (double z = -L0/2 ; z < L0/2; z += dx){
        fprintf (fpx, "%g %g %g %g %g \n",
                 x, z, 1-masque(x,z),  interpolate (flux, x, 0, z)/c00/bi, interpolate (c, x, 0, z)/c00);}
        fprintf (fpx,"\n");}
   fclose(fpx);
    
    FILE *  fpc = fopen("Fcut.txt", "w");
    for (double x = -L0/2 ; x < L0/2; x += dx){
        for (double z = -L0/2 ; z < L0/2; z += dx){
            fprintf (fpc, "%g %g %g %g  \n",
                     x,  interpolate (flux, x, 0, 0)/c00/bi, z, interpolate (flux, 0, 0, z)/c00/bi);}
        fprintf (fpc,"\n");}
    fclose(fpc);
    fprintf(stdout," end\n");
}
/**
## Run
 Then compile and run:
 
~~~bash
 rm sag3D; qcc  -g -O2 -DTRASH=1 -Wall  sag3D.c -o sag3D ; ./sag3D
~~~

or better 

~~~bash
 make sag3D.tst;make sag3D/plots    
 make sag3D.c.html ; open sag3D.c.html 
~~~ 
 
## Results

plots of flux

~~~gnuplot  3D plot of flux
 set hidden3d
 sp 'Fcxz.txt' u 1:2:4 w l
~~~

 concentration mask ...

~~~gnuplot  multiplot
 set multiplot;set hidden3d;
 set size 0.6,0.6
 set origin .0,.5 ;;set view 30,30
 sp [-5:5][:][0:]'Fcxz.txt'  u 1:2:4 t'flux' w l
 set origin .5,.5;set view 90,0
 sp [:][][0:]'Fcxz.txt'   u 1:2:4t'flux' w l
 set origin .5,.0
 set view 0,0;set contour;set nosurface;
 sp [:][:][0:]'Fcxz.txt' u 1:2:3 t'mask'w l
 set origin .0,.0 ;
 set view 90,270
 set size 0.4,0.4
 p'Fcut.txt'  u 1:2 t'cut x' w l,''u 3:4 t'cut z'w l
 unset multiplot
~~~


~~~gnuplot  with colors
 reset
 set pm3d; set palette rgbformulae 22,13,-31;unset surface;
 set ticslevel 0;
 unset border;
 unset xtics;
 unset ytics;
 unset ztics;
 unset colorbox;
 #set xrange[-3:3];set yrange[-3:3];
 set view 0,0
sp'Fcxz.txt' u 1:2:(($4)>1 ? $4:  1) not
~~~
 
## Bibliography
 
* N. Dupuis, J. Décobert, P.-Y. Lagrée, N. Lagay, D. Carpentier, F. Alexandre (2008):
 "Demonstration of planar thick InP layers by selective MOVPE".
 Journal of Crystal Growth
 issue 23, 15 November 2008, Pages 4795-4798
 
* N. Dupuis, J. Décobert, P.-Y. Lagrée , N. Lagay, C. Cuisin, F. Poingt, C. Kazmierski, A. Ramdane, A. Ougazzaden (2008):
 "Mask pattern interference in AlGaInAs MOVPE Selective Area Growth : experimental and modeling analysis".
 Journal of Applied Physics 103, 113113 (2008)
 
 Version 1: april 2015

ready for new site 09/05/19
*/