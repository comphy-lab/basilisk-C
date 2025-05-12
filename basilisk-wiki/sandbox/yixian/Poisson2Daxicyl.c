/**We solve an  axisymmetric Poisson’s Problem from Wang’s thesis.

The longitudinal coordinate ($z$-axis) is *x*, and the radial coordinate
($\rho$- or $r$-axis) is *y*. Note that *y* (and so *Y0*) cannot be negative.
 
$$
\frac {\partial_y (y \epsilon \partial_y \phi)}{y} + \partial_x(\epsilon \partial_x \phi) =  - \frac{sin(y)}{e^{x} y}
$$
 
with Boundary Condition :
$$ \phi(0,y) =cos(y) $$
$$ \phi(1,y) =e^{-1}cos(y) $$
$$ \phi(x,1) =e^{-x}cos(1) $$
The exact solution of the above axisymmetric problem is given by
$$ \phi = e^{-x} cos y $$
 */


#include "grid/multigrid.h"
#include "axi.h"
#include "run.h"
#include "poisson.h"
#define LEVEL 6
mgstats mg;
scalar phi[];
face vector epsilon[];
face vector a[], alpha[];
face vector gradphi[];
double dx,dy;
/**
The source term  */
#define RHOE(x,y) sin(y)/(y*exp(x))

/**
The exact solution  */
#define PHIE(x,y) exp(-x)*cos(y)

/**
The boundary condition  */
phi[right] = dirichlet(exp(-1)*cos(y));
phi[left] = dirichlet(cos(y));
phi[top] = dirichlet(exp(-x)*cos(1.));


/**
Main Program  */
int main()
{
L0 = 1.;
N = 1 << LEVEL;
dx = 1./N;
dy = 1./N;
run();
}

event init (i = 0) {
scalar rhoe[], rhs[];
foreach() {
    rhoe[] = RHOE(x,y);
    rhs[] = -rhoe[]*cm[];
    phi[] = 0.;
}
foreach_face()
    epsilon.x[] = fm.x[];
boundary ((scalar *){rhoe, phi, epsilon});
mg = poisson (phi, rhs, epsilon);
}
/**
 save field
 */
event interface (t = 0 ) {
char s[80];
sprintf (s, "field-%g.txt", t);
FILE * fp = fopen (s, "w");
output_field ({phi}, fp, linear = true);
fclose (fp);
    
/**
analytical solution */
FILE *  fpx = fopen("Fce.txt", "w");
for (double x = 0 ; x <=1; x += dx){
    for (double y = 0; y <= 1; y += dy){
            fprintf (fpx, "%g %g %g \n", x, y,  PHIE (x,y));}
            fprintf (fpx,"\n");}
fclose(fpx);
    
/**
 Error */
FILE *  fpe = fopen("Fcerror.txt", "w");
static double cerror=0., ce=0., cs=0.;
for (double x = 0 ; x <=1; x += dx){
    for (double y = 0; y <= 1; y += dy){
            ce = fabs(PHIE (x,y));
            cs = fabs( interpolate (phi, x, y));
            cerror+=(ce-cs)*(ce-cs);
        }}
fprintf (fpe, " %g \n", fabs(cerror));
fclose(fpe);
    
}

/**
# Run

to run

~~~bash
qcc -g -O2 -Wall -o Poisson Poisson2Daxicyl.c -lm
./Poisson > out
~~~
Plot of the numerical and analytical solution
~~~gnuplot 
sp 'field-0.txt', 'Fce.txt'
~~~

# Bibliography

Kaichun Wang, BEM simulation for glass parisons, Ph.D. thesis, Eindhoven Technical
University, the Netherlands, 2002.
*/