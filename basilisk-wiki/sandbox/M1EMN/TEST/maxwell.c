/**
# Wave equation with dispersion due to a moving charge
 
 Explicit resolution of Maxwell in 1D
 $$\frac{\partial B}{\partial t}= - \frac{\partial E}{\partial x}$$
 $$\frac{\partial E}{\partial t}= - \frac{\partial B}{\partial x}
 +  \rho u_p$$
 + plus Newton Law
 
 
*/
//#include "grid/cartesian1D.h"
#include "grid/multigrid1D.h"
#include "run.h"

scalar E[],B[];
double xp,up,Mp,En,f,tmax;
scalar proba[];
double dt, le;
FILE * out;
/**
 Outflow boundary conditions. */
B[left] = 0;
E[left] = neumann(0);
proba[left] = neumann(0);
B[right] = 0 ;
E[right] = neumann(0);
proba[left] = neumann(0);

/**
 density of charge */
double rhoe(double x1, double size)
{
    double r=0;
    if (fabs( x1) < size/2.)
        r= ( 1.+cos( 2.*pi*( x1 ) / size ) ) /2.;
    return r;
}
/**
 Parameters, siez of domain, origin, number of points, time step
 */
int main() {
    L0 = 6.;
    X0 = -3;
    N = 128*2*2;
    DT = .0001;
    tmax = 10000;
/**
 initial position velocity mass and size
 */
 
    xp=0;
    up=.2;
    Mp=2;
    le=.5;
    
    FILE *g = fopen("F.IN", "w");
  //  fprintf(g,"xp=%lf          \n",xp);
  //  fprintf(g,"up=%lf          \n",up);
    fprintf(g,"Mp=%lf          \n",Mp);
    fprintf(g,"le=%lf          \n",le);
    fclose(g);
/**
     some Results are written in a file   */
        char name[80];
        sprintf (name, "out");
        out = fopen (name, "w");
        run();
}
/**
 The initial field is zero */
event init (t = 0) {
    foreach()
    E[] = 0. + 1*exp(-x*x/.10);
    foreach()
    B[] = 0*exp(-x*x);
    foreach()
    proba[] = 0;
    boundary ({E,B,proba});
}
/**
 Output fields. */
event printdata (t += 1; t <= tmax) {
    double pr=1e-99;
    foreach()
    pr +=proba[]*proba[]*Delta;
    
    fprintf (stdout,"p[][-1:2] '-' u 1:2 t'E' ,'' t'B' w l,'' w l, '' w l\n");
    foreach()
    fprintf (stdout, "%g %g  \n", x, E[] );
    fprintf(stdout,"e \n");
    foreach()
    fprintf (stdout, "%g %g  \n", x, B[] );
    fprintf(stdout,"e \n");
    foreach()
    fprintf (stdout, "%g %g  \n", x, rhoe(x-xp,le) );
    fprintf(stdout,"e \n");
    foreach()
    fprintf (stdout, "%g %g  \n", x,  proba[]/sqrt(pr));
    fprintf(stdout,"e \n");
}
/** 
 Save Energy
*/
event printE (t += .1 ) {
    En=0;
    double pr=1e-99;
    foreach()
    pr +=proba[]*proba[]*Delta;
    
    foreach()
    En += (sq(E[])/2+sq(B[])/2) * Delta;
    
    //    En += Mp * sq(up)/2;
    fprintf (out, "%g %g %g %g %g \n", t, xp, up , En , Mp * sq(up)/2);
    
    FILE *g = fopen("proba.out", "w");
    foreach()
    fprintf (g, "%g %g  \n", x,  proba[]/sqrt(pr));
    fclose(g);
}
/**
 change the mass on the fly for the fun
*/
event readata (t += 1){
    FILE *g = fopen("F.IN", "r");
    fscanf(g,"xp=%lf          \n",&xp);
    fscanf(g,"up=%lf          \n",&up);
    fscanf(g,"Mp=%lf          \n",&Mp);
    fscanf(g,"le=%lf          \n",&le);
    fclose(g);
}

/**
 Integration. */
event integration (i++) {
    double dt = DT;
  
    dt = dtnext (t, dt);
    /**
     Variation of B
     $$\frac{\partial B}{\partial t}= - \frac{\partial E}{\partial x}$$ */
    foreach()
    B[] += - dt* (E[0] - E[-1])/Delta;
     boundary ({B});
    /**
     Variation of E
     $$\frac{\partial E}{\partial t}= - \frac{\partial B }{\partial x} + \rho u_p$$ */
    foreach()
    E[] += - dt* (B[1] - B[0])/Delta - dt * rhoe(x-xp,le)*up;
    boundary ({E});
    /** 
     compuation of the force
     $f = \int \rho E dx$
    */
    f=0;
    foreach()
    f += rhoe(x-xp,le)*E[]*Delta;
 
    
    /** dispalcmeent of the mass:
     $M_p\frac{d u_p}{dt} = f $ change positio
     $\frac{d x_p}{dt} = u_p$
     */
    up += f * dt / Mp;
    xp += up * dt;
    /**
     refelxion */
    
    if(xp>=(L0/2-le/2)){
        xp = L0 -  le -xp;
        up = -up;
    }
    if(xp<=-L0/2+le/2){
        xp = - L0 +  le - xp;
        up = -up;
    }
 /**
  update probality */
    foreach()
    proba[] +=   rhoe(x-xp,le)*Delta;
    boundary ({proba});

    
}
/**
#Run
 
~~~bash
 qcc -g -O3 -o maxwell maxwell.c -lm
 ./maxwell
~~~
 
#Results
 
 cumalative proba
 
~~~gnuplot proba,
 set xlabel "x"
 p'proba.out'
~~~
 

 Energy
 
 
~~~gnuplot energy as function of time,
 set xlabel "t"
  p'out' u 1:4 t'ElecMagn' w lp,'' u 1:($5+$4) t'total' w lp
~~~
 
 
 
 
 Paris  03/16
 
*/
