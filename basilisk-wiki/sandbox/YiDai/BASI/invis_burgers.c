/**
# Solving 1D inviscid burgers equaiton 

We solve 1D inviscid burgers equation:
$$\partial_{t}U + U \partial_{x}U = 0$$

using Bell-Collela-Glaz [bcg](http://basilisk.fr/src/bcg.h#advection) advection scheme. See viscid case [here](http://basilisk.fr/sandbox/YiDai/BASI/burgers.c)
*/
/**
~~~gnuplot
reset
set xlabel "x"
set ylabel "T"
p[][] 'out' u ($1):($2) t 'case0' w l,
~~~
*/
#include "grid/cartesian1D.h"
#include "run.h"
#include "timestep.h"
#include "bcg.h"

vector u[];
face vector uf[];


int main(){
    L0 = 1;
    // X0 = -L0/2.;
    N = 1 << 9;
    DT = sq(L0/N)/4;
    run();
}

u.n[left] = dirichlet(0.);
u.n[right] = dirichlet(0.);

event init(t = 0){
    foreach(){
        u.x[] = sin(2*pi*x);
    }
    boundary ((scalar *){u});
    trash ({uf});
    foreach_face()
        uf.x[] = fm.x[]*face_value (u.x, 0);
}
// event printdata (t = 0; t <= 1000 * DT; t += 100 * DT) {
event printdata (t = 0; t <= 0.3; t += 0.05) {

// event printdata (t = 0.3) {
    foreach()
        printf ("%g %g %g\n", x, u.x[], t);
    printf ("\n\n");
}

event advection_term (i++,last){
    double dt = DT;
    dt = dtnext (dt);
    advection ((scalar *){u}, uf, dt);
}