#include "grid/octree.h"
#include "embed.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "tension.h"

#include "utils.h"


/**
This test is trying to simulate a two phase flow with embed boundary and surface tension.

The initial flow is a drop in the air, with a wall above the drop.

In this simulation, we want to simulate a drop of liquid in a gaz, with surface tension, a wall above and an adaptive mesh. However, it seems that there is an incompatibility and the velocity field diverge (1e+100) once we start adapting the mesh.

If we just have a drop in a normal gaz field, the simulation works fine (wich sounds normal, but it is important to check the obvious)
*/

int LEVEL = 6;
double uemax = 0.01;

int main(){
    /**
    We put gentle condition for the 2 fluids*/
    rho1 = 1.;
    rho2 = 1./10.;
    mu1 = 1.;
    mu2 = 1./5.;
    f.sigma = 1.;
    origin (-L0/2., 0., -L0/2.);
    
    /**
    The grid  is not very refined*/
    init_grid(64);

    run();
}

/**
The boundary condition is non slip boundary condition. The condition are set after the main, like in the sphere.c example*/

u.n[embed] = dirichlet(0);
u.t[embed] = dirichlet(0); 
u.r[embed] = dirichlet(0);

event init(t = 0){
    /**
    Initialisation of the drop in the middle of the simulation*/
    fraction(f, -(sq(x)+sq(y)+sq(z)-sq(0.1)));
    
    /**
    We add a wall in the simulation, above the drop. Cs = 0 for y>0.31415, and Cs = 1 below.
    
    If we comment this line, the simulation works just fine.*/
    solid (cs, fs, (0.31415-y));
    dump("init");
}


/**
Event were the bug occurs

If we turn on the adaptation, the velocity field start to diverge (1e+100). We can observe it with the max value of u.x.*/
event logfile(i++){
    stats ux = statsf(u.x);
    stats gx = statsf(g.x);
    fprintf(stderr, "i = %d, t = %g, maxUx = %g, maxGx = %g\n",i,t,ux.max, gx.max);
    adapt_wavelet ({f,u}, (double[]){0.01, uemax, uemax, uemax}, LEVEL, 2);
}


/**
End of the simulation.*/
event end(i = 10){
    dump("end");
}