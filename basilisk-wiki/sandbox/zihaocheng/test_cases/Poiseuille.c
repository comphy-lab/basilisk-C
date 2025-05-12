/**
# Introduction
This is a simple test on the 2D planar Poiseuille flow. p_n stands for the number of solid particles. In this case we only consider single-phase flow:
*/

#define lmin (1) 
#define lmax ((int)log2(L0))
#define p_n (0)

/**
Source files of my fluid solver (lattice Boltzmann method) and fluid-particle interaction (immersed boundary method) are on the way :) We define the domain size L0, kinetic viscosity nu, total time steps and boundary types:
left: inlet;
right outlet;
top and bottom: no-slip
*/

#include "grid/quadtree.h"
#include "LBM.h"

int main()
{
    L0 = 64;

    periodic(right);

    origin(-L0/2,-L0/2);

    init_grid(1 << (lmax-0));

    nu = 0.02;

    step_max = 20000;

    bd.left.type   = 3;
    bd.right.type  = 3;
    bd.top.type    = 1;
    bd.bottom.type = 1;

    run();
} 

/**
# Initialization
We write some key parameters in the log file. We also initialize the velocity, density field, and distribution functions. Once we restore the flow filed from the previous simulation, we generate some plots: 
*/

event start(i = 0)
{
    if (!restore (file = "restart"))
    {
        fprintf(stderr, "Max refinement level        = %12d\n", lmax);
        fprintf(stderr, "Total time steps            = %12d\n", step_max);
        fprintf(stderr, "Max vecolity                = %12.6f\n", uref);
        fprintf(stderr, "Viscosity                   = %12.6f\n", nu);
        fprintf(stderr, "Channel length              = %12.6f\n", L0);

        foreach()
        {
            u.x[] = 0.;
            u.y[] = 0.;
            rho[] = 1.;
        }

        for (int k = 0; k < q; k++)
        {
            Calculate_feq(k, rho, u, f[k]);
        }
    }
    else
    {
        clear();
        box();
        squares("u.x", map = cool_warm, max = uref);
        save("xx-flow.jpg");
    }
}

/**
# L2-error
We monitor the L2-error of ux: 
*/

event output_error(i += 100)
{
    vector ua[];
    double error = 0.;
    int cell_tn = 0;
    foreach(reduction(+:error) reduction(+:cell_tn))
    {
        if (x != left_bd && x != right_bd && y != top_bd && y != bottom_bd)
        {
            ua.x[] = 4.0 * uref * ((y - bottom_bd) * (top_bd - bottom_bd) - sq(y - bottom_bd)) / sq(top_bd - bottom_bd);
            ua.y[] = 0.;
            error += sq(u.x[] - ua.x[]);
            cell_tn++;
        }
    }
    error = sqrt(error / cell_tn);
    fprintf(stderr,"\n=========Error at %10d step is: %15.12f=========\n", i, error);
  
}

event final(i = step_max)
{
    dump("restart");
    fprintf(stderr,"\n=========End Simulation=========\n");
}

/**

![Horizontal velocity contour](http://www.basilisk.fr/sandbox/zihaocheng/test_cases/post_processing/xx-flow.png)(width="800" height="800")

~~~gnuplot L2 error of the horizontal velocity velocity
set xlabel 'L0'
set ylabel 'L_2 norm'
set logscale
set xrange [32:512]
set yrange [1.e-8:1.e-4]
set format y "10^{%L}" [_$_]
set xtics 32,2,512
set grid
plot '../ref_data/2_nd' u 1:2 w l lw 1 lc rgb "blue" title "Second-order convergence",\
'../my_data/Poiseuille_ux' u 1:2 w p ps 1.5 pt 6 lw 3 lc rgb "red" title "Present"
~~~

*/


