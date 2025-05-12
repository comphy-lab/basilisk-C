/**
###Introduction

We simulate the flow past a bounded cylinder at Reynolds number 20 and 100. We consider the following geometry: 

* top wall (to the cylinder center): 2.1D
* bottom wall: 2D
* inlet: 2D
* outlet: 20D

Flow is characterized by:

$$Re = \frac{u_{avg}D}{\nu}, u_{avg} = \frac{2}{3}u_{max}$$
*/

#define lmin (1) 
#define lmax ((int)log2(L0))
#define p_n (1)
#define top_bd    (4.2*pr[0] + pc[0].y)
#define bottom_bd (-4.*pr[0] + pc[0].y)
#define right_bd  (40.*pr[0] + pc[0].x)
#define IBM_stencil (1)

#include "grid/quadtree.h"
#include "lag_to_eul.h"
#include "LBM.h"
#include "eul_to_lag.h"

int main()
{
    pr[0] = 10.;

    L0 = (int)(51.2*pr[0]);

    pc[0].x = pr[0] * 4. + 0.5 - L0 / 2.;
    pc[0].y = 0.5;

    origin(-L0/2,-L0/2);

    init_grid(1 << (lmax-0));

    nu = 0.02;

    step_max = 10000;

    bd.left.type   = 3;
    bd.right.type  = 3;
    bd.top.type    = 1;
    bd.bottom.type = 1;

    run();
} 

/**
We use mask() function to create a rectangular domain, somehow restore() is not compatible with mask() so restart part is not included in this case. 
*/

event start(i = 0)
{
    fprintf(stderr, "Max refinement level        = %12d\n", lmax);
    fprintf(stderr, "Total time steps            = %12d\n", step_max);
    fprintf(stderr, "Particle radius             = %12.6f\n", pr[0]);
    fprintf(stderr, "Max vecolity                = %12.6f\n", uref);
    fprintf(stderr, "Viscosity                   = %12.6f\n", nu);
    fprintf(stderr, "Particle Re                 = %12.6f\n", Re_l);

    mask(y < (bottom_bd - 0.5) ? bottom : y > (top_bd + 0.5) ? top : none);
    mask(x > (right_bd + 0.5) ? right : none);

    foreach ()
    {
        u.x[] = 4.0 * uref * ((y - bottom_bd) * (top_bd - bottom_bd) - sq(y - bottom_bd)) / sq(top_bd - bottom_bd);
        u.y[] = 0.;
        rho[] = 1.;
    }

    for (int k = 0; k < q; k++)
    {
        Calculate_feq(k, rho, u, f[k]);
    }
}

/**
We visualize the vorticity contour and generate a movie. 
*/
event image_movie(i += 40)
{
    vorticity(u,omega);
    stats omega_s = statsf(omega);
    p_shape (cc, ff, particle);

    clear();
    view(fov = 5, width = 1800, height = 600, tx = -(pc[0].x + 16.5 * pr[0]) / L0);
    box();
    draw_vof("cc", "ff", filled = -1, fc = {1, 1, 1});
    draw_vof("cc", lw = 3, lc = {0, 0, 0});
    squares("omega", map = cool_warm, max = omega_s.max/4., min = omega_s.min/4.);
    save("vorticity.mp4");
}

/**
We record the drag and lift coefficient at each time step. 
*/
event force_record(i++, i <= step_max)
{
    write_force(particle,(2.*uref/3.));
    static FILE *fp = fopen("force", "w");
    fprintf(fp, "%d %g %g\n", i, particle[0].f_tot.x, particle[0].f_tot.y);
    fflush(fp);
}    

/**
It is important to free the memory of cache which is created for storing the interpolation stencil of each Lagrangian node. 
*/
event final(i = step_max)
{
    free_stencil_cache(particle);
    fprintf(stderr,"\n=========End Simulation=========\n");
}

/**
###Simulation results at Re = 100

~~~gnuplot Time evolution of the drag coefficient
set xlabel 't^*'
set ylabel 'C_D'
set xrange [10:50]
set yrange [3:3.3]
plot '../my_data/Cylinder_bounded_force' every ::32000::160000 u ($1/3200):($2) w l lw 3 lc rgb "blue" title "Drag",\
3.22 w p pt 2 lc rgb "red" t "Schafer et al., 1996, Re=100, C_{Dmax,lower}",\
3.24 w p pt 3 lc rgb "dark-green" t "Schafer et al., 1996, Re=100, C_{Dmax,upper}"
~~~

~~~gnuplot Time evolution of the lift coefficient
set xlabel 't^*'
set ylabel 'C_L'
set xrange [10:50]
set yrange[-1.5:2]
plot '../my_data/Cylinder_bounded_force' every ::32000::160000 u ($1/3200):($3) w l lw 3 lc rgb "blue" title "Lift",\
0.99 w p pt 2 lc rgb "red" t "Schafer et al., 1996, Re=100, C_{Lmax,lower}",\
1.01 w p pt 3 lc rgb "dark-green" t "Schafer et al., 1996, Re=100, C_{Lmax,upper}"
~~~

![Vorticity contour](http://www.basilisk.fr/sandbox/zihaocheng/test_cases/post_processing/Cylinder_bounded_vorticity.mp4)(width = "800")
*/

/**
###Reference 

~~~bib
@InBook{Schaefer1996,
  author    = {Sch{\"a}fer, M. and Turek, S. and Durst, F. and Krause, E. and Rannacher, R.},
  pages     = {547--566},
  publisher = {Vieweg+Teubner Verlag},
  title     = {Benchmark Computations of Laminar Flow Around a Cylinder},
  year      = {1996},
  address   = {Wiesbaden},
  isbn      = {978-3-322-89849-4},
  booktitle = {Flow Simulation with High-Performance Computers \uppercase\expandafter{\romannumeral2}: {DFG} Priority Research Programme Results 1993--1995},
  doi       = {10.1007/978-3-322-89849-4_39},
  file      = {:Flow_Past_Cylinder/dfg_benchmark_results.pdf:PDF},
  url       = {https://doi.org/10.1007/978-3-322-89849-4_39},
}
~~~

*/
