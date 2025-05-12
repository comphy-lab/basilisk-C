/**
###Introduction

We simulate a 2D cylinder oscillating vertically in a fluid at rest. Motion of the cylinder is precribed by:

$$U_y(t) = -A\cos(2\pi ft)$$

The oscillating amplitude $A$ and frequency $f$ are related to the  Keulegan-Carpenter number:

$$KC = \frac{U_{max}}{Df} = \frac{2\pi A}{D}$$

The flow is characterized by:
$$Re = \frac{U_{max}D}{\nu}$$

We consider two flow configurations:

* Re = 100 & KC = 5: L0 = 12.8D
* Re = 200 & KC = 10: L0 = 25.6D
*/

#define lmin (1) 
#define lmax ((int)log2(L0))
#define p_n (1)
#define IBM_stencil (1)
#define KC 10
#define freq (uref / (KC * 2. * pr[0]))
#define p_displacement_y(i) (- KC * pr[0] * sin(2.*M_PI*freq*i) / M_PI)
#define p_acceleration_y(i) (2. * M_PI * freq * uref * sin(2.*M_PI*freq*i))

#include "grid/quadtree.h"
#include "lag_to_eul.h"
#include "LBM.h"
#include "eul_to_lag.h"

int main()
{
    pr[0] = 20.;

    L0 = (int)(51.2*pr[0]);

    pc[0].x = 0.;
    pc[0].y = 0.;

    origin(-L0/2,-L0/2);

    init_grid(1 << (lmax-0));

    nu = 0.02;

    step_max = 100000;

    bd.left.type   = 2;
    bd.right.type  = 2;
    bd.top.type    = 2;
    bd.bottom.type = 2;

    run();
} 

/**
We initialize the quadtree grid based on the solid object's geometry. 
*/
event start(i = 0)
{
    if (!restore (file = "restart"))
    {
        fprintf(stderr, "Max refinement level        = %12d\n", lmax);
        fprintf(stderr, "Total time steps            = %12d\n", step_max);
        fprintf(stderr, "Particle radius             = %12.6f\n", pr[0]);
        fprintf(stderr, "Max vecolity                = %12.6f\n", uref);
        fprintf(stderr, "Viscosity                   = %12.6f\n", nu);
        fprintf(stderr, "Particle Re                 = %12.6f\n", Re_l);
        fprintf(stderr, "Period                      = %12d\n", (int)(1/freq));

        astats ss;
        int ic = 0;
        do
        {
            ic++;
            p_shape(cc, ff, particle);
            ss = adapt_wavelet({cc}, (double[]){1.e-30},
                               maxlevel = (lmax), minlevel = (lmin));
        } while ((ss.nf || ss.nc) && ic < 100);
        Rearrange_indices(particle);

        foreach()
        {
            u.x[] = 0.;
            if (sq(x - pc[0].x) + sq(y - pc[0].y) > sq(pr[0]))
                u.y[] = 0.;
            else
                u.y[] = -uref;
            rho[] = 1.;
        }

        for (int k = 0; k < q; k++)
        {
            Calculate_feq(k, rho, u, f[k]);
        }
    }
/**
We restart the simulation for the post-processing purpose. 
*/
    else
    {
        restore_particle(particle,1000);
        Rearrange_indices(particle);
        pc[0].y = p_displacement_y(1000);
        pa[0].y = p_acceleration_y(1000);
    }
}

/**
We dynamically refine and coarsen grid at each time step. 
*/
event adapt(i++)
{
    p_shape(cc, ff, particle);
    adapt_wavelet({cc, u}, (double[]){1.e-3, (1.e-3)*uref, (1.e-3)*uref}, maxlevel = (lmax), minlevel = (lmin));
}

/**
We record the total number of grid cells and force coefficients of each time step. 
*/
event force_record(i++, i <= step_max)
{
    write_force(particle,uref);
    static FILE *fp = fopen("force", "w");
    fprintf(fp, "%g %ld %g %g\n", i*freq, grid->tn, particle[0].f_tot.x, particle[0].f_tot.y);
    fflush(fp);
}    

/**
We update the cylinder's position and acceleration. 
*/
event update_particle(i++, i <= step_max)
{
    pc[0].y = p_displacement_y(i);
    pa[0].y = p_acceleration_y(i);
}

/**
Generate movie. 
*/
event image_movie(i += 50)
{
    vorticity(u,omega);
    p_shape(cc, ff, particle);

    clear();
    view(fov = 6, width = 800, height = 1000);
    box();
    cells();
    draw_vof("cc", "ff", filled = -1, fc = {1, 1, 1});
    draw_vof("cc", lw = 3, lc = {0, 0, 0});
    squares("omega", map = cool_warm, linear = 1);
    save("vorticity.mp4");
}

/**
We write the flow field and cylinder's information to the dump files. 
*/
event my_dump(i += 1000)
{
    char name[80];
    sprintf(name, "dump_%d", i);
    dump(name);
    dump_particle(particle, i);
}

/**
Free the memory of caches. 
*/
event final(i = step_max)
{
    free_stencil_cache(particle);
    dump("final");
    dump_particle(particle, i);
    fprintf(stderr,"\n=========End Simulation=========\n");
}

/**
###Quantitative results at Re = 200 & KC =10 (with R = 40)

~~~gnuplot Time evolution of the force coefficients
set xlabel 't/T'
set ylabel 'C_D, C_L'
set yrange[-4:4]
set xrange[0:40]
plot '../my_data/Oscillating_Cylinder_force' u ($1/16000):($3) w l lw 3 lc rgb "blue" title "Drag",\
'../my_data/Oscillating_Cylinder_force' u ($1/16000):($2) w l lw 3 lc rgb "red" title "Lift"
~~~

~~~gnuplot Drag coefficient over one cycle
set xlabel 't/T'
set ylabel 'C_D'
set xrange [0:1]
set yrange [-4:4]
set grid
plot '../my_data/Oscillating_Cylinder_force' every ::303850::319850 u ($1-303850)/16000:3 w l lw 3 lc rgb "blue" title "Drag",\
'../ref_data/Dutsch1998-fig15-cycle14.csv' u 1:($2/2) w p ps 1.5 pt 6 lw 2 lc rgb "red" t 'fig. 15, Dutsch et al., 1998, Re=200 and KC=10'
~~~

~~~gnuplot Lift coefficient over one cycle
set xlabel 't/T'
set ylabel 'C_L'
set xrange [0:1]
set yrange [-4:4]
set grid
plot '../my_data/Oscillating_Cylinder_force' every ::311900::327900 u ($1-311900)/16000:2 w l lw 3 lc rgb "blue" title "Lift",\
'../ref_data/Dutsch1998-fig16-cycle14.csv' u 1:($2/2) w p ps 1.5 pt 6 lw 2 lc rgb "red" t 'fig. 16, Dutsch et al., 1998, Re=200 and KC=10'
~~~

*/

/**
###Qualitative results at Re = 200 & KC =10 (with R = 20)

![Vorticity contour](http://www.basilisk.fr/sandbox/zihaocheng/test_cases/post_processing/Oscillating_cylinder_vorticity.mp4)(width = "600")
*/

/**
###Reference 

~~~bib
@Article{Dutsch1998,
  author    = {D{\"u}tsch, H. and Durst, F. and Becker, S. and Lienhart, H.},
  journal   = {Journal of Fluid Mechanics},
  title     = {Low-{R}eynolds-number flow around an oscillating circular cylinder at low {K}eulegan-{C}arpenter numbers},
  year      = {1998},
  pages     = {249-271},
  volume    = {360},
  doi       = {10.1017/S002211209800860X},
  file      = {:Oscillating_Cylinder/s002211209800860x.pdf:PDF},
  publisher = {Cambridge University Press},
  url       = {https://www.cambridge.org/core/journals/journal-of-fluid-mechanics/article/abs/lowreynoldsnumber-flow-around-an-oscillating-circular-cylinder-at-low-keulegancarpenter-numbers/4FA2C65919CD7E90E9C9FA3D2076AF47},
}
~~~

*/