/**

###Introduction

We simulate the migration of a neutrally bouyant cylinder (in 2D) in a planar Couette flow and Poiseuille flow respectively. In case of planar Couette flow (linear shear flow), top and bottom walls are moving towards opposite directions at an equal velocity $U_w/2$. The fluid bulk Reynolds number and particle Reynolds number are:
$$Re_f = \frac{U_w H}{\nu}, Re_p = \frac{\gamma R^2}{\nu}$$
in which $H$ represents the channel height, $\gamma = U_w/H$ is the shear rate. We consider $Re_f=40$ and $Re_p=0.625$ and $H=8R$ specifically, and the initial position of the cylinder center $y_c = -H/4$:

*/

#define lmin (1) 
#define lmax ((int)log2(L0))
#define p_n (1)
#define IBM_stencil (1)
#define Re (40.)
#define channel_width (4.*pr[0])
#define top_bd    (channel_width+0.5)
#define bottom_bd (-channel_width+0.5)
#define vg (uref/channel_width)
#define Re_p (vg*sq(pr[0])/nu)
#define free_motion 1

#include "grid/quadtree.h"
#include "lag_to_eul.h"
#include "LBM.h"
#include "eul_to_lag.h"

/**
Differing from the other cases, we set $\nu=1/30$ acoording to [Feng et al, 2004](*Feng2004) (which also uses the feedback IB-LBM to solve for the FSI problems)

*/
int main()
{
    pr[0] = 16.;

    L0 = (int)(128*pr[0]);

    pc[0].x = -L0/2.+4.*pr[0]+0.5;
    pc[0].y = -channel_width/2.+0.5;

    origin(-L0/2,-L0/2);

    init_grid(1 << (lmax-0));

    nu = 1./30.;
    uref = Re*nu/channel_width/2./2.;

    step_max = 200000;

    run();
} 

/**
We apply a slip velocity (difference between particle and fluid at the initial state).
*/
event start(i = 0)
{
    fprintf(stderr, "Max refinement level        = %12d\n", lmax);
    fprintf(stderr, "Total time steps            = %12d\n", step_max);
    fprintf(stderr, "Particle radius             = %12.6f\n", pr[0]);
    fprintf(stderr, "Max vecolity                = %12.6f\n", uref);
    fprintf(stderr, "Viscosity                   = %12.6f\n", nu);
    fprintf(stderr, "Flow Re                     = %12.6f\n", Re);
    fprintf(stderr, "Particle Re                 = %12.6f\n", Re_p);
    fprintf(stderr, "Characteristic time scale   = %12.6f\n", channel_width/uref);

    mask (y < (bottom_bd-0.5) ? bottom : y > (top_bd+0.5) ? top : none);

    foreach ()
    {
        if (sq(x - pc[0].x) + sq(y - pc[0].y) > sq(pr[0]))
            u.x[] = -vg * y;
        else
            u.x[] = 0.;
        u.y[] = 0.;
        rho[] = 1.;
    }

    for (int k = 0; k < q; k++)
    {
        Calculate_feq(k, rho, u, f[k]);
    }
}

/**
We output the dimensionless displacement, velocity, force and torque. Be careful of the reference quantity of each variable for non-dimensionalization.  
*/
event force_record(i++, i <= step_max)
{
    write_force(particle, uref);
    static FILE *fp = fopen("force", "w");
    fprintf(fp, "%d %g %g %g %g %g %g %g %g %g %g\n", i, i/(channel_width/uref), 
                (particle[0].center.pos.x+L0/2-0.5)/(channel_width*2.), (particle[0].center.pos.y+channel_width-0.5)/(channel_width*2.), particle[0].angle.z,
                particle[0].center.vel.x/(uref*2.), particle[0].center.vel.y/(uref*2.), particle[0].center.agl_vel.z,
                particle[0].F_tot.x, particle[0].F_tot.y, particle[0].T_tot.z);
    fflush(fp);
}    

event final(i = step_max)
{
    free_stencil_cache(particle);
    fprintf(stderr,"\n=========End Simulation=========\n");
}

/**
###Quantitative results

We plot the cylinder's vertical displacement, and also the time evolution of each velocity component. Our simulation results are compared with those reported in literature as well. 

~~~gnuplot Time evolution of the cylinder center postion
set xlabel 't^*'
set ylabel 'y/H'
set xrange [0:55]
set yrange [0.2:0.7]
plot '../ref_data/Feng2004-fig4' u ($1):($2) w l lw 8 lc rgb "light-green" title 'Feng et al.,2004',\
'../ref_data/Feng1994-fig2' u ($1):($2) w p ps 1.6 pt 6 lw 1.8 lc rgb "orange-red" title 'Feng et al.,1994',\
'../my_data/Neutrally_bouyant_Couette' u ($2):($4) w l lw 2 lc rgb "blue" title 'Present'
~~~

~~~gnuplot Time evolution of the cylinder velocity
set xlabel 't^*'
set ylabel 'Horizontal velocity u_c/U_w'
set y2label 'Vertical velocity v_c/U_w'
set xrange [0:55]
set yrange [0:0.3]
set y2range [-0.005:0.01]
set y2tics -0.005,0.005,0.01
plot '../ref_data/Feng2004-fig5-ux' u ($1):($2) w l lw 8 lc rgb "light-blue" t 'Feng et al,2004:u' axes x1y1,\
'../ref_data/Feng2004-fig5-uy' u ($1):($2) w l lw 8 lc rgb "light-pink" t 'Feng et al,2004:v' axes x1y2,\
'../my_data/Neutrally_bouyant_Couette' u ($2):($5) every::100 w l lw 2 lc rgb "blue" title 'Present: u' axes x1y1,\
'../my_data/Neutrally_bouyant_Couette' u ($2):($6) every::00 w l lw 2 lc rgb "red" title 'Present: v' axes x1y2
~~~

*/


/**
###Reference

~~~bib
@article{Feng1994, title={Direct simulation of initial value problems for the motion of solid bodies in a {N}ewtonian fluid. {P}art 2. {C}ouette and {P}oiseuille flows}, volume={277}, DOI={10.1017/S0022112094002764}, journal={Journal of Fluid Mechanics}, publisher={Cambridge University Press}, author={Feng, J. and Hu, H. H. and Joseph, D. D.}, year={1994}, pages={271–301}}

@Article{Feng2004,
  author   = {Z. G. Feng and E. E. Michaelides},
  journal  = {Journal of Computational Physics},
  title    = {The immersed boundary-lattice {B}oltzmann method for solving fluid-particles interaction problems},
  year     = {2004},
  issn     = {0021-9991},
  number   = {2},
  pages    = {602-628},
  volume   = {195},
  doi      = {https://doi.org/10.1016/j.jcp.2003.10.013},
  file     = {:LBM/1-s2.0-S0021999103005758-main.pdf:PDF},
  url      = {https://www.sciencedirect.com/science/article/pii/S0021999103005758},
}
~~~

*/

