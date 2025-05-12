/**
# The canopy in nighttime inversion
This is the testing case for [canopy scheme](canopy2d.h). For illustration purpose we present it in 2D. The dissipation term needs to be included in 3D case. 

*/


#include "navier-stokes/centered.h"
#include "tracer.h"
#include "diffusion.h"
#include "view.h"
#include "../Antoonvh/profile5b.h"

#define CP 1005.      // C_p for air
#define gCONST 9.81   // Gravitational constant
#define TREF 273.     // Kelvin
#define INVERSION 0.4 // Kelvin per meter
#define karman 0.4    // von Karman constant

#define roughY0u 0.1               // roughness wind length
#define WIND(s) 0.5                // log Wind profile
#define QFLX (-0.001 * (t / 5000)) // 0 (0.001 = 30wm2)
#define BSURF (1.5 * b[] - 0.5 * b[0, 1])
#define GFLX (-Lambda * (BSURF - bd))
double Lambda = 0.005, bd = 0.; // Grass coupling
#define STRAT(s) (gCONST / TREF * INVERSION * (s < 100 ? s : 100))

scalar b[];
scalar *tracers = {b};

int maxlevel, minlevel;
double eps;
double Re = 2000;

#define CANOPY 1
#if CANOPY
#include "canopy2d.h"
#endif

face vector av[];
face vector muc[];
double TEND = 300.;

int main()
{
    minlevel = 4;
    maxlevel = 8;
    a = av;
    L0 = 400;
    X0 = Y0 = Z0 = 0;
    N = 64;
    mu = muc;
    foreach_dimension()
    {
        u.x.refine = refine_linear; // Momentum conserved     }
    }
    b.gradient = minmod2; // Flux limiter
    run();
}

event properties(i++)
{
    foreach_face()
        muc.x[] = fm.x[] / Re;
    boundary((scalar *){muc});
}

event init(t = 0)
{
    b[bottom] = BSURF;
    b[top] = dirichlet(STRAT(y));
    u.n[bottom] = dirichlet(0.);
    u.t[bottom] = dirichlet(0.);
    u.n[top] = dirichlet(0.);
    u.t[top] = dirichlet(WIND(y));

    periodic(left);
    eps = .1;
    while (adapt_wavelet((scalar *){u, b}, (double[]){eps, eps, eps, 0.35 * 9.81 / 273}, maxlevel, minlevel).nf)
    {
        foreach ()
        {
            b[] = STRAT(y);
            u.x[] = WIND(y);
            TV[] = 275;
        }
    }
}

/* Gravity forcing */
event acceleration(i++)
{
    foreach_face(y)
    {
        av.y[] = (b[] + b[0, -1]) / 2.;
    }
#if CANOPY
    foreach_face()
    {
        if (y <= Zh)
            av.x[] = av.x[] - Cd * PAD(y) * fabs(u.x[]) * u.x[];
    }
#endif
}

mgstats mgb;
/* Diffusion */
event tracer_diffusion(i++)
{
    scalar r[];
    foreach ()
    {
        r[] = 0;
        if (y < Delta)
            r[] = (QFLX + GFLX) / sq(Delta); // div needed as normalization
#if CANOPY
        if (y <= Zh && y >= Delta)
        {
            r[] = r[] + H[] / Cp_a * rho_a;
        }
#endif
    }
    double flx = 0, bt = 0;
    double fctr = CP * TREF / gCONST;
    foreach_boundary(bottom, reduction(+
                                       : flx) reduction(+
                                                        : bt))
    {
        flx = flx + (QFLX + GFLX) * sq(Delta);
        bt = bt + BSURF * sq(Delta);
    }
    bt = bt / sq(L0);
    flx = flx / sq(L0);
    fprintf(stderr, "soil=%g %g %g %d\n", t, fctr * flx, fctr * bt / CP, i);
    mgb = diffusion(b, dt, mu, r = r);
}

event adapt(i++)
{
    adapt_wavelet((scalar *){u, b}, (double[]){eps, eps, eps, .35 * 9.81 / 273}, maxlevel, minlevel);
}

event profiler(t += 100.)
{
    char fname[99];
    sprintf(fname, "proft=%dT", (int)(t / 100.));
    scalar Ta[];
    foreach ()
    {
        Ta[] = b[] * T_ref / Gconst + T_ref;
    }
    profile((scalar *){u, Ta, TV}, fname);
}

event movie(t += 1.)
{
    output_ppm(b, file = "b.mp4", n = 400, linear = true, min = 0, max = 1);
}

event end(t = TEND)
{
}

/**
## Results
We show the mixing buoyancy layer near the canopy layer
![buoyancy field](canopy_test/b.mp4)

We show the air temperature profiles over time
~~~gnuplot
set xr [270:285]
set yr [0:50]
set xlabel 'Buoyancy' font ",15"
set ylabel 'height' font ",15"
set key box top left
set grid
plot 'proft=1T' u 4:1 w l lw 3 t 't = 1T' , \
'proft=2T' u 4:1 w l lw 3 t 't = 200T',	    \
'proft=3T' u 4:1 w l lw 3 t 't = 300T',     \
'proft=0T' u 4:1 w l lw 3 t 't = initial',
~~~

Let's check the wind speed
~~~gnuplot
set xr [0:0.6]
set yr [0:50]
set xlabel 'wind [m s-1]' font ",15"
set ylabel 'height' font ",15"
set key box top left
set grid
plot 'proft=1T' u 2:1 w l lw 3 t 't = 1T' , \
'proft=2T' u 2:1 w l lw 3 t 't = 200T',	    \
'proft=3T' u 2:1 w l lw 3 t 't = 300T',     \
'proft=0T' u 2:1 w l lw 3 t 't = initial',
~~~

*/
