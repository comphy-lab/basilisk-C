/** 
St. Venant using Wave Propagation Algorithm 
*/

#include "grid/quadtree.h"

#define USE_WAVEPROP 1

#if USE_WAVEPROP
# include "rpn2_swe.h"
#else
# include "saint-venant.h"
scalar * scalars = {h};  /* Needed for conservation check */
#endif

#define MAXLEVEL 8
#define MINLEVEL 5

static double sum0[3];
static double threshold = 4e-3;

int main() 
{
    origin(-0.5,-0.5);

    //CFL = 0.9;  /* Use defaults set by solver. */

    periodic(right);
    periodic(top);

    N = 1 << MINLEVEL;
    run();
}

event init (i = 0) 
{
#if USE_WAVEPROP
    dt_fixed = false;
    dt_initial = 2.5e-3;   /* to get CFL= 0.9 on first step */

    conservation_law = true;
    order = 2;
#endif    

#if TREE
    astats s;
    int k = 0;
    do {
        foreach()       
        {
            h[] = 0.1 + 1.*exp(-200.*(x*x + y*y));
        }
        boundary (all);      
        s = adapt_wavelet ({h}, (double []){threshold}, maxlevel = MAXLEVEL);
        fprintf(stdout,"Initial conditions : k = %d\n",k);
        if (k > MAXLEVEL-MINLEVEL)
        {
            fprintf(stdout,"Initial conditions : Maximum refinement exceeded\n");
            break;
        }
        k++;
    } while (s.nf > 0);
#else
    foreach()       
    {
        h[] = 0.1 + 1.*exp(-200.*(x*x + y*y));
    }
#endif    
    boundary(scalars);
}


event logfile (i++)
{
#if USE_WAVEPROP    
    vector u[];
    foreach()
    {
        u.x[] = hu.x[]/h[];   
    }
#endif    
    fprintf (stderr, "%g %g\n", t, normf(u.x).max);    
}

event movie (i += 10)
{
  output_ppm (h, linear = true, file = "h.mp4", n = 512);
  scalar l[];
  foreach()
      l[] = level;
  output_ppm (l, file = "level.mp4", min = MINLEVEL, max = MAXLEVEL, n = 512);
}

event conservation_check (i += 160)
{
    double sum[3] = {0,0,0};
    foreach()
    {
        int m = 0;
        for(scalar s in scalars)
        {
            sum[m] += Delta*Delta*s[];
            m++;
        }
    }

    /* Save initial mass for comparison */
    if (i == 0)
    {
        for(int m = 0; m < 3; m++)
            sum0[m] = sum[m];
    }

    fprintf(stdout,"Conservation check: Step %3d : t = %16.8e\n",i,t);
    int m = 0;
    for(scalar s in scalars)
    {
        fprintf(stdout,"sum[%d] = %24.16f %24.16e\n",m, sum[m], fabs(sum0[m]-sum[m]));
        m++;
    }
    fprintf(stdout,"\n");
}


event end (t = 4.0)
{
  printf ("i = %d t = %g\n", i, t);
}


#if TREE
event adapt (i++) 
{
  adapt_wavelet ({h}, (double []){threshold}, maxlevel = MAXLEVEL);
}
#endif

/**
Here are the results.

![Animation of free surface height](swe/h.mp4)

![Animation of the level of refinement](swe/level.mp4)

~~~gnuplot Maximum velocity as a function of time
set xlabel 'Time'
set ylabel 'Maximum velocity'
unset key
set yrange [0 to 0.6]
plot 'log' w l
~~~
*/




