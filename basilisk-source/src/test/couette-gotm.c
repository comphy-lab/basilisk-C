/**
# Turbulent Couette flow

See [the GOTM web site](https://gotm.net/cases/couette) and
section 12.1.1 of [the GOTM manual](#gotm).

"The Couette scenario is the most basic of all GOTM scenarios. It
represents a shallow (10 m deep), unstratified layer of fluid above a
flat bottom that is driven by a constant surface stress in the
x-direction. Earth's rotation is ignored. This flow is often referred
to as turbulent Couette flow. After the onset of the surface stress, a
thin turbulent near-surface layer is generated that rapidly entrains
into the non-turbulent deeper parts of the water column. The solution
at the end of the simulation, when the problem has become fully
stationary, is shown in the figure below."

## Results

~~~gnuplot Stationary velocity profile
set term SVG size 600,300
set xlabel 'u (m/s)'
set ylabel 'z (m)'
plot [0:1][-10:0]'log' u 2:1 w l t '' lw 2
~~~

~~~gnuplot Turbulent diffusivity
set xlabel 'ν_t (m^2/s)'
plot [0:0.05][-10:0]'log' u 3:1 w l t '' lw 2
~~~

## References

~~~bib
@manual{gotm, 
  title = {{GOTM} Source Code and Test Case Documentation}, 
  author = {Lars Umlauf and Hans Burchard and Karsten Bolding}, 
  edition = {version 4},
  year = {2018}, 
  pdf = {https://gotm.net/manual/stable/pdf/a4.pdf}
}
~~~
*/

#include "grid/multigrid1D.h"
#include "layered/hydro.h"
#include "layered/gotm.h"

int main()
{
  G = 9.81;
  N = 1;
  nl = 100;
  DT = 20;
  size (100e3);
  periodic (right);

  /**
  We use the k-$\epsilon$ model. */
  
  turbulence_turb_method      = turbulence_first_order;
  turbulence_tke_method       = turbulence_tke_keps;
  turbulence_len_scale_method = turbulence_generic_eq;  

  /**
  The bottom roughness leng scale (of GOTM) needs to be adjusted. */
  
  meanflow_z0s_min = 0.003;
  meanflow_h0b = 0.1;

  /**
  The surface wind stress is constant. */

  airsea_tau[] = { 1.027 };
  
  run();
}

/**
## Initial conditions

Ten metre deep, constant layer thicknesses. */

event init (i = 0)
{
  foreach() {
    zb[] = -10.;
    foreach_layer()
      h[] = 10./nl;
  }
}

/**
## Outputs

24 hours is enough to reach a stationary profile. The turbulent
diffusivity $\nu_t$ computed by GOTM is accessed through the [C
interface](/src/gotm/common.h) of the corresponding Fortran field. */

event profiles (t = 24*3600)
{
  foreach() {
    double z = zb[];
    foreach_layer()
      fprintf (stderr, "%g %g %g\n", z + h[]/2., u.x[],
	       (turbulence_num.a[point.l] + turbulence_num.a[point.l + 1])/2.),
      z += h[];
  }
  fprintf (stderr, "\n");
}
