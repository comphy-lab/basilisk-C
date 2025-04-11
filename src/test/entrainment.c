/**
# Entrainment case from GOTM

See [the GOTM web site](https://gotm.net/cases/entrainment/) and
section 12.1.4 of [the GOTM manual](#gotm).

"The entrainment scenario is ideally suited to benchmark the model
performance in stress-driven entrainment situations against available
experiments (see [Umlauf and Burchard, 2005](#umlauf2005)) The results
shown in the figure below illustrate that after the onset of the
constant surface stress in the x-direction, a thin near-surface layer
is accelerated (top figure), and slowly entrains into the stratified,
non-turbulent interior region. Shear-driven turbulence in this region
is mirrored in the large turbulent diffusivities shown in the third
figure, which generate a nearly well-mixed surface layer that is
separated from the interior by a pycnocline of gradually increasing
strength (second figure)."

## Results

~~~gnuplot Evolution of the horizontal velocity (m/s)
set term PNG size 800,400
set output 'u.png'
set pm3d map interpolate 4,4
# the matlab 'parula' colormap
set palette defined (0 '#352a87',1 '#0363e1',2 '#1485d4',3 '#06a7c6',\
           4 '#38b99e',5 '#92bf73',6 '#d9ba56',7 '#fcce2e',8 '#f9fb0e')
set cbrange [0:0.45]
set ylabel 'z (m)'
set xlabel 'time (hours)'
set xrange [0:24]
set yrange [-35:0]
splot 'kepsilon' u ($1/3600.):2:3
~~~

~~~gnuplot Evolution of the squared buoyancy frequency $N^2$ (s^-2^) (x 10^4^)
set output 'N2.png'
set cbrange [*:*]
splot 'kepsilon' u ($1/3600.):2:($5*1e4)
~~~

~~~gnuplot Evolution of the (log of the) turbulent diffusivity $\nu_t^h$ (m^2^/s)
set output 'nut.png'
set cbrange [-6:-2.5]
splot 'kepsilon' u ($1/3600.):2:(log10($6))
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

@article{umlauf2005,
  title = {Second-order turbulence closure models for geophysical 
            boundary layers. A review of recent work},
  journal = {Continental Shelf Research},
  volume = {25},
  number = {7},
  pages = {795 - 827},
  year = {2005},
  doi = {https://doi.org/10.1016/j.csr.2004.08.004},
  author = {Lars Umlauf and Hans Burchard}
}

@article{price1979, 
  title = {On the scaling of stress-driven entrainment experiments}, 
  volume = {90}, 
  DOI = {10.1017/S0022112079002366}, 
  number = {3}, 
  journal = {Journal of Fluid Mechanics}, 
  publisher = {Cambridge University Press}, 
  author = {Price, James F.}, 
  year = {1979}, 
  pages = {509–529}
}
~~~
*/

#include "grid/multigrid.h"
#include "layered/hydro.h"
#include "layered/gotm.h"

/**
We add the temperature field for each layer. */

event defaults (i = 0)
{
  T = new scalar[nl];
}

event cleanup (t = end) {
  delete ({T});
}

/**
The code runs with three different turbulence models (which give very
similar results). */

char * casename;
FILE * casefp = NULL;

int main()
{
  periodic (right);
  size (10e3);
  G = 9.81;
  N = 1;
  nl = 100;

  /**
  The surface wind stress is constant. */
  
  airsea_tau[] = { 0.1027 };

  meanflow_maxitz0b = 1;	// default = 10

  /**
  ## Turbulence mode setup

  ### The generic mixing length model
  */
  
  // generic.xml
  turbulence_alpha = 0.256;	// default = 0
  turbulence_turb_method = 3;	// default = 2
  turbulence_len_scale_method = 10;	// default = 8
  turbulence_stab_method = 1;	// default = 3
  turbulence_compute_kappa = 1;	// default = 0
  turbulence_compute_c3 = 1;	// default = 0
  
  turbulence_gen_m = 1;	// default = 1.5
  turbulence_gen_n = -0.67;	// default = -1
  turbulence_cpsi1 = 1;	// default = 1.44
  turbulence_cpsi2 = 1.22;	// default = 1.92
  turbulence_sig_kpsi = 0.8;	// default = 1
  turbulence_sig_psi = 1.07;	// default = 1.3

#if 1 // necessary for all models
  turbulence_k_min = 1e-10;	// default = 1e-08
  turbulence_eps_min = 1e-14;	// default = 1e-12
  turbulence_kb_min = 1e-10;	// default = 1e-08
  turbulence_epsb_min = 1e-14;	// default = 1e-12

  turbulence_scnd_method = 2;	// default = 0
  turbulence_kb_method = 1;	// default = 0
  turbulence_epsb_method = 1;	// default = 0
  turbulence_scnd_coeff = 5;	// default = 0
  turbulence_cc1 = 5;	// default = 0
  turbulence_cc2 = 0.8;	// default = 0
  turbulence_cc3 = 1.968;	// default = 0
  turbulence_cc4 = 1.136;	// default = 0
  turbulence_cc6 = 0.4;	// default = 0
  turbulence_ct1 = 5.95;	// default = 0
  turbulence_ct2 = 0.6;	// default = 0
  turbulence_ct3 = 1;	// default = 0
  turbulence_ct5 = 0.3333;	// default = 0
  turbulence_ctt = 0.72;	// default = 0
  turbulence_nuhiw = 1e-05;	// default = 5e-05
#endif
  
  casename = "generic";
  casefp = fopen (casename, "w");
  run();

  /**
  ### The k-$\epsilon$ model
  */
  
  // kepsilon.xml
  turbulence_alpha = 0;	        // default = 0  
  turbulence_turb_method = 3;	// default = 2
  turbulence_len_scale_method = 8; // default = 8
  turbulence_stab_method = 1;	// default = 3
  turbulence_compute_kappa = 1;	// default = 0
  turbulence_compute_c3 = 1;	// default = 0

  turbulence_gen_m = 1.5;	// default = 1.5
  turbulence_gen_n = -1;	// default = -1
  //  turbulence_gen_p = 3;	        // default = 3
  turbulence_cpsi1 = 1.44;	// default = 1.44
  turbulence_cpsi2 = 1.92;	// default = 1.92
  turbulence_sig_k = 1;
  turbulence_sig_kpsi = 1;	// default = 1
  turbulence_sig_psi = 1.3;	// default = 1.3
  
  turbulence_cpsi3plus = 0;	// default = 1
  turbulence_gen_d = -1.087;	// default = -1.2
  turbulence_gen_alpha = -4.97;	// default = -2
  turbulence_gen_l = 0.09;	// default = 0.2

  casefp = fopen ("kepsilon", "w");
  casename = "kepsilon";
  casefp = fopen (casename, "w");
  run();

  /**
  ### The k-$\omega$ model
  */
  
  // komega.xml
  turbulence_alpha = 0.256;	// default = 0
  turbulence_turb_method = 3;	// default = 2
  turbulence_len_scale_method = 10;	// default = 8
  turbulence_stab_method = 1;	// default = 3
  turbulence_compute_kappa = 1;	// default = 0
  turbulence_compute_c3 = 1;	// default = 0
  turbulence_gen_m = 0.5;	// default = 1.5
  turbulence_gen_p = -1;	// default = 3
  turbulence_cpsi1 = 0.55555;	// default = 1.44
  turbulence_cpsi2 = 0.83333;	// default = 1.92
  turbulence_sig_kpsi = 2;	// default = 1
  turbulence_sig_psi = 2;	// default = 1.3
  turbulence_gen_alpha = -2.53;	// default = -2
  turbulence_gen_l = 0.25;	// default = 0.2

  casename = "komega";
  casefp = fopen (casename, "w");
  run();
}

/**
## Initial conditions

A constant temperature stratification with buoyancy frequency $N^2 =
10^{-4}$ s$^{-2}$, starting at 20$^\circ$C and with a constant
background salinity of $20$ psu. 

The depth is 50 metres and the layers are equal thickness. */
 
event init (i = 0)
{
  turbulence_report_model();
  foreach() {
    zb[] = -50.;
    foreach_layer()
      h[] = 50./nl;
  }
  constant_NNT (20, 0, 1e-4, T);
}

/**
The timestep is set to 60 seconds. */
 
event timestep (t += 60);

/**
## Outputs */
 
void profile (FILE * fp)
{
  foreach() {
    double z = zb[];
    foreach_layer()
      fprintf (fp, "%g %g %g %g %g %g\n", t, z + h[]/2.,
	       u.x[], T[],
	       meanflow_nn.a[point.l + 1], turbulence_nuh.a[point.l + 1]),
      z += h[];
  }
  fprintf (fp, "\n");  
}

event profiles (t += 240)
  profile (casefp);

event end (t = 24*3600) {
  fprintf (stderr, "# %s\n", casename);
  profile (stderr);
}
