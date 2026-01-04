/**
## Modified wavenumber analysis of 12-th order accurate CCD scheme

We analyse the combined-compact finite difference (ccd) scheme
presented [here](order12_ccd.c). We use a numerical approximation of the modified wavenumber. 

~~~gnuplot Modified wavenumber analysis
set xr [0:3.15]
set yr [0:3.15]
set xlabel 'Delta k [-]'
set ylabel 'Delta k^* [-]'
set grid
set size square
set key top left
plot 'out' w l lw 4 t '2nd order central', sin(x) w l lw 3 t 'analytical 2nd order', x w l lw 2 t 'perfect scheme', '' u 1:3 w l lw 2 t '4th order pade', '' u 1:4 w l lw 2 t '12th order ccd'
~~~
 */
#include "grid/multigrid1D.h"
#include "solve.h"

int main() {
  periodic (left);
  N = 64; {
    init_grid(N);
    scalar s[], c[], ds[], dds[], ds_pade[];
    for (int n = 100; n >= 1 ; n--) {
      double k = n*2*pi/L0;
      foreach() {
	ds_pade[] = 0;
	ds[] = 0;
	dds[] = 0;
	s[] = sin(k*x);
	c[] = cos(k*x);
      }
      
      double k_star = 0, k_star_pade = 0, k_star_ccd = 0;
      int l = 0;
      solve (ds_pade, 1./4.*(ds_pade[-1] + ds_pade[1]) + ds_pade[],
	     3./4.*(s[1] - s[-1])/Delta);
      for (int j = 0; j < 50; j++) {
	foreach() {
	  ds[] = (+445./2592.*(s[2] - s[-2])/Delta
		  - 100./81.*(s[1] - s[-1])/Delta
		  + 23./432.*(ds[2] + ds[-2])
		  - 4./9.*(ds[1] + ds[-1])
		  + 1./216.*(dds[2] - dds[-2])*Delta
		  - 4./27.*(dds[1] - dds[-1])*Delta); 
	}
	foreach() {
	  dds[] = (-65/324.*(s[2] + s[-2])/sq(Delta)
		   + 320./81.*(s[1] + s[-1])/sq(Delta)
		   - 15./2.*(s[])/sq(Delta)
		   - 25/432.*(ds[2] - ds[-2])/Delta
		   + 40./27.*(ds[1] - ds[-1])/Delta
		   - 1./216.*(dds[2] + dds[-2])
		   + 8./27.*(dds[1] + dds[-1])); 
	}
      }
      foreach(reduction(+:l) reduction(+:k_star)) {
	k_star += fabs(c[]) > 0 ? (s[1] - s[-1])/(2*Delta) / c[] : 0;
	k_star_pade += fabs(c[]) > 0 ? ds_pade[]/c[] : 0;
	k_star_ccd += fabs(c[]) > 0 ? -ds[]/c[] : 0;
	l += fabs(c[]) > 0 ? 1 : 0;
      }
      printf ("%g %g %g %g\n", L0/N*k, L0/N*k_star/l, L0/N*k_star_pade/l, L0/N*k_star_ccd/l);
    }
  }
}