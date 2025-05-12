/**
[Return to my homepage](http://basilisk.fr/sandbox/nlemoine/README)

# MacDonald test problem #3 -- stationary, periodic subcritical flow  (MacDonald *et al.*, 1997)

## Flow equation

For a 1D stationary flow, the continuity equation in the Saint-Venant system boils down to $q = q_x = \textrm{cst} = q_0$ and the momentum equation becomes:

$$0\quad =\quad  -\partial_x\left[\frac{q^2}{h} + \frac{1}{2}g h^2\right]- gh\ \partial_x z_b-\frac{\tau}{\rho}$$
In the Manning-Strickler model, bed shear stress is given by:
$$\boldsymbol{\tau}\quad=\quad\rho\,g\,n^2 h^{-\frac{1}{3}}\left|\mathbf{u} \right|\mathbf{u}\quad=\quad\rho\,g\,n^2 h^{-\frac{7}{3}}\left|\mathbf{q} \right|\mathbf{q}$$
so that
$$0\quad  =\quad  +\frac{q_0^2}{h^2}\partial_x h- gh\ \partial_x h - gh\ \partial_x z_b - n^2 g\, h^{-\frac{7}{3}}q_0^2$$
Hence
$$\partial_x z_b\quad=\quad \underbrace{\left[\frac{q_0^2}{g h^3}-1\right]}_{\textrm{Fr}^2-1}\partial_x h\quad-\quad n^2 q_0^2\,h^{-\frac{10}{3}}$$
The rationale behind [MacDonald's test case](https://doi.org/10.1061/(ASCE)0733-9429(1997)123:11(1041)) is that of an inverse problem: we determine the bed elevation profile corresponding to a desired water depth profile $h(x)$ and uniform flow rate $q_0$, i.e., the right-hand side of the latter equation is fully known. The method can be straightforwardly extended to the case of a variable Manning coefficient $n(x)$.

For a periodic geometry,  the integration of $\partial_x z_b$ can be easily performed using Fourier transform. Indeed, the Fourier spectrum of the detrended topography $z_b^\ast$ is a function of wavenumber $k$ (in $\textrm{rad}\cdot\textrm{m}^{-1}$) which satisfies:

$$\mathcal{F}[z_b^\ast](k) = \begin{cases}
\displaystyle\frac{1}{ik}\mathcal{F}[\partial_x z_b](k) & \textrm{for }k\neq 0 \\   
0 & \textrm{for }k=0\textrm{ (zero-mean elevation)}
\end{cases}$$
  
The inverse transform yields $z_b^\ast(x)$. The average slope $I=I_x$, hence the linear trend which must be added to $z_b^\ast$, is simply the spatial mean of $\partial_x z_b$ over the length $\lambda$ of a pattern:

$$I\quad=\quad -\frac{1}{\lambda}\int_\lambda (\partial_x z_b) dx\quad \propto \quad\mathcal{F}[\partial_x z_b](k=0)$$

We use the [detrended bathymetry](manning-tilt.h) $z_b^\ast(x)$ as an input to Basilisk and we set $\ \texttt{tilt.x} = I$, together with a periodic boundary condition on all fields.

## Implementation

This test case makes use of Fourier transforms with the [FFTW](https://www.fftw.org/) library.
*/

#include <fftw3.h>
#include <complex.h>
#pragma autolink -lfftw3 
#include "grid/multigrid1D.h"
#include "saint-venant.h"
#include "nlemoine/SWE/manning-tilt.h"

#define LEVEL 7
#define tfin 900.
#define sec_per_day 86400.
#define M_PI acos(-1.0)

#define Lpattern 1000.
#define npatterns 3
#define q0 2.0
#define hmean (9./8.)
#define hrange (1./4.)
#define nharmo 10	// nharmo < N/(2*npatterns)
/**
### Functions returning analytical depth and Manning coefficient for any given position
*/

double get_analytical_depth(double x, int derivative)
{
   double wavenum = 2.0*M_PI/Lpattern;
   double res;
   
   if(!derivative)
     res = hmean + hrange*sin(wavenum*x);
   else
     res = wavenum*hrange*cos(wavenum*x);

   return res;
}

double get_manning(double x)
{
//   return 0.03-0.01*sin(2.0*M_PI*x/Lpattern);
   return 0.03;
}

double AMPLITUDE[nharmo], PHASESHIFT[nharmo];

int init_harmo_zb ()
{
   fftw_complex *in, *out;
   fftw_plan p;
   in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
   out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);

   double xm,hm,nm,dh_dx,Sf,Fr2,dzb_dx;

   /* initialize input array */

   for(int m=0;m<N;m++)
   {
      xm = L0*m/((double)N);
      hm = get_analytical_depth(xm,0);
      dh_dx = get_analytical_depth(xm,1);
      nm = get_manning(xm);
      Sf = sq(q0*nm)*pow(hm,-10./3.); // friction slope
      Fr2 = sq(q0)*pow(hm,-3.)/G; // squared Froude
      if(Fr2>1)
        printf("error : analytical solution is locally supercritical (Fr>1)\n");

      dzb_dx = (Fr2-1.)*dh_dx - Sf ;

      in[m][0] = dzb_dx;
      in[m][1] = 0.;
   }

   /* compute forward Fourier transform of dzb_dx */

   p = fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
   fftw_execute(p);
   fftw_destroy_plan(p);

   /* retrieve average slope (real part of 1st element out[0] */

   tilt.x = (-1./N)*out[0][0];
   printf("tilt.x = %g\n",tilt.x);
   in[0][0] = 0.;
   in[0][1] = 0.;

   /* integrate in Fourier space */
   
   double k0 = (2.*M_PI/L0);

   for(int m=1;m<N;m++)
   {
      double wavenum = m<(N/2) ? k0*m : k0*(m-N);
      in[m][0] = (1./wavenum)*out[m][1];
      in[m][1] = -(1./wavenum)*out[m][0];
   }

/*   // compute inverse Fourier transform to get detrended bathymetry

   p = fftw_plan_dft_1d(N, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);
   fftw_execute(p);
   fftw_destroy_plan(p);

   // write result

   if(pid()==0){
     FILE * fp = fopen("FFT.txt","w");
     for(int m=0;m<N;m++)
       fprintf(fp,"%g %g %g %g\n",in[m][0],in[m][1],(1./N)*out[m][0],(1./N)*out[m][1]);
     fclose(fp);
   }*/

   // WARNING : fftw3 does not behave like Scilab, we have to scale the output by 1/N

/**
In order to be able to evaluation bed elevation at any position, we compute the amplitude $A_m$ and phase shift $\phi_m$ for a given number of harmonics so that
$$z_b^\ast(x) = \sum_{m=1}^{N_{\textrm{harmo}}}A_m \cos\left(m\,k_\textrm{pattern}\,x +\phi_m\right)$$
where $k_\textrm{pattern} = \displaystyle\frac{2\pi}{\lambda_\textrm{pattern}}$ is the pattern's wave number.
*/
   for(int m=0;m<nharmo;m++)
   { 
      double complex zz = in[npatterns*m][0] + I*in[npatterns*m][1];
      AMPLITUDE[m] = 2.0*cabs(zz)/N;
      PHASESHIFT[m] = cabs(zz)>0. ? cimag(clog(zz)) : 0.;
   }

   /* clean-up */

   fftw_free(in);
   fftw_free(out);

   return(0);
}
/**
Now we can define the function returning bed elevation for any position
*/

double get_zb(double x)
{
   double res = 0., kpattern = 2.*M_PI/Lpattern;

   for(int m=1;m<nharmo;m++)
     res += AMPLITUDE[m]*cos(kpattern*m*x+PHASESHIFT[m]);

   return res;
}

/**
### Main
*/

int main (int argc, char * argv[])
{
  G = 9.81;
  N = 1 << LEVEL;
  L0 = npatterns*Lpattern;
  size (L0);

  periodic(right);
  run();
  return(0);
}

/**
In order to ensure convergence towards the analytical solution, we have to initialize with the correct volume in the domain. Since $z_b^\ast$ has zero-mean, we initialize with $h_\textrm{ini}(x) = \overline{h} - z_b^\ast(x)$ i.e. $\eta_\textrm{\,ini}(x) = \overline{h}$. It is a ``lake-at-rest'' condition in the detrended frame i.e. a flowline of constant slope $\texttt{tilt.x}$ (see [manning-tilt.h](../manning-tilt.h)).
*/

event init (i=0)
{
   (void) init_harmo_zb ();

   FILE * fp = fopen("zb.txt","w");
   foreach()
   {   
      zb[] = get_zb(x);
      h[] = hmean-zb[];
      eta[] = hmean;
      u.x[] = 0.;
      nmanning[] = get_manning(x);
      fprintf(fp,"%g %g\n",x,zb[]);
   }
   fclose(fp);

   boundary(all);

   DT = 0.5;
}

event snapshot (t=0 ; t+=60.) {

    char name[100];
    FILE * fp;
    sprintf(name,"profile-%g.dat",t);
    fp=fopen(name,"w");
    foreach() {
      fprintf(fp,"%g %g %g %g %g %g %g %g\n"
	      ,x,zb[],h[],get_analytical_depth(x,0),u.x[],zb[]-x*tilt.x,h[]+zb[]-x*tilt.x,get_analytical_depth(x,0)+zb[]-x*tilt.x);
    }
    fclose(fp);
}

event stop (t = tfin)
{
    // Write final profile in log file
    foreach()
      fprintf(stderr,"%g %g %g %g %g %g %g %g\n"
	      ,x,zb[],h[],get_analytical_depth(x,0),u.x[],zb[]-  x*tilt.x,h[]+zb[]-x*tilt.x,get_analytical_depth(x,0)+zb[]-x*tilt.x);
}

/**
~~~gnuplot
   reset
   set xlabel "x (m)"
   set ylabel "elevation (m)"
   ydim = 800 
   xdim = 480
   plot 'profile-0.dat' u 1:6 w lines lc rgb "black" lw 3 title 'zb', \
        'profile-0.dat' u 1:7 w point pointtype 6 ps 0.8 lc rgb "red" title 'initial', \
        'profile-900.dat' u 1:7 w point pointtype 6 ps 0.8 lc rgb "blue" title 't = 15 min', \
        'profile-900.dat' u 1:8 w lines lc rgb "green" lw 2 title 'analytical'
        
~~~
*/