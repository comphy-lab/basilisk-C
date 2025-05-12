/**
# Tempeature-gradient sharpening in a canopy.
 
The evolution of a temperature profile is studied with the minimal(?)
ingredients for stable layer formation. This is the [flux version of related axample](canopy.c).  Consider a canopy (top at $z = z_c$) with a fixed temperature flux `Flx`. The evolution of the temperature ($T$) profile is described with an effective diffusivity $K$,

$$\frac{\partial T}{\partial t} = \frac{\partial}{\partial z}K\frac{\partial T}{\partial z }.$$

$K$ is *modelled* by considering 3 effects:   

* A constant Background diffusivity due to shear, top down mixing etc ($K_b$).   
* A Height-dependend diffusivity correction for the canopy details.  
* A stability correction.  

In a dimenonless setting, we choose $z_c = 1, \Delta T = 1$ and $K_b =
1$. The canopy correction defines $K_c$:

$$K_c = K_b \left( 1 + 4\frac{C_k}{z_c^2} \left(z\left(z - z_c\right)\right)\right),$$  
*/
#define Kc (Kb*(1 + 4*Ck*(x*(x - 1)))) 
/**
Which gives a maximum of $K_b$ at top and bottom and a minimum of
$(1-C_k)K_b$ halfway the canopy. This is to model the effects of the
turbulent air induced by the surface of the Earth and the canopy
crown.

The stability correction gives $K_s$,

$$K_s = K_c e^{-a\frac{\partial T}{\partial z}},$$

with $a$ an inverse temperature-gradient scale.

Finally, to omit the most prominent issue with this description,
there is an minimum value for $K$; $K_{\mathrm{min}}$.

$$K = \mathrm{max} \left(K_s, K_{\mathrm{min}}\right).$$
*/
#define K(grad) (max(Kc*exp(-a*(grad)), Kmin))
/**
## Numerical set-up

The system is set up and solved for 
 */
#include "grid/multigrid1D.h" // a 1D grid
#include "diffusion.h"         
#include "run.h"              // Timeloop
/**
The values of the physical parameters are chosen to be *just*
collapse.
 */
double Kb = 1, Kmin = 0.025, Ck = 0.4, a = 0.6;
double Flx = 0.368; // > (1 - Ck)/(exp(1)*a) 
double tend = 15;

scalar T[], Th[];
T[left]  = neumann (-Flx/K(face_gradient_x(Th,0)));  //bottom
T[right] = neumann ( Flx/K(face_gradient_x(Th,0)));  //top

int main() {
  N  = 512;
  TOLERANCE = 1e-4;
  DT = 0.001;
  run();
}

event init (t = 0) {
  foreach()
    T[] = x;                // Linearly stratified
  boundary ({T});
}

event diff (i++) {
  dt = dtnext (DT);
  face vector D[];
  foreach_cell_all()
    Th[] = T[];
  foreach_face() 
    D.x[] = K (face_gradient_x (T, 0));
  boundary ({D.x});
  diffusion (T, dt, D);
}
/**
## Output

![A layer emerges](canopyflx/mov.mp4)

The movie is generated with `gnuplot` and `ffmpeg`.
*/
FILE * gnuplotPipe;

event init (t = 0) {
  gnuplotPipe = popen ("gnuplot", "w");
  fprintf(gnuplotPipe,
	  "set term pngcairo\n"
	  "set xr [0: 1]\n"
	  "set yr [0: 1]\n"
	  "set key top left\n"
	  "set grid\n"
	  "set title 'Temperature and its Diffusivity'\n"
	  "set xlabel 'T/ΔT, K/K_b'\n"
	  "set ylabel 'z/z_c'\n");
}

int frame = 0;
event movie (t += 0.1) {
  fprintf (gnuplotPipe,
	   "set output 'plot%d.png'\n"
	   "plot 0.6*(x-0.5) + 0.5 w l lw 1 t 'Max. gradient', ",
	   frame);
  fprintf (gnuplotPipe, "'-' w l lw 5 t 'T', '' w l lw 4 t 'K'\n");
  foreach()
    fprintf (gnuplotPipe, "%g %g\n", T[], x);
  fprintf (gnuplotPipe, "e\n");
  foreach_face()
    fprintf (gnuplotPipe, "%g %g\n", K (face_gradient_x(T, 0)), x);
  fprintf (gnuplotPipe, "e\n");
  frame++;
}

event stop (t = tend) {
  system ("rm mov.mp4");
  system ("ffmpeg -r 25 -f image2 -i plot%d.png -c:v libx264 \
          -vf format=yuv420p -y mov.mp4");
  system ("rm plot*");
  return 1;
}

/**
## A criterion for layer formation

A steady state solution is characterized by a constant flux:

$$K\frac{\partial T}{\partial z} = H$$

For convinience we write,

$$y = \frac{\partial T}{\partial z}$$

then, 

$$K_c e^{-\alpha y}y = H,$$

There exists a maximum for H when $y = \alpha^{-1}$. At this value for
$y$ the flux decreases with increasing gradients. It is well known
that such systems Collapse (van de Wiel et al., van de Wiel et
al.). Thus we have a constraint $y < \frac{1}{\alpha}$.

$$H_{\mathrm{max}} = \frac{K_c}{e\alpha},$$

$$H_{\mathrm{max}} = \frac{K_b \left(1-4C_k(z^2 - z)\right)}{e\alpha}.$$

The most sustainable heatflux is the minimal value of
$H_{\mathrm{max}}(z)$, which is at height $z = 0.5z_c$,

$$\frac{H_{\mathrm{max, s}}}{K_b} = \frac{1 - C_k}{e\nalpha}.$$
*/
