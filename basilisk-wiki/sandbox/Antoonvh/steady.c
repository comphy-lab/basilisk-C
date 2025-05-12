/**
# Steady-state respecting drift or diffusion

We consider an erchodic drift-and-diffusion system that is
characterized by a bimodality in the long-term histogram. Here we show
methods to find the drift *or* diffusion fields as a function of
the steady-state PDF and the diffusion *or* drift field,
respectively.

## Examples:

We compute the Drift coefficent as a function of a prescribed steady
distribution and diffusion field.

![Evolution of the PDF towards the steady stady based on a prescribed
 steady state and  diffusion field](steady/mov0.mp4)

Next, we also compute a diffusion field from a prescribed steady state
and drift field.

![Evolution of the PDF towards the steady stady based on a prescribed
 steady state and  drift field](steady/mov.mp4)
 */
#include "grid/multigrid1D.h"
#include "fpe.h"
#include "utils.h"

#define PDFS ((exp(-sq(x - 1.)) + exp(-sq(x + 1.)))/(sqrt(pi)/5.))

void find_D1 (scalar rhos, scalar D2, face vector D1) {
  foreach_face () {
    D1.x[] = 0.5*(D2[] + D2[-1])/2*((log(D2[]) - log(D2[-1]))/Delta +
				    (log(rhos[0]) - log(rhos[-1]))/Delta);
  }
}

void find_D2 (scalar rhos, scalar D2, face vector D1) {
  scalar  D1c[];
  foreach() {
    D1c[] = (D1.x[1] + D1.x[])/2.;
    D2[] = log(D2[]);
  }
  for (int j = 0; j < 99; j++) {
    foreach() 
      D2[] = D2[]/2. + (D2[-1] + Delta*(2*D1c[]/exp(D2[]) -
					(log(rhos[]) - log(rhos[-1]))/(Delta)))/2;
    foreach() 
      if (x > 0)
	D2[] = interpolate (D2, -x);
  }
  foreach()
    D2[] = exp (D2[]);
}

int mode;

int main() {
  L0 = 10;
  X0 = -5;
  N = 128;
  DT = 0.001;
  run();
  mode = 1;
  run();
}
 
  scalar pdfs[];
event init (t = 0) {
  foreach()
    pdfs[] = PDFS;
  pdfs[left] = dirichlet (PDFS);
  pdfs[right] = dirichlet (PDFS);
  foreach()
    rho[] = exp(-sq(x - 2));
  stats b = statsf(rho);
  foreach()
    rho[] /= (b.sum/b.volume);

  if (mode == 0) {
    foreach() 
      D2[] = 1 + x/10;
    find_D1 (pdfs, D2, D1);
  } else {
    foreach()
      D2[] = 1.;
    foreach_face()
      D1.x[] = -x/2;
    find_D2 (pdfs, D2, D1);
  }
}

event stop (t = 15) {
  return 1;
}

/**
## Movie making
*/
FILE * gnuplotPipe = NULL;
int frame = 0;

event init (t = 0) {
  frame = 0;
  if (gnuplotPipe == NULL) {
    gnuplotPipe = popen ("gnuplot", "w");
    fprintf(gnuplotPipe,
            "set term pngcairo\n"
            "set xr [-5: 5]\n"
            "set yr [-5.2: 5.2]\n"
            "set key box top left\n"
            "set grid\n"
            "set title 'PDF evolution'\n"
            "set xlabel 'x'\n"
	    "set key bottom left\n"
            "set ylabel 'pdf, D_1, D_2'\n");
  }
}

event movie(t += 0.1){
  fprintf(gnuplotPipe, "set output 'plot%d.png'\n", frame);
  fprintf(gnuplotPipe, "plot \
          '-' w l lw 5 t 'pdf',			\
          '-' w l lw 2 t 'Steady solution' ,	\
	  '-' w l lw 2 t 'Diffusion',			\
	  '-' w l lw 2 t 'Drift'\n");
  for (scalar s in {rho, pdfs, D2}) {
    foreach()
      fprintf(gnuplotPipe, "%g\t%g\n",x, s[]);
    fprintf(gnuplotPipe, "e\n");
  }
  foreach_face()
    fprintf(gnuplotPipe, "%g\t%g\n",x, D1.x[]);
  fprintf(gnuplotPipe, "e\n");
  
  
  frame++;
}

event finalize_movie (t = end) {
  system("ffmpeg -r 25 -f image2 -i plot%d.png -c:v libx264 -vf format=yuv420p -y mov.mp4");
  system("rm plot*");
  if (mode == 0)
    system ("mv mov.mp4 mov0.mp4");
  pclose (gnuplotPipe);
  gnuplotPipe = NULL;
}
