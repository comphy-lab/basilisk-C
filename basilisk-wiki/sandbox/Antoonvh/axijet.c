/**
# The Creation of a laminar vortex ring

Instead of leading a [turbulent jet](smokey_puff.c), a vortex ring can
be part of a laminar flow. According to the literature, for a
stroke-to-radius ratio of the injecting piston $S/d = 4$, a nice-enough
vortex ring is generated. On this page we study the grid sensitivity
of the formation process. Also we aim to find values of three
different formulations of the reynolds number:

The first is the one we control in the experiment and is based on the
injection velocity ($U_j$), the orifice radius ($R$) and the fluids
viscosity ($\nu$),

 $$Re_1 = \frac{U_jR}{\nu} = 2000.$$

Alternatively, in their seminal work, T.T. Lim and T.B. Nickels(1992) defined a Reynolds number based on the
size of the vortex ring ($D$) and translation velocity ($U_t$). Noting
that it quite ambiguous what is meant with "vortex size". Anyhow, 

$$Re_{\text{lim}} = \frac{U_tD}{\nu},$$

Finally, it seems that in the literature the ambiguity is elimited by
using the circulation of the vortex structure ($\Gamma$), which is
defined as the integral of the velocity vector projected along an
infinite line that passes trough the centre of the ring.

$$Re_{\Gamma} = \frac{\Gamma}{\nu}$$
*/

#include "axi.h"
#include "navier-stokes/centered.h"
#include "fractions.h"
#define JET (tanh((1-y)*5) * (y<1))

scalar f[];
double ti = 5;
double ue; 
int maxlevel = 11 ;
double ov = 2500;
int j = 0;
int m = 0;
const face vector muc[] = {1./ov, 1./ov};
FILE * fpRe;
double xx[2], yy[2], circ[2];

u.t[left] = dirichlet(0.);
u.n[left] = dirichlet(JET * min(ti - t, 1.) * min(t, 1.) * (t <= ti)); 
u.n[top] = neumann(0.);
p[top] = dirichlet(0.);
/**
## Runs

Sixs runs are performed for thee different refinement criteria ($\zeta$) and two maximum levels of refinement (ML);

1. ML = 9, $\zeta = 0.02 U_j$  
2. ML = 9, $\zeta = 0.015 U_j$  
3. ML = 9, $\zeta = 0.01 U_j$  
4. ML = 10, $\zeta = 0.02 U_j$ 
5. ML = 10, $\zeta = 0.015 U_j$  
6. ML = 10, $\zeta = 0.01 U_j$
*/

int main(){
  fpRe = fopen("Reynoldsnumbers", "w");
  L0 = 16;
  init_grid(32);
  foreach_dimension()
    u.x.refine = refine_linear;
  for (ue = 0.01; ue >= 0.001; ue /= 2.){
    circ[0] = 0;
    circ[1] = 0;
    j++;
    printf("running j = %d, ue = %g, ml = %d\n", j, ue, maxlevel);
    run();
  }
}

event init(t = 0){
  mu = muc;
  refine(x < (X0 + L0/10) && y < 1.1 && level < maxlevel);
  f.refine = f.prolongation = fraction_refine;
  fraction(f, 1.-y);
  boundary({f});
  restriction({f});
}

event adapt(i++)
  adapt_wavelet((scalar*){u}, (double[]){ue, ue}, maxlevel);


/**
## Output

Here we discuss the output

### Movie

First a movie is generated for the fifth experiment.
*/
event movie(t += 0.1; t <= 5*ti){
  scalar omega[];
  vorticity(u, omega);
  char pname[99];
  sprintf(pname,"ppm2mp4 j%d.mp4", j);
  static FILE * fp = popen (pname, "w");
  output_ppm(omega, fp, n = 512, min = -5, max = 5, box = {{0, 0},{8, 4}});
}
/**
![Experiment number 5](axijet/j5.mp4)

### Reynolds numbers

We calculate $Re_{\text{lim}}$ and $Re_{\Gamma}$. We define $D$
as the distance between the centres of vorticity of the sections of
the ring with oposing swirl direction. Their values correspond to
those at $t = 4\times t_i$, where $t_i$ is the duration of the
injection. The results are written to a file.
 */

event Re(t=3*ti; t+= ti){
  double xp = 0, yp = 0, w = 0;
  scalar omega[];
  vorticity(u, omega);
  foreach(){
    w += fabs(omega[]) * sq(Delta);
    xp += fabs(omega[]) * sq(Delta) * x;
    yp += fabs(omega[]) * sq(Delta) * y;
    circ[m] += omega[] * sq(Delta);
  }
  xx[m] = xp / w;
  yy[m] = yp / w;
  m++;
  if (m == 2){
    double Regamma = ((circ[0] + circ[1])/2.)*ov;
    double U = (abs(xx[1] - xx[0])) / (ti);
    double D = (yy[1] + yy[0]);
    double ReUDv  = U*D*ov; 
    fprintf(fpRe, "%d\t%g\t%g\n", j, Regamma, ReUDv);
    fflush(fpRe);
    m = 0;
  }
}

/**
and plotted here

~~~gnuplot
set terminal pngcairo enhanced font 'Verdana,12'
set xr [0.5:6.5]
set yr [0:4000]
set xlabel 'Run number'
set ylabel 'Re'
set key outside
plot 'Reynoldsnumbers' u 1:3 pt 2 t 'Re lim' , \
     'Reynoldsnumbers' u 1:2 pt 3 t 'Re Gamma'
~~~

It appears for once that a higher resolution mesh corresponds to
a lower Reynolds numbers. Let us investigate it a bit more...

### Vorticity structure 

Three iso lines of the vorticiy field are outputted for all six runs .
 */
event vortex_ring(t = 5*ti){
  scalar f[], omega[];
  f.refine = f.prolongation = fraction_refine;
  boundary((scalar*){u});
  vorticity(u, omega);
  boundary({omega});
  vertex scalar Phi[], vort[];
  foreach_vertex()
    vort[] = interpolate(omega, x, y);
  double iso[3] = {1, 2 , 4};
  for (int n = 0; n < 3; n++){
    foreach_vertex()
      Phi[] = vort[] - iso[n];
    fractions(Phi, f);
    boundary({f});
    char fname[99];
    sprintf(fname, "facets%d%g", j, iso[n]);
    FILE * fp = fopen (fname, "w");
    output_facets(f, fp);
    fclose(fp);
  }
}
/**
A few are plotted,

~~~gnuplot
set yr[0.6:1.9]
set xr[5.8:7.1] 
set xlabel 'z/R'
set ylabel 'r/R'
set size ratio -1
plot 'facets21' w l lw 2 lc 1  t 'ex. 2, Omega = 1' , \
'facets22' w l lw 3 lc 1  t 'ex. 2, Omega = 2' ,      \
'facets24' w l lw 4 lc 1  t 'ex. 2, Omega = 4' ,      \
'facets41' w l lw 2 lc 3  t 'ex. 4, Omega = 1' ,      \
'facets42' w l lw 3 lc 3  t 'ex. 4, Omega = 2' ,      \
'facets44' w l lw 4 lc 3  t 'ex. 4, Omega = 4' ,      \
'facets61' w l lw 2 lc 5  t 'ex. 6, Omega = 1' ,      \
'facets62' w l lw 3 lc 5  t 'ex. 6, Omega = 2' ,      \
'facets64' w l lw 4 lc 5  t 'ex. 6, Omega = 4'
~~~

Apart from a small (spatial and/ or temporal) shift, the structures
found between ex. #4 and #6 are arguably similar, wheareas ex #2 shows a
too-large structure. This could be due to excessive numerical
diffusion. So the 10 levels of refinement is a nesessity.
 */
