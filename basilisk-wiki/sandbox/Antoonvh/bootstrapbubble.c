/**
   ![Bootstrapping is not often encountered is reality. Image courtesy of Theodor Hosemann [via Wordsmith](http://wordsmith.org/words/bootstrap.html).](http://wordsmith.org/words/images/bootstrap.jpg)
   
   
#A bootstrapping bubble 
   
   In the example of a [bouncing droplet](bounce.c), some strange
   behaviour was observed regarding the advection of a droplet that is
   subject to surface tension. It appeared that for specific values of
   the dimensionless ratios, the kinetic energy of the system
   increases as the droplet advects trough the domain. This is
   unexpected as there was only supposed to be a viscous drainage of
   energy due to the finite Reynoldsnumber.
  
  
   
## Set-up 
   
   We follow the set-up of the aforementioned bouncing droplet
   example, except that on this page we are not interested in the
   bounce event itself and hence use periodic boundaries instead of
   the (hydrophobic) no-slip bottom.
   
   The dimensionless parameters are:
  
 $$ Re = \frac{\rho_a UR}{\mu_a} \approx 30, $$ 

 $$ We = \frac{\rho_w U^2R}{\sigma} \approx 0.2, $$

 $$ \Pi_1 = \frac{\mu_w }{\mu_a} \approx 100, $$

 $$ \Pi_1 = \frac{\rho_w }{\rho_a} \approx 1000 \rightarrow 100, $$
*/
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "tension.h"
#include "view.h"

int maxlevel = 8; //The presented results are converged-ish with
		  //respect to those obtained when using 9 levels of
		  //refinement (not shown)
int j;

int main(){
  mu2 = 1./30.; 
  mu1 = 100./30.;
  rho2 = 1.; 
  rho1 = 100;
  f.sigma = 500.;
  foreach_dimension()
    periodic(left);
  init_grid (1<<7);
  L0 = 10;
  X0 = Y0 = -L0/2;
  /**
## Three cases 
     
     For our analysis of the issue we run the set-up for three cases:
     First a (naive) default run is performed, second refinement is
     enfored arround the droplet's interface, and finally we re-run de
     default case with a decreased value for the surface tension
     parameter ($\sigma$). Run 1 to 3 are tagged with value $j=0$ to
     $2$, respectively.
  */
  j = 0;
  run();
  j = 1;
  run();
  j = 2;
  f.sigma=250;
  run();
}

event init (i = 0){
  refine (sq(x) + sq(y) < sq(1. + 0.5) && sq(x) + sq(y) > sq(1.-0.25) && level < maxlevel);
  fraction(f, sq(1.) - sq(x) - sq(y));
  foreach()
    u.y[] = -f[];
}
/**
## Force refinement on the interface 
   
   It is not well known (?) that the default (*smart*) adaptation
   rules tend to apply a coarser resolution when an interface is
   alligned with the grid (see [here](http://basilisk.fr/sandbox/Antoonvh/circle.c)). For a near circular bubble, I would *guess*
   that this makes no sense when evaluating the curvature
   (?). Therefore, during the second run (i.e. $j=1$), the algorithm
   is forced to adapt all interfacial cells to the maximum level of
   refinement, by using a dummy field (`ff`) that uses the default
   bilinear prolongation technique.
*/
event adapt (i++){
  if (j == 1){
    scalar ff[];
    foreach()
      ff[] = f[];
    boundary({ff});
    adapt_wavelet((scalar *){u,f,ff}, (double []){0.02, 0.02, 0.001, 0.001}, maxlevel);
  }  else {
    adapt_wavelet((scalar *){u,f}, (double []){0.02, 0.02, 0.001}, maxlevel);
  }
}
/**
## Output
   
   We monitor the evolution of the kinetic energy and write the result
   to seperate files for the different runs.
*/
event diagnose (i += 25){
  double e = 0;
  double e1 = 0;
  double e2 = 0;
  foreach(reduction(+:e) reduction(+:e1) reduction(+:e2)){
    e +=  0.5*sq(Delta)*(sq(u.x[]) + sq(u.y[])) * ((rho1*f[]) + (rho2*(1-f[])));
    e1 += 0.5*sq(Delta)*(sq(u.x[]) + sq(u.y[])) * ((rho1*f[]));
    e2 += 0.5*sq(Delta)*(sq(u.x[]) + sq(u.y[])) * (rho2*(1. - f[]));
  }
  char fname[100];
  sprintf(fname,"dataj=%d", j);
  static FILE * fp = fopen (fname,"w");
  fprintf (fp, "%g\t%g\t%g\t%g\n", t, e, e1, e2);
  fflush (fp);
}

/**
   For the all-important visual reference, some movies are rendered
   that display the bubble and the used computational mesh.
*/

event bviewer (t += 0.1; t <= 5){
  clear();
  view (fov = 20, width = 720, height = 720);
  squares ("f", min = 0, max = 1);
  cells();
  draw_vof ("f", lc = {1,0,1}, lw = 3.);
  char fname[99];
  sprintf (fname, "%d.mp4", j);
  save (fname);
}
/**
## Results
   
   One may view the difference in the gridding between the default
   case (top) and the radially refined case (bottom).
   
   ![Reference case](bootstrapbubble/0.mp4)
   ![Radially refined case](bootstrapbubble/1.mp4)
   
   The used strategy of enforcing refinment seems to have worked as
   intented. A plot of the evolution of the kinetic energy as was
   diagnosed for the three different runs is presented below:
   
   ~~~gnuplot
   set xlabel 'time'
   set ylabel 'Energy'
   set key box top left
   set size square
   plot 'dataj=0' u 1:2 w l lw 3 t 'Reference (default)' ,	\
   'dataj=1' u 1:2 w l lw 3 t 'Radially refined',		\
   'dataj=2' u 1:2 w l lw 3 t 'Reduced surface tension'
   ~~~
   
   From run 1 and 3 it appears that the surface-tension-force term is
   responsible for a spurrious source of kinetic energy. The ansantz
   is that this is due to improper reconstruction of the interface at
   the locations of the coarser level interfacial cells. These are not
   balanced arround the droplet because it does not have an exact
   circular shape and also the flow arround the droplet warrants an
   non-radially symmetric refinement. Fortunately this ansantz can be
   tested using the results from the second run. It appears from the
   plot that when refinement is enforced at the droplets interface,
   the spurrious source of energy is atleast much less prominent.
   
   I am *guessing* that this issue arrises from the fact that the
   wavelet-based error estimation for a volume-fraction field as
   defined in `two-phase.h`/`vof.h` is designed to only concern the
   errors for advection of an interface. However, when the additional
   module `tension.h` is included, the estimation of the curvature
   gives rise to an addional source of numerical errors. Without being
   burdened with any knowedge of the details, I expect that the latter
   procedure is different enough from evaluating the advection scheme
   such that the default volume-of-fluid-fiield wavelet-estimated
   error is not able to identify the challinging regions for curvature
   reconstruction correctly.
   
   Some addional analysis is presented [here](http://basilisk.fr/sandbox/Antoonvh/circle.c).
*/
