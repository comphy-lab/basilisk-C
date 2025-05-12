/**
# Another dipolar vortex-wall collision 

Kramer et al. (2007) studied the collision of a dipolar vortex with a
no-slip wall. They were inspired by the set-up of Orlandi (1990), but
For some reason however, Kramer et al. decided to used two Gaussian
vortices to model a dipolar vortex rather than the Lamb-Chaplygin
vortex, [that can also collides with walls on it's
trajectory](lamb-dipole.c).

## The novelty of this page

The reason we study this specific case is because Kramer et
al. present nice benchmark results in their study. Most notably, the
maximum enstrophy ($Z$) during the simulation run is of interest
here. In order to reproduce their results we do a convergence study
where we iteratively decrease the allowed tolerance on the wavelet
estimated error in the representaiton of the velocity component fields
($\zeta_u$). Eventough it is formally not required to limit the
maximum resolution, in practice it helps solver performance a lot to
do so. It is infact a mandatory argument for the `adapt_wavelet()`
function. However, it is not always clear what the maximum resolution
should be and from a user's perspective it is subjectively much more
natural to think in terms of an error threshold to reduce numerical
errors rather than directly in terms grid resolution (or is
it?).

Traditionally, static grid codes have been bamzoozling their users by
making it seem normal that users should figure out a suitable mesh
resolution for their simulations. Eventough it is very hard to do
consistently, Nobody really complains, and researchers and engineers
consider it part of their job to come up with some estimate,
typically based on experience and emperism.

On this page, a method to diagnose the solution and adjust the maximum
level of refinement accordingly is presented. The methods rely on the
wavelet estimated error in the representation of the solution field.

## Set-up 

The case is set-up after Kramer et al. (2017) with a Reynolds number ($Re$)
of $Re=1250$. This number is defined using the RMS velocity($U_{RMS}$)
of the initialized flow field and the domain half width ($W$, :S), $$Re
= \frac{U_{RMS}W}{\nu},$$

with $\nu$ the fluid's viscosity.
*/
#include "navier-stokes/centered.h"

#define  sq_dist ((sq(x-xp)+sq(y-yp))/sq(r_0))

double frac = 0.35;
double r_0 = 0.1;
double x_0 = 1.;
double y_0 = 0.;
double omg_0 = 301.94;
double init_refine_level = 11;
double ue;
double tend = 1;
double Re = 1250;
int maxlevel = 10;
int startlog = 5;
FILE * fpgif, * fplev;
FILE * fpm;
void add_gaussian_vortex(scalar s,double xp, double yp, double sign){
  foreach()
    s[] += sign * omg_0*(1-sq_dist)*exp(-sq_dist);      
}
/**
A second-order-accurate domain integrator function for the second order moment of a scalar field *s* is defined below. 
*/
double sndom(scalar s){
  double som=0;
  boundary({s});
  foreach(reduction(+:som)){
    double cv = 1.;
    foreach_dimension()
      cv*=Delta;
    som+=cv*sq(s[]);
    foreach_dimension()
      som+=cv*sq((s[1] - s[-1]))/48.;
  }
  return som;
}

u.t[bottom]=dirichlet(0.);
u.t[top]=dirichlet(0.);

int main(){
  periodic(left);
  init_grid(1<<5);
  L0=2;
  Y0=-1.;
  /**
  The case is run for $\zeta_u=[0.2, 0.1, 0.05]$. 
  */
  for (ue = 0.2;ue>0.04; ue/=2.){
    char fname[99];
    sprintf(fname, "ppm2mp4 omega%g.mp4",ue);
    fpgif = popen (fname, "w");
    sprintf(fname, "ppm2mp4 level%g.mp4",ue);
    fplev = popen (fname, "w");
    sprintf(fname,"log_ue=%g",ue);
    fpm = fopen(fname,"w");
    run();
    fclose(fpm);
  }
}

event init(i=0){
  CFL=0.7;
  u.x.refine=u.y.refine=refine_linear;
  TOLERANCE=10-5;
  scalar omega[], psi[], lev[];
  psi[bottom]=dirichlet(0.);
  refine((sq(x-x_0)+sq(y-y_0))<sq(0.8) && level < init_refine_level-2);
  refine((sq(x-x_0)+sq(y-y_0))<sq(0.6) && level < init_refine_level-1);
  refine((sq(x-x_0)+sq(y-y_0))<sq(0.4) && level < init_refine_level);
  add_gaussian_vortex(omega, x_0 + r_0, y_0, 1);
  add_gaussian_vortex(omega, x_0 - r_0, y_0, -1);
  boundary({omega});
  poisson(psi,omega);
  boundary({psi});
  double U=0;
  foreach(reduction(+:U)){
    u.x[]=(psi[0,-1]-psi[0,1])/(2*Delta);
    u.y[]=(psi[1]-psi[-1])/(2*Delta);
    U+=sq(Delta)*(sq(u.x[])+sq(u.y[])); 
  }
  U/=sq(L0);
  U=sqrt(U);
  double muz = U/Re;
  const face vector muc[]={muz,muz};
  mu=muc;
  boundary((scalar*){u,omega,psi});
}

event adapt(i++){
  adapt_wavelet((scalar*){u},(double[]){ue,ue},maxlevel);
}

event movie(t+=0.005; t<=tend){
  scalar omega[];
  scalar lev[];
  vorticity(u,omega);
  boundary({omega});
  foreach()
    lev[]=level;
  output_ppm(omega, fpgif, n = 512,min = -100, max = 100);
  output_ppm(lev, fplev, n = 512,min = 3, max = depth());
}

event log_event(i=startlog; t<=tend; i+=10){
  scalar omega[],chi[];
  vorticity(u,omega);
  double Z=sndom(omega);
  Z/=2.;
  double eZ=0;
  foreach_dimension(){
    wavelet(u.y, chi);
    foreach(reduction(+:eZ))
      eZ += sq(chi[]);
  }
  eZ/=2.;
  long int n = grid->tn;
  printf("%g\t%d\t%g\t%g\t%d\t%ld\n",t, i, Z, eZ, maxlevel, n);
  long int gc[20];
  memset(gc,0,sizeof(long int)*(20));
  foreach()
    gc[level]++;
  if (pid()== 0){
    if (i == startlog){
      fprintf(fpm,"t\ti\tZ\teZ\tml\tn\trtime");
      for (int j = 1 ; j<=19; j++)
	fprintf(fpm,"\tn%d",j);
      fprintf(fpm,"\n");
    }
    fprintf(fpm,"%g\t%d\t%g\t%g\t%d\t%ld\t%g",t, i, Z, eZ, maxlevel, n, perf.t);
    for (int j = 1 ; j<=19; j++)
      fprintf(fpm,"\t%ld",gc[j]);
    fprintf(fpm,"\n");
  }
  fflush(fpm);
}
/**
## Objectively determine a suitable maximum level of refinement
Since we are interested in the enstrophy ($Z$) that is defined as follows:

$$Z=\int\int \omega^2 \mathrm{d}A,$$
where,
$$\omega = \nabla \times \mathbf{u}. $$

When we follow the anzats that any cell-centered value of $\mathbf{u}$ (e.g. $u^i$) is accurately determined down to its estimated error:
$$u^i \equiv u^i \pm \chi_{u^i},$$
we can estimate an error ($\epsilon_Z$) in the
determination of the enstrophy based on how well the numerical solution
($\mathbf{u}$) is spatially discretized.

$$\epsilon_Z=\int\int \left(\frac{\chi_{u.x}}{\Delta}\right)^2 + \left(\frac{\chi_{u.y}}{\Delta}\right)^2 \mathrm{d}.A$$

If the grid cells at the maximum resolution contribute *relatively*
little to the $\epsilon_Z$ value, the maximumlevel can be reduced
without too much additional error, and avoid over refinement. Alternatively, if the grid cells at
the maximum level of refinement that have an estimated error ($\chi$)
larger than the error threshold ($\chi>\zeta_u$) contribute a lot to the
total error, the maximum level needs to increased to have a consistent
gridding, and prevent under refinement.

This means that if the algorithm decides to increase or decrease the maximum level of refinement, there will be sudden jumps in estimated error and number of grid cells, eventough the solution evolves smoothly. Notice that this is rooted in the fact that the resolution can only vary by a factor of two when using a tree grid, and is therefore consistent with the limitations of the used methods.
*/
event adapt_maxlevel(i=10; i++){
  scalar chi[];
  double eZ=0;
  double eZmlr = 0;
  double eZmlc = 0;
  double acc = 2.; // Accuracy order of the prolongation opperation
  double fac = sq(pow(acc,2.))/(pow(2.,(double)(dimension))); // The fraction the error at the maximum level will reduce with due to a possible increase in the maximum level.  
  foreach_dimension(){
    wavelet(u.y, chi);
    foreach(reduction(+:eZ) reduction(+:eZmlr) reduction(+:eZmlc)){
      eZ += sq(chi[]);
      if (level == maxlevel){
	eZmlc += sq(chi[]);
	if(chi[]>ue)
	  eZmlr += sq(chi[]);
      }
    }
  }
  /**
  In some exotic(?) senario it could be possible that both conditions are statisfied. Then, the maxlevel is increased.
  */
  if (eZmlr > (frac * eZ))
    maxlevel++;
  else if (eZmlc < (frac * eZ /(fac*3./2.) ))
    maxlevel--;
}


/**
## Results

![The vorticity dynamics during the collision for $\zeta_u=0.05$](kramer/omega0.05.mp4)

~~~gnuplot Enstrophy
set xlabel 'time'
set ylabel 'Z'
plot 'log_ue=0.2' u 1:3 w l lw 2 t 'ue = 0.2' ,\
     'log_ue=0.1' u 1:3 w l lw 2 t  'ue = 0.1' ,\
     'log_ue=0.05'u 1:3 w l lw 2 t  'ue = 0.05'#,\
     #'log_ue=0.025'u 1:3 w l lw 2 t  'ue = 0.025'
~~~

![The grid resolution $\zeta_u=0.05$](kramer/level0.05.mp4)

~~~gnuplot Maximum level
set yr [7.5 : 12.5]
set xlabel 'time'
set ylabel 'Maximum Level'
plot 'log_ue=0.2' u 1:5 w l lw 2 t 'ue = 0.2' ,\
     'log_ue=0.1' u 1:5 w l lw 2 t  'ue = 0.1',\
     'log_ue=0.05'u 1:5 w l lw 2 t  'ue = 0.05'#,\
     #'log_ue=0.025'u 1:5 w l lw 2 t  'ue = 0.025'
~~~

~~~gnuplot time
set yr [10 : 2500]
set logscale y 2
set xlabel 'Physical time'
set ylabel 'wall-clock time'
set key top left
plot 'log_ue=0.2' u 1:7 w l lw 2 t 'ue = 0.2' ,\
     'log_ue=0.1' u 1:7 w l lw 2 t  'ue = 0.1',\
     'log_ue=0.05'u 1:7 w l lw 2 t  'ue = 0.05'#,\
     #'log_ue=0.025'u 1:5 w l lw 2 t  'ue = 0.025'
~~~

The cluster is about 4 times slower tham my laptop. Therefore we cannot show here that the $\zeta_u=0.05$ results can be considered to be converged against the results obtained from running with $\zeta_u=0.025$.
*/
