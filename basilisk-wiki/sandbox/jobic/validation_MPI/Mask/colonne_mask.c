
//#include "grid/multigrid.h"
#include "grid/quadtree.h"
#include "granular.h"
#if QUADTREE
#include "output_vtu_foreach.h"
#endif
// Domain extent : Lx/D_g
#define LDOMAIN 120.
#define maskr 0.75

// heap definition 2W is the size of the gate of the silo
double  H0,R0,D,W,tmax,Q,Wmin,DW;
//film size
double  Lz;
// Maximum refinement, ne pas oublier de changer W, la taille de l'orifice
#define LEVEL 8
///6 here to speed up for web site but 7 better
char s[80];
FILE * fpf,*fwq;

double grav=9.81, rho0=1.;


double redim_time(double t_adim) {
  return t_adim*sqrt(Dg/grav);
}

double redim_len(double x_adim) {
  return x_adim*Dg;
}

double adim_len(double x_dim) {
  return x_dim/Dg;
}

double adim_time(double t_dim) {
  return t_dim/sqrt(Dg/grav);
}

double adim_visco(double eta_dim) {
  return eta_dim/(Dg*sqrt(Dg*grav));
}

double redim_debit(double Q_adim) {
  return Q_adim*Dg*sqrt(Dg*grav);
}

void mytimer_print (FILE * pfile, timer t, int i, size_t tnc) {
  timing s = timer_timing (t, i, tnc, NULL);
  fprintf (pfile,
	   "\n# " GRIDNAME 
	   ", %d steps, %g CPU, %.4g real, %.3g points.step/s, %d var\n",
	   i, s.cpu, s.real, s.speed, (int) (datasize/sizeof(double)));
@if _MPI
  fprintf (pfile,
	   "# %d procs, MPI: min %.2g (%.2g%%) "
	   "avg %.2g (%.2g%%) max %.2g (%.2g%%)\n",
	   npe(),
	   s.min, 100.*s.min/s.real,
	   s.avg, 100.*s.avg/s.real,
	   s.max, 100.*s.max/s.real);
@endif
}

p[top] = dirichlet(0);
u.n[top] = neumann(0);
u.t[bottom] =  fabs(x-(maskr + 1)*LDOMAIN/2)<= W ? neumann(0):  dirichlet(0);
u.n[bottom] =  fabs(x-(maskr + 1)*LDOMAIN/2)<= W ? neumann(0):  dirichlet(0);
p[bottom]   =  fabs(x-(maskr + 1)*LDOMAIN/2)<= W ? dirichlet(0): neumann(0); 
u.n[right] = dirichlet(0);
u.n[left] = dirichlet(0);
u.t[right] = dirichlet(0);
u.t[left] = dirichlet(0);
f[left]= neumann(0);
/**
the three cases in the main */
int main() {
  char nomfichier[250];
  double dt_dim;
  
  L0 = LDOMAIN;
  // number of grid points
  N = 1 << LEVEL;
  // maximum timestep
  dt_dim=0.01;
  DT = dt_dim;
  TOLERANCE = 1e-4; 
  NITERMIN = 5;
// Initial conditions a=.5
  H0=0.9*adim_len(4.);
  R0=500.000;
  // size of the hole
  Wmin = 0.157;
  DW = LDOMAIN/N;
  
  const face vector g[] = {0.,-1.};
  a = g;
  
  fwq = fopen ("outWQ", "w");
  fclose(fwq);
  Lz = LDOMAIN;
 
  //mus=0.4;
  /* LEVEL 6 : 8 */
  /* LEVEL 7 : 16 */
  /* LEVEL 8 : 32 */
  /* LEVEL 9 : 64 */
  W = 8*DW;
  sprintf(nomfichier, "./results_N%1d_H0%g_W%f_dt%f.txt",N,redim_len(H0),redim_len(2*W),dt);
  fpf = fopen (nomfichier, "w");
  Q = 0; 
  tmax = adim_time(7.);
  run();
  mytimer_print (fpf,perf.gt, iter, perf.tnc);
  fclose (fpf);
}
/**
initial heap, a rectangle
*/
event init (t = 0) {
  mask (x < maskr * L0 ? left : none);

  scalar phi[];
  foreach_vertex()
    phi[] = H0 - y;
  fractions (phi, f);
/**
lithostatic pressure, with a zero pressure near the hole
to help 
*/
   foreach()
     p[] = (fabs(x-(maskr + 1)*LDOMAIN/2)<= W && fabs(y)<= adim_len(0.1)) ?  0 : max(H0 - y,0) ;    
  #if QUADTREE
  char name[80];
  sprintf(name, "./init_sol_at_%g_for_%g",t,redim_len(2*W));
  output_vtu((scalar *) {f,p,phi}, (vector *) {u}, name);
  #endif
}

/**
convergence outputs
*/
void mg_print (mgstats mg)
{
  if (mg.i > 0 && mg.resa > 0.)
    fprintf (stderr, "#   %d %g %g %g\n", mg.i, mg.resb, mg.resa,  
	  mg.resb > 0 ? exp (log (mg.resb/mg.resa)/mg.i) : 0); 
}
/**
convergence stats
*/
event logfile (i++) {
  stats s = statsf (f);
  fprintf (stderr, "%g %d %g %g %g %g\n", t, i, dt, s.sum, s.min, s.max - 1.);
  mg_print (mgp);
  mg_print (mgpf);
  mg_print (mgu);
  fflush (stderr);
}
/**
save some interfaces
*/
//event interface (t = {0, 1. , 2., 3., 4.}) {
 event interface ( t = 0; t += 1 ; t <= tmax) {
}
/**
Rate of flowing materials across the hole
 $$-\frac{dV}{dt} \text{ with } V = \int f dv $$
*/
event debit (t += 0.05 ) {
  static double Vold,V=1,Qinst=0;
  Vold=V; 
  V=0;
  foreach(reduction(+:V))
    V = V + f[]* Delta * Delta;
  Qinst = -(V-Vold)/.05;  
  if(Qinst > Q) Q = Qinst;
  if(t>=.1) {
    fprintf (stdout,"%lf %lf %lf %lf %lf %lf %lf %lf\n",t,V,redim_len(V/LDOMAIN),redim_len(2*W),Q,Qinst,redim_time(t),redim_debit(Qinst)); 
    fflush (stdout);
    fprintf (fpf,"%lf %lf %lf %lf %lf %lf %lf %lf\n",t,V,redim_len(V/LDOMAIN),redim_len(2*W),Q,Qinst,redim_time(t),redim_debit(Qinst)); 
    fflush (fpf);
  }
}  
/**
film output
*/
#if 0
event movie (t += 0.05) {
  static FILE * fp1 = popen ("ppm2mpeg > level.mpg", "w");
  scalar l[];
  foreach()
    l[] = level;
  output_ppm (l, fp1, min = 0, max = LEVEL, 
	      n = 2048, box = {{0,0},{Lz,Lz}});

  foreach()
    l[] = f[]*(1+sqrt(sq(u.x[]) + sq(u.y[])));
  boundary ({l});
  static FILE * fp2 = popen ("ppm2mpeg > velo.mpg", "w");
  output_ppm (l, fp2, min = 0, max = 2., linear = true, 
	      n = 2048, box = {{0,0},{Lz,Lz}});

  static FILE * fp3 = popen ("ppm2mpeg > pressure.mpg", "w");
  foreach()
    l[] = f[]*p[];
  output_ppm (l, fp3, min = 0, linear = true,
	      n = 2048, box = {{0,0},{Lz,Lz}});
}
#else


event pictures (t+=10) {
  scalar InG[];
  foreach()
     InG[] = f[]*In[];
  output_ppm (InG, file = "f.png", min=0, max = 2,  spread = 2, n = 256, linear = true,
  box = {{0,0},{Lz,Lz}});

  #if QUADTREE  
  char name[80];
  sprintf(name, "./sol_at_%g_for_%g",t,redim_len(2*W));
  output_vtu((scalar *) {f,p,InG,eta}, (vector *) {u}, name);
  #endif
}

#endif

#if 0
event movie (t += 0.05) {
    scalar l[];
    foreach()
    l[] = level;
    output_ppm (l, file = "level.mp4", min = 0, max = LEVEL,
                n = 2048, box = {{0,-1},{Lz,Lz}});
    
    foreach()
    l[] = f[]*(1+sqrt(sq(u.x[]) + sq(u.y[])));
    boundary ({l});
    output_ppm (l, file = "velo.mp4", min = 0, max = 2., linear = true,
                n = 2048, box = {{0,-1},{Lz,Lz}});
    foreach()
    l[] = f[]*p[];
    output_ppm (l, file = "pressure.mp4", min = 0,max = 2., linear = true,
                n = 2048, box = {{0,-1},{Lz,Lz}});
}/**

event pictures (t=3) {
    output_ppm (f, file = "f.png", min=0, max = 2,  spread = 2, n = 512, linear = true,
                box = {{0,0},{2,2}});
    output_ppm (p, file = "p.png", min=0, max = 1,  spread = 2, n = 512, linear = true,
                box = {{0,0},{2,2}});
}
  */
#endif
/**
If gfsview is installed on your system you can use this to visualise the simulation as it runs.
 If you manage to install `Bview`it is better !
*/
#if 0
event gfsview (i += 10) {
  static FILE * fp = popen ("gfsview2D -s  ", "w");
  output_gfs (fp, t = t);
}
#endif
/**
if "grid/multigrid.h" is not included, then this is the quadtree with possibility of adaptation
*/
#if QUADTREE
/*
event adapt (i++) {
  adapt_wavelet ({f,u.x,u.y}, (double[]){5e-3,0.001,0.001}, LEVEL, LEVEL-2,
     list = {p,u,pf,uf,g,f});
}
*/
#endif

