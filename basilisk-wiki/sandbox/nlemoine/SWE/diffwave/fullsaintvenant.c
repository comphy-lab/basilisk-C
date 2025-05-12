#include "terrain.h"
#include "saint-venant.h"
#include "utils.h"
#include "output.h"
#include "nlemoine/SWE/manning-tilt.h"

#define ETABF 10.6
#define n_ch 0.03
#define n_fp 0.06

#define sec_per_day 86400.
#define M_PI acos(-1.0)

#define ETAE     1e-2 // error on water depth (1 cm)
#define HE       1e-2 // error on water depth (1 cm)
#define UE       1e-2 // 0.01 m/s

#define LEVEL 7
#define MINLEVEL 5

double stage_hydrograph (double tt)
{
   return( ETABF-2.0*cos(2.0*M_PI*tt/sec_per_day) );
}

double segment_Q (coord segment[2], scalar h, vector u)  // vient de "layered/hydro.h"
{
  coord m = {segment[0].y - segment[1].y, segment[1].x - segment[0].x};
  normalize (&m);
  double tot = 0.;

  foreach_segment (segment, p) {
    double dl = 0.;
    foreach_dimension() {
      double dp = (p[1].x - p[0].x)*Delta/Delta_x*(fm.y[] + fm.y[0,1])/2.;
      dl += sq(dp);
    }
    dl = sqrt (dl);    
    for (int i = 0; i < 2; i++) {
      coord a = p[i];
      tot += dl/2.*
	interpolate_linear (point, (struct _interpolate)
			    {h, a.x, a.y, 0.})*
	(m.x*interpolate_linear (point, (struct _interpolate)
				 {u.x, a.x, a.y, 0.}) +
	 m.y*interpolate_linear (point, (struct _interpolate)
				 {u.y, a.x, a.y, 0.}));
    }
  }
  // reduction
#if _MPI
  MPI_Allreduce (MPI_IN_PLACE, &tot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif
  return tot;
}

double segment_eta (coord segment[2], scalar h, scalar eta)
{
  coord m = {segment[0].y - segment[1].y, segment[1].x - segment[0].x};
  normalize (&m);
  double tot = 0., wtot = 0.;

  foreach_segment (segment, p) {
    double dl = 0.;
    foreach_dimension() {
      double dp = (p[1].x - p[0].x)*Delta/Delta_x*(fm.y[] + fm.y[0,1])/2.;
      dl += sq(dp);
    }
    dl = sqrt (dl);    
    for (int i = 0; i < 2; i++) {
      coord a = p[i];
      tot += dl/2.*
	interpolate_linear (point, (struct _interpolate)
			    {h, a.x, a.y, 0.})*
	interpolate_linear (point, (struct _interpolate)
				 {eta, a.x, a.y, 0.}) ;

      wtot += dl/2.* interpolate_linear (point, (struct _interpolate)
			    {h, a.x, a.y, 0.});
    }
  }
  // reduction
#if _MPI
  MPI_Allreduce (MPI_IN_PLACE, &tot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce (MPI_IN_PLACE, &wtot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif
  return (tot/wtot);
}

FILE * fpcontrol;

int main (int argc, char * argv[])
{

  size (2048.);
  origin (-424.,-2332.);
  init_grid (1 << LEVEL);

  // Get .kdt file
  system("wget https://dropsu.sorbonne-universite.fr/s/sT5GkqPqzsaMYab/download");
  system("mv download garonne_local.kdt");
  // Get .pts file
  system("wget https://dropsu.sorbonne-universite.fr/s/5TWiFJaAXGKcb6L/download");
  system("mv download garonne_local.pts");
  // Get .sum file
  system("wget https://dropsu.sorbonne-universite.fr/s/sz5qF6fc7MwkPzz/download");
  system("mv download garonne_local.sum");
  
  G = 9.81;
  tilt.x = 0.;
  tilt.y = 1.103e-3;

  run();

  if(pid()==0.)
   fclose(fpcontrol);

}

event initialize (i=0)
{
    terrain(zb,"garonne_local", NULL);
    scalar zmin = zb.dmin;
  
    output_ppm (zb, min = 3, max = 20, file = "topo6.png",
		n = 512, linear = true);
    
   foreach()
   {
     double etaini = stage_hydrograph(0.); 
     eta[] = zb[]>etaini ? zb[] : etaini;
     h[] = eta[] > zb[] ? eta[]-zb[] : 0.;
     eta[] = zb[] + h[];
     u.x[] = 0.;
     u.y[] = 0.;
     double zzmin = zmin[]<nodata ? zmin[] : zb[];
     nmanning[] = zzmin<ETABF ? n_ch : n_fp ;
   } 

  if(pid()==0)
   fpcontrol = fopen("outlet_info.txt","w");

}

event update_BC(i++)
{
    // Conditions aux limites
  u.n[top] = neumann(0.);
//  u.t[top] = dirichlet(0.);
  u.t[top] = neumann(0.);
  h[top] = neumann(0.);
  zb[top] = neumann(0.);
  eta[top] = neumann(0.);
  nmanning[top] = neumann(0.);

  u.n[bottom] = neumann(0.);
  u.t[bottom] = dirichlet(0.);
  zb[bottom] = neumann(0.);
  nmanning[bottom] = neumann(0.);

  foreach_boundary(bottom)
  { 
    eta[] = fmax(zb[],stage_hydrograph(t));
    h[] = eta[]-zb[];
  }

  eta[bottom] = neumann(0.);
  h[bottom] = neumann(0.);
 
  boundary(all);
}

scalar l[];

event snapshot (t+=300. ; t<=86400.)
{
  scalar m[], etam[];
  foreach() {
    m[] = eta[]*(h[] > dry) - zb[];
    etam[] = h[] < dry ? 0. : (zb[] > 0. ? eta[] - zb[] : eta[]);
  }
  output_ppm (h, mask = m, min = 0, max = 10, 
	      n = 2048, linear = true, file = "h.mp4");
}

event info_conservation (t+=30. ; t<=86400.)
{ 
    coord Section[2];
    Section[0] = (coord) { X0 , Y0+L0-L0/N };
    Section[1] = (coord) { X0+L0 , Y0+L0-L0/N };
  
    // "Measure" discharge at top (outlet) section
    double Q = segment_Q (Section,h,u);
    
    // "Measure" mean elevation at top (outlet) section
    double etam = segment_eta (Section,h,eta);
  
    // Total water volume in domain  
    double Vol = 0.;
    foreach()
      Vol += h[] * Delta * Delta ;
  
    if(pid()==0){
      fprintf(fpcontrol,"%g %g %g %g\n",t/3600.,Q,etam,Vol);
      fflush(fpcontrol);
    }
}

// Adaptivity

scalar absu[];

int adapt() {
#if TREE
  scalar eta[];
  foreach()
    eta[] = h[] > dry ? h[] + zb[] : 0;
  boundary ({eta});

  foreach()
    absu[]= h[] > dry ? sqrt(u.y[]*u.y[]+u.x[]*u.x[]) : 0.;

  astats s = adapt_wavelet ({h,eta,absu}, (double[]){HE,ETAE,UE},
			    LEVEL, MINLEVEL);

  scalar zmin = zb.dmin;
  foreach(){
    double zzmin = zmin[]<nodata ? zmin[] : zb[];
    nmanning[] = zzmin<ETABF ? n_ch : n_fp;
    l[] = level;    
  }

  boundary(all);
  
//  fprintf (ferr, "# refined %d cells, coarsened %d cells\n", s.nf, s.nc);
  return s.nf;

#else // Cartesian
  return 0;
#endif
}

event do_adapt (i++) adapt();

/**
# Results 

![Simulated water depth (flow direction is from bottom to top)](fullsaintvenant/h.mp4)(width=50% )
![outlet_info.txt](fullsaintvenant/outlet_info.txt)
*/