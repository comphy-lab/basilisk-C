/**
# Tsunami runup onto a complex three-dimensional beach

This example is a classical validation test case for tsunami models. It was proposed at the ["Third international workshop on long-wave runup models"](http://isec.nacse.org/workshop/2004_cornell/bmark2.html). It is based on experimental data obtained in a wave tank in order to understand the extreme runups observed near the village of Monai during the 1993 Okushiri tsunami.*/

#include "terrain.h"
#include "utils.h"
//#define SAINT_VENANT 1
#if SAINT_VENANT
# include "saint-venant.h"
#else
# include "green-naghdi.h"
#endif
#include "output.h"


#define MAXLEVEL 9
#define MINLEVEL 6
/**
Need to import input water elevation on the left boundary which can be downloaded from 
[here](http://isec.nacse.org/workshop/2004_cornell/data/Benchmark_2_input.txt)
*/
// time and water elevation arrays, number of points
double *tp=NULL; 
double *qp=NULL; 
int pt_num = 0; 

void read_hydrograph (const char * name)
{
  FILE * fp;
  if ((fp = fopen (name , "r")) == NULL) {
    fprintf (stderr,"cannot open hydrograph data file.\n name = %s ",name);
    assert (false);
  }
  pt_num = 0;
  if (fscanf(fp, "%d", &pt_num) == EOF) {
      fclose (fp);
      return;
  }
  tp = malloc(pt_num*sizeof(double));
  qp = malloc(pt_num*sizeof(double));
  
  double time,q;
  for(int i=0; i<pt_num; ++i)
  {
      if (fscanf (fp, "%lf \t %lf \n", &time, &q) == EOF) break;
      tp[i] = time; qp[i] = q;
  }
  fclose (fp);
}


// Plain linear interpolation. I'm sure this can be better.
double interpolate_input(double _t)
{
    if(!tp) {read_hydrograph("../Benchmark_2_input.txt");}
    int idx = 0;
    double time = tp[idx];
    while( time<_t )
    {
        idx++;
        time = tp[idx];
    }
    if(idx > 0){idx--;}
    double res;
    if( idx < (pt_num-1) )
    { res = qp[idx] + (qp[idx+1]-qp[idx])*(_t-tp[idx])/(tp[idx+1]-tp[idx]); }
    else res = qp[pt_num-1];

    return res;
}

/**
Initialization
*/
int main()
{
    dry = 1e-04;
    G = 9.81;
    N = 1 << MAXLEVEL;
    size(5.448);
    origin(0,0);
    init_grid (1 << MINLEVEL);
    run();
}

u.n[top] = dirichlet(0);
u.t[top] = dirichlet(0);
u.n[bottom] = dirichlet(0);
u.t[bottom] = dirichlet(0);
u.n[right] = dirichlet(0);
u.t[right] = dirichlet(0);


/**
Need to impose correct boundary conditions on left side. 
*/
u.t[left] = 0.;
h[left] = dirichlet(interpolate_input(t)-zb[]);
u.n[left] = -2*(sqrt (G*max(h[],0.)) - sqrt(G*max((interpolate_input(t)) - zb[], 0.)))+u.x[];

/**
Define gauges and print to files every time step
*/
Gauge gauges[] = {
    {"myinput", 1e-3, 1.7, "Left Boundary"},
    {"myp5", 4.521, 1.196, "P5"},
    {"myp7", 4.521, 1.696, "P7"},
    {"myp9", 4.521, 2.196, "P9"},
    {NULL}
};

event gauges1 (i++) output_gauges (gauges, {eta,h,zb});


/**
Top boundary is at 3.402. Terrain is generated with 
xyz2kdt -v output_name < input_file where the input_file is downloaded from [here](http://isec.nacse.org/workshop/2004_cornell/data/Benchmark_2_Bathymetry.txt).
*/

event init (t = 0) {

    mask (y > 3.402 ? top : -1);
    terrain (zb, "../new_terrain", NULL);

    conserve_elevation();
    foreach() {
      zb [] = x > 5.448 ? 0.13535 : zb[];
      h[] = max(0., - zb[]);
    }
    boundary ({h});
}

event init(i++)
{
    foreach() {
        zb [] = x > 5.448 ? 0.13535 : zb[];
    }
    // testing input hydrograph interpolation
    static FILE* fpinput = fopen ("input_wl.txt", "w");
    fprintf (fpinput, "%g %g\n", t,
             interpolate_input (t)
            );
}

/**
Quadratic friction
*/
event friction(i++) {

  foreach() {
        u.x[] = h[]>dry ?  u.x[]/(1+dt*1e-03*norm(u)/h[]) : 0.;
        u.y[] = h[]>dry ?  u.y[]/(1+dt*1e-03*norm(u)/h[]) : 0.;
        h[]   = h[]>dry ?  h[] : 0.;
    }
}


#if QUADTREE
event adapt (i++) {
  scalar eta[];
  foreach()
    eta[] = h[] > dry ? h[] + zb[] : 0;
  boundary ({eta});

  astats s = adapt_wavelet ({eta}, (double[]){2e-4}, MAXLEVEL, MINLEVEL);
  fprintf (stderr, "# refined %d cells, coarsened %d cells\n", s.nf, s.nc);
}
#endif

/**
Export for gfsview - just for fun
*/
void write_to_gfs(char* fname)
{
    FILE* fp = fopen (fname, "w");
    output_gfs (fp, {zb,h});
    fclose(fp);
}

event outputfile (t += 2)
{
    char fname[50];
    sprintf(fname,"grid-%f.gfs",t);
    write_to_gfs(fname);
}
/**
Need to end simulation at some point and clean up the memory.
*/
event end (t = 22) {
  printf ("i = %d t = %g\n", i, t);
  if(tp) {free(tp);}
  if(qp) {free(qp);}
}


