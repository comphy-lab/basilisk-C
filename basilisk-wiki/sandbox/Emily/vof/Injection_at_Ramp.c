/**
# Wave generation by fluidised flow. 

This is developed from Lily Battershill's code for the Navier-Stokes VOF case of the [Bougouin, 2020](#references) fluidized granular flow tsunami generation experiment except just doing 2 phase flow rather than 3 phase. The rotated version was giving a lot of trouble so this is a version that is not rotated

This particular version is set up to implement a PELE experiment from Gert Lube and Goeffrey Roberts - see O:\MAU25501\Working\Modelling\PELE

Angle of slope 5deg  0.087266 Radians
Water depth 0.5 m
Length of flume 8.4m

08/10/2025

*/
#include "grid/quadtree.h"
#include "navier-stokes/centered.h"

#define THETA0 0.087266
#define DEPTH0 0.5 //depth of water
#define FLUME_LENGTH 8.4
#define NINPUT 11 // size of input file

#define MAXLEVEL 12

#include "two-phase.h"
#include "waveprobes.h"
//#include "reduced.h" //Uncomment this line if using G vector for gravity
#include "tag.h"
#include "utils.h"

double ue = 1e-3,H0=0.,UF=0.;
FILE * fp;
//# Hacked to force with starting input
double in_time[NINPUT] = { 0.,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5 };
double in_h[NINPUT] = { 0.,0.03996,0.01015,0.03172,0.01522,0.04504,0.03616,0.03806,0.02918,0.03679,0. };
double in_vel[NINPUT] = { 0.,1.92,1.935516153,1.822921891,1.579404068,1.160774194,1.120440025,1.084493185,0.901428016,0.578008094,0. };
double hin0=0.,hin1=0.,tin0=0.,tin1=0.,velin0=0.,velin1=0.,tinBC=0.,hinBC=0.,velinBC=0.;
int inI=0;

/** Define parameters */
#define io (0.01) // Width of initial jet
#define ramp (DEPTH0/sin(THETA0)) // 0.5 IS THE WATER DEPTH (M) (Length of the ramp)
#define runout (2.67)            // Length of edge in solid bottom of the flume
#define MUR 54                   // Water/Air mu ratio?
#define imagerate 0.01
#define endlim 8.		// Time limit? 60 seconds?
#define TRIANGLE (-tan(THETA0)*(x+io)) // Solid bottom of the flume

int main() {
  L0 = FLUME_LENGTH;       // Edge length
  origin (0, -DEPTH0);
  rho1 = 998; //water
  rho2 = 1.27; //air
  mu1 = 1.002e-2; //water
  mu2 = 1.002e-2/MUR; //air
  init_grid (1 << 7);
  fp=fopen("reporting.txt", "w");
  (const) face vector av[] = {0, -9.81, 0.};
  a = av;
  //G.y=-9.81;  //Uncomment this line if using G vector for gravity and comment two lines above
  CFL = 0.3;
  DT=0.005;
  run();
}


/** Set boundary conditions */

// time varying BC 
u.n[left] = dirichlet((y < H0) && (y >= 0.) ? UF*cos(THETA0) : 0 );
u.t[left] = dirichlet((y < H0) && (y >= 0.) ? -UF*sin(THETA0) : 0 );
f[left] = dirichlet((y < H0) ? 1 : 0 );


/** Initilaise domain */
event init0 (i = 0)
{

  refine((((y < 0.1) && (y > - 0.1)) ||((y < 0.1 - tan(THETA0)*(x+io)) && (y > - 0.1 - tan(THETA0)*(x+io)))) && level < MAXLEVEL);          //refining around water surface and ramp
  scalar f1[],f2[];
  foreach_vertex(){
    f1[] = io - x;
    f2[] = H0 - y;
    f2[] = min(f1[],f2[]);
    f1[] = - y;
    f1[] = max(f1[],f2[]);
  }
  fractions (f1, f);
  foreach(){
    foreach_dimension(){
      u.x[]=0.;
    }
  }
  boundary({f, u});

}


event updateBC(i++){
  if (i==0){
    hin0=in_h[0];
    tin0=in_time[0];
    velin0=in_vel[0];
    hin1=in_h[inI+1];
    tin1=in_time[inI+1];
    velin1=in_vel[inI+1];
    inI++;
    tinBC=(t-tin0)/(tin1-tin0);
    H0=hin0+tinBC*hin1;
    UF=velin0+tinBC*velin1;
    
    fprintf(stdout,"BC: %g %d %g %g %g %g %g %g\n",t,inI,tin0,tin1,hin0,hin1,velin0,velin1);
    fflush(stdout);
  }
  if ( t >= tin1 && inI < NINPUT ) {
    hin0=hin1;
    tin0=tin1;
    velin0=velin1;
    hin1=in_h[inI+1];
    tin1=in_time[inI+1];
    velin1=in_vel[inI+1];
    inI++;
    fprintf(stdout,"%g %d %g %g %g %g %g %g\n",t,inI,tin0,tin1,hin0,hin1,velin0,velin1);
    fflush(stdout);
  }
  else if (inI == NINPUT){
    hin0=0.;
    hin1=0.;
    velin0=0.;
    velin1=0.;
    tin1=1.e10;
  }
  tinBC=(t-tin0)/(tin1-tin0);
  H0=hin0+tinBC*hin1;
  UF=velin0+tinBC*velin1;
}
/** Slope implementation */
#include "fracface.h"
scalar tri[];
event make_slope(i++) {
  fraction (tri, y <= TRIANGLE);    //Bottom of flume (solid)
  foreach()
    foreach_dimension()
    u.x[] -= u.x[]*tri[];
  face vector tri_face[];
  face_fraction (tri, tri_face);
  boundary((scalar *){tri_face});
  foreach_face()
    uf.x[] -= tri_face.x[]*uf.x[];
  boundary ((scalar *){u,uf});
}

/** Log files */
#include "view.h"
event logfile (i++){
  fprintf (stderr, "%d %g %d %d\n", i, t, mgp.i, mgu.i);
}
event dump(i+=1000){
  dump();
}


event pictures_t ( t+=imagerate ; t<=endlim)
{
  char three[80];
  sprintf (three, "f-%g.asc",t);
  FILE * fp1 = fopen (three, "w");
  foreach()
    {
      fprintf(fp1,"%g %g %g %g %g\n",x,y,f[],u.x[],u.y[]);
    }
  fclose (fp1);
  view (fov = 5, tx = -.5, ty = -0.,
	bg = {1,1,1},
	width = 800, height = 200);
  draw_vof ("tri", filled = 1, fc = {-1,1,1});
  draw_vof ("f", filled = -1, fc = {1,1,1});
  squares ("u.x", min = -0.2, max = 0.2, linear = true);
  isoline ("u.x", 0., lc = {1,1,1}, lw = 2);
  sprintf (three, "u.x.%g.png",t);
  save (three);
  
  view (fov = 5, tx = -.5, ty = 0.,
	bg = {1,1,1},
	width = 800, height = 200);
  draw_vof ("tri", filled = 1, fc = {-1,1,1});
  draw_vof ("f", filled = -1, fc = {1,1,1});
  squares ("u.y", min = -0.1, max = 0.1, linear = true);
  isoline ("u.y", 0., lc = {1,1,1}, lw = 2);
  sprintf (three, "u.y.%g.png",t);
  save (three);
}
/** Output files */
int time_step = 0.;
double shift_control = 1.;
event output_files  (t += imagerate; t < endlim) {
  // Images (to postprocess in python)
  char three[80];
  // Velocity profiles
  sprintf (three, "vprof-%d",time_step);
  FILE * fp1 = fopen (three, "w");
  //double step = (L0)/(1 << (lmax));
  for (double y = 0.; y <= 0.06; y += 8/(pow(2,MAXLEVEL))){
    fprintf (fp1, "%g %g %g %g\n", y,interpolate (u.x, 0.02, y),interpolate (u.y, 0.02, y),interpolate(f, 0.02,y));} //Outputs y, u.x, u.y, f
  fclose (fp1);
  time_step += 1.;
}


/** Output VOF facets */
int add = 0;
event interface (t+= 0.004) {
  char fname[99];
  sprintf(fname, "water_elevation%d.txt", add);
  FILE * fp = fopen(fname, "w");
  output_facets(f,fp);
  fclose(fp);
  add += 1.;  
}


/** Adaption */
event adapt (i++)
{
  //Remove small droplets
  remove_droplets(f,-5);
  adapt_wavelet((scalar *){u,f,tri},(double[]){ue,ue,0.01,0.01},MAXLEVEL-1);
}

event reporting(i+=10){
  double umaxA=0.,umaxW=0.,sumf=0.,xlocA=0.,ylocA=0.,xlocW=0.,ylocW=0.;
  foreach(serial){
    if ( (1.-f[])*(u.x[] * u.x[] + u.y[] * u.y[]) > umaxA * umaxA ) {
      umaxA=(1.-f[])*sqrt( u.x[] * u.x[] + u.y[] * u.y[] );
      xlocA=x;
      ylocA=y;
    }
    if ( f[]*(u.x[] * u.x[] + u.y[] * u.y[]) > umaxW * umaxW ) {
      umaxW=f[]*sqrt( u.x[] * u.x[] + u.y[] * u.y[] );
      xlocW=x;
      ylocW=y;
    }
    sumf += f[] * Delta * Delta;
  }
  fprintf(fp,"%g %g %g %g %g %g %g %g %g\n",t,dt,sumf,umaxA,xlocA,ylocA,umaxW,xlocW,ylocW);
  fflush(fp);
}
event ending(t=endlim){
  fclose(fp);
}

