/**
##  Copy of axisymmetric toroidal bubble field
The same code as the [toroidal](http://basilisk.fr/sandbox/yonghui/vortex/toroidalbubble.c)
code ,
we want to output a small region near the biggest bubble 
and then "copy" it to 
[another](http://basilisk.fr/sandbox/yonghui/vortex/planarvort.c)
2D simulation.
*/
#include "axi.h"  
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "tension.h"
#include "tag.h"  
#include "view.h"
/**
Parameters
*/
#define LVM 8        // max refine level 
#define thick 0.015625 //
#define R0 0.125       // solid nozzle's inner layer radius
#define Re 5000.       // Dimensionless Re = Re_sim * R_0
#define We 100.        // Dimensionless We = We_sim * R_0
#define mur 0.02       //ratio viscosity water-air
#define tramp 0.1*R0   //start time 
#define deltag sqrt(R0/Re)*R0 //vorticity layer thickness based on Re
#define ainj 4.6*R0    // x position of bubble's right interface 
#define xinj 4.8*R0    // x position of nozzle's right interface 
// Wb: initial bubble size, take a try with 0.6, 1.2, 2.4, 4, etc) 
#define Wb 4.*R0  
#define tend 3.0  
FILE * fp1, * fp2, * fp3;
scalar omega[];
/**
## BC
*/
u.n[top] = dirichlet(0.0);
u.n[right] = neumann(0.0);
p[right]   = dirichlet(0.0);

u.n[left] = y <= R0 ? dirichlet(erf(t/tramp)*erf((-y+R0)/(deltag * sqrt(t/R0 + 1.e-6) ) )) : dirichlet(0.0);
u.t[left] = dirichlet(0.);
f[left] = dirichlet(0.);

bid solid;
u.n[solid] = dirichlet(0.0);
u.t[solid] = dirichlet(0.0);
/**
## main 
*/
int main() {
  size(2.);
  init_grid (64);
  //origin(-1./2.,0);
  rho1 = 1.e-3, mu1 = mur*R0/Re;
  rho2 = 1., mu2 = R0/Re, f.sigma = R0/We;
  TOLERANCE = 1e-5;
  //DT = 0.001;
  fp1 = fopen("bub1","w");
  fp2 = fopen("bub2","w");
  fp3 = fopen("bub3","w");
  run();
}
/**
## init
*/
event init (t = 0) {
  refine ( y < 0.6  && x < 1.5 && level < LVM);
  mask ( (x < xinj && y > R0 && y + 0.5*x< R0 + 0.5*xinj && y < R0 + 2.*thick) ? solid :none);
  fraction (f, (x < ainj && x > ainj - Wb && y < R0*1.01) ? 1 : 0  );
}

#if 1 //ELIMINATE_AIR
event bubblecount( i += 20){
  scalar m[];
  foreach(){
    m[] = f[] > 1e-4;}
  int n = tag (m);
  double v[n];     //volume /(2 PI) !!!!!
  double sft[n];   //surface
  double cnc[n];   //cell number count
  coord b[n];
  for (int j = 0; j < n; j++){
    v[j] = b[j].x = b[j].y =  sft[j]= cnc[j] = 0.;
  }

  foreach_leaf(){
    if (m[] > 0) {
      int j = m[] - 1;
      v[j] += dv()*f[];
      sft[j] += sq(Delta)*f[];
      cnc[j] += f[];
      coord p = {x,y,z};
      foreach_dimension()
        b[j].x += dv()*f[]*p.x;
    }
  }
  // less than $5 \%$ initial surface
  double sfmax = 0.05*R0*Wb;		


  for (int j = 0; j < n; j++) {
    fprintf (fp1, "%d %g %d %g %g %g %g\n", i, t, j, v[j],sft[j], b[j].x,b[j].y );

    fflush (fp1);  

    double cellv=0.;
    foreach(){ 
      if (m[] > 0) {
        int j = m[]-1;
        if(sft[j] < sfmax){
          f[] = 0.;
          fprintf (stdout, "## CLEAN %d %g %g %g %g %g\n", i, t, x, y, sft[j], cnc[j]);
          cellv += f[]*sq(Delta);
        }
      }
    }
    fflush (stdout);  
  }

}
#endif

#if 1 //ELIMINATE_WATER
event dropcount( i += 20){
  scalar mm[];
  foreach(){
    mm[] = f[] < 0.99;}
  int n = tag (mm);
  double wsft[n]; //surface
  coord b[n];
  for (int j = 0; j < n; j++){
    b[j].x = b[j].y =  wsft[j] = 0.;
  }

  foreach_leaf(){
    if (mm[] > 0) {
      int j = mm[] - 1;
      wsft[j] += sq(Delta)*(1.-f[]);
      coord p = {x,y,z};
      foreach_dimension()
        b[j].x += dv()*(1.-f[])*p.x;
    }
  }
  // less than $5 \%$ initial surface
  double sfmax = 0.05*R0*Wb;		

  foreach(){ 
    if (mm[] > 0) {
      int j = mm[]-1;
      if(wsft[j] < sfmax){
        f[] = 1.;
      }
    }
  }
  fflush (stdout);  
}
#endif

/**
In this event we calculate the center of the biggest bubble and saved the field nearby.
*/

event locbubble (t += 0.1){
  scalar m[];
  foreach(){
    m[] = f[] > 1e-4;}
  int n = tag (m);
  double sft[n]; //surface
  double uxc[n]; //surface
  double uyc[n]; //surface
  coord b[n];
  for (int j = 0; j < n; j++){
    b[j].x = b[j].y =uxc[j]=uyc[j]=sft[j]= 0.;
  }

  foreach_leaf(){
    if (m[] > 0) {
      int j = m[] - 1;
      sft[j] += sq(Delta)*f[];
      uxc[j] += .q(Delta)*u.x[];
      uyc[j] += sq(Delta)*u.y[];
      coord p = {x,y,z};
      foreach_dimension()
        b[j].x += sq(Delta)*f[]*p.x;
    }
  }

  int nmax=0;
  if (n > 0){
    for (int j = 0; j < n; j++){
      nmax = sft[nmax] <= sft[j] ? j : nmax;
    }}

  double bubcx= b[nmax].x/sft[nmax];
  double bubcy= b[nmax].y/sft[nmax];
  double uxx = uxc[nmax];
  double uyy = uyc[nmax];
  fprintf (fp3, "%d %g %g %g %g %g \n", i, t, bubcx,bubcy,uxx,uyy);
  fflush (fp3);  

  char datavor[2000]; //time
  sprintf(datavor, "data_t%03g", t*100);
  FILE * fpf3 = fopen (datavor, "w");
  /**
We want to locate the cell (rcx,rcy) closest to the bubble surface center (bubcx,bubcy).*/
  double rcx = bubcx;
  double rcy = bubcy;
  double rrr = 0.05;
  foreach(){
    if (sq(y-bubcy)+sq(x-bubcx) <= sq(rrr) ){
      rcx = x;
      rcy = y;
      rrr = sqrt( sq(y-bubcy)+sq(x-bubcx) );
    }
  }

  fprintf(stderr,"#FOUND %g %g %g %g %g",rrr,rcx,bubcx,rcy,bubcy);
  /**
We saved the field nearby if the distance $d \leq 0.5$, 
We moved the origin towrds the centered cell!!
*/
  foreach(){
    if (sq(y-rcy) + sq(x-rcx) < sq(0.5)) {
      fprintf(fpf3,"%g %g %g %g %g %g\n",x-rcx+Delta/2.,y-rcy+Delta/2., \
              omega[],f[],u.x[]-uxx,u.y[]-uyy);
    }
  }

  char picname[800];
  sprintf (picname, "Re%gWe%grb%02gt%03g.png",Re,We,Wb*800, t*100);
  view (fov = 6.,tx = -0.5, bg = {1,1,1}, width = 2048, height = 600, samples = 1);
  clear();
  squares("f",max=1.,min=0.);
  draw_vof("f", lw = 1.5);
  mirror (n = {0., 1.}, alpha = 0.) {
    squares("omega",max=30.,min=-30.);
    draw_vof("f", lw = 3.);
  }
  save (picname);
}
