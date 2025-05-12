
/**
    This file contains routines used for the python interface. most of
    these routines can't be in the .i file because they need to be
    compiled with qcc.
*/


#include "grid/multigrid.h"
#include "pg.h"

scalar * dbl = NULL; 
scalar * field = NULL; 
scalar * bl_lt = NULL; 

scalar psibt_lt[]; 

/**
   Continuation parameters and fields used by the bifurcation
   algorithm
 */

int contpar = 0;
double forcing_mag = 1.0;
scalar b_surf0[]; 

/**
   dtconv is the pseudo convection time scale. should be small
   compared to the horizontal diffusion time scale but not too small
   to avoid convergence issues.
 */ 

double dtconv = 1e-1; 


double d = 0.05;
#define g(x)    ( 1-exp(-sq(x)/(2*sq(d))) + 0.1)
#define gp(x)   ( (x)/sq(d)*exp(-sq(x)/(2*sq(d))))

double h (double x, double y){
//  return ( g(x-0.0)*g(1.0-x)*g(y-0.0)*g(1.0-y) );
  return ( 1.0);
}

double hx (double x, double y){
  return ( 0.0);
}

double hy (double x, double y){
  return ( 0.0);
}

double k (double x, double y, double s){
  return (1e-2);
}

/* double h (double x, double y){ */
/*   return ( g(x-0.0)*g(1.0-x)*g(y-0.0)*g(1.0-y) ); */
/* } */

/* double hx (double x, double y){ */
/*   return ( (gp(x-0.0)*g(1.0-x)-g(x-0.0)*gp(1.0-x))*g(y-0.0)*g(1.0-y)); */
/* } */

/* double hy (double x, double y){ */
/*   return ( g(x-0.0)*g(1.0-x)*(gp(y-0.0)*g(1.0-y)-g(y-0.0)*gp(1.0-y))); */
/* } */


/**
   wind stress and wind stress derivative
*/

double tau0 = 1e-2;

/* double taux   (double x, double y){ return (tau0*sin(2*(y-ys)*pi));} */
/* double taux_y (double x, double y){ return (2*pi*tau0*cos(2*(y-ys)*pi));} */
double taux   (double x, double y){ return (0.);}
double taux_y (double x, double y){ return (0.);}
double tauy   (double x, double y){ return (0.);}
double tauy_x (double x, double y){ return (0.);}


/**
    Explicit vertical diffusion routine for the bifuraction solver */

trace
void vdiff_explicit  (scalar * bl, scalar * dbl) {
  vertbc(bl);
  double vflux[nl+1];
  foreach(){
    vflux[0] = 0.; // no flux bottom
    for (int l = 1; l < nl+1 ; l++) {    
      scalar b0 = bl[l];
      scalar b1 = bl[l+1];
      scalar kfs = kfsl[l];
      scalar db = dbl[l];

      vflux[l] = kfs[]*(1+sq(a)*sq(sf[l])*(sq(hpc.x[])+sq(hpc.y[])))/hc[]
        *(b1[] - b0[])/ds;
      db[] += (vflux[l] - vflux[l-1])/(ds*hc[]);
    }
  }
}
/**
    Since we use a static restoring for the convection algorithm,
    there is no explicit tendency. we compute an artificial tendency
    with a time scale dtconv. */


void convection_tend(scalar * bl, scalar * dbl, double dt)
{

  foreach() {
    for (int l = 1; l < nl+1 ; l++) {		
      scalar b = bl[l];
      scalar b_sav = pl[l-1];
      b_sav[] = 1.0*b[];
    }
  }

  convection(bl);

  foreach() {
    for (int l = 1; l < nl+1 ; l++) {
      scalar b = bl[l];
      scalar b_sav = pl[l-1];
      scalar db = dbl[l];
      db[] += (b[]-b_sav[])/dt;
      b[] = b_sav[];
    }
  }
}

void convection_tend_lt(scalar * bl, scalar * dbl, scalar *bl_lt, double dt)
{

  foreach() {
    for (int l = 1; l < nl+1 ; l++) {		
      scalar b = bl[l];
      scalar b_sav = pl[l-1];
      b_sav[] = 1.0*b[];
    }
  }
      
  convection(bl);

  foreach() {
    for (int l = 1; l < nl+1 ; l++) {		
      scalar b = bl[l];
      scalar b_lt = bl_lt[l];
      scalar b_sav = pl[l-1];
      scalar db = dbl[l];
      db[] -= (b[] + b_lt[])/dt;
      b[] = b_sav[] + b_lt[];
    }
  }

  convection(bl);

  foreach() {
    for (int l = 1; l < nl+1 ; l++) {		
      scalar b = bl[l];
      scalar b_sav = pl[l-1];
      scalar db = dbl[l];
      db[] += (b[])/dt;
      b[] = b_sav[];
    }
  }

}

void adjust_kv(scalar *bl){

  foreach(){
    for (int l = 2; l < nl ; l++) {
      scalar b0 = bl[l-1];
      scalar b1 = bl[l];
      scalar kfs = kfsl[l-1];

      if (b1[] - b0[] < 0)
	kfs[] = 100*k(x,y,sf[l-1]);
      else
	kfs[] = k(x,y,sf[l-1]);
    }
  }
}

void forcing(scalar * bl, scalar *dbl)
{

  scalar b = bl[nl];
  scalar db = dbl[nl];
  foreach() 
    db[] += (b_surf[] - b[])/tau_surf;
}

void pyset_vars()
{
  set_vars();


/**
   constants
*/
  a = 1.; 
  r = 0.02; 
  tau_surf = 5e-2;
  ys = 0.2;
  origin (0.0, ys);

  assert (dbl    == NULL);
  assert (bl_lt  == NULL);



  for (int l = 0; l < nl+2; l++) {
    scalar db = new scalar;
    dbl = list_append (dbl, db);
    scalar b_lt = new scalar;
    bl_lt = list_append (bl_lt, b_lt);
  }

  foreach() {
    for (scalar db in dbl)
      db[] = 0.0;
  }


  psibt_lt[right]  = dirichlet(0);
  psibt_lt[left]   = dirichlet(0);
  psibt_lt[top]    = dirichlet(0);
  psibt_lt[bottom] = dirichlet(0);


  foreach() {
//    b_surf[] = 1e-3*y;
    b_surf[] = 2*cos(pi*(y-ys));
    b_surf0[] = b_surf[];
}


  // write topo
  char name[80];
  sprintf (name, "topo.dat");
  FILE * fp = fopen (name, "w");
  foreach()
    fprintf (fp, "%g %g %g\n", x, y, hc[]);
  fflush(fp);
  fclose(fp);

  boundary (all);

}

/**
   Python interface routines (should be in .i file but I need foreach)
 */

void pyset_field (int ifield, double * val1, int len1){
  if (ifield == 1)      // 1: buoyancy
    field = bl;
  else if (ifield == 2) //2: buoyancy tendancy
    field = dbl;
  else if (ifield == 3) //3: buoyancy linear tangent
    field = bl_lt;
  else
    printf("pyset_field, ifield = %i not supported\n",ifield);
  
  int i = 0;
  for (int l = 0; l < nl; l++) {
    scalar b = field[l+1];
    foreach() {
      b[] =  val1[i];
      i++;
    }
  }
  boundary(field);
  vertbc(field);
} 

void pyget_field (int ifield, double * val2, int len2){
  if (ifield == 1)      // 1: buoyancy
    field = bl;
  else if (ifield == 2) //2: buoyancy tendancy
    field = dbl;
  else if (ifield == 3) //3: buoyancy linear tangent
    field = bl_lt;
  else
    printf("pyset_field, ifield = %i not supported\n",ifield);
  
  int i = 0;
  for (int l = 0; l < nl; l++) {
    scalar b = field[l+1];
    foreach() {
      val2[i] = b[];
      i++;
    }
  }
} 

void pystep ( double * val1, int len1,
              double * val2, int len2){
  
  foreach()
    for (scalar s in dbl)
      s[] = 0.;

  int ifield = 1;
  pyset_field ( ifield, val1, len1);

  double dtloc = 1.0;
//  adjust_kv(bl);
  dtloc = velocity(bl, ul, wl, psibt, dtloc );
  advection(bl, ul, wl, dbl);
  hdiffusion(bl, dbl);
  forcing(bl, dbl);
  vdiff_explicit(bl, dbl);
  convection_tend(bl, dbl, dtconv);

  ifield = 2;  
  pyget_field ( ifield, val2, len2);
} 

void pystep_lt ( double * val1, int len1,
		 double * val2, int len2,
		 double * val3, int len3){
  
  int ifield = 1;
  pyset_field ( ifield, val1, len1);
  ifield = 3;
  pyset_field ( ifield, val3, len3);

  foreach()
    for (scalar s in dbl)
      s[] = 0.;

 // adjust_kv(bl);

  double dtloc = 1.0;
  dtloc = velocity(bl, ul, wl, psibt, dtloc );
  advection(bl_lt, ul, wl, dbl);

  scalar b_surf_sav[];
  foreach() {
    b_surf_sav[] = b_surf[];
    b_surf[] = 0.;
  }

  dtloc = velocity(bl_lt, ul, wl, psibt_lt, dtloc );
  advection(bl, ul, wl, dbl);
  hdiffusion(bl_lt, dbl);
  forcing(bl_lt, dbl);
  vdiff_explicit(bl_lt, dbl);
  convection_tend_lt(bl, dbl, bl_lt, dtconv);

  ifield = 2;  
  pyget_field ( ifield, val2, len2);

  foreach() 
    b_surf[] = b_surf_sav[];

} 

/**
   Set and adjust the continuation parameter via the python interface
 */
void pyset_contpar(int pycontpar) { 
  contpar = pycontpar; 
}

void pyadjust_contpar(double contpar_val){
  if (contpar == 1) {
    forcing_mag = contpar_val;
    foreach() 
      b_surf[] = forcing_mag*b_surf0[];
  }
}

void pyinit_vertgrid(int pynl){ nl = pynl;}

void pytrash_vars(){
  trash_vars();

  free (bl_lt), bl_lt = NULL;
  free (dbl), dbl = NULL;
}
