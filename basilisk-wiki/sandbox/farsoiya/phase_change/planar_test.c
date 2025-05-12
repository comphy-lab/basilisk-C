/**
# Phase change planar test
*/
#include "centered.h"
#include "two-phase.h"
#include "navier-stokes/conserving.h"
#include "tension.h"
#include "reduced.h"
#include "henry2.h"

scalar c[], * stracers = {c};
double bubble_radius = 1.;
double box_size = 5;
double conc_liq1 = 0, conc_gas1 = 1.;

scalar mdot[];

double mwt = 1.; //Molecular weight

int MAXLEVEL;

int main (int argc, char **argv)
{
  size (box_size);
	
  MAXLEVEL = 8;
  N = 1 << MAXLEVEL;


  rho1 = 1.;
  rho2 = 1.;
	f.sigma = 0;
  c.alpha = 1;
  TOLERANCE = 1e-4;
	//~ DT = 1.e-4;
  //this is constrained because of pure single component bubble
  conc_gas1 = rho2/mwt ; 
  c.D1 = 1.;
  c.D2 = 1.;

    run();
  
}

event vof (i++)
{
	
  scalar ff[];
  foreach(){
    ff[] = f[];
  }
  foreach(){

    //~ if (interfacial(point, f)) {
    if (f[] > 1e-6 && f[] <1. - 1e-6) {

	
      coord n = interface_normal (point, f);  
      double alpha = plane_alpha (f[], n);
  
		
      double deltaD = mdot[]*dt*sqrt(sq(n.x)+sq(n.y));	
      //~ double deltaD = mdot2[]*dt;	
	
      ff[] = plane_volume (n, alpha + deltaD);

      //~ ff[] = f[] + deltaD;


      //~ if ( ff[] >= 1.){
      if ( alpha + deltaD >= 0.5){
	
		
	if (fabs(n.y) >= fabs(n.x)){
	  if (n.y <= 0.){
	    if (f[0,-1] <= 1.e-6){
	      ff[0, -1] =  plane_volume (n,   - 1. + ( alpha + deltaD) );
	      //~ ff[0, -1] = ff[] - 1.;
						
	    }
	    else if (n.x > 0 ){
	      //~ ff[1, -1] = ff[] - 1.; // plane_volume (n,   0.5 - (alpha + mdot2[]*dt - 0.5) );
	      ff[1, -1] = plane_volume (n,   - 1. + ( alpha + deltaD) );
	    }
	    else {
	      //~ ff[-1, -1] = ff[] - 1.; //plane_volume (n,   0.5 - (alpha + mdot2[]*dt - 0.5) );
	      ff[-1, -1] = plane_volume (n,   - 1. + ( alpha + deltaD) );
						
	    }
	  }
	  else{
	    if (f[0,1] < 1.e-6){
	      //~ ff[0, 1] = ff[] - 1.; //plane_volume (n,   0.5 - (alpha + mdot2[]*dt - 0.5) );
	      ff[0, 1] = plane_volume (n,  - 1. + ( alpha + deltaD) );
	    }
	    else if (n.x > 0 ){
	      //~ ff[1,1] = ff[] - 1.; //plane_volume (n,   0.5 - (alpha + mdot2[]*dt - 0.5) );
	      ff[1,1] = plane_volume (n, -  1. + ( alpha + deltaD) );
	      //~ if (i==22){
	      //~ double jk = ff[];
	      //~ double jk1 = deltaD;
	      //~ double jk2 = f[1,1];
	      //~ printf("some");
	      //~ }
	    }
	    else {
	      //~ ff[-1, 1] = ff[] - 1.; //plane_volume (n,   0.5 - (alpha + mdot2[]*dt - 0.5) );
						
	      ff[-1, 1] = plane_volume (n,  - 1. + ( alpha + deltaD) );
	    }
	  }
	  //~ cc[] = 1.;
	  //~ else {
	  //~ if (deltaD<0)
	  //~ cc[0, 1] = plane_volume (n, - 0.5 + alpha - deltaD*sqrt(sq(n.x)+sq(n.y)) );
	  //~ else
	  //~ cc[0, -1] = plane_volume (n,  0.5 - alpha - deltaD*sqrt(sq(n.x)+sq(n.y)) );
	}
	else {
					
	  if (n.x <= 0.){
	    if (f[-1,0] < 1.e-6){
	      //~ ff[-1, 0] = ff[] - 1.; //plane_volume (n,   0.5 - (alpha + mdot2[]*dt - 0.5) );
	      ff[-1, 0] = plane_volume (n, -  1. + ( alpha + deltaD) );
	    }
	    else if (n.y > 0 ){
	      //~ ff[-1, 1] = ff[] - 1.; //plane_volume (n,   0.5 - (alpha + mdot2[]*dt - 0.5) );
	      ff[-1, 1] = plane_volume (n, -  1. + ( alpha + deltaD) );
	    }
	    else {
	      //~ ff[-1, -1] = ff[] - 1.; //plane_volume (n,   0.5 - (alpha + mdot2[]*dt - 0.5) );
						
	      ff[-1, -1] = plane_volume (n, -  1. + ( alpha + deltaD) );
	      //~ printf("\n deb = %lf %lf %lf %lf %lf %lf %lf", alpha, n.x, n.y, ff[], mdot2[]*dt, sqrt(sq(n.x)+sq(n.y)), ff[-1,-1]);
	    }
	  }
	  else{
	    if (f[1,0] < 1.e-6){
	      //~ ff[1, 0] = ff[] - 1.; //plane_volume (n,   0.5 - (alpha + mdot2[]*dt - 0.5) );
	      ff[1, 0] = plane_volume (n, -  1. + ( alpha + deltaD) );
	    }
	    else if (n.y > 0 ){
	      //~ ff[1, 1] = ff[] - 1.; //plane_volume (n,   0.5 - (alpha + mdot2[]*dt - 0.5) );
	      ff[1, 1] = plane_volume (n,  - 1. + ( alpha + deltaD) );
						
	    }
	    else {
	      //~ ff[1, -1] = ff[] - 1.; //plane_volume (n,   0.5 - (alpha + mdot2[]*dt - 0.5) );
	      ff[1, -1] = plane_volume (n,  - 1. + ( alpha + deltaD) );
	    }
	  }
	}
				

      }
	
    }
  }
  boundary({ff});
	
  foreach(){
    f[] = ff[];
  }
  boundary({f});
	
	
}


#define co ((coord){1., cm[]})
event stability (i++) {
	
    
  
  foreach(){
    mdot[] = 0;
    if (interfacial(point, f)) {
      //~ if (f[] > 1e-6 && f[] <1. - 1e-6) {
    
			
      foreach_dimension(){

	double kc = co.x*c.D1*c.D2/(c.D1*(1. - f[]) + f[]*c.D2);
				
	//	double normal_dir = ( f[1] - f[-1] ) < 0 ? -1. : 1. ;
				
	//	int flux_dir = ( c[1] - c[-1] ) <= 0 ? 0 : 1 ;
				
	// Guess mdot worked fine with static planar test and agrees interface with moving planar test 
	mdot[]  +=    -  ( kc*(c[1] - c[])/(Delta/2. + f[]*Delta) )/Delta * mwt/0.5; 
				
	// central difference mdot, didn't agreed well for planar test
	//~ mdot[]  +=    - ( f[1] - f[-1] )/Delta/2. * kc*(c[1] - c[-1])/Delta/2. * mwt/0.5; 
				
	// as given in Maes 2020, doesn't work well, division by f[] makes mdot very large for less filled cells
	//~ mdot[]  +=    -   kc*(c[1] - c[-1])/Delta/2. * mwt * ( f[1] - f[-1] )/Delta/2./f[] ; 
				

      }
      //~ if (i!=0)
      //~ mdot[]  =      (sqrt(c.D1/M_PI/t) )/Delta*mwt; //working ( f[1] - f[-1] )/Delta/2. *
				
      
      mdot[] *= mwt/rho2 ; 
		if (mdot[] < 0.){
	printf("\nerror -ve mdot, %lf %lf %lf",x, y, mdot[]);
	exit(0);
      }	
      // time step restriction the volfrac change per cell should not exceed 0.5
      if (fabs(mdot[])>0.){
	double dt = 0.5/mdot[];
	if (dt < dtmax)
	  dtmax = dt;
      }
    }
		
		
  }
	
}


event init (t = 0)
{

  //~ fraction (f,  -(sq(bubble_radius) - sq(x - box_size*0.5) - sq(y- box_size*0.5)));
  //~ fraction (f,  -(sq(bubble_radius) - sq(x - box_size*0.5) - sq(y)));
  //~ fraction (f,   (2.4 + x - y));
  fraction (f,  - (2.49 - y));

  foreach()
    c[] = conc_liq1*f[] + conc_gas1*(1. - f[]);
	
  boundary({c});
}

event loger(i++){
	
  printf("\n%d %lf %lf", i, t, dt);
	
  scalar pos[];
  position (f, pos, {0,1});
  
  scalar derv[];
  foreach(){
    derv[] = 0;
	
    derv[] =  - mwt/rho2 * ( c.D1*(c[0,1] - c[])/(Delta/2. + f[]*Delta)  );
    //~ derv[] +=  - mwt/rho2 * ( c.D1*(c[1] - c[-1])/Delta/2.  );

			
	
  }
  

  char name[80];
  sprintf (name, "wave-%d", N);
  static FILE * fp = fopen (name, "w");
  if(i>10)
    fprintf (fp, "%g %g %g %g %g %g %g\n", t, 2.49 - 2.*conc_gas1*sqrt(t*c.D1/M_PI), statsf(pos).max,
	     statsf(derv).max, mwt/rho2*sqrt(c.D1/M_PI/t), interpolate(c, 2.5, statsf(pos).max + 0.5 ), erfc(0.5/sqrt(c.D1*t)/2. ) );
  fflush (fp);
  
 
}	

event end(t = 0.2){
	
  printf("\n analytical depth = %lf\n",2.49 - 2.*conc_gas1*sqrt(t*c.D1/M_PI));
	
  foreach(){
    if (f[] > 1e-6 && f[] < 1. - 1e-6) {
      printf("\n numerical depth = %lf\n",y); fflush(stdout);
      exit(0);
    }
  }
}


/**
# Plots

<center>
<table>
<tr>
<td>![](planar_test/grad_moving.png){ width="400px" }</td>
<td>![](planar_test/conc_planar_moving.png){ width="400px" }</td>
<td>![](planar_test/interface_moving.png){ width="400px" }</td>
</tr>
<tr>
<td><center>Gradient at interface</center></td> 
<td><center>Concentration at x -x0 = 0.5</center></td> 
<td><center>Interface position</center></td> 
</tr>
</table>
</center>

~~~pythonplot

import numpy as np
import matplotlib.pyplot as plt
	
plt.figure()


t128, analytical128, num128, mdot128, anmdot128,cnum128, can128 = np.loadtxt('wave-128',delimiter=' ',unpack=True)
t256, analytical256, num256, mdot256, anmdot256,cnum256, can256 = np.loadtxt('wave-256',delimiter=' ',unpack=True)

plt.plot(t256[1:-1:1], anmdot256[1:-1:1],'k+',label='Analytical')

plt.plot(t128, mdot128,'k--',label='Basilisk L7')

plt.plot(t256, mdot256,'k',label='Basilisk L8')


plt.xlabel(r'$t$',fontsize=20)
plt.ylabel(r'$D\nabla c$',rotation=0,fontsize=20)
plt.tight_layout()
plt.legend()
plt.savefig('grad_moving.png')

plt.figure()

plt.plot(t128, can128,'k+',label='Analytical')

plt.plot(t128, cnum128,'k--',label='Basilisk L7')

plt.plot(t256, cnum256,'k',label='Basilisk L8')
plt.legend()

plt.xlabel(r'$t$',fontsize=20)
plt.ylabel(r'$c(x-x0 = 0.5)$',rotation=90,fontsize=20)
plt.tight_layout()

plt.savefig('conc_planar_moving.png')

plt.figure()

plt.plot(t256, analytical256,'k+',label='Analytical')

plt.plot(t128, num128,'k--',label='Basilisk L7')

plt.plot(t256, num256,'k',label='Basilisk L8')
plt.legend()
plt.xlabel(r'$t$',fontsize=20)
plt.ylabel(r'$x0$',rotation=90,fontsize=20)
plt.tight_layout()

plt.savefig('interface_moving.png')

plt.figure()
plt.savefig(' ')


~~~
*/
