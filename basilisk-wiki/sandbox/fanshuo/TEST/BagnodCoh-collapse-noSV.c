#include "grid/cartesian1D.h"
#include "run.h"

scalar h[];
scalar Y[];
scalar Q[];
scalar Fr[];
scalar hp[],hpc[],nu[],cu[];
scalar dQ[]; 

double theta,dt,S,tmax,inct=0.0001;
double I0, deltamu, tauc, mu0 ,G, dg;
char s[80];
/**
Main with definition of parameters
*/
int main() {
  L0 = 6.;
  X0 = -1;
	theta = 0.8;
	N = 128*4;
  DT = (L0/N)*(L0/N)/10 ;
  tmax =10; 
	I0 = 0.3;
	dg = 0.04;
	G = 1.;
	deltamu = 0.26;
	mu0 = 0.38;
	tauc = 0.05;
 /**
 a loop for the three values of Bingham parameter
 */ 
 for (S = 0.6 ; S <= 0.6 ; S += 0.1){ //  B = 1.25;  //.5 1.25 2
  sprintf (s, "Cohe-x-tauc%.3f-theta%.2f.txt", tauc,theta);
  FILE * fp = fopen (s, "w"); 
  fclose(fp);
  sprintf (s, "Cohe-shape-tauc%.3f-theta%.2f.txt", tauc,theta);
  FILE * fs = fopen (s, "w");  
  fclose(fs);
  inct = DT;
  run();
 }
}
/** 
initial elevation: a "rectangular column" of surface 0.5:
*/
event init (t = 0) {
 foreach(){
    h[] = 0.5*(-1*pow(x,2)+ 1)*(fabs(x)<1) ; 
    Q[]=0;
    }
	h[left] = neumann(0);
	h[right] = neumann(0);
	Q[left] = neumann(0);
	Q[right] = neumann(0);
  boundary ({h,Q});
  }
/** 
print data
*/
event output (t += 1; t <= tmax) {
	printf("tauc=%.2f,theta =%.3f, t=%f\n",tauc,theta,t);
//event output (t += inct; t < tmax) {   
  double  xf=0,xe=0;
  inct = inct*2;
  inct = min(1,inct);
/**
tracking the front and the end of the heap
*/   
  foreach(){
   xf = h[] > 1e-4 ?  max(xf,x) :  xf ;
   xe = h[] > 1e-4 ?  min(xe,x) :  xe ;
	Fr[] = h[]>0.3 ? Q[]/h[]/sqrt(h[]*cos(theta)):0;
  } 
 /**
save them
 */  
  sprintf (s, "Cohe-shape-tauc%.2f-theta%.2f.txt", tauc,theta);
  FILE * f = fopen (s, "a");  
	fclose(f);
  sprintf (s, "Cohe-x-tauc%.2f-theta%.2f.txt", tauc,theta);
  FILE * fp = fopen (s, "w");   
  foreach()
    fprintf (fp, "%g %g %g %g %g\n", fmin((x-xf),0), h[],Y[],Fr[] ,xe-xf);   
  fclose(fp);
}
/**
save the hight the flux and the yield surface as a function of time
*/ 
//event printdata (t = {0, 0.0625, 0.25, 1, 4, 100 }  ){
event printdata (t+=1;t<=tmax){

  sprintf (s, "Cohe-shape-tauc%.2f-theta%.2f.txt", tauc,theta);
  FILE * fp = fopen (s, "a");  
  foreach(){
		Fr[] = h[]>0.3 ? Q[]/h[]/sqrt(h[]*cos(theta)):0;
    fprintf (fp, "%g %g %g %g %g %g\n", x, h[], Q[], Y[],Fr[], t);}
  fprintf (fp, "\n\n");
  fclose(fp);
}
/** 
integration 


Definition of the flux 
*/

double FQ(double h,double Y)
{
	return I0*sqrt(cos(theta))/(deltamu*dg)*(6./15.*pow(h,2.5)+4./15.*pow((h-Y),2.5)-2./3.*(h-Y)*pow(h,1.5));
  //return 1./(6*mu)*(3*h - Y)*Y*Y ;
  //return (2*(2* pow(h - Y ,2.5) + pow(h,1.5)*(-2*h + 5*Y)))/15.;
  //return (2*( h* pow(h,1.5)))/5.;
  //
}
/**
definition of the derivative of the flux 

$Y=h- B/|S-\partial_xh|$ is such that $dY/dh$ is almost 1
 $$ c = \frac{d}{dh} (1/6*(3*h - Y(h))*(Y(h))^2)   \simeq 1./6*Y*Y* (3 - 1) +  1./3*(3*h - Y)*Y* 1$$
 */
double FQh(double h,double Y)
{
    //return 1./6*Y*Y* (3 - 1) +  1./3*(3*h - Y)*Y* 1;
		return I0*sqrt(cos(theta))/(deltamu*dg)*Y*sqrt(h); 

}



event integration (i++) {
  double dt = DT;
/**
finding the good next time step
*/
  dt = dtnext (dt);
/**
  $O(\Delta)$ down stream derivative
*/
  foreach(){
    hp[] =  ( h[0,0] - h[-1,0] )/Delta;
}
  boundary ({hp});
/**
Centered derivative for $h'$ used in the yield criteria for the yield surface $Y$
 
 Maybe not the best choice...
*/  
  foreach()
    hpc[] =  ( h[1,0] - h[-1,0] )/2/Delta;
  boundary ({hp});
/**
Yield surface
*/
  foreach(){
    Y[] =  max( h[] - tauc/fabs(sin(theta) -  hpc[]*cos(theta)-cos(theta)*mu0)  ,0);
}
  boundary ({Y});
/**
  the flux is taken with the mean value with the next cell, $\nu$ is an intermediate variable, a kind of viscosity (that is why we center)
*/   
  foreach(){
     nu[] =  (FQ(h[1,0],Y[1,0]) + FQ(h[0,0],Y[0,0]))/2; 
	}
  boundary ({nu}); 
/**
  To avoid oscillations in case of very large slope S the advective part of the flux is corrected with $-c (h_i -h_{i-1})$, where $c=\partial Q_c/\partial h$
*/  
  foreach()
    cu[] =  (tan(theta)-mu0)*(FQh(h[-1,0],Y[-1,0]) + FQh(h[0,0],Y[0,0]))/2;
    //cu[] =  0.;
  boundary ({cu});
  
  foreach(){
    Q[] = nu[]*(tan(theta)-hp[]-mu0)- cu[]* (h[0,0]-h[-1,0]);
		}
   //Q[] = -cos(theta)*nu[] * hp[] + nu[] * sin(theta) - cu[] * (h[0,0]-h[-1,0]) ; 
  //  Q[] = -cos(theta)*nu[] * hp[] + nu[] * tan(theta);
  boundary ({Q});  
/**
   derivative $O(\Delta)$ up stream like for heat equation
*/ 
  foreach(){
    dQ[] =  ( Q[1,0] - Q[0,0] )/Delta;
}
  boundary ({dQ});
/** 
update $h$  
*/
  foreach(){
    h[] +=  - dt*dQ[];
  }
  boundary ({h});  
}

