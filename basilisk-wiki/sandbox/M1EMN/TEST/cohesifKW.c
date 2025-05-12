/**

# Problem 
Kinetic Wave approximation of Shallow Water Savage Hutter description 
$$\frac{\partial h}{\partial t}+\frac{\partial Q}{\partial x}=0$$

with $Q$ function of ....



# Code 
*/

#include "grid/cartesian1D.h"
#include "run.h"


scalar h[];
scalar Y[];
scalar Q[];
scalar hp[],hpc[],nu[];
scalar dQ[]; 
double dt,SmmuS,B,tmax;
char s[80];

double  xf=0,xe=0;
double a=1;
double g=1;
double d=0.033;

int main() {
  L0 = 20;
  X0=0;
  SmmuS = 0.1;
  N = 512;
  DT = .00005;
  tmax =60;


 
  for (B = 0. ; B <= 0.15 ; B += 0.05){
  xf =0;  
  sprintf (s, "xKW-%.2f.txt", B);
  FILE * fp = fopen (s, "w"); 
  fclose(fp);
  sprintf (s, "shapeKW-%.2f.txt", B);
  FILE * fs = fopen (s, "w");  
  fclose(fs);
  run();
 }
    
}

h[left]=neumann(0);
Y[left]=neumann(0);
Q[left]=neumann(0);


event init (t = 0) {

  foreach(){
    h[] = (fabs(x)<.62) + 0.00001 ; // (0.55 -   .64 +)
    Q[]=0;
    }
  boundary ({h,Q});
  }

event printdata (t += 1; t <=   tmax ) {
 
  foreach(){
   xf = h[] > 1e-4 ?  max(xf,x) :  xf ;
   xe = h[] > 1e-4 ?  min(xe,x) :  xe ;
  }

  sprintf (s, "frontKW-%.2f.txt", B);
  FILE * f = fopen (s, "w");  
  foreach()
    fprintf (f, "%g %g %g \n", fmin((x-xf),0), h[], xe-xf);   
  fclose(f);
  double u01,h01;
  h01= interpolate(h,.1);
  u01=interpolate(Q,.1)/h01;
  fprintf (stderr, "%g %g %g %g %g \n", t, xf, xe,u01,h01);
  sprintf (s, "xKW-%.2f.txt", B);
  FILE * fp = fopen (s, "a");   
   fprintf (fp, "%g %g  \n", t, xf); 
  fclose(fp);
 
  sprintf (s, "shapeKW-%.2f.txt", B);
  fp = fopen (s, "a");
  foreach()
    fprintf (fp, "%g %g %g %g %g %g \n", x, h[], Q[], Y[], t, xf);
  fprintf (fp, "\n");
  fclose(fp);
}



double FQ(double h,double Y)
{
return (3*pow(h,2.5)/5 -pow(h,1.5)*(h-Y)/3 - 4*pow(h-Y,2.5)/15);
}

event integration (i++) {
  double dt = DT;
  dt = dtnext (dt);

  foreach()
    hp[] =  ( h[0,0] - h[-1,0] )/Delta;
  boundary ({hp});

  foreach()
    hpc[] =  ( h[1,0] - h[-1,0] )/2/Delta;
  boundary ({hpc});


  foreach()
   Y[] =  max( h[]*(1 - B/fabs(SmmuS -  hpc[] )) ,0);
  boundary ({Y});


  foreach()
     nu[] =   (FQ(h[-1,0],Y[-1,0]) + FQ(h[0,0],Y[0,0]))/2; 
  boundary ({nu});
    
  foreach()
    Q[] = - 2*sqrt(g)/(3*a*d)*nu[] * (-SmmuS + hp[]) ; //
  boundary ({Q});

  foreach()
    dQ[] =  ( Q[1,0] - Q[0,0] )/Delta;
  boundary ({dQ});

  foreach(){
    h[] +=  - dt*dQ[];
  }
  boundary ({h});  
}




/**
# Results 

Compare Kinetic Wave approximation in case of no cohesion with Savage Hutter 

we have to adjust a little bit the initial values and boundary condition to try to obtain an agremeent 

~~~gnuplot
  set xlabel "x-xf"
  set key left bottom
  p[-20:]'frontKW-0.00.txt' w l ,'../cohesifSH/frontSH-0.000.txt' w l

~~~

next step introduce cohesion and see the $Y$ height:



~~~gnuplot
  set xlabel "x-xf"
  set key left bottom
  p[0:15]'shapeKW-0.05.txt' w l ,''u 1:4 t'Y' w l

~~~


Compare front for different values of cohesion

~~~gnuplot
  set xlabel "x-xf"
  set key left bottom
  p[-20:]'frontKW-0.00.txt' w l,'frontKW-0.05.txt' w l,\
  'frontKW-0.10.txt' w l,'frontKW-0.15.txt' w l,

~~~

Compare velocity of the front for different values of cohesion

~~~gnuplot
 reset
 set key left bottom
 set xlabel "t"
 set ylabel "x f"
 p'xSH-0.000.txt',.2363*x+1.6 t'line',\
 '../cohesifKW/xKW-0.00.txt',\
 'xKW-0.05.txt','xKW-0.10.txt'
~~~




and check the factor two ratio .... 



# Links



* [http://basilisk.fr/sandbox/M1EMN/TEST/cohesifSH.c]()

*/


