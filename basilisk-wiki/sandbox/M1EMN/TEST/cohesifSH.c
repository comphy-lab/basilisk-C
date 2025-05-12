/**


solve Savage Hutter front with cohesion



 
 ![animation of the front](cohesifSH/animate.gif)
 
*/

#include "grid/multigrid1D.h"
#include "saint-venant.h"
#define pourtoutes(item, array) \
for(int keep = 1, \
count = 0,\
size = sizeof (array) / sizeof *(array); \
keep && count != size; \
keep = !keep, count++) \
for(item = (array) + count; keep; keep = !keep)

double a=1;
double d=0.033;
double  xf=0,xe=0;


u.n[left] = dirichlet(.232); //.28
h[left] = neumann(0);
u.n[right] = neumann(0);
h[right] = neumann(0);

double Ltas,Htas,tmax,xfront,tauY;
char s[80];
int main()
{
  X0 = 0.;
  L0 = 30;
  G = 1.; 
  N = 1024;
  Htas = .5;
  Ltas = .5;
  tmax=80;
  DT=0.002;
  FILE * f = fopen ("tauYxmaxSH.txt", "w");
  sprintf (s, "shapeSH-%.3f.txt", tauY);
  FILE * fs = fopen (s, "w");
  fclose(fs);

  //double values[] ={ 0, 0.01, 0.02, 0.05,  0.1 , 0.2 , 0.3, 0.5 , 1 };
  double values[] = {0.0, 0.005, 0.010 };
  pourtoutes(double *v, values) {
  tauY=*v;
    xf=0;  
  printf("# tauY %lf\n",tauY);
    
  run();
  fprintf (f, "%g %g  \n", tauY, xfront );
  }
}

event init (i = 0){
  double h0;
   h0=Htas;
   foreach(){
    zb[] = -0.4*x;
    h[] = (x < Ltas) ? h0 : 0;
    u.x[]= 0;}
}



event coulomb_friction (i++) {
  double mu=0,U=0;
  foreach() {
    U=norm(u);
    /** ici on rajoute  mu=.3 + a*d*U/pow(h[],1.5);  */
    mu=0.3 + 5./2*a*d*U/pow(h[],1.5);
      if(U>0){
          foreach_dimension()
          u.x[] = max(U -dt *( mu * G + tauY/h[]),0)*u.x[]/U;}
  }
  boundary ({u.x});
}




event outputfront (t += 1 ) {
    sprintf (s, "xSH-%.3f.txt", tauY);
    static FILE * ff = fopen(s, "w");

    foreach(){
        xf = h[] > 1e-20 ?  max(xf,x) :  xf ;
        xe = h[] > 1e-20 ?  min(xe,x) :  xe ;
    }
     fprintf (ff, "%g %g %g \n", t, xf , xe);
     xfront = xf;
}


event output (t += .2; t < tmax) {
    static FILE * fp = popen ("gnuplot -persist 2> /dev/null", "w");
#ifdef gnuX
    
#else
    if(t==0) fprintf (fp,"set term gif animate;set output 'animate.gif';set size ratio .5\n");
#endif
    fprintf (fp,"set title ' collapse tauY=%.3lf--- t= %.2lf '\n"
             "p[0:10][-.5:2]  '-' u 1:($2) t'free surface' w lp lt 3,0 not w l linec -1\n",tauY,t);
  foreach()
      fprintf (fp, "%g %g \n", x, h[]);
  fprintf (fp,"e\n\n");
}



event outputlog (t += 1; t < tmax) {
   foreach()
     fprintf (stderr, "%g %g %g \n", x, h[], t);
   fprintf (stderr, "\n");
   char s[80];
   sprintf (s, "shapeSH-%.3f.txt", tauY);
   FILE * fp = fopen (s, "a");
    foreach()
      fprintf (fp, "%g %g %g %g %g \n",x,h[],u.x[],zb[],t);
      fprintf (fp, "\n");
    fclose(fp);
    
    sprintf (s, "frontSH-%.3f.txt", tauY);
    FILE * f = fopen (s, "w");
    foreach()
      fprintf (f, "%g %g %g \n", fmin((x-xf),0), h[],u.x[]);
    fclose(f);
}




/**
 
 
 vérification de la relation entre U et h
 $$U=\frac{2(\alpha-\mu_s) h^{3/2}}{5 a d}$$
 
~~~gnuplot
 set key left bottom
 set xlabel "x"
 p[-15:][:]'frontSH-0.000.txt' u 1:3 w l,'' u 1:(.1*($2**1.5)/(2.5*.033)) w l
~~~
 
 The front shouw be superposed with Louise analytical solution :
 
 $$X(H) =9 −3f_1 H+54 f_1(\tau −9 uv) H+....$$ 
 
 
compare Kinetik Wave and Savage-Hutter
 
~~~gnuplot
 set key left bottom
 p[-20:]'../cohesifKW/frontKW-0.00.txt' w l ,'../cohesifSH/frontSH-0.000.txt' w l
~~~  
 
 
Advance of front: $x_f \sim U t$, as we push, it is the same 

~~~gnuplot
 reset
 set key left bottom
 set xlabel "t"
 p'xSH-0.000.txt',.2363*x+1.6 t'line','../cohesifKW/xKW-0.00.txt'
~~~

but when change-ing cohesion, the depth changes

~~~gnuplot
 reset
 set key left bottom
 set xlabel "t"
 p'xSH-0.000.txt',.2363*x+1.6 t'line','xSH-0.005.txt','xSH-0.010.txt'
~~~
 

~~~gnuplot
 reset
 set key left bottom
 set xlabel "x"
 p'frontSH-0.000.txt','frontSH-0.005.txt','frontSH-0.010.txt'
~~~
 
 
 
~~~gnuplot
 set key left bottom
 p[-20:]'frontSH-0.001.txt' w l,'../cohesifKW/frontKW-0.05.txt' w l
~~~  


 
 
 
 
# Links



* [http://basilisk.fr/sandbox/M1EMN/Exemples/cohesive_savagehutter.c]()

* [http://basilisk.fr/sandbox/M1EMN/TEST/cohesifKW.c]()
*/
