/**

# Resolution de Poiseuille
 
with 'Chinease Method':  solve NS with penalization for the wall. For, $r>1$ impose at each time step $u=v=0$
*/
#include "axi.h"
#include "navier-stokes/centered.h" 
#include "two-phase.h"
 
double dR,xst,Re;
int LEVEL;

double stenose( double x){
  double xb,xl;
  xb=xst-.5;
  xl=1.;
  return (fabs((x-xb)/xl-.5)<.5?  (1+cos(pi+2*pi*(x-xb)/xl))/2.   : 0 );
}

int main(int argc, char const *argv[]) {

  dR = 0.0;
  xst = 1.5;
  Re = 100.;  // Reynolds
 /*
  if (argc >= 2){   dR = atof(argv[1]);}
  if (argc >= 3){   xst = atof (argv[2]);}
  if (argc >= 4){   Re = atof (argv[3]);  } */

  L0 = 1.;         // la boite originale est 10x10

  TOLERANCE = 1e-6;
  // two-phase gives f=1 to 1 and f=0 to phase 2  
  rho1 = 1,rho2 = 1.;
  mu1 = 1./Re ;
  mu2 = 1.;

  
  LEVEL=6;    // Nombre de oints par cote  delta x = L0/(2**LEVEL)
  N = 1 << LEVEL;
  system("rm puv.dat");
  run();
  //2.736 real
    
  L0=2;
  LEVEL=7;
  N = 1 << LEVEL;
  system("rm puv.dat");
  run();
  //44.75 real
    
  LEVEL=8;
  N = 1 << LEVEL;
  system("rm puv.dat");
  run();
  //694 real
  
  // un comment 
  /*
  LEVEL=9;
  N = 1 << LEVEL;
  system("rm puv.dat");
  run();
*/
}

/**
 Boundary conditions, see below
 */
u.t[top] = dirichlet(0); 
u.n[top] = dirichlet(0); 
u.n[left] = ( y <= 1 ? dirichlet(2*(1-sq(y))):dirichlet(0.)); 

p[left] =  neumann(8./Re) ; 
pf[left]   = neumann(8./Re) ;
u.n[right] =  ( y <= 1 ? neumann(0.) : dirichlet(0.) ) ;
p[right]   = dirichlet(0.);
pf[right]  = dirichlet(0.);

event init (t = 0) {
  
fraction (f, (1. - dR *stenose(x))- y  );
/** Poiseuille init
 
 $$0 = -\frac{\partial p}{\partial x} + \frac{1}{r} \frac{\partial }{\partial r}  r \frac{\partial u}{\partial r}$$
 
 gives $u = 2 (1-r^2)$ for a flux equal to $\int u r dr = 1/2$ and pressure gradient $-8/Re$ (see BC)

*/
   foreach() {
     u.x[]=2.*(1-y*y)*(f[]);
     u.y[]=0.0;
     p[]=8.*(L0-x)/Re*(f[]);
    }

/**

 save the shape of the stenosis, here flat!
*/
  FILE * fp2 = fopen ("interface.csv", "w");
  output_facets (f, fp2);
  fprintf(fp2,"\n");
  fclose (fp2);
}

// at each iteration do zero everywhere in phase 2: 'chinease method'
event bound (i++,last){ 
  fraction (f, (1. - dR *stenose(x))- y  );

  foreach() {
    if (f[]>1e-6 && f[]<1.-1e-6){
     u.x[]=u.y[]=p[]=0.0;
    } else {
    u.x[] = u.x[]*(f[]);
    u.y[] = u.y[]*(f[]);
     p[] = p[]*(f[]);
   }
  }

}
/**
 old value of the velocity is saved
 for time convergence tests
*/
scalar un[];
event init_un (i = 0) {
    foreach()
    un[] = u.x[];
}

event conv (t += 1 ) {
    double du = change (u.x, un);
    fprintf(stdout,"t= %g %g \n",t, du);
    if (i > 0 && du < 1.0e-4)
        return 1; /* stop */
}

#if gfv
event gfsview (i += 10){
    static FILE * fp = popen ("gfsview2D -s pois.gfv", "w");
    output_gfs (fp);
}
#endif

event sauve(t += 1; t <= 200)
{
  char filename[80];
  sprintf (filename, "puv.dat");
  // remove(filename);
  FILE * fp = fopen (filename, "a");
  output_field ({p,u}, fp, linear = true);
  fprintf(fp,"\n\n");
  fclose (fp);
}
  
event print (t=end) {
/** for 'bview'
*/
  scalar press[];
  foreach(){
    press[]=p[];
  }
  dump (file = "snapshot");

  char fichier[80];
  sprintf(fichier, "pfinalL0=%3.2fN=%1d.dat", L0, N);
  FILE  *fp = fopen (fichier, "w");
  output_field ({p,u}, fp, linear = true);
  fclose (fp);

}

/**



*/
 
/**
if no '#include "grid/multigrid.h"'' then adapt using a diffusive trick on 'f'
*/
#if QUADTREE
 
event adapt(i+=1){
 scalar K1[],K2[]; 
 foreach()
   K1[]=(f[0,1]+f[0,-1]+f[1,0]+f[-1,0])/4;
 boundary({K1});
 
 for(int k=1;k<3;k++)  
 {
 foreach()
   K2[]=(K1[0,1]+K1[0,-1]+K1[1,0]+K1[-1,0])/4;
 boundary({K2});

 foreach()
   K1[]=K2[]; 
}

 foreach()
   K1[]=(K2[0,1]+K2[0,-1]+K2[1,0]+K2[-1,0])/4*noise();;
 boundary({K1});

 adapt_wavelet({K1,f},(double[]){0.001,0.01}, maxlevel = LEVEL, minlevel = 2);
}
#endif


/**

# Run

pour lancer 

~~~bash

qcc -O2  -Dgfv=1 -o poiseuilleCM poiseuilleCM.c -lm
./poiseuilleCM


make poiseuilleCM.c.tst
make poiseuilleCM.c/plots
make poiseuilleCM.c.c.html

~~~


# Results
 
## Effective Poiseuille

 For a domain of size $L0=1$, no approximation, we recover the Poiseuille solutuon:
 
~~~gnuplot pressure
 reset
 set xlabel "x"
 set ylabel "p"
 p'pfinalL0=1.00N=64.dat' u 1:($2==0?$3:NaN) ,-8*(x-1)/100
~~~
 
 the velocity is the parabola:

~~~gnuplot  velocity
 reset
 set dummy r
 set xlabel "r"
 set ylabel "u(x,r)"
 p[0:1.]'pfinalL0=1.00N=64.dat'u 2:4 w l,(r<1? 2*(1-r*r) : NaN)
 
 ~~~

 
## Poiseuille CM
 
For a domain of size $L0=2$, half of the domain is penlaized,
 
~~~gnuplot
 reset
 set xlabel "x"
 set ylabel "p"
 p'pfinalL0=2.00N=128.dat' u 1:($2==0?$3:NaN) t'N=128',\
  'pfinalL0=2.00N=256.dat' u 1:($2==0?$3:NaN) t'N=256',\
  'pfinalL0=2.00N=512.dat' u 1:($2==0?$3:NaN) t'N=512',\
 -8*(x-2)/100 t'Poiseueille'
~~~
 
 
~~~gnuplot  velocity
 reset
 set dummy r
 set xlabel "r"
 set ylabel "u(x,r)"
 p[0:1.]'pfinalL0=2.00N=128.dat'u 2:4 w l,\
        'pfinalL0=2.00N=256.dat'u 2:4 w l,\
    (r<1? 2*(1-r*r) : NaN)
 
 ~~~

   
in file 'puv.dat' are all the blocks, take only the last one with 'stats'
 
~~~gnuplot  pressure
 reset
 set xlabel "x"
 set ylabel "p"
stats 'puv.dat'  u 1
p'puv.dat' u 1:($2==0?$3:NaN) i (STATS_blocks-2),-8*(x-2)/100,0
 
~~~

 ~~~gnuplot  central velocity
 reset
 set xlabel "x"
 set ylabel "u(x,r=1)"
 p'puv.dat'  u 1:($2==0?$4:NaN) i (STATS_blocks-2) w d,2
 
 ~~~
 

# bibliography

* any textbook

* [The RNS/Prandtl equations and their link with other asymptotic descriptions: Application to the wall shear stress scaling in a constricted pipe
Pierre-Yves Lagrée, Sylvie Lorthois](http://www.lmm.jussieu.fr/~lagree/TEXTES/PDF/lagreelorthois05.pdf)

*/













