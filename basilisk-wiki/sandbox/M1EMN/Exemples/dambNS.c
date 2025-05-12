/**
# Dam Break

We solve by `navier-stokes/centered.h` the dam-break and compare with Ritter's solution from 1D Saint-Venant.

~~~gnuplot  configuration
set size ratio 4/10. 
set arrow from 4,0 to 4,1 front
set arrow from 4,1 to 4,0 front
set arrow from 0.1,.1 to 9.9,0.1 front
set arrow from 9.1,.1 to 0.1,0.1 front
set label  "L0" at 6,.25 front
set label  "depth h=1" at 2.5,.5 front
set label  "water" at 1.2,.9 front
set label  "air" at 1.2,1.05 front
set xlabel "x"
set ylabel "h"
p [0:10][0:2.5]0 not, (x<=5) w filledcurves x1 linec 3 t'free surface' 
~~~

# Code

viscosity MU is `1/Re`, the value is indicative... The `LEVEL` as well, a finer grid will create instabilities...
*/
#include "navier-stokes/centered.h"
#include "vof.h"
#define RHOF 1.e-3
#define MU 1./1000
#define LEVEL 8

scalar f[],*interfaces = {f};
face vector alphav[];
face vector muv[];

int main() {
    L0 = 10.;
    DT = 1e-2;
    init_grid (1 << LEVEL);
    u.n[bottom] = dirichlet(0);
    u.t[bottom] = neumann(0);
 
    run();  
}
/**

Initial condition, a heap $h=1$ for $x<L_0/2$
*/
event init (t = 0) {
    mask (y >  L0/4 ? top :
          none);
    const face vector g[] = {0,-1};
    a = g;
    alpha = alphav;
    scalar phi[];
    foreach_vertex(){
      phi[] =  min(1 - y, L0/2 - x);
    // 4studium; with jump 
    // phi[] =  min((.9*(x<L0/2)+.1) - y, L0 - x);
}
    fractions (phi, f);
    
    foreach(){  	
        u.x[] =    f[] * (1e-8);
        u.y[] =    0;     
    }
}
/**
 total density
 */
#define rho(f) ((f) + RHOF*(1. - (f)))
#define muc(f) ((f)*MU + MU/100.*(1. - (f)))

event properties (i++) {
/**
  We set a  ad hoc viscosity field, and ad hoc density */
  mu = muv;
  	
    foreach_face() {
        double fm = (f[] + f[-1])/2.;
        alphav.x[] = 1./rho(fm);
        muv.x[] = muc(fm);
    }
    boundary ((scalar *){muv,alphav});
}
/**
 convergence outputs
 */
void mg_print (mgstats mg)
{
    if (mg.i > 0 && mg.resa > 0.)
        fprintf (stderr, "#   %d %g %g %g\n", mg.i, mg.resb, mg.resa,
                 exp (log (mg.resb/mg.resa)/mg.i));
}
/**
 convergence outputs
*/
event logfile (i++) {
    stats s = statsf (f);
    fprintf (stderr, "%g %d dt=%g %g %g %g\n", t, i, dt, s.sum, s.min, s.max - 1.);
    mg_print (mgp);
    mg_print (mgpf);
    mg_print (mgu);
    fflush (stderr);
}
/**
for plots: the interface every $t$
*/
event interface (t +=1;t <= 3){
     output_facets (f);
    fprintf(stdout,"\n"); 
    fprintf(stderr,"------~~~~~~-----  \n" );
}

event movie (t += 0.05) {
    scalar l[];
    foreach()
    l[] = -1 + f[]* (.25+    sqrt(sq(u.x[]) + sq(u.y[])));
    boundary ({l}); 
    output_ppm (l, min = -1,  max= 3.25,
      linear = true,
                n = 1024, box = {{0,-1},{10,2.5}},
                file = "dambNS.mp4"
                );                   
    foreach()
    l[] = level;
    boundary ({l}); 
    output_ppm (l, n = 1024,  min = 0.,   max= 12.,
      box = {{0,-1},{10,2.5}},
      file = "level.mp4")  ;         
}
#if 0
event movie (t += 0.05) {
    scalar l[];
    foreach()
    l[] = -1 + f[]* (.25+    sqrt(sq(u.x[]) + sq(u.y[])));
    boundary ({l});
    static FILE * fp2 = popen ("ppm2mpeg > dambNS.mpg", "w");
    output_ppm (l, fp2 , min = -1,  max= 3.25,
      linear = true,
                n = 2048, box = {{0,-1},{10,2.5}}
                );          
    
    foreach()
    l[] = level;
    boundary ({l});
    static FILE * fp1 = popen ("ppm2mpeg > level.mpg", "w");
    output_ppm (l, fp1,   box = {{0,0},{10,2.5}})  ;          
}
#endif
event pictures (t=0.05) {
    scalar l[];
    foreach()
    l[] = f[]* (   sqrt(sq(u.x[]) + sq(u.y[])));; 
    boundary ({l});
     output_ppm (l, file = "dambNS.png", min = 0,   max= 2,
     // linear = true,
                //n = 1024 ,
                 box = {{0,0},{10,2.5}}
                );
    foreach()
    l[] = level; 
    boundary ({l});
     output_ppm (l, file = "level.png", min = 0,   max= 12,
     // linear = true,
                //n = 1024 ,
                 box = {{0,0},{10,2.5}}
                );
             
}

#if QUADTREE
event adapt(i++){
 scalar g[];
 foreach()
   g[]=f[]*(1+noise())/2;
 boundary({g});
 //adapt_wavelet ({g,f,u.x,u.y}, (double[]){0.001,.01,0.01,0.01}, minlevel = 5,maxlevel = LEVEL);
 adapt_wavelet({g,f},(double[]){0.001,0.01},minlevel = 5,maxlevel = LEVEL);
} 
#endif

/**

# Run

to run

~~~bash
 qcc -g -O2 -Wall  -o dambNS dambNS.c -lm
 ./dambNS > out
~~~

or with `makefile`

~~~bash
make dambNS.tst;make dambNS/plots;make dambNS.c.html
~~~


# Results

Plot of the interface compared with Ritter solution from Saint-Venant, of course it is different as SV supposes 
small slope for $h$. 
 
~~~gnuplot plot vs Shallow Water Solution
reset
h(x,t)=(((x-5)<-t)+((x-5)>-t)*(2./3*(1-(x-5)/(2*t)))**2)*(((x-5)<2*t))
set xlabel "x"
set ylabel "h(x,t)"
p[0:10][0:2]'out' t 'comp t=0,1,2,3 'w l, h(x,1)w l linec -1 not ,h(x,2)w l linec -1 not ,h(x,1)w l linec -1 not ,h(x,3) w l linec -1 t'Ritter'

~~~
 
 ![Animation, hight colored by velocity](dambNS/dambNS.mp4)   


 ![level](dambNS/level.png)
 ![Animation  level ](dambNS/level.mp4)  


# Links
 
see the non viscous dam break with [Basilisk](http://basilisk.fr/sandbox/M1EMN/Exemples/damb.c)  

# Bibliography

* [PYL](http://www.lmm.jussieu.fr/~lagree/COURS/MFEnv/MFEnv.pdf)
* P. K. Stansby, A. Chegini And T. C. D. Barnes The initial stages of dam-break flow J. Fluid Mech. (1998), vol. 374, pp. 407–424. 
*  Marcela A. Cruchaga, Diego J. Celentano Tayfun E. Tezduyar, Collapse of a liquid column: numerical simulation and experimental validation
Comput. Mech. (2007) 39: 453–476 DOI 10.1007/s00466-006-0043-z


Paris 10/18
*/
 
    


