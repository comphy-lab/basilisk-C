/**
# blood flow in curved pipe

Idea from Antoonvh [tube](http://basilisk.fr/sandbox/Antoonvh/tube.c)

Some exemples of the functions used in this code.
[RRC](biomeca/2dworrc.c) \&
[Womersley](biomeca/2dwoflow.c) 

We want to test if the new "[embed](http://basilisk.fr/src/embed.h)" function can better treat with the non-slip wall conditions (and the sharp velocity/pressure differences at the inlet outlet surface).

We want to simulate the blood circulation in a torus shape pipe (a simplified model for the aortic arche) with major radius $R_t$ and minor radius $r_t$. 

We propose a sinusoidal shape input flux and a KV model as output condition to simulate the effect of heart beat and following artery.*/
#include "grid/octree.h"
#include "embed.h"
#include "navier-stokes/centered.h"
#include "view.h"
face vector av[], muc[];
/**
   The torus is defined using these variables and macros:
*/
double minor_rt = 1., major_rt = 3.;
#define R      (sqrt(sq(x) + sq(y)))
#define SINPSI (y/R)
#define COSPSI (x/R)
#define TUBE   (R <= 0.01 ? major_rt				\
		: (sqrt(sq(x - major_rt*COSPSI) +		\
			sq(y - major_rt*SINPSI) + sq(z))))
/**
 Macros for KV model*/
#define alpha_w 10.
#define  ORR 0.1
#define  ORC 0.1
#define  OCC 10
#define tend 50.
double	Pout,Pold,Qold,PPPP,DEBIT;


/**
We set the inlet and outlet conditions on the bottum.


Z - front (z =L0) back (z=0) 

Y - top bottom 

X - left right
*/

u.n[embed] = dirichlet (0);
u.t[embed] = dirichlet (0);
#if (dimension == 3)
u.r[embed] = dirichlet (0);
#endif
u.n[bottom] = ( x >= 0)? dirichlet(5.*cs[]*cos(t)):neumann(0);
u.t[bottom] = ( x >= 0)? dirichlet(0):neumann(0);
#if (dimension == 3)
u.r[bottom] = ( x >= 0)? dirichlet(0):neumann(0);
#endif

//p[bottom] = ( y >= 0)? dirichlet(0):neumann(0) ;
//pf[bottom] = ( y >= 0)? dirichlet(0):neumann(0) ;
p[top] = neumann(0) ;
pf[top] = neumann(0) ;
int j;
int main() {
  L0 = 10;
  X0 = Z0 = -L0/2.;
  N = 64;
  mu = muc;
  a = av;
  run();
}
event init (t = 0) {
  vertex scalar phi[];
  foreach_vertex() 
    phi[] = minor_rt - TUBE;
  fractions (phi, cs, fs);
  view (phi = 0.1, psi = 0., theta = 0.8);
  box();
  draw_vof ("cs", edges = true, lw = 0.5);
  save ("vof.png");//		f0.refine = f0.prolongation = fraction_refine;		
  Pout = DEBIT = Pold = Qold = 0.;
  j=0;
}


/**
KV model bc */
event mybc (i++){
  double RRe =0.8;
  double alpha = (RRe + ORR) / ORC + 1.;
  double beta =  1. - RRe/(RRe + ORR +ORC);
  double QQQ = beta * (1.-exp( -alpha*t/((RRe+ORR)*OCC)))* DEBIT;
  double term1 = ORC * OCC * Pold / dt;
  double term2 = (ORR + ORC) * DEBIT ;
  double term3 = ORR * ORC * OCC * (DEBIT - Qold) / dt;
  double deno = 1. + ORC * OCC / dt ;

  Pout =  (term1 + term2 + term3) / deno;

  p[bottom] = ( y >= 0)? dirichlet(Pout):neumann(0) ;
  pf[bottom] = ( y >= 0)? dirichlet(Pout):neumann(0) ;

  Pold = Pout;
  Qold = DEBIT;
}

/**
movie for fun*/
event movies (t += 0.05) {
  scalar omega[];
  vorticity(u,omega);
  view(width=1280,height=640);
  box();
  cells();
  squares ("u.y",max=5.,min=-5.);
  //squares ("omega");
  save ("cs.mp4");	
}



/**
We note the velocity of axial direction at the middle of the model (x,z=0)*/
event veloprofil (t += tend*0.04)
{
  if (t > tend*0.8){
    j += 1;
    char name2[80];
    sprintf (name2, "test-%d", j);	
    FILE *fp2;
    fp2 = fopen (name2, "w");
    for (double uyy = 1.5; uyy <= 4.5; uyy += 0.01)
      fprintf(fp2,"%g %g %g\n", t, uyy,interpolate(u.x , 0., uyy,0.));		
  }
}

/**
Dynamic adaption, with maximum resolution corresponds to a $64^3$.
*/

event adapt (i++) {
  fprintf(stderr, "%d %g %g\n", i, t,interpolate(u.x, 0., 3., 0.));
  adapt_wavelet ((scalar*) {cs,u},(double[]){1e-4,0.01, 0.01, 0.01}, 7);
}

event stop (t = tend);


/**
# Results

![Animation of axial velocity and mesh adaption](torusrrc/cs.mp4)


~~~gnuplot axial velocity 
plot 'test-1' u 2:3 w l t'time1',\
'test-2' u 2:3 w l t'time2',\
'test-3' u 2:3 w l t'time3',\
'test-4' u 2:3 w l t'time4',\
'test-5' u 2:3 w l t'time5'
~~~

*/


