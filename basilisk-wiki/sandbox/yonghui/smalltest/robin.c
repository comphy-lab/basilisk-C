/**
# Robin condition test 
See how to use the patch from [Thomas](http://basilisk.fr/sandbox/tfullana/patch_robinBC).


My defintion of Robin condition is a little bit different
$$ a \boldsymbol{u} + b \frac{\partial{\boldsymbol{u}}}{\partial \boldsymbol{n}} = c 
$$
Write in discretized form at boundary
$$
v[ ghost ] = \frac{2 c \Delta + (2b - a\Delta) v [\quad]}{2b + a\Delta}
$$
The code I used:

~~~
double aax = bbx = 1.;

@define robin(aax,bbx,ccx) ((2.*(ccx)*Delta/(2.*(bbx)+(aax)*Delta)) +((2.*(bbx)-(aax)*Delta)/(2.*(bbx)+(aax)*Delta))*val(_s,0,0,0))

@define robin_homogeneous() (((2.*(bbx)-(aax)*Delta)/(2.*(bbx)+(aax)*Delta))*val(_s,0,0,0))
~~~

Until now you have to define the values of $aax$ and $bbx$ in the code, in order to make sure robin_homogeneous() works.
*/



#include "grid/multigrid.h"
#include "axi.h"
#include "navier-stokes/centered.h"
#include "navier-stokes/swirl.h"
#include "view.h"
#define ROBIN 1

/**
## Boundary conditions */
u.n[top] = dirichlet(0);
u.t[top] = dirichlet(0);
#if ROBIN
w[top] = robin(aax,bbx,fc3);
#else
w[top]   = dirichlet(aax);
#endif

/**
## Setup */
FILE *fp;
double fa3,fb,fc3,AAA;

int main()
{
  fp=fopen("file1","w");
  N = 64;
  periodic(right);
  const face vector muc[] = {.1,.1};
  stokes = true;
  mu = muc;
  DT = 0.1;
  run();
}

event init (i = 0) {
  // a b c
  aax = 2.;
  bbx = 2.;
  fc3 =	2.;
}

event logfile (i += 10;i<400)
{
  fprintf(stderr,"%d %g\n",i,t);
  if (i > 380) {
    double wwt,ppt;
    for (double xxx=0.01; xxx <= L0 ;xxx += 0.01){
      wwt = interpolate(w,L0/2.,xxx);
      ppt = interpolate(p,L0/2.,xxx);
      fprintf(fp,"%g %g %g\n",xxx,wwt,ppt);
    }
  }
}