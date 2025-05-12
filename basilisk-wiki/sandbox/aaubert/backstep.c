#include "embed.h"
#include "navier-stokes/centered.h"
// #include "navier-stokes/perfs.h"
#include "tracer.h"
#include "view.h"

scalar f[];
scalar * tracers = {f};
double Reynolds =800.;
int maxlevel;
face vector muv[];

double hs = 0.5,hi=0.5, U0 = 1.;       //size of the step
double Li=3.;                          //size of the inlet
double Lo=30.;                         //size of the outlet

scalar un[];
int fin;

int main()
{
  double L0 = Li+Lo;
  size(L0);
  origin (-Li, -hs-0.5);
  N = 512;
  init_grid(N);
  mu = muv;
  fin=0;
  TOLERANCE=1e-5 [*];
  run();
}

/** 
The viscosity is based on the Reynolds number */

event properties (i++)
{
  foreach_face()
    muv.x[] =fm.x[]*2*hi*U0/Reynolds; 
}

/**
The fluid is injected on the left boundary with a Poiseuille flow of mean velcoity $U0$. An outflow condition is used on the right boundary. */

u.n[left]  = cs[]? dirichlet(6*U0*y/hi*(1-y/hi)):dirichlet(0.);

p[right]   = dirichlet(0.);
pf[right]  = dirichlet(0.);

u.n[bottom]=dirichlet(0.);
u.t[bottom]=dirichlet(0.);

/**
No-slip condition on the solid */

u.n[embed] = dirichlet(0.);
u.t[embed] = dirichlet(0.);


event init (t = 0)
{

  solid (cs, fs, intersection(intersection(union (y+1e-14, x+1e-14),hi-y+1e-14),y+0.5+1e-14));

  /**
  We set the initial velocity field. */
  
  foreach()
    u.x[] = cs[] ? U0 : 0.;

  foreach()
    un[]=u.x[];

}

/**
We check the number of iterations of the Poisson and viscous
problems and if we have reached a stationnary flow */

event logfile (i++) {
  
  double avg=normf(u.x).avg;
  double du=change(u.x,un)/(avg+1e-30);
  if (i>1&&du<1e-3) {
    fin=1;
  }
  
  fprintf (stderr, "%g %d %g %g %d %d\n", Lo, i, t,du, mgp.i, mgu.i);

}

/**
Determination of the recirculation length define as the position where the velocity on the bottom wall goes to zero.
For this, we do a dichotomy algorithm just above the bottom wall*/

event recirculation (t+=1.; t<=400.) {
  char name[1024];
  sprintf(name,"recirculation_length_Li=%g_Lo=%g_Re=%g.txt",Li,Lo,Reynolds);
  static FILE * recirculation=fopen(name,"w");
  double xa=1.;
  double xb=10.;
  double ua=interpolate(u.x,xa,-0.45);
  double ub=interpolate(u.x,xb,-0.45);
  double xm=0.;
  double um=0.;
  int compteur=0;
  if (ua*ub>0) {
    fprintf(recirculation,"%g,%g,%d\n",t,0.,compteur);
    fflush(recirculation);
  }
  else {
    while (xb-xa>1e-3) {
      xm=(xa+xb)/2;
      um=interpolate(u.x,xm,-0.45);
      compteur+=1;
      if (ua*um>0) {
	xa=xm;
	ua=um;
      }
      else {
	xb=xm;
	ub=um;
      }
    }
    fprintf(recirculation,"%g,%g,%d\n",t,xm,compteur);
    fflush(recirculation);
  }

}

/**
We produce animations of the vorticity and horizontal velocity field as well as a movie with the evolution of the isolines of horizontal velocity. */

event movies (t+=1; t<=400.)
{
  scalar omega[], m[];
  scalar llevel[];
  vorticity (u, omega);
  foreach() {
    m[] = cs[] - 0.5;
    llevel[]=level;
  }
  output_ppm (omega, file = "vort.mp4",linear = true, mask = m,box = {{-Li,-hs},{Lo,hi}});
  char name[1024];
  sprintf(name,"velocity_x_Li=%g_Lo=%g_Re=%g.mp4",Li,Lo,Reynolds);
  
  scalar m2[];
  foreach() {
    if (cs[]>0) {
      m2[]=u.x[];
    }
    else {
      m2[]=nodata;
    }
  }
  view(tx=-0.4,ty=0.01,sy=2.,fov=5,width=800,height=200,camera="front");
  clear();
  isoline("m2",n=21);
  squares("m2");//,max=U0,min=-U0);
  squares("cs",map=gray,spread=-1);
  box();
  char name4[1024];
  sprintf(name4,"isoline_Re=%g.mp4",Reynolds);
  save(name4);
}

/**
We get the final state of the vorticity and the horizontal velocity fields and the position of the isolines. */

event images (t=400.)
{
  scalar omega[], m[];
  scalar llevel[];
  vorticity (u, omega);
  foreach() {
    m[] = cs[] - 0.5;
    llevel[]=level;
  }
  output_ppm (omega, file = "vort.png",linear = true, mask = m,box = {{-Li,-hs},{Lo,hi}});
  char name[1024];
  sprintf(name,"velocity_x_Li=%g_Lo=%g_Re=%g.png",Li,Lo,Reynolds);

  scalar m2[];
  foreach() {
    if (cs[]>0) {
      m2[]=u.x[];
    }
    else {
      m2[]=nodata;
    }
  }

  view(tx=-0.4,ty=0.01,sy=2.,fov=5,width=800,height=200,camera="front");
  clear();
  isoline("m2",n=21);
  squares("m2");//,max=U0,min=-U0);
  squares("cs",map=gray,spread=-1);
  box();
  char name4[1024];
  sprintf(name4,"isoline_Re=%g.png",Reynolds);
  save(name4);

/**
We output also the  horizontal profile of the velocity at the oultet to check that we recover a Poiseuille flow. */
  
  double xf=Lo-1.;
  int number2=50;
  double yp;
  char name2[1024];
  sprintf(name2,"velocity_end_Li=%g_Lo=%g_Re=%g.txt",Li,Lo,Reynolds);
  FILE * velocity=fopen(name2,"w");
  for (int kk=0;kk<=number2;kk++) {
    yp=-hs+(hi+hs)*kk/number2;
    res=interpolate(u.x,xf,yp);
    fprintf(velocity,"%g,%g\n",yp,res);
    fflush(velocity);
  }

/**
We get data for comparison with Ertruk (2008) */

  vector gu[],gv[];
  gradients ({u.x,u.y},{gu,gv});
  double res2,res3,res4,res5,res6,res7;
  double x6=6*hs;
  int number4=50;
  double yp6;
  char name5[1024];
  sprintf(name5,"data_x=6h_Re=%g_Li=%g_Lo=%g.txt",Reynolds,Li,Lo);
  FILE * data6=fopen(name5,"w");
  for (int kk=0;kk<=number4;kk++) {
    yp6=-hs+(hi+hs)*kk/number4;
    res=interpolate(u.x,x6,yp6);
    res2=interpolate(u.y,x6,yp6);
    res3=interpolate(omega,x6,yp6);
    res4=interpolate(gu.x,x6,yp6);
    res5=interpolate(gu.y,x6,yp6);
    res6=interpolate(gv.x,x6,yp6);
    res7=interpolate(gv.y,x6,yp6);
    fprintf(data6,"%g,%g,%g,%g,%g,%g,%g,%g\n",yp6,res,res2,res3,res4,res5,res6,res7);
    fflush(data6);
  }

  double x14=14*hs;
  int number5=50;
  double yp14;
  char name6[1024];
  sprintf(name6,"data_x=6h_Re=%g_Li=%g_Lo=%g.txt",Reynolds,Li,Lo);
  FILE * data14=fopen(name6,"w");
  for (int kk=0;kk<=number5;kk++) {
    yp14=-hs+(hi+hs)*kk/number5;
    res=interpolate(u.x,x14,yp14);
    res2=interpolate(u.y,x14,yp14);
    res3=interpolate(omega,x14,yp14);
    res4=interpolate(gu.x,x14,yp14);
    res5=interpolate(gu.y,x14,yp14);
    res6=interpolate(gv.x,x14,yp14);
    res7=interpolate(gv.y,x14,yp14);
    fprintf(data14,"%g,%g,%g,%g,%g,%g,%g,%g\n",yp14,res,res2,res3,res4,res5,res6,res7);
    fflush(data14);
  }

  double x30=30*hs;
  int number6=50;
  double yp30;
  char name7[1024];
  sprintf(name7,"data_x=6h_Re=%g_Li=%g_Lo=%g.txt",Reynolds,Li,Lo);
  FILE * data30=fopen(name7,"w");
  for (int kk=0;kk<=number6;kk++) {
    yp30=-hs+(hi+hs)*kk/number6;
    res=interpolate(u.x,x30,yp30);
    res2=interpolate(u.y,x30,yp30);
    res3=interpolate(omega,x30,yp30);
    res4=interpolate(gu.x,x30,yp30);
    res5=interpolate(gu.y,x30,yp30);
    res6=interpolate(gv.x,x30,yp30);
    res7=interpolate(gv.y,x30,yp30);
    fprintf(data30,"%g,%g,%g,%g,%g,%g,%g,%g\n",yp30,res,res2,res3,res4,res5,res6,res7);
    fflush(data30);
  }


  //return 1; //end the code
}
  
  
/**
We adapt according to the error on the embedded geometry and velocity. */
#if 1
event adapt (i++) {
  adapt_wavelet ({cs,u}, (double[]){1e-2,3e-2,3e-2}, maxlevel, 5);
}
#endif

/**
<video width="368" height="320" controls>
<source src="isoline_Re800.mp4" type="video/mp4">
</video> 
*/