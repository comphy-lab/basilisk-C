#include "embed.h"
#include "navier-stokes/centered.h"
// #include "navier-stokes/perfs.h"
#include "LESSmagorinsky.h"
#include "tracer.h"
#include "view.h"

scalar f[];
scalar * tracers = {f};
double Reynolds =80000.;
int maxlevel=10;


double hs = 0.5,hi=0.5, U0 = 1.;
double Li=3.;
double Lo=30.;   //Lo=9.

double molvis;  //molecular viscosity
double Cs;      //Smagorinsky constant

scalar un[];
int fin;

int main()
{

  double L0 = Li+Lo;
  size(L0);
  origin (-Li, -hs-0.5);
  N = 512;
  init_grid(N);
  fin=0;
  TOLERANCE=1e-5 [*];
  run();

}


/**
The fluid is injected on the left boundary with velocity $U0$. An outflow
condition is used on the right boundary. */

u.n[left]  = cs[]? dirichlet(6*U0*y/hi*(1-y/hi)):dirichlet(0.);
//p[left]    = neumann(0.);
//pf[left]   = neumann(0.);
//f[left]    = dirichlet(y < 0);

//u.n[right] = neumann(0.);
p[right]   = dirichlet(0.);
pf[right]  = dirichlet(0.);

u.n[bottom]=dirichlet(0.);
u.t[bottom]=dirichlet(0.);

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

  molvis=2*hi*U0/Reynolds;
  Cs=0.1;

}

/**
We check the number of iterations of the Poisson and viscous
problems. */

event logfile (i++) {
  
  double avg=normf(u.x).avg;
  double du=change(u.x,un)/(avg+1e-30);
  if (i>1&&du<1e-3) {
    fin=1;
  }
  
  fprintf (stderr, "%g %d %d %g %g %d %d\n", Lo,maxlevel, i, t,du, mgp.i, mgu.i);

}

event recirculation (t+=1.; t<=200.) {
  char name[1024];
  sprintf(name,"recirculation_length_Li=%g_Lo=%g_level=%d.txt",Li,Lo,maxlevel);
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
We produce animations of the vorticity and tracer fields... */

event movies (t+=1; t<=200.)
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
  sprintf(name,"velocity_x_Li=%g_Lo=%g_level=%d.mp4",Li,Lo,maxlevel);
  output_ppm (u.x, file = name,linear = true, mask = m,box = {{-Li,-hs},{Lo,hi}});
 output_ppm (u.y, file = "velocity_y.mp4",linear = true, mask = m,box = {{-Li,-hs},{Lo,hi}});
 output_ppm (p, file = "pressure.mp4",linear = true, mask = m,box = {{-Li,-hs},{Lo,hi}});
  output_ppm (llevel, file = "level.mp4",linear = false, mask = m,box = {{-Li,-hs},{Lo,hi}});

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
  isoline("u.x",n=21);
  squares("m2",max=U0,min=-U0);
  squares("cs",map=gray,spread=-1);
  box();
  char name2[1024];
  sprintf(name2,"isoline_Re=%g.mp4",Reynolds);
  save(name2);

  scalar mueddy[];
  double sum;
  foreach() {
    if (cs[]>0) {
      sum=0.;
      foreach_neighbor() {
	sum+=mu.x[];
      }
      mueddy[]=sum/4.;
    }
    else {
      mueddy[]=nodata;
    }
  }
  view(tx=-0.4,ty=0.01,sy=2.,fov=5,width=800,height=200,camera="front");
  clear();
  //isoline("u.x",n=21);
  squares("mueddy");//,max=U0,min=-U0);
  squares("cs",map=gray,spread=-1);
  box();
  char name3[1024];
  sprintf(name3,"mueddy_Re=%g.mp4",Reynolds);
  save(name3);
}

event images (t=200.)
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
  sprintf(name,"velocity_x_Li=%g_Lo=%g_level=%d.png",Li,Lo,maxlevel);
  output_ppm (u.x, file = name,linear = true, mask = m,box = {{-Li,-hs},{Lo,hi}},min=0,max=3/2*U0);
 output_ppm (u.y, file = "velocity_y.png",linear = true, mask = m,box = {{-Li,-hs},{Lo,hi}});
 output_ppm (p, file = "pressure.png",linear = false, mask = m,box = {{-Li,-hs},{Lo,hi}});
  output_ppm (llevel, file = "level.png",linear = false, mask = m,box = {{-Li,-hs},{Lo,hi}});
  output_ppm (fm.x, file = "test.png",linear = true, mask = m,box = {{-Li,-hs},{Lo,hi}});

  FILE * result=fopen("velocity_x_Re=80000.txt","w");
  int number=50;
  double xp;
  double res;
  for (int kk=0;kk<number;kk++) {
    xp=kk*10.0/number;
    Point point=locate(xp,-0.45);
    res=interpolate_linear(point,u.x,xp,-0.45);
    fprintf(result,"%g,%g\n",xp,res);
    fflush(result);
  }

  double xf=Lo-1.;
  int number2=50;
  double yp;
  char name2[1024];
  sprintf(name2,"velocity_end_Li=%g_Lo=%g_level=%d.txt",Li,Lo,maxlevel);
  FILE * velocity=fopen(name2,"w");
  for (int kk=0;kk<=number2;kk++) {
    yp=-hs+(hi+hs)*kk/number2;
    res=interpolate(u.x,xf,yp);
    fprintf(velocity,"%g,%g\n",yp,res);
    fflush(velocity);
  }

  double xi=-Li+1.5;
  int number3=50;
  double ypi;
  char name3[1024];
  sprintf(name3,"velocity_init_Li=%g_Lo=%g_level=%d.txt",Li,Lo,maxlevel);
  FILE * velocity2=fopen(name3,"w");
  for (int kk=0;kk<=number3;kk++) {
    ypi=hi*kk/number3;
    res=interpolate(u.x,xi,ypi);
    fprintf(velocity2,"%g,%g\n",ypi,res);
    fflush(velocity2);
  }

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
  isoline("u.x",n=21);
  squares("m2",max=U0,min=-U0);
  squares("cs",map=gray,spread=-1);
  box();
  char name4[1024];
  sprintf(name4,"isoline_Re=%g.png",Reynolds);
  save(name4);

  scalar mueddy[];
  double sum;
  foreach() {
    if (cs[]>0) {
      sum=0.;
      foreach_neighbor() {
	sum+=mu.x[];
      }
      mueddy[]=sum/4.;
    }
    else {
      mueddy[]=nodata;
    }
  }
  view(tx=-0.4,ty=0.01,sy=2.,fov=5,width=800,height=200,camera="front");
  clear();
  //isoline("u.x",n=21);
  squares("mueddy");//,max=U0,min=-U0);
  squares("cs",map=gray,spread=-1);
  box();
  char name5[1024];
  sprintf(name5,"mueddy_Re=%g.png",Reynolds);
  save(name5);

  //return 1; //end the code
}
  
  
/**
We adapt according to the error on the embedded geometry and velocity. */
#if 1
event adapt (i++) {
  adapt_wavelet ({cs,u}, (double[]){1e-2,3e-2,3e-2}, maxlevel, 3);
}
#endif

