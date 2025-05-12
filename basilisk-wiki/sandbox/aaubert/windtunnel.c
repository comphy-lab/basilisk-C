/**
Test case for he k-epsilon model. It corresponds to a flow of velocity $U_1$ in a wind tunnel. The flow is allow to mix with ambient air through the bottom boundary. This creates a mixing layer between the fast and the slow air. For more details about the geometry, look at [Liepmann and Laufer](https://ntrs.nasa.gov/api/citations/19930081851/downloads/19930081851.pdf) and some result for different turbulent model [here](https://ntrs.nasa.gov/api/citations/19970017828/downloads/19970017828.pdf) */

#include "navier-stokes/centered.h"

#define WALL 0
#define WALL_FUNCTION 0
#define LINEARISATION 0
#define BVP4 0

#include "RANSkepsilon.h"
#include "navier-stokes/perfs.h"
#include "view.h"


double U1=18.;  //velocity of the air in the wind tunnel

//properties of air
double mu_air=1.85e-5;
double rho_air=1.204;

int maxlevel=12;

scalar rhov[];

scalar un[];

vector h2[];

int main() {
  for (maxlevel=8;maxlevel<9;maxlevel++) {
    size(1.);
    DT=0.1;
    origin(0.,-L0/2.);
    N=128;
    init_grid(N);
    TOLERANCE=1e-5 [*];
    rho=rhov;
    Deltamin=L0/pow(2.,maxlevel);
    theta=1.;         //we use limiter to avoid creation of region where k or epsilon become negative in the advection term
    k.gradient=minmod2;
    epsilon.gradient=minmod2;
    u.x.gradient=minmod2;
    u.y.gradient=minmod2;
    p.gradient=minmod2;

    run();
  }
}


event init(t=0) {
  //CFL=0.1;
  //ct3=0.;
  if (!restore(file="dump")) {  //we can restore from a previous simulation
    foreach() {
        molvis[]=mu_air/rho_air;
        u.x[]=U1*(0.5+0.5*tanh(10.*(y/(x+0.1)+0.007)));
        u.y[]=0.;
        un[]=u.x[];
        k[]=3./2.*sq(1./100.*U1);   //initial turbulence of 1%
        epsilon[]=rho_air*C_mu*sq(k[])/(1.*mu_air);  //mu_t=1 mu_air at the inlet
        mutke[]=rho_air*C_mu*sq(k[])/epsilon[];
        rhov[]=rho_air;
    }
  }
  else {
    N=512;
  }
  
  u.n[left]=dirichlet(U1*(0.5+0.5*tanh(10.*(y/(0.1)+0.007))));
  u.t[left]=dirichlet(0.);
  k[left]=dirichlet(3./2.*sq(1./100.*U1));
  epsilon[left]=dirichlet(rho_air*C_mu*sq(3./2.*sq(1./100.*U1))/(1.*mu_air));
  mutke[left]=dirichlet(1.*mu_air);

  p[right]=dirichlet(0.);
  pf[right]=dirichlet(0.);
  u.t[right]=neumann(0.);
  u.n[right]=neumann(0.);
  k[right]=neumann(0.);
  epsilon[right]=neumann(0.);
  mutke[right]=neumann(0.);


}

event logfile(i++) {
    
  stats smu=statsf(mutke);
  stats sk=statsf(k);
  stats seps=statsf(epsilon);
  double du=change(u.x,un);
  double uu=interpolate(u.x,0.5,0.001);
  stats ssource=statsf(sourcek);
  double xmax,ymax;
  double valuemax=mu_air/rho_air;
  double xmin,ymin;
  double valuemin=mu_air/rho_air;
  foreach() {
    if (mutke[]>valuemax) {
      xmax=x;
      ymax=y;
      valuemax=mutke[];
    }
    if (mutke[]<valuemin) {
      xmin=x;
      ymin=y;
      valuemin=mutke[];
    }
  }
 
  fprintf(stderr,"%g %d %d %g %g %d %d %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g\n",rho_air*U1*L0/mu_air,maxlevel,i,t,dt,mgp.i,mgu.i,smu.min,smu.max,du,uu,ssource.min,ssource.max,sk.min,sk.max,seps.min,seps.max,xmax,ymax,valuemax,xmin,ymin,valuemin);
  fflush(stderr);
  
  //We add dissipation at the right boundary to limit the return of the viscosity
  foreach() {
    if (x>L0-5.*L0/pow(2.,maxlevel)) {
      k[]=(1.*k[-1,1] + 2.*k[0,1] + 1.*k[1,1] + \
           2.*k[-1,0] + 4.*k[0,0] + 2.*k[1,0] + \
           1.*k[-1,-1] + 2.*k[0,-1] + 1.*k[1,-1])/16.;
      epsilon[]=(1.*epsilon[-1,1] + 2.*epsilon[0,1] + 1.*epsilon[1,1] + \
                 2.*epsilon[-1,0] + 4.*epsilon[0,0] + 2.*epsilon[1,0] + \
                 1.*epsilon[-1,-1] + 2.*epsilon[0,-1] + 1.*epsilon[1,-1])/16.;
    }
  }
  foreach() {
    mutke[]=rho[]*C_mu*sq(k[])/epsilon[];
  }
  foreach_face() {
    mut.x[]=fm.x[]*(face_value(rho,0)*face_value(molvis,0)+face_value(mutke,0));
  }
}



#if 1
event movies (t+=0.05;t<=5.) {

  scalar omega[];
  vorticity(u,omega);



  output_ppm (omega, file = "vort.mp4",linear = true);
  output_ppm (p, file = "pressure.mp4",linear = true);
  output_ppm (u.x, file = "velocity_x.mp4",max=U1,min=0.,linear = true);
  output_ppm (u.y, file = "velocity_y.mp4",max=1.,min=-1.,linear = true);
  output_ppm (mutke, file = "viscosity.mp4",linear = true);
  output_ppm (k,file="k.mp4",min=0.,linear=false);
  output_ppm (epsilon,file="epsilon.mp4",min=0.,linear=false);

}
#endif

event images(t=5.) {

  double xp1=0.1*L0;
  char name1[1024];
  sprintf(name1,"velocity_profile_x=%g_level=%d.txt",xp1,maxlevel);
  FILE * velocity1=fopen(name1,"w");
  int number1=pow(2,maxlevel);
  double yp1,res1,res1k,res1eps,res1mutke;
  for (int kk=0;kk<=number1;kk++) {
    yp1=-L0/2.+kk*1.0/(number1*1.0)*L0;
    foreach_point(xp1,yp1) {
      yp1=y;
      xp1=x;
      res1=u.x[];
      res1k=k[];
      res1eps=epsilon[];
      res1mutke=mutke[];
    }
    fprintf(velocity1,"%g,%g,%g,%g,%g,%g\n",xp1,yp1,res1,res1k,res1eps,res1mutke);
    fflush(velocity1);
  }

  double xp2=0.3*L0;
  char name2[1024];
  sprintf(name2,"velocity_profile_x=%g_level=%d.txt",xp2,maxlevel);
  FILE * velocity2=fopen(name2,"w");
  int number2=pow(2,maxlevel);
  double yp2,res2,res2k,res2eps,res2mutke;
  for (int kk=0;kk<=number2;kk++) {
    yp2=-L0/2.+kk*1.0/(number2*1.0)*L0;
    foreach_point(xp2,yp2) {
      yp2=y;
      xp2=x;
      res2=u.x[];
      res2k=k[];
      res2eps=epsilon[];
      res2mutke=mutke[];
    }
    fprintf(velocity2,"%g,%g,%g,%g,%g,%g\n",xp2,yp2,res2,res2k,res2eps,res2mutke);
    fflush(velocity2);
  }

  double xp3=0.5*L0;
  char name3[1024];
  sprintf(name3,"velocity_profile_x=%g_level=%d.txt",xp3,maxlevel);
  FILE * velocity3=fopen(name3,"w");
  int number3=pow(2,maxlevel);
  double yp3,res3,res3k,res3eps,res3mutke;
  for (int kk=0;kk<=number3;kk++) {
    yp3=-L0/2.+kk*1.0/(number3*1.0)*L0;
    foreach_point(xp3,yp3) {
      yp3=y;
      xp3=x;
      res3=u.x[];
      res3k=k[];
      res3eps=epsilon[];
      res3mutke=mutke[];
    }
    fprintf(velocity3,"%g,%g,%g,%g,%g,%g\n",xp3,yp3,res3,res3k,res3eps,res3mutke);
    fflush(velocity3);
  }

  double xp4=0.7*L0;
  char name4[1024];
  sprintf(name4,"velocity_profile_x=%g_level=%d.txt",xp4,maxlevel);
  FILE * velocity4=fopen(name4,"w");
  int number4=pow(2,maxlevel);
  double yp4,res4,res4k,res4eps,res4mutke;
  for (int kk=0;kk<=number4;kk++) {
    yp4=-L0/2.+kk*1.0/(number4*1.0)*L0;
    foreach_point(xp4,yp4) {
      yp4=y;
      xp4=x;
      res4=u.x[];
      res4k=k[];
      res4eps=epsilon[];
      res4mutke=mutke[];
    }
    fprintf(velocity4,"%g,%g,%g,%g,%g,%g\n",xp4,yp4,res4,res4k,res4eps,res4mutke);
    fflush(velocity4);
  }

  double xp5=0.9*L0;
  char name5[1024];
  sprintf(name5,"velocity_profile_x=%g_level=%d.txt",xp5,maxlevel);
  FILE * velocity5=fopen(name5,"w");
  int number5=pow(2,maxlevel);
  double yp5,res5,res5k,res5eps,res5mutke;
  for (int kk=0;kk<=number5;kk++) {
    yp5=-L0/2.+kk*1.0/(number5*1.0)*L0;
    foreach_point(xp5,yp5) {
      yp5=y;
      xp5=x;
      res5=u.x[];
      res5k=k[];
      res5eps=epsilon[];
      res5mutke=mutke[];
    }
    fprintf(velocity5,"%g,%g,%g,%g,%g,%g\n",xp5,yp5,res5,res5k,res5eps,res5mutke);
    fflush(velocity5);
  }

  char name6[1024];
  sprintf(name6,"spreading_rate_level=%d.txt",maxlevel);
  FILE * spreading6=fopen(name6,"w");
  double Xlist[5]={0.1*L0,0.3*L0,0.5*L0,0.7*L0,0.9*L0};
  double xp6;
  double ya0,yb0,ym0;
  double resm0;

  double ya9,yb9,ym9;
  double resm9;

  double ya1,yb1,ym1;
  double resm1;
  for (int ii=0;ii<5;ii++) {
    xp6=Xlist[ii];

    ya0=Y0;

    yb0=Y0+L0;

    while (fabs(yb0-ya0)>1e-5) {
      ym0=(ya0+yb0)/2.;
      resm0=interpolate(u.x,xp6,ym0);
      if (resm0<0.5*U1) {
        ya0=ym0;
      }
      else {
        yb0=ym0;
      }
    }

    ya9=Y0;

    yb9=Y0+L0;

    while (fabs(yb9-ya9)>1e-5) {
      ym9=(ya9+yb9)/2.;
      resm9=interpolate(u.x,xp6,ym9);
      if (resm9<0.9*U1) {
        ya9=ym9;
      }
      else {
        yb9=ym9;
      }
    }

    ya1=Y0;

    yb1=Y0+L0;

    while (fabs(yb1-ya1)>1e-5) {
      ym1=(ya1+yb1)/2.;
      resm1=interpolate(u.x,xp6,ym1);
      if (resm1<0.1*U1) {
        ya1=ym1;
      }
      else {
        yb1=ym1;
      }
    }

    fprintf(spreading6,"%g,%g,%g,%g\n",xp6,ym0,ym1,ym9);
    fflush(spreading6);
  }
}


#if 1
event adapt(i++) {
  adapt_wavelet({u,k,epsilon},(double[]){1e-5,1e-5,1e-5,1e-5},maxlevel,5);
}
#endif

#if 0
event snapshot(t+=5.) {
  char name[80];
  sprintf(name,"dump-t=%g_level=%d",t,maxlevel);
  dump (file=name);
}
#endif



/** 
![Animation of the horizontal velocity](windtunnel/velocity_x.mp4)(width=50%, loop)

We see that we recover the self similar solution for the mixing layer
~~~gnuplot Velocity
set datafile separator(',')
set yrange [-0.25:0.20]
set ylabel "y/x$"
set xlabel "U/U_1"
plot "velocity_profile_x=0.1_level=8.txt" using ($3/18.):($2/$1) title "x=0.1", \
"velocity_profile_x=0.3_level=8.txt" using ($3/18.):($2/$1) title "x=0.3", \
"velocity_profile_x=0.5_level=8.txt" using ($3/18.):($2/$1) title "x=0.5", \
"velocity_profile_x=0.7_level=8.txt" using ($3/18.):($2/$1) title "x=0.7", \
"velocity_profile_x=0.9_level=8.txt" using ($3/18.):($2/$1) title "x=0.9"
~~~

![Animation of the turbulent viscosity](windtunnel/viscosity.mp4)(width=50%, loop)

The turbulent viscosity is to high in the slow layer due to a return at the right boundary. This can probably be corrected with a better outflow boundary condition.
~~~gnuplot Turbulent viscosity
set datafile separator(',')
set xrange [-0.25:0.20]
set yrange [*:*]
set xlabel "y/x"
set ylabel "mu_t/mu"
plot "velocity_profile_x=0.1_level=8.txt" using 2:($6/1.85e-6) title "x=0.1", \
"velocity_profile_x=0.3_level=8.txt" using 2:($6/1.85e-6) title "x=0.3", \
"velocity_profile_x=0.5_level=8.txt" using 2:($6/1.85e-6) title "x=0.5", \
"velocity_profile_x=0.7_level=8.txt" using 2:($6/1.85e-6) title "x=0.7", \
"velocity_profile_x=0.9_level=8.txt" using 2:($6/1.85e-6) title "x=0.9"
~~~

The spreading of the mixing layer is linear, characteristic of the self similar variable
~~~gnuplot Spreading
set datafile separator(',')
set xrange [0.:1.]
set xlabel "x"
set ylabel ""
plot "spreading_rate_level=8.txt" using 1:2 with line title "U/U1=0.5", \
"spreading_rate_level=8.txt" using 1:3 with line title "U/U1=0.1", \
"spreading_rate_level=8.txt" using 1:4 with line title "U/U1=0.9"
~~~

*/