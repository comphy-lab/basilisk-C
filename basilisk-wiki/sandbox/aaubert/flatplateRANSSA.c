/**
Test case of a flat plate to validate the implementation of the Spalart-Allmaras model */

#include "embed.h"
#include "navier-stokes/centered.h"
#include "distance.h"

#define WALL 1             //wall is present : need to define a distance_to_wall function
double distance_to_wall(double x,double y) {
  if (x>1.) {
    return y;
  }
  else {
    return sqrt(pow(x-1.,2.)+pow(y,2.));
  }
}

#define WALL_FUNCTION 1
#define LINEARISATION 1
#define BVP4 0  //fixme : not working yet

//scalar cs[];

#include "RANSSpalartAllmaras.h"
//#include "navier-stokes/double-projection.h"
#include "navier-stokes/perfs.h"
#include "view.h"


double Reynolds=5000000.;
int maxlevel=12;

double U0=1.;
double t0=30.;
face vector muv[];
coord * listcoord;


scalar un[];
scalar pavgg3[];

int main() {
  for (maxlevel=10;maxlevel<11;maxlevel++) {
    size(3.);
    N=1024;
    init_grid(N);
    TOLERANCE=1e-5 [*];
    theta=1.;
    Deltamin=L0/pow(2.,maxlevel);
    //u.x.gradient=minmod2;
    //u.y.gradient=minmod2;
    //p.gradient=minmod2;
    //muhat.gradient=minmod2;
    run();
  }
}



event init(t=0) {
  //molvis=U0/Reynolds;
  //ct3=0.;
  leading=1.;
  trailing=3.;
  if (!restore(file="dump")) {  //we can restore from a previous simulation
    foreach() {
      cs[]=1.;
      molvis[]=U0/Reynolds;
      u.x[]=U0;
      un[]=U0;
      muhat[]=3*molvis[];
      muhatSA[]=3.*molvis[];
    }
  }
  else {
    N=512;
  }
  u.n[left]=dirichlet(U0);
  //muhat[left]=dirichlet(3.*molvis);

  p[right]=dirichlet(0.);
  pf[right]=dirichlet(0.);
  u.n[right]=neumann(0.);
  u.t[right]=neumann(0.);

  u.t[top]=dirichlet(U0);
  u.n[top]=neumann(0.);
  muhat[top]=dirichlet(3.*U0/Reynolds);
  muhatSA[top]=dirichlet(3.*U0/Reynolds);

  u.t[bottom]=x<leading?neumann(0.):dirichlet(0.);
  u.n[bottom]=dirichlet(0.);
  muhat[bottom]=x<leading?neumann(0.):dirichlet(0.);
  muhatSA[bottom]=x<leading?neumann(0.):dirichlet(0.);
  p[bottom]=neumann(0.);
}


#if WALL_FUNCTION
event boundary_condition(i++) {
  if (i>100) {
    u.t[bottom]=(x<leading)?neumann(0.):dirichlet(wall_condition_velocity(x,y,true,muhat));
    u.n[bottom]=(x<leading)?dirichlet(0.):dirichlet(wall_condition_velocity(x,y,false,muhat));
    muhat[bottom]=(x<leading)?neumann(0.):dirichlet(wall_condition_viscosity(x,y,muhat));
    muhatSA[bottom]=(x<leading)?neumann(0.):dirichlet(wall_condition_viscosity_SA(x,y,muhat));
  }
}
#endif

event logfile(i++) {
  stats smuhat=statsf(muhat);
  double xmax=0.;
  double ymax=0.;
  double changemax=0.;
  double xmaxplate=0.;
  double ymaxplate=0.;
  double changemaxplate=0.;
  foreach() {
    if (fabs(u.x[]-un[])>changemax) {
      xmax=x;
      ymax=y;
      changemax=fabs(u.x[]-un[]);
    }
    if (y<1.) {
      if (fabs(u.x[]-un[])>changemaxplate) {
	xmaxplate=x;
	ymaxplate=y;
	changemaxplate=fabs(u.x[]-un[]);
      }
    }
  }
  double du=change(u.x,un);
  double uu=interpolate(u.x,2.,0.001);
  double pp=interpolate(p,1.01514,0.00146484);
  double uu2=interpolate(u.x,1.01514,0.00146484);
  double nuhat2=interpolate(muhat,1.01514,0.00146484);
  if (i==1000) {
    foreach() {
      pavgg3[]=p[];
    }
  }
  else {
    foreach() {
      pavgg3[]=0.9*pavgg3[]+0.1*p[];
    }
  }
  double pp3;
  if (i<1000) {
    pp3=0.;
  }
  else {
    pp3=interpolate(pavgg3,1.01514,0.00146484);
  }
  double maxnuhat=0.;
  double xnuhat=0.;
  double ynuhat=0.;
  foreach() {
    if (muhat[]>maxnuhat) {
      xnuhat=x;
      ynuhat=y;
      maxnuhat=muhat[];
    }
  }
  double maxnuhatSA=0.;
  double xnuhatSA=0.;
  double ynuhatSA=0.;
  foreach() {
    if (muhatSA[]>maxnuhatSA) {
      xnuhatSA=x;
      ynuhatSA=y;
      maxnuhatSA=muhatSA[];
    }
  }
  //double pp1=interpolate(pavg,1.01514,0.00146484);
  //double pp3=interpolate(pavg3,1.01514,0.00146484);
  double utauex=0.;
  foreach_point(1.01514,0.00146484) {
    utauex=muhat.utau[];
  }
  fprintf(stderr,"%g %d %d %g %g %d %d %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g\n",Reynolds,maxlevel,i,t,dt,mgp.i,mgu.i,smuhat.min,smuhat.max,du,uu,xmax,ymax,changemax,pp,uu2,nuhat2,pp3,maxnuhat,xnuhat,ynuhat,utauex,maxnuhatSA,xnuhatSA,ynuhatSA,changemaxplate,xmaxplate,ymaxplate);
  fflush(stderr);
}

#if 1
event movies (t+=0.1;t<=5.) {

  scalar omega[];
  vorticity(u,omega);

  output_ppm (omega, file = "vort.mp4",linear = true);
  output_ppm (p, file = "pressure.mp4",linear = true);
  output_ppm (u.x, file = "velocity_x.mp4",max=U0,min=0.,linear = true);
  output_ppm (u.y, file = "velocity_y.mp4",max=U0,min=0.,linear = true);
  output_ppm (muhat, file = "viscosity.mp4",linear = true);

  view(tx=-0.5,ty=-0.,height=200,width=1000,fov=5,camera="front");
  clear();
  //isoline("u.x",n=21);
  squares("u.x",max=U0,min=-U0);
  box();
  char name2[1024];
  sprintf(name2,"isoline_Re=%g_level=%d.mp4",Reynolds,maxlevel);
  save(name2);

  view(tx=-0.5,ty=-0.,height=200,width=1000,fov=5,camera="front");
  clear();
  //isoline("muhat",n=21);
  squares("muhat");
  box();
  char name3[1024];
  sprintf(name3,"isoviscosity_Re=%g_level=%d.mp4",Reynolds,maxlevel);
  save(name3);

  view(tx=-0.5,ty=-0.,height=200,width=1000,fov=5,camera="front");
  clear();
  //isoline("p",n=21);
  squares("p");
  box();
  char name4[1024];
  sprintf(name4,"isopressure_Re=%g_level=%d.mp4",Reynolds,maxlevel);
  save(name4);
}
#endif

double center_gradient_2(double s0,double s1,double s2) {
  return s1-s0;
}

#if 1
event images (t=5.) {

  scalar omega[];
  vorticity(u,omega);

  output_ppm (omega, file = "vort.png",linear = true);
  output_ppm (u.x, file = "velocity_x.png",max=U0,min=0.,linear = true);
  output_ppm (u.y, file = "velocity_y.png",max=U0,min=0.,linear = true);
  output_ppm (muhat, file = "viscosity.png",linear = true);

  view(tx=-0.5,ty=-0.,height=200,width=1000,fov=5,camera="front");
  clear();
  //isoline("u.x",n=21);
  squares("u.x",max=U0,min=-U0);
  box();
  char namevid[1024];
  sprintf(namevid,"isoline_Re=%g_level=%d.png",Reynolds,maxlevel);
  save(namevid);

  view(tx=-0.5,ty=-0.,height=200,width=1000,fov=5,camera="front");
  clear();
  //isoline("u.x",n=21);
  squares("p",max=U0,min=-U0);
  box();
  char namevid2[1024];
  sprintf(namevid2,"isopressure_Re=%g_level=%d.png",Reynolds,maxlevel);
  save(namevid2);

  face vector gu[],gv[];
  u.x.gradient=center_gradient_2;
  u.y.gradient=center_gradient_2;
  gradients ({u.x,u.y},{gu,gv});


  double xp1=1.0;
  double xp1real;
  foreach_point(1.+xp1,0.) {
    xp1real=x;
  }
#if WALL_FUNCTION
  double ut1=interpolate(u.x,xp1real,d_IP);
  double utauprec1;
  double utau1=0.1;
  double err1=1.;
  int n_iter1=0;
  while ((n_iter1<50)&&(err1>1e-3)) {
    utauprec1=utau1;
    utau1=utau1-gwallSA(utau1,d_IP,ut1,U0/Reynolds)/gwallSAprime(utau1,d_IP,ut1,U0/Reynolds);
    err1=fabs(utau1-utauprec1)/(utau1+1e-10);
    n_iter1+=1;
  }
#else
  double tau_1=U0/Reynolds*interpolate(gu.y,xp1real,0);
  double utau1=sqrt(tau_1);
#endif
  char name1[1024];
  sprintf(name1,"velocity_profile_x=%g_Re=%g_level=%d.txt",xp1,Reynolds,maxlevel);
  FILE * velocity1=fopen(name1,"w");
  int number1=100;
  double yp1;
  for (int kk=0;kk<=number1;kk++) {
    yp1=exp(kk*log(0.5)/number1+(1-kk*1.0/number1)*log(1e-9));
    foreach_point(1+xp1,yp1) {
      fprintf(velocity1,"%g,%g,%g,%g\n",x-1.,y*utau1/(U0/Reynolds),u.x[]/utau1,utau1);
      fflush(velocity1);
    }
  }

  double xp2=1.0;
  char name2[1024];
  sprintf(name2,"ratio_turbulent_viscosity_x=%g_Re=%g_level=%d.txt",xp2,Reynolds,maxlevel);
  FILE * ratio_visc_2=fopen(name2,"w");
  int number2=100;
  double yp2,muhat2,res2;
  for (int kk=0;kk<=number2;kk++) {
    yp2=kk*0.025/number2;
    foreach_point(1.+xp2,yp2) {
      muhat2=muhat[];
      xp2=x-1.;
      yp2=y;
    }
    res2=muhat2*fv1(muhat2/(U0/Reynolds));
    fprintf(ratio_visc_2,"%g,%g,%g\n",xp2,yp2,res2/(U0/Reynolds));
    fflush(ratio_visc_2);
  }

  char name3[1024];
  sprintf(name3,"skin_friction_Re=%g_level=%d.txt",Reynolds,maxlevel);
  FILE * skin_friction3=fopen(name3,"w");
  int number3=100;
  double xp3,res3,muhat3;
  for (int kk=0;kk<=number3;kk++) {
    xp3=kk*2.0/number3;
    foreach_point(1.+xp3,0.) {
      xp3=x-1.;
    }
#if WALL_FUNCTION
    double ut3=interpolate(u.x,1.+xp3,d_IP);
    double utauprec3;
    double utau3=0.1;
    double err3=1.;
    int n_iter3=0;
    while ((n_iter3<50)&&(err3>1e-3)) {
      utauprec3=utau3;
      utau3=utau3-gwallSA(utau3,d_IP,ut3,U0/Reynolds)/gwallSAprime(utau3,d_IP,ut3,U0/Reynolds);
      err3=fabs(utau3-utauprec3)/(utau3+1e-10);
      n_iter3+=1;
    }
    res3=sq(utau3);
#else
    muhat3=interpolate(muhat,1.+xp3,0.);
    res3=(U0/Reynolds+muhat3*fv1(muhat3/(U0/Reynolds)))*interpolate(gu.y,1.+xp3,0);
#endif
    fprintf(skin_friction3,"%g,%g\n",xp3,res3*2.);
    fflush(skin_friction3);
  }

  char name4[1024];
  sprintf(name4,"turbulent_viscosity_Re=%g_level=%d.txt",Reynolds,maxlevel);
  FILE * turbulent_viscosity4=fopen(name4,"w");
  int number4x=100;
  int number4y=50;
  double xp4,yp4,muhat4,res4;
  for (int kkx=0;kkx<=number4x;kkx++) {
    for (int kky=0;kky<=number4y;kky++) {
      xp4=kkx*3.0/number4x;
      yp4=kky*0.05/number4y;
      muhat4=interpolate(muhat,xp4,yp4);
      res4=muhat4*fv1(muhat4/(U0/Reynolds));
      fprintf(turbulent_viscosity4,"%g,%g,%g\n",xp4,yp4,res4);
      fflush(turbulent_viscosity4);
    }
  }

}
#endif

#if 1
event adapt(i++) {
  scalar nn[];
  foreach() {
    double d=distance_to_wall(x,y);
    if (d<2.*d_IP) {
      nn[]=noise();
    }
  }
  adapt_wavelet({u,muhat,nn},(double[]){1e-3,1e-3,1e-3,1e-5},maxlevel,5);
}
#endif

#if 0
event snapshot(t+=1.) {
  char name[80];
  sprintf(name,"dump-t=%g_level=%d",t,maxlevel);
  dump (file=name);
}
#endif

/**
We can plot the velocity profile and compare it with the profile obtain using CFL3D (dash line).


~~~pythonplot Velocity profile
import numpy as np
import matplotlib.pyplot as plt

level=10

Xpos=[1]
Ypos=[]
Upos=[]
Yt=[]
Ut=[]
for x in Xpos:

    source=open("velocity_profile_x="+str(x)+"_Re=5e+06_level="+str(level)+".txt",'r')

    Y=[]
    U=[]
    for ligne in source:
        L=ligne.strip().split(',')
        if len(Y)<=100:
            Y.append(float(L[1]))
            U.append(float(L[2]))

    source.close()
    Ypos.append(Y)
    Upos.append(U)

    Y=[]
    U=[]
    source=open("../data/Flat_plate/Data_y+u+_NASA_CFL3D.txt",'r')
    
    for ligne in source:
        L=ligne.strip().split(' ')
        if float(L[0])>5:
            break
        Y.append(float(L[0]))
        U.append(float(L[1]))

    source.close()
    Yt.append(Y)
    Ut.append(U)

plt.figure()
for i in range(len(Xpos)):
   plt.scatter(Ypos[i],Upos[i],label="level="+str(level),marker='+')
   plt.plot(10**np.array(Yt[i]),Ut[i],label="NASA CFL3D",linestyle="--",color='k',linewidth=2)
plt.xlabel("$y^+$",fontsize=15)
plt.ylabel("$u^+$",fontsize=15)
plt.xscale('log')
plt.title("Velocity profile for level="+str(level),fontsize=20)
plt.legend()
plt.savefig('velocity.png')
~~~

We can also look at the turbulent viscosity along the plate.

~~~pythonplot Turbulent viscosity
level=10

source=open("ratio_turbulent_viscosity_x=1_Re=5e+06_level="+str(level)+".txt",'r')
Y=[]
Nuhat=[]
for ligne in source:
   L=ligne.strip().split(',')
   if len(Y)<=100:
        Y.append(float(L[1]))
        Nuhat.append(float(L[2]))
source.close()

Yt=[]
Nuhatt=[]
source=open("../data/Flat_plate/Data_ratio_viscosity_turbulent_NASA_CFL3D.txt",'r')
for ligne in source:
    L=ligne.strip().split(' ')
    if float(L[1])>0.025:
       break
    Yt.append(float(L[1]))
    Nuhatt.append(float(L[2]))

source.close()

plt.figure()
plt.scatter(Nuhat,Y,label="level="+str(level),marker='+')
plt.plot(Nuhatt,Yt,label="NASA CFL3D",linestyle='--',color="k",linewidth=2)
plt.xlabel("Ratio turbulent viscosity",fontsize=15)
plt.ylabel("$y$",fontsize=15)
plt.title("Profile of turbulent viscosity at $x=1$",fontsize=20)
plt.legend()
plt.ylim([0.,0.025])
plt.savefig('Turbulent_viscosity.png')
~~~

Finally, we can see the friction coefficient along the flat plate.

~~~pythonplot Skin Friction coefficient

level=10
source=open("skin_friction_Re=5e+06_level="+str(level)+".txt",'r')

X=[]
Cf=[]
for ligne in source:
    L=ligne.strip().split(',')
    X.append(float(L[0]))
    Cf.append(float(L[1]))

source.close()


Xt=[]
Cft=[]
source=open("../data/Flat_plate/Data_skin_friction_coefficient_NASA_CFL3D.txt",'r')

for ligne in source:
    L=ligne.strip().split(' ')
    Xt.append(float(L[0]))
    i=1
    while L[i]=='':
       i+=1
    Cft.append(float(L[i]))

source.close()

plt.figure()
plt.scatter(X,Cf,label="level="+str(level),marker='+')
plt.plot(Xt,Cft,label="NASA CFL3D",linestyle='--',color="k",linewidth=2)
plt.ylim([0.002,0.006])
plt.xlabel("$x$",fontsize=15)
#plt.ylabel("Skin friction",fontsize=15)
plt.title("Skin friction coefficient along the plate",fontsize=20)
plt.legend()
plt.savefig('Skin friction.png')
~~~

*/