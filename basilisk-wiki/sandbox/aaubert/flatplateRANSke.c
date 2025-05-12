//#include "embed.h"
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
#define BVP4 0

scalar cs[];

#include "RANSkepsilon.h"
//#include "navier-stokes/double-projection.h"
#include "navier-stokes/perfs.h"
#include "view.h"


double Reynolds=5000000.;
int maxlevel=12;

double U0=1.;
double rho_air=1.;
double t0=30.;
face vector muv[];
coord * listcoord;
scalar rhov[];


scalar un[];
scalar pavgg3[];

int main() {
  for (maxlevel=10;maxlevel<11;maxlevel++) {
    size(3.);
    N=1024;
    init_grid(N);
    DT=0.1;
    TOLERANCE=1e-5 [*];
    rho=rhov;
    theta=1.;          //we use gradient limiter
    Deltamin=L0/pow(2.,maxlevel);
    u.x.gradient=minmod2;
    u.y.gradient=minmod2;
    p.gradient=minmod2;
    k.gradient=minmod2;
    epsilon.gradient=minmod2;

    run();
  }
}



event init(t=0) {
  //molvis=U0/Reynolds;
  //ct3=0.;
  leading=1.;
  if (!restore(file="dump")) {  //we can restore from a previous simulation
    foreach() {
      cs[]=1.;
      molvis[]=U0/Reynolds;
      u.x[]=U0;
      un[]=U0;
      k[]=3./2.*sq(1./100.*U0);  //turbulent intensity of 1% at the inlet
      epsilon[]=rho_air*C_mu*sq(k[])/(10.*molvis[]*rho_air);  //mu_t=10mu at the inlet
      mutke[]=rho_air*C_mu*sq(k[])/epsilon[];
      rhov[]=rho_air;
    }
  }
  else {
    N=512;
  }
  u.n[left]=dirichlet(U0);
  k[left]=dirichlet(3./2.*sq(1./100.*U0));
  epsilon[left]=dirichlet(rho_air*C_mu*sq(3./2.*sq(1./100.*U0))/(10.*U0/Reynolds*rho_air));

  p[right]=dirichlet(0.);
  pf[right]=dirichlet(0.);
  u.n[right]=neumann(0.);
  u.t[right]=neumann(0.);
  k[right]=neumann(0.);
  epsilon[right]=neumann(0.);

  u.t[top]=dirichlet(U0);
  u.n[top]=neumann(0.);
  k[top]=neumann(0.);
  epsilon[top]=neumann(0.);

  u.t[bottom]=x<leading?neumann(0.):dirichlet(0.);
  u.n[bottom]=dirichlet(0.);
  k[bottom]=x<leading?neumann(0.):dirichlet(3./2.*sq(1./100.*U0));
  epsilon[bottom]=x<leading?neumann(0.):dirichlet(rho_air*C_mu*sq(3./2.*sq(1./100.*U0))/(10.*U0/Reynolds*rho_air));
  p[bottom]=neumann(0.);
}


#if WALL_FUNCTION
event boundary_condition(i++) {
  if (i>100) {
    u.t[bottom]=(x<leading)?neumann(0.):dirichlet(wall_condition_velocity(x,y,true,mutke));
    u.n[bottom]=(x<leading)?dirichlet(0.):dirichlet(wall_condition_velocity(x,y,false,mutke));
    k[bottom]=(x<leading)?neumann(0.):dirichlet(wall_condition_k(x,y,mutke));
    epsilon[bottom]=(x<leading)?neumann(0.):dirichlet(wall_condition_epsilon(x,y,mutke));
  }
}
#endif

event logfile(i++) {
  stats smutke=statsf(mutke);
  stats sk=statsf(k);
  stats seps=statsf(epsilon);
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
  double nuhat2=interpolate(mutke,1.01514,0.00146484);
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
    if (mutke[]>maxnuhat) {
      xnuhat=x;
      ynuhat=y;
      maxnuhat=mutke[];
    }
  }
  double maxnuhatSA=0.;
  double xnuhatSA=0.;
  double ynuhatSA=0.;
  foreach() {
    if (mutke[]>maxnuhatSA) {
      xnuhatSA=x;
      ynuhatSA=y;
      maxnuhatSA=mutke[];
    }
  }
  //double pp1=interpolate(pavg,1.01514,0.00146484);
  //double pp3=interpolate(pavg3,1.01514,0.00146484);
  double utauex=0.;
  foreach_point(1.01514,0.00146484) {
    utauex=u.x[];
  }
  fprintf(stderr,"%g %d %d %g %g %d %d %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g\n",Reynolds,maxlevel,i,t,dt,mgp.i,mgu.i,smutke.min,smutke.max,du,uu,sk.min,sk.max,seps.min,seps.max,xmax,ymax,changemax,pp,uu2,nuhat2,pp3,maxnuhat,xnuhat,ynuhat,utauex,maxnuhatSA,xnuhatSA,ynuhatSA,changemaxplate,xmaxplate,ymaxplate);
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
  output_ppm (mutke, file = "viscosity.mp4",linear = true);
  output_ppm (k, file = "k.mp4",linear = false);
  output_ppm (epsilon, file = "epsilon.mp4",linear = false);
  output_ppm (production,file="production.mp4",linear=false);

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
  //isoline("mutke",n=21);
  squares("mutke");
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
  output_ppm (mutke, file = "viscosity.png",linear = true);

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
      muhat2=rho[]*mutke[];
      xp2=x-1.;
      yp2=y;
    }
    res2=muhat2;
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
    muhat3=interpolate(mutke,1.+xp3,0.);
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
      muhat4=interpolate(mutke,xp4,yp4);
      res4=muhat4;
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
  adapt_wavelet({u,mutke,nn},(double[]){1e-3,1e-3,1e-3,1e-5},maxlevel,5);
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
We compare this result with the one obtain using the Spalart-Allmaras model in Basilisk and in the NASA solver CFL3D

~~~pythonplot Velocity profile
import numpy as np
import matplotlib.pyplot as plt

level=10

Xpos=[1]
Ypos=[]
Upos=[]
Ysa=[]
Usa=[]
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
    
    source=open("../flatplateRANSSA/velocity_profile_x="+str(x)+"_Re=5e+06_level="+str(level)+".txt",'r')

    Y=[]
    U=[]
    for ligne in source:
        L=ligne.strip().split(',')
        if len(Y)<=100:
            Y.append(float(L[1]))
            U.append(float(L[2]))

    source.close()
    Ysa.append(Y)
    Usa.append(U)

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
   plt.scatter(Ypos[i],Upos[i],label="level="+str(level),marker='+',color='red')
   plt.scatter(Ysa[i],Usa[i],label="level="+str(level)+" SA",marker='o',color='blue')
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

source=open("../flatplateRANSSA/ratio_turbulent_viscosity_x=1_Re=5e+06_level="+str(level)+".txt",'r')
Ysa=[]
Nuhatsa=[]
for ligne in source:
   L=ligne.strip().split(',')
   if len(Ysa)<=100:
        Ysa.append(float(L[1]))
        Nuhatsa.append(float(L[2]))
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
plt.scatter(Nuhat,Y,label="level="+str(level),marker='+',color="red")
plt.scatter(Nuhatsa,Ysa,label="level="+str(level)+" SA",marker='o',color="blue")
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

source=open("../flatplateRANSSA/skin_friction_Re=5e+06_level="+str(level)+".txt",'r')

Xsa=[]
Cfsa=[]
for ligne in source:
    L=ligne.strip().split(',')
    Xsa.append(float(L[0]))
    Cfsa.append(float(L[1]))

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
plt.scatter(X,Cf,label="level="+str(level),marker='+',color='red')
plt.scatter(Xsa,Cfsa,label="level="+str(level)+" SA",marker='o',color='blue')
plt.plot(Xt,Cft,label="NASA CFL3D",linestyle='--',color="k",linewidth=2)
plt.ylim([0.002,0.006])
plt.xlabel("$x$",fontsize=15)
#plt.ylabel("Skin friction",fontsize=15)
plt.title("Skin friction coefficient along the plate",fontsize=20)
plt.legend()
plt.savefig('Skin friction.png')
~~~

*/