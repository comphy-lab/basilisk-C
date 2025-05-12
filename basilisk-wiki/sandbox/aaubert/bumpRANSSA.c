#include "embed.h"
#include "navier-stokes/centered.h"
#include "redistance.h"

#define WALL 1             //wall is present : need to define a distance_to_wall function

scalar d_wall[];

double distance_to_wall(double x1,double y1) {
  double value;
  if (y1<0.) {
    return y1;
  }
  foreach_point(x1,y1) {
    if (x<0.) {
      value=min(d_wall[],sqrt(sq(x-0.)+sq(y)));
    }
    else if (x>1.5) {
      value=min(d_wall[],sqrt(sq(x-1.5)+sq(y)));
    }
    else {
      value=min(d_wall[],y);
    }
  }
  return value;
}

#define WALL_FUNCTION 1
#define LINEARISATION 1
#define BVP4 0  //fixme : not working yet

#include "RANSSpalartAllmaras.h"
#include "navier-stokes/double-projection.h"
#include "navier-stokes/perfs.h"
#include "view.h"


double Reynolds=3000000.;
int maxlevel=12;

double U0=1.;
face vector muv[];

scalar un[];

double def_solid(double x, double y) {
  if (x<0.3) {
    return y;//-x+0.3+y;
  }
  if (x>1.2) {
    return y;//-1.2+x+y;
  }
  else {
    double l=0.9;
    double a=0.05;
    return -a*pow(sin(pi*x/l-pi/3.),4.)+y-1e-10;
  }
}

int main() {
  for (maxlevel=10;maxlevel<11;maxlevel++) {
    size(5.);
    DT=0.1;
    origin(-1.75,0.);
    N=512;
    init_grid(N);
    TOLERANCE=1e-7 [*];
    Deltamin=L0/pow(2.,maxlevel);
    //theta=1.3;
    //u.x.gradient=minmod2;
    //u.y.gradient=minmod2;
    //p.gradient=minmod2;
    //muhatSA.gradient=minmod2;
    run();
  }
}


event init(t=0) {
  //CFL=0.1;
  //ct3=0.;
  if (!restore(file="dump")) {  //we can restore from a previous simulation
    refine ((y<0.06) && (x>-0.1) && (x<1.6) && (level<maxlevel));
    solid (cs,fs,def_solid(x,y));
    foreach() {
      d_wall[]=def_solid(x,y);
    }
    redistance(d_wall);
    foreach() {
      if (x<0.1) {
        cs[]=1.;    //inconsistancy in the computation of cs
      }
      molvis[]=U0/Reynolds;
      u.x[]=cs[]? U0 :0.;
      un[]=U0;
      muhat[]=cs[]?3.*molvis[]:0.;
      muhatSA[]=cs[]?3.*molvis[]:0.;
      muhatSAdiss[]=cs[]?3.*molvis[]:0.;
    }
  }
  else {
    N=512;
  }

  foreach() {
    molvis[]=U0/Reynolds;
  }
  leading=0.;
  trailing=1.5;
  
  u.n[left]=dirichlet(U0);
  muhat[left]=dirichlet(3.*U0/Reynolds);
  muhatSA[left]=dirichlet(3.*U0/Reynolds);
  
  p[right]=dirichlet(0.);
  pf[right]=dirichlet(0.);
  u.n[right]=neumann(0.);
  u.t[right]=neumann(0.);
  
  p[top]=neumann(0.);
  u.t[top]=neumann(0.);
  u.n[top]=dirichlet(0.);
  muhat[top]=neumann(0.);
  muhatSA[top]=neumann(0.);
  
  u.t[bottom]=((x<0.)||(x>1.5))?neumann(0.):dirichlet(0.);
  u.n[bottom]=dirichlet(0.);
  muhat[bottom]=((x<0.)||(x>1.5))?neumann(0.):dirichlet(0.);
  muhatSA[bottom]=((x<0.)||(x>1.5))?neumann(0.):dirichlet(0.);
  p[bottom]=neumann(0.);
  
  
  u.n[embed]=dirichlet(0.);
  u.t[embed]=dirichlet(0.);
  muhat[embed]=dirichlet(0.);
  muhatSA[embed]=dirichlet(0.);
  p[embed]=neumann(0.);

}

#if WALL_FUNCTION
event boundary_condition(i++) {
  if (i>100) {
    u.t[bottom]=((x<0.)||(x>1.5))?neumann(0.):dirichlet(wall_condition_velocity(x,y,true,muhat));
    muhat[bottom]=((x<0.)||(x>1.5))?neumann(0.):dirichlet(wall_condition_viscosity(x,y,muhat));
    muhatSA[bottom]=((x<0.)||(x>1.5))?neumann(0.):dirichlet(wall_condition_viscosity_SA(x,y,muhat));
    muhatSAdiss[bottom]=((x<0.)||(x>1.5))?neumann(0.):dirichlet(wall_condition_viscosity_SA_diss(x,y,muhat));
  }
}
#endif


event logfile(i++) {
  stats smuhat=statsf(muhat);
  double xmax,ymax;
  double changemax=0.;
  foreach() {
    if (fabs(u.x[]-un[])>changemax) {
      changemax=fabs(u.x[]-un[]);
      xmax=x;
      ymax=y;
    }
  }
  double du=change(u.x,un);
  double uu=interpolate(u.x,2.,0.001);
  stats ssource=statsf(source);
 
  fprintf(stderr,"%g %d %d %g %g %d %d %g %g %g %g %g %g %g %g %g %g %g %g\n",Reynolds,maxlevel,i,t,dt,mgp.i,mgu.i,smuhat.min,smuhat.max,du,uu,ssource.min,ssource.max,statssolver.niteravg,statssolver.norm2avg,statssolver.norm2max,changemax,xmax,ymax);
  fflush(stderr);


}

#if 1
event movies (t+=0.1;t<=5.) {

  scalar omega[];
  vorticity(u,omega);

  scalar m2[];
  foreach() {
    m2[]=cs[]-0.5;
  }

  output_ppm (omega, file = "vort.mp4",linear = true,mask=m2);
  output_ppm (p, file = "pressure.mp4",linear = true,mask=m2);
  output_ppm (u.x, file = "velocity_x.mp4",max=U0,min=0.,linear = true,mask=m2);
  output_ppm (u.y, file = "velocity_y.mp4",max=U0,min=0.,linear = true,mask=m2);
  output_ppm (muhat, file = "viscosity.mp4",linear = true,mask=m2);

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

  if ((i>1500)&&(i<1700)) {
    view(tx=-0.5,ty=-0.,height=200,width=1000,fov=5,camera="front");
    clear();
    //isoline("muhat",n=21);
    squares("muhat");
    box();
    char name5[1024];
    sprintf(name5,"test_isoviscosity_Re=%g_level=%d.mp4",Reynolds,maxlevel);
    save(name5);

    view(tx=-0.5,ty=-0.,height=200,width=1000,fov=5,camera="front");
    clear();
    //isoline("muhat",n=21);
    squares("u.x");
    box();
    char name6[1024];
    sprintf(name6,"test_isoline_Re=%g_level=%d.mp4",Reynolds,maxlevel);
    save(name6);

    view(tx=-0.5,ty=-0.,height=200,width=1000,fov=5,camera="front");
    clear();
    //isoline("muhat",n=21);
    squares("p");
    box();
    char name7[1024];
    sprintf(name7,"test_isopressure_Re=%g_level=%d.mp4",Reynolds,maxlevel);
    save(name7);

    output_ppm (u.x, file = "test_velocity_x.mp4",max=U0,min=0.,linear = true,mask=m2);
    output_ppm (u.y, file = "test_velocity_y.mp4",max=U0,min=0.,linear = true,mask=m2);
    output_ppm (muhat, file = "test_viscosity.mp4",linear = true,mask=m2);
  }
}
#endif

double center_gradient_2(double s0,double s1,double s2) {
  return s1-s0;
}

event images (t=5.) {

  scalar omega[];
  vorticity(u,omega);

  scalar m2[];
  foreach() {
    m2[]=cs[]-0.5;
  }

  double molvis2=U0/Reynolds;
  double ddelta=2.*Deltamin;

  output_ppm (omega, file = "vort.png",linear = true,mask=m2);
  output_ppm (u.x, file = "velocity_x.png",max=U0,min=0.,linear = true,mask=m2);
  output_ppm (u.y, file = "velocity_y.png",max=U0,min=0.,linear = true,mask=m2);
  output_ppm (muhat, file = "viscosity.png",linear = true,mask=m2);

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

  double xp1=0.75;
  double a=0.05;
  double l=0.9;
  double b=1.;
  double yp01=a*pow(sin(pi*xp1/l-pi/3.),4.);
  char name1[1024];
  sprintf(name1,"velocity_profile_x=%g_Re=%g_level=%d.txt",xp1,Reynolds,maxlevel);
  FILE * velocity1=fopen(name1,"w");
  int number1=100;
  double yp1,res1;
  for (int kk=0;kk<=number1;kk++) {
    yp1=yp01+b*kk*0.04/number1;
    foreach_point(xp1,yp1) {
      yp1=y;
      xp1=x;
    }
    res1=interpolate(u.x,xp1,yp1);
    fprintf(velocity1,"%g,%g,%g\n",xp1,yp1-yp01,res1);
    fflush(velocity1);
  }

  double xp2=1.20148;
  char name2[1024];
  sprintf(name2,"velocity_profile_x=%g_Re=%g_level=%d.txt",xp2,Reynolds,maxlevel);
  FILE * velocity2=fopen(name2,"w");
  int number2=100;
  double yp2,res2;
  for (int kk=0;kk<=number2;kk++) {
    yp2=b*kk*0.04/number2;
    foreach_point(xp2,yp2) {
      yp2=y;
      xp2=x;
    }
    res2=interpolate(u.x,xp2,yp2);
    fprintf(velocity2,"%g,%g,%g\n",xp2,yp2,res2);
    fflush(velocity2);
  }


  double xp3=0.75;
  char name3[1024];
  sprintf(name3,"ratio_turbulent_viscosity_x=%g_Re=%g_level=%d.txt",xp3,Reynolds,maxlevel);
  FILE * ratio_visc_3=fopen(name3,"w");
  int number3=100;
  double yp3,muhat3,res3,d3;
  for (int kk=0;kk<=number3;kk++) {
    yp3=0.05+b*kk*0.015/number3;
    foreach_point(xp3,yp3) {
      yp3=y;
      xp3=x;
      muhat3=muhat[]/molvis[];
    }
    d3=distance_to_wall(xp3,yp3);
    res3=muhat3*fv1(muhat3);
    fprintf(ratio_visc_3,"%g,%g,%g\n",xp3,yp3,res3);
    fflush(ratio_visc_3);
  }

  char name4[1024];
  sprintf(name4,"skin_friction_Re=%g_level=%d.txt",Reynolds,maxlevel);
  FILE * skin_friction4=fopen(name4,"w");
  int number4=100;
  double xp4,yp4,res4,d4;
  for (int kk=0;kk<=number4;kk++) {
    xp4=b*kk*1.5/number4;
#if WALL_FUNCTION
    if ((xp4<=0.3)||(xp4>=1.2)) {
      yp4=0.;
    }
    else {
      yp4=a*pow(sin(pi*xp4/l-pi/3.),4.);
    }
    foreach_point(xp4,yp4) {
      double d4=distance_to_wall(x,y);
      double utau4=0.1;

      double dx=(distance_to_wall(x+Delta,y)-distance_to_wall(x-Delta,y))/(2.*Delta);
      double dy=(distance_to_wall(x,y+Delta)-distance_to_wall(x,y-Delta))/(2.*Delta);
	  
      double nx=dx/sqrt(sq(dx)+sq(dy)+1e-10);
      double ny=dy/sqrt(sq(dx)+sq(dy)+1e-10);
	  
      double alpha4=ddelta-d4;
      double x1=x+alpha4*nx;
      double y1=y+alpha4*ny;
      
      double u1x=interpolate(u.x,x1,y1);
      double u1y=interpolate(u.y,x1,y1);
      
      double ut1=ny*u1x-nx*u1y;

      if (ut1<0.) {
	ut1=-ut1;
	ny=-ny;
	nx=-nx;
      }
      
      double utauprec;
      double err=1.;
      int n_iter=0;
      while ((n_iter<50)&&(err>1e-3)) {
	utauprec=utau4;
	utau4=utau4-gwallSA(utau4,ddelta,ut1,molvis[])/gwallSAprime(utau4,ddelta,ut1,molvis[]);
	err=fabs(utau4-utauprec)/(fabs(utau4)+1e-10);
	n_iter+=1;
      }
      fprintf(skin_friction4,"%g,%g\n",x,2.*sq(utau4));
      fflush(skin_friction4);
    }
#else
    if ((xp4<=0.3)||(xp4>=1.2)) {
      yp4=0.;
      double muhat4;
      muhat0=interpolate(muhat,xp4,yp4);
      d4=distance_to_wall(xp4,yp4);
      res4=(molvis2+muhat4*fv1(muhat4/molvis2))*interpolate(gu.y,xp4,yp4);
    }
    else {
      yp4=a*pow(sin(pi*xp4/l-pi/3.),4.);
      coord b,n;
      coord dudn;
      double muhat4;
      double muturbulent=0.;
      foreach_point(xp4,yp4) {
	embed_geometry(point,&b,&n);
	dudn=embed_gradient(point,u,b,n);
	d4=distance_to_wall(x,y);
	//muturbulent=muhat[]*fv1(i,muhat[]/molvis,ddelta/d6);
	muhat4=embed_interpolate(point,muhat,b);
	muturbulent=muhat4*fv1(muhat4/molvis[]);
      }
      res4=(molvis2+muturbulent)*(n.y*dudn.x-n.x*dudn.y);
    } 
    fprintf(skin_friction4,"%g,%g\n",xp4,res4*2.);
    fflush(skin_friction4);
#endif
  }

  char name5[1024];
  sprintf(name5,"turbulent_viscosity_Re=%g_level=%d.txt",Reynolds,maxlevel);
  FILE * turbulent_viscosity5=fopen(name5,"w");
  int number5x=100;
  int number5y=50;
  double xp5,yp5,muhat5,res5,d5;
  for (int kkx=0;kkx<=number5x;kkx++) {
    for (int kky=0;kky<=number5y;kky++) {
      xp5=b*kkx*1.5/number5x;
      yp5=b*kky*0.08/number5y;
      muhat5=interpolate(muhat,xp5,yp5);
      d5=distance_to_wall(xp5,yp5);
      res5=muhat5*fv1(muhat5/molvis2);
      if ((xp5<=0.3)||(xp5>=1.2)) {
	fprintf(turbulent_viscosity5,"%g,%g,%g\n",xp5,yp5,res5);
	fflush(turbulent_viscosity5);
      }
      else if (yp5>=a*pow(sin(pi*xp5/l-pi/3.),4.)) {
	fprintf(turbulent_viscosity5,"%g,%g,%g\n",xp5,yp5,res5);
	fflush(turbulent_viscosity5);
      }
      else {
	fprintf(turbulent_viscosity5,"%g,%g,%g\n",xp5,yp5,-1.);
	fflush(turbulent_viscosity5);
      }
    }
  }
  
  char name6[1024];
  sprintf(name6,"skin_pressure_Re=%g_level=%d.txt",Reynolds,maxlevel);
  FILE * skin_pressure6=fopen(name6,"w");
  int number6=1100;
  double xp6,yp6,res6;
  for (int kk=0;kk<=number6;kk++) {
    xp6=b*kk*1.5/number6;
    if ((xp6<=0.3)||(xp6>=1.2)) {
      yp6=0.;
      foreach_point(xp6,yp6) {
	res6=p[];
	xp6=x;
      }
    }
    else {
      yp6=a*pow(sin(pi*xp6/l-pi/3.),4.);
      foreach_point(xp6,yp6) {
      	coord n,b;
      	embed_geometry(point,&b,&n);
      	res6=embed_interpolate(point,p,b);
	xp6=x;
      }
    }
    fprintf(skin_pressure6,"%g,%g,%g\n",xp6,yp6,res6);
    fflush(skin_pressure6);
  }

  char name7[1024];
  sprintf(name7,"utau_Re=%g_level=%d.txt",Reynolds,maxlevel);
  FILE * utaufile7=fopen(name7,"w");
  int number7=100;
  foreach() {
    if ((x>=0.)&&(x<=1.5)) {
      if (cs[]<1.) {
	double d7=distance_to_wall(x,y);
	if ((d7>=0)&&(d7<Delta)) {
	  double utau7=0.1;

	  double dx=(distance_to_wall(x+Delta,y)-distance_to_wall(x-Delta,y))/(2.*Delta);
	  double dy=(distance_to_wall(x,y+Delta)-distance_to_wall(x,y-Delta))/(2.*Delta);
	  
	  double nx=dx/sqrt(sq(dx)+sq(dy)+1e-10);
	  double ny=dy/sqrt(sq(dx)+sq(dy)+1e-10);
	  
	  double alpha7=ddelta-d7;
	  double x1=x+alpha7*nx;
	  double y1=y+alpha7*ny;
	  
	  double u1x=interpolate(u.x,x1,y1);
	  double u1y=interpolate(u.y,x1,y1);

	  double ut1=ny*u1x-nx*u1y;
	  
	  if (ut1>0.) {
	    ut1=-ut1;
	    nx=-nx;
	    ny=-ny;
	  }
	  double utauprec;
	  double err=1.;
	  int n_iter=0;
	  while ((n_iter<50)&&(err>1e-3)) {
	    utauprec=utau7;
	    utau7=utau7-gwallSA(utau7,ddelta,ut1,molvis[])/gwallSAprime(utau7,ddelta,ut1,molvis[]);
	    err=fabs(utau7-utauprec)/(fabs(utau7)+1e-10);
	    n_iter+=1;
	  }
	  fprintf(utaufile7,"%g,%g,%g,%g\n",x,utau7,d7*utau7/molvis[],ddelta*utau7/molvis[]);
	  fflush(utaufile7);
	}
      }
    }
  }

  char name8[1024];
  sprintf(name8,"max_turbulent_viscosity_Re=%g_level=%d.txt",Reynolds,maxlevel);
  FILE * turbulent_viscosity8=fopen(name8,"w");
  int number8x=100;
  int number8y=50;
  double xp8,yp8,muhat8,res8,d8;
  double yp80=0.;
  double max8=0.;
  for (int kkx=0;kkx<=number8x;kkx++) {
    xp8=b*kkx*1.5/number8x;
    max8=0.;
    if ((xp8<=0.3)||(xp8>=1.2)) {
      yp80=0.;
    }
    else {
      yp80=a*pow(sin(pi*xp8/l-pi/3.),4.);
    }
    for (int kky=0;kky<=number8y;kky++) {
      yp8=yp80+b*kky*0.08/number8y;
      muhat8=interpolate(muhat,xp8,yp8);
      d8=distance_to_wall(xp8,yp8);
      res8=muhat8/molvis2*fv1(muhat8/molvis2);
      if ((d8>Deltamin)&&(res8>max8)) {
	max8=res8;
      }
    }
    fprintf(turbulent_viscosity8,"%g,%g\n",xp8,max8);
    fflush(turbulent_viscosity8);
  }
}

#if 1
event adapt(i++) {
  adapt_wavelet({u,muhat},(double[]){1e-5,1e-5,1e-5},maxlevel,5);
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
We can plot the velocity profile and compare it with the profile obtain using CFL3D (dash line).


~~~pythonplot Velocity profile
import numpy as np
import matplotlib.pyplot as plt

Xpos=[0.75,1.20148]
Ypos10=[]
Upos10=[]
Ypos11=[]
Upos11=[]
Yt=[]
Ut=[]
for x in Xpos:

    source=open("velocity_profile_x="+str(x)+"_Re=3e+06_level=10.txt",'r')
    Y=[]
    U=[]
    for ligne in source:
        L=ligne.strip().split(',')
        if len(Y)<=100:
            Y.append(float(L[1]))
            U.append(float(L[2]))
    source.close()
    Ypos10.append(Y)
    Upos10.append(U)
    
    source=open("velocity_profile_x="+str(x)+"_Re=3e+06_level=11.txt",'r')
    Y=[]
    U=[]
    for ligne in source:
        L=ligne.strip().split(',')
        if len(Y)<=100:
            Y.append(float(L[1]))
            U.append(float(L[2]))
    source.close()
    Ypos11.append(Y)
    Upos11.append(U)

    Y=[]
    U=[]
    source=open("../data/Bump/Data_ratio_velocity_profile_x="+str(x)+"_NASA_CFL3D.txt",'r')
    
    for ligne in source:
        L=ligne.strip().split(' ')
        if float(L[1])>0.04:
            break
        Y.append(float(L[1]))
        U.append(float(L[2]))

    source.close()
    Yt.append(Y)
    Ut.append(U)

Color=['blue','darkblue']
plt.figure()
for i in range(len(Xpos)):
   if i==0:
      plt.scatter(Upos10[i],Ypos10[i],color=Color[i],marker='v',label="level=10")
      plt.scatter(Upos11[i],Ypos11[i],color=Color[i],marker='.',label="level=11")
   else:   
      plt.scatter(Upos10[i],Ypos10[i],color=Color[i],marker='v')
      plt.scatter(Upos11[i],Ypos11[i],color=Color[i],marker='.')
   plt.plot(Ut[i],Yt[i],color=Color[i],linestyle="--",label="x="+str(Xpos[i]))
plt.xlabel("Vertical position",fontsize=15)
plt.ylabel("Velocity",fontsize=15)
plt.title("Velocity profile",fontsize=20)
plt.legend()
plt.savefig('velocity.png')
~~~

We can also look at the turbulent viscosity at the top of the bump.

~~~pythonplot Turbulent viscosity
Level=[10,11]

Yl=[]
Nuhatl=[]

level=Level[0]

source=open("ratio_turbulent_viscosity_x=0.75_Re=3e+06_level="+str(level)+".txt",'r')
Y=[]
Nuhat=[]
for ligne in source:
   L=ligne.strip().split(',')
   if len(Y)<=100:
        Y.append(float(L[1]))
        Nuhat.append(float(L[2]))
source.close()
Yl.append(Y)
Nuhatl.append(Nuhat)

for level in Level[1:]:
   source=open("../data/Bump/ratio_turbulent_viscosity_x=0.75_Re=3e+06_level="+str(level)+".txt",'r')
   Y=[]
   Nuhat=[]
   for ligne in source:
      L=ligne.strip().split(',')
      if len(Y)<=100:
           Y.append(float(L[1]))
           Nuhat.append(float(L[2]))
   source.close()
   Yl.append(Y)
   Nuhatl.append(Nuhat)

Yt=[]
Nuhatt=[]
source=open("../data/Bump/Data_ratio_viscosity_turbulent_NASA_CFL3D.txt",'r')
for ligne in source:
    L=ligne.strip().split(' ')
    if float(L[1])>0.065:
       break
    Yt.append(float(L[1]))
    Nuhatt.append(float(L[2]))

source.close()

plt.figure()
for i in range(len(Level)):
   plt.scatter(Nuhatl[i],Yl[i],label="level="+str(Level[i]),marker='+')
plt.plot(Nuhatt,Yt,label="NASA CFL3D",linestyle='--',color="k",linewidth=2)
plt.xlabel("Ratio turbulent viscosity",fontsize=15)
plt.ylabel("$y$",fontsize=15)
plt.title("Profile of turbulent viscosity at $x=0.75$",fontsize=20)
plt.legend()
plt.ylim([0.050,0.065])
plt.savefig('Turbulent_viscosity.png')
~~~

Finally, we can see the friction coefficient and the pressure coefficient along the bump.

~~~pythonplot Skin Friction coefficient
Level=[10,11]

Xl=[]
Cfl=[]

level=Level[0]
source=open("skin_friction_Re=3e+06_level="+str(level)+".txt",'r')
X=[]
Cf=[]
for ligne in source:
    L=ligne.strip().split(',')
    X.append(float(L[0]))
    Cf.append(float(L[1]))
source.close()
Xl.append(X)
Cfl.append(Cf)

for level in Level[1:]:
   source=open("../data/Bump/skin_friction_Re=3e+06_level="+str(level)+".txt",'r')
   X=[]
   Cf=[]
   for ligne in source:
       L=ligne.strip().split(',')
       X.append(float(L[0]))
       Cf.append(float(L[1]))
   source.close()
   Xl.append(X)
   Cfl.append(Cf)

Xt=[]
Cft=[]
source=open("../data/Bump/Data_skin_friction_coefficient_NASA_CFL3D.txt",'r')

for ligne in source:
    L=ligne.strip().split(' ')
    Xt.append(float(L[0]))
    i=1
    while L[i]=='':
       i+=1
    Cft.append(float(L[i]))

source.close()

plt.figure()
for i in range(len(Level)):
   plt.scatter(Xl[i],Cfl[i],label="level="+str(Level[i]),marker='+')
plt.plot(Xt,Cft,label="NASA CFL3D",linestyle='--',color="k",linewidth=2)
plt.xlabel("$x$",fontsize=15)
plt.ylabel("Skin friction",fontsize=15)
plt.ylim([0,0.009])
plt.title("Skin friction coefficient along the plate",fontsize=20)
plt.legend()
plt.savefig('Skin friction.png')
~~~

~~~pythonplot Pressure coefficient
Level=[10,11]

Xl=[]
Cpl=[]

level=Level[0]
source=open("skin_pressure_Re=3e+06_level="+str(level)+".txt",'r')
X=[]
Cp=[]
for ligne in source:
    L=ligne.strip().split(',')
    if len(X)<1100 and float(L[0])<2:
        X.append(float(L[0]))
        Cp.append(-2*float(L[2]))
source.close()
Xl.append(X)
Cpl.append(Cp)

for level in Level[1:]:
   source=open("../data/Bump/skin_pressure_Re=3e+06_level="+str(level)+".txt",'r')
   X=[]
   Cp=[]
   for ligne in source:
       L=ligne.strip().split(',')
       if len(X)<1100 and float(L[0])<2:
           X.append(float(L[0]))
           Cp.append(-2*float(L[2]))
   source.close()
   Xl.append(X)
   Cpl.append(Cp)

Xt=[]
Cpt=[]
source=open("../data/Bump/Data_skin_pressure_coefficient_NASA_CFL3D.txt",'r')

for ligne in source:
    L=ligne.strip().split(' ')
    Xt.append(float(L[0]))
    i=1
    while L[i]=='':
        i+=1
    Cpt.append(-float(L[i]))

source.close()

plt.figure()
for i in range(len(Level)):
   plt.scatter(Xl[i],Cpl[i],label="level="+str(Level[i]),marker='+')
plt.plot(Xt,Cpt,label="NASA CFL3D",linestyle='--',color="k",linewidth=2)
plt.xlabel("$x$",fontsize=15)
plt.ylabel("Skin pressure",fontsize=15)
plt.title("Skin pressure coefficient along the plate",fontsize=20)
plt.legend()
plt.savefig('Pressure coefficient.png')
~~~

*/