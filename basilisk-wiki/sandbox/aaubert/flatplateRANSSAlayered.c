#include "grid/multigrid1D.h"
#include "layered_perso/hydro.h"
#include "layered_perso/nh.h"
#include "layered_perso/remap.h"


#define IMPLICIT_SA 0
#define WALL 1             //wall is present : need to define a distance_to_wall function
double distance_to_wall(double x,double y) {
  if (x>1.) {
    return y;
  }
  else {
    return sqrt(pow(x-1.,2.)+pow(y,2.));
  }
}
#include "layered_perso/RANSSpalartAllmaraslayered.h"

#include "layered/perfs.h"


double Reynolds=5000000.;
int maxlevel=11;

double Tc=1.;
double U0=1.;
double U00=1.;
double t0=30.;
double molvis;
face vector muv[];
coord * listcoord;

double ddelta;

scalar un[];

int main() {
  for (nl=40;nl<=40;nl+=30) {  //fixme:2 simulations in a row does not work
    maxlevel=9;
    size(3. [1]);
    N=1024;
    init_grid(N);
    //mu=muv;
    G=0.;
    //theta_H=1.;
    //nl=51;
    ddelta=5.*L0/pow(2.,maxlevel);
    theta=1.3;
    u.x.gradient=minmod2;
    w.gradient=minmod2;
    phi.gradient=minmod2;
    muhat.gradient=minmod2;
    un=new scalar[nl];
    run();
  }
}

event init(t=0) {
  //geometric_beta(rmin=0.2,top=false);
  geometric_beta_perso(H=1,rmin=0.2,nc=10);
  foreach() {
    foreach_dimension() {
      //lambda_b.x[]=x<1.?1.:0.;
    }
    zb[]=0.;
  }
  //CFL=1.0;
  molvis=U0/Reynolds;
  //ct3=0.;
  if (!restore(file="dump-t=25_nl=50hhh")) {
    foreach() {
      foreach_layer() {
	h[]=L0/(3.*nl);
      }
    }
    vertical_remapping(h,tracers);
    foreach() {
      double zlayer=0.;
      foreach_layer() {
	zlayer+=h[]/2.;
	u.x[]=U0;
	muhat[]=3.*molvis;
	zlayer+=h[]/2.;
      }
    }
  }
  else {
    N=1024;
  }
  u.n[left]=dirichlet(U0);
  w[left]=dirichlet(0.);
  muhat[left]=neumann(0.);
  phi[right]=dirichlet(0.);
  u.n[right]=neumann(0.);
  w[right]=neumann(0.);
  muhat[right]=neumann(0.);

}


event cleanup(t=end,last) {
  delete ({un});
}

event image(t=5) {
  char name1[1024];
  sprintf(name1,"test_velocity_distribution_level=%d.txt",maxlevel);
  FILE * test1=fopen(name1,"w");
  double xp1,yp1,res1;
  foreach() {
    xp1=x;
    yp1=h[0,0];
    res1=u.x[0,0];
    fprintf(test1,"%g,%g,%g\n",xp1,yp1,res1);
    fflush(test1);
  }

  char name2[1024];
  sprintf(name2,"test_velocity_vertical_level=%d.txt",maxlevel);
  FILE * test2=fopen(name2,"w");
  double xp2,yp2,res2;
  foreach_point(1.0049) {
    xp2=x;
    yp2=0.;
    foreach_layer() {
      yp2+=h[]/2.;
      res2=u.x[];
      fprintf(test2,"%g,%g,%g\n",xp2,yp2,res2);
      fflush(test2);
      yp2+=h[]/2.;
    }
  }

  double utau=0.;
  char name3[1024];
  sprintf(name3,"velocity_profile_x=1_Re=%g_nlayer=%d.txt",Reynolds,nl);
  FILE * test3=fopen(name3,"w");
  double xp3,yp3,res3;
  foreach_point(2.) {
    xp3=x;
    yp3=zb[];
    foreach_layer() {
      if (point.l==0) {
	utau=sqrt(molvis*u.x[]*2./h[]);
      }
      yp3+=h[]/2.;
      res3=u.x[];
      fprintf(test3,"%g,%g,%g\n",xp3,yp3*utau/molvis,res3/utau);
      fflush(test3);
      yp3+=h[]/2.;
    }
  }

  char name4[1024];
  sprintf(name4,"ratio_turbulent_viscosity_x=1_Re=%g_nlayer=%d.txt",Reynolds,nl);
  FILE * test4=fopen(name4,"w");
  double xp4,yp4,res4;
  foreach_point(2.) {
    xp4=x;
    yp4=zb[];
    foreach_layer() {
      yp4+=h[]/2.;
      res4=muhat[]*fv1(muhat[]/molvis)/molvis;
      fprintf(test4,"%g,%g,%g\n",xp4,yp4,res4);
      fflush(test4);
      yp4+=h[]/2.;
    }
  }

  char name5[1024];
  sprintf(name5,"turbulent_viscosity_Re=%g_nlayer=%d_layered.txt",Reynolds,nl);
  FILE * turbulent_viscosity5=fopen(name5,"w");
  double xp5,yp5,muhat5,res5;
  foreach() {
    xp5=x;
    yp5=zb[];
    foreach_layer() {
      yp5+=h[]/2.;
      muhat5=muhat[];
      res5=muhat5*fv1(muhat5/molvis)/molvis;
      fprintf(turbulent_viscosity5,"%g,%g,%g\n",xp5,yp5,res5);
      fflush(turbulent_viscosity5);
      yp5+=h[]/2.;
    }
  }

  char name6[1024];
  sprintf(name6,"vertical_velocity_profile_x=1_Re=%g_nlayer=%d.txt",Reynolds,nl);
  FILE * test6=fopen(name6,"w");
  double xp6,yp6,res6;
  foreach_point(2.) {
    xp6=x;
    yp6=zb[];
    foreach_layer() {
      yp6+=h[]/2.;
      res6=w[];
      fprintf(test6,"%g,%g,%g\n",xp6,yp6,res6);
      fflush(test6);
      yp6+=h[]/2.;
    }
  }

  char name7[1024];
  sprintf(name7,"skin_friction_Re=%g_nlayer=%d.txt",Reynolds,nl);
  FILE * test7=fopen(name7,"w");
  double xp7,yp7,res7;
  foreach() {
    if (x>=1.) {
      xp7=x;
      yp7=h[0,0,0]/2.;
      res7=molvis*u.x[0,0,0]/yp7;
      fprintf(test7,"%g,%g,%g\n",xp7,yp7,2.*res7);
      fflush(test7);
    }
  }

  char name8[1024];
  sprintf(name8,"horizontal_velocity_Re=%g_nlayer=%d_layered.txt",Reynolds,nl);
  FILE * velocity8=fopen(name8,"w");
  double xp8,yp8,ux8,res8;
  foreach() {
    xp8=x;
    yp8=zb[];
    foreach_layer() {
      yp8+=h[]/2.;
      ux8=u.x[];
      res8=ux8;
      fprintf(velocity8,"%g,%g,%g\n",xp8,yp8,res8);
      fflush(velocity8);
      yp8+=h[]/2.;
    }
  }

  char name9[1024];
  sprintf(name9,"vertical_velocity_Re=%g_nlayer=%d_layered.txt",Reynolds,nl);
  FILE * velocity9=fopen(name9,"w");
  double xp9,yp9,ux9,res9;
  foreach() {
    xp9=x;
    yp9=zb[];
    foreach_layer() {
      yp9+=h[]/2.;
      ux9=w[];
      res9=ux9;
      fprintf(velocity9,"%g,%g,%g\n",xp9,yp9,res9);
      fflush(velocity9);
      yp9+=h[]/2.;
    }
  }
  
}

#if 0
event adapt(i++) {
  //adapt_wavelet({u},(double[]){1e-3,1e-3},maxlevel,5);
  adapt_wavelet({u,muhat},(double[]){1e-3,1e-3},maxlevel,5);
}
#endif

#if 0
event snapshot(t+=5.) {
  char name[80];
  sprintf(name,"dump-t=%g_nlayer=%d",t,nl);
  dump (file=name);
}
#endif

/**
~~~pythonplot Skin friction coefficient
import numpy as np
import matplotlib.pyplot as plt

Nl=[40,50]

Xn=[]
Yn=[]
Cfn=[]

for nl in Nl:
    X=[]
    Y=[]
    Cf=[]
    source=open("skin_friction_Re=5e+06_nlayer="+str(nl)+".txt",'r')
    
    for ligne in source:
        L=ligne.strip().split(',')
        X.append(float(L[0])-1)
        Y.append(float(L[1]))
        Cf.append(float(L[2]))
        
    source.close()
    Xn.append(X)
    Yn.append(Y)
    Cfn.append(Cf)

Xt=[]
Cft=[]
source=open("../data/Flat_plate/Data_skin_friction_coefficient_NASA_CFL3D.txt",'r')

for ligne in source:
   L=ligne.strip().split(' ')
   Xt.append(float(L[0]))
   i=1
   while L[i]=="":
       i+=1
   Cft.append(float(L[i]))

source.close()

plt.figure()
for i in range(len(Nl)):
    plt.scatter(Xn[i],Cfn[i],label="nl="+str(Nl[i]),marker="+")
plt.plot(Xt,Cft,linestyle="--",label="CFL3D",color='k')
plt.legend()
plt.ylim([0,0.01])
plt.savefig('Skin_friction_coefficient.png')
~~~

~~~pythonplot Velocity profile at x=1
Xn=[]
Yn=[]
Un=[]

for nl in Nl:
   X=[]
   Y=[]
   U=[]
   source=open("velocity_profile_x=1_Re=5e+06_nlayer="+str(nl)+".txt",'r')

   for ligne in source:
       L=ligne.strip().split(',')
       X.append(float(L[0]))
       Y.append(float(L[1]))
       U.append(float(L[2]))

   source.close()
   Xn.append(X)
   Yn.append(Y)
   Un.append(U)

Yt=[]
Ut=[]
source=open("../data/Flat_plate/Data_y+u+_NASA_CFL3D.txt",'r')

for ligne in source:
   L=ligne.strip().split(' ')
   if float(L[0])>5.23:
       break
   Yt.append(np.exp(float(L[0])*np.log(10)))
   Ut.append(float(L[3]))

source.close()

plt.figure()
for i in range(len(Nl)):
   plt.scatter(Yn[i],Un[i],label="$y+=$"+str(Yn[i][0]),marker="+")
plt.plot(Yt,Ut,linestyle="--",label="CFL3D",color='k')
plt.legend()
plt.xscale('log')
plt.savefig('Velocity_profile.png')
~~~

~~~pythonplot Turbulent viscosity at x=1

Xn=[]
Yn=[]
Mun=[]

for nl in Nl:
   X=[]
   Y=[]
   Mu=[]
   source=open("ratio_turbulent_viscosity_x=1_Re=5e+06_nlayer="+str(nl)+".txt",'r')

   for ligne in source:
       L=ligne.strip().split(',')
       if float(L[1])<=0.025:
          X.append(float(L[0]))
          Y.append(float(L[1]))
          Mu.append(float(L[2]))

   source.close()
   Xn.append(X)
   Yn.append(Y)
   Mun.append(Mu)

Yt=[]
Mut=[]
source=open("../data/Flat_plate/Data_ratio_viscosity_turbulent_NASA_CFL3D.txt",'r')

for ligne in source:
    L=ligne.strip().split(' ')
    if float(L[1])>0.025:
        break
    Yt.append(float(L[1]))
    Mut.append(float(L[2]))

source.close()

plt.figure()
for i in range(len(Nl)):
   plt.scatter(Mun[i],Yn[i],label="nl="+str(Nl[i]),marker="+")
plt.plot(Mut,Yt,linestyle="--",label="CFL3D",color='k')
plt.legend()
plt.ylim([0,0.025])
plt.savefig('Turbulent_viscosity.png')
~~~
*/