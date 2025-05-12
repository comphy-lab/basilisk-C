/**
## Variations on the second GABLS intercomparison case
This page presents how the mixing length $l$ controls the slope in the observed temperature profile. For a descroption of the set-up see [the default set-up](GABLS2.c)    
*/
#include "grid/bitree.h"
#include "diffusion.h"
#include "run.h"

#define T1bottom (((((-25*cos(0.22*((t/3600)+16) + 0.2)-10) *(((t/3600)+16)<=17.4))+((((t/3600)+16)<=30)*(((t/3600)+16)>17.4)*((-0.54*((t/3600)+16))+15.2))+((((t/3600)+16)<=41.9)*(((t/3600)+16)>30)*(-7-(25*cos(0.21*((t/3600)+16)+1.8))))+((((t/3600)+16)<=53.3)*(((t/3600)+16)>41.9)*(-0.37*((t/3600)+16)+18.0))+((((t/3600)+16)<=65.6)*(((t/3600)+16)>53.3)*(-4-25*cos(0.22*((t/3600)+16)+2.5)))+((((t/3600)+16)>65.6)*(4.4)))+273.15))
#define maxlevel 9
// Holtslag and Boville with Stable F(Ri) according to vdW
#define fris(Ri) (sq((1-(Ri/0.20)))*(Ri<0.20)) //vdW 2017 Critical Ri
//#define fris(x) (1/(1+(10*x*(1+8*x))))
#define friu(Ri) (sqrt(1-(18.*Ri)))            // Holtslag en Boville 1992
#define friubm(Ri,y) ((1-((10*Ri)/(1+75*y*sqrt((x+zo/zo)*fabs(Ri))))))  // Louis 1982
#define friubh(Ri,y) ((1-((15*Ri)/(1+75*y*sqrt((x+zo/zo)*fabs(Ri))))))  // louis 1982

/**
Below some constants are decleared.  
*/

double Tref=283.15;   // Reference temperture
double Ugeo=-8.;      // u-component of the Geostrophic wind
double Vgeo = 3.;     // v-component of the Geostrophic wind
double zo=0.01;       // Roughness length for heat and momentum
double He= 65.6;      //Total hours of simulation

double Rv= 461.5;     // Gas constant for water vapor
double Rd = 287.;     // Gas constant for dry air
double L= 2500000.;   // Latent heat release for evaporization
double cpd= 1004;     // Density of air.

double exner;         
int nn;
double Up[100];
double Cm,Ch,Cq,qtsat,TK1;

scalar u[],v[],T1[],qt[];

mgstats mgb;
double eu=0.25;
double eb=0.5;
/**
Here the maxim mixing length is declared as a variable `lz` so that it can be varied. 
*/
double lz = 70;
int j;

int main(){
  /**
  ## 4 runs. 
  We run four different runs; 
  
  1. Default run
  2. Using a decreased max. mixing length $lz$
  3. Using an increased max. mixing length $lz$
  4. Default run, without adaptivity. 
  */
  init_grid(1<<(maxlevel));
  L0=4096;
  X0=0;
  j = 0;
  run();
  j = 1;
  lz=35;
  run();
  j=2;
  lz=140;
  run();
  j=3;
  lz=70;
  run();
}	

u[left]=dirichlet(0.);
v[left]=dirichlet(0.);
T1[left]=dirichlet(T1bottom);
T1[right]=neumann((20./3000.));

event init(i=0){
  exner = pow(972./1000.,Rd/cpd);
  TOLERANCE=10E-8;
  DT=1;
  u.refine=refine_linear;
  v.refine=refine_linear;
  T1.refine=refine_linear;
  foreach(){
    u[]=Ugeo;
    v[]=Vgeo;
    T1[]=(((x<=200))*(288-(2*x/200)))+
      ((x>200)*(x<=850)*286)+
      ((x>850)*(x<=900)*(286+2*(x-850)/50))+
      ((x>900)*(x<=1000)*(288+4*(x-900)/100))+
      ((x>1000)*(292+20*((x-1000)/3000)));
    qt[]=0.0025;
  }
  if (j!=3){
    while(adapt_wavelet({u,v,T1},(double[]){eu,eu,eb},maxlevel,4,{u,v,T1,qt}).nc){
      foreach(){
	u[]=Ugeo;
	v[]=Vgeo;
	T1[]=(((x<=200))*(288-(2*x/200)))+
	((x>200)*(x<=850)*286)+
	  ((x>850)*(x<=900)*(286+2*(x-850)/50))+
	  ((x>900)*(x<=1000)*(288+4*(x-900)/100))+
	  ((x>1000)*(292+20*((x-1000)/3000)));
	qt[]=(0.0025*(x<=900))+
	  ((x>900)*(x<=1000)*(0.0025-(0.002*(x-900)/100)))+
	  ((x>1000)*(x<=2000)*(0.0005+(0.0025*(x-1000)/1000)))+
	  ((x>2000)*(x<=3500)*(0.003-(0.001*(x-2000)/1500)))+
	  ((x>3500)*(0.002-(0.0005*(x-3500)/500)));
      }
    }
  }
}

event Diffusion(i++){
  boundary({u,v,T1,qt});
  nn=0;
  scalar rx[],ry[],rT1[],T1f[],rqt[],b[],w[];
  face vector kh[],sqd[],Ri[],fRi[];
  double CN,U,bbottom,es,qsl;
  kh.x.refine=no_restriction;
  sqd.x.refine=no_restriction;
  Ri.x.refine=no_restriction;
  fRi.x.refine=no_restriction;
  foreach()
    b[]=((T1[]-Tref)*(9.81/Tref))*(1+(0.608*qt[]));
  foreach(){
    w[]=-0.005*(((x/1000)*(x<=1000))+(x>1000))*(((t/3600)+16-24)>=24);
    rx[]=0.000139*(v[]-Vgeo);
    rx[]-=w[]*((u[1]-u[-1])/(2*Delta));
    rT1[]=-w[]*((T1[1]-T1[-1])/(2*Delta));
    ry[]=0.000139*(Ugeo-u[]);
    ry[]-=w[]*((v[1]-v[-1])/(2*Delta));
    rqt[]=-w[]*((qt[1]-qt[-1])/(2*Delta));
    if (x<Delta){ //Compute Surface fluxes
      bbottom=((T1bottom-Tref)*(9.81/Tref))*(1+(0.608*qt[]));
      TK1= T1bottom*exner;
      es=610.78*exp(17.27*(TK1-273.16)/(TK1-35.86));
      qsl= (Rd/Rv)*(es/(97200.-(1-((1-Rd/Rv)*es))));
      qtsat=qsl*((1+((sq(L)/(Rv*cpd*sq(TK1)))*qt[]))/(1+((sq(L)/(Rv*cpd*sq(TK1)))*qsl)));
      CN=sq(0.4/log((x)/zo));
      if (b[]>bbottom){ //Stable
	Cm=CN*fris(((x-zo)*(b[]-(bbottom))/(sq(u[])+sq(v[]))));
	Ch=Cm;
	Cq=Cm*0.025;
      }
      else{ //Unstable
	CN = sq(0.4/log((x)/zo));
	Cm = CN*friubm((x-zo)*(b[]-(bbottom))/(sq(u[])+sq(v[])),CN);
	Ch = CN*friubh((x-zo)*(b[]-(bbottom))/(sq(u[])+sq(v[])),CN);
	Cq = Ch*0.025;
      }
      U=sqrt(sq(u[])+sq(v[]));
      rx[]-=(u[]*Cm*U)/Delta;
      ry[]-=(v[]*Cm*U)/Delta;
      rT1[]-=((T1[]-T1bottom)*Cm*U)/Delta;
      rqt[]-=((qt[]-qtsat)*Cq*U)/Delta;
    }
  }
   
  boundary(all);
  foreach_face()//Compute turbulent diffusivities
    {
      sqd.x[]=(sq((u[]-u[-1])/(Delta))+sq((v[]-v[-1])/(Delta)));
      Ri.x[]= ((b[]-b[-1])/(Delta))/(sqd.x[]+0.00001);
      if (Ri.x[]<0)
	fRi.x[]=friu(Ri.x[]);
      else
	fRi.x[]=fris(Ri.x[]);
      kh.x[]=sq(min(0.4*x,lz))*(sqrt(sqd.x[]))*fRi.x[];
    }
  boundary({kh.x});
  dt=dtnext(DT);
  
  // Time integrate and log the total number of Multigrid Cycles
  mgb=diffusion(u,dt,kh,rx);
  nn+=mgb.i;
  mgb=diffusion(v,dt,kh,ry);
  nn+=mgb.i;
  mgb=diffusion(T1,dt,kh,rT1);
  nn+=mgb.i;
  mgb=diffusion(qt,dt,kh,rqt);
  nn+=mgb.i;
}
/**
## output
We output the profiles at the specified time.
*/
event outprof (t = 79200){
  char fname[99];
  sprintf(fname, "%d", j);
  FILE * fp = fopen (fname, "w");
  foreach()
    fprintf(fp, "%g\t%g\n", x,T1[]);
  fclose(fp);
}

/**
# Stop
The simulation is stopped right after we have outputted the desired profiles. 
*/
event stop (t = 79202){
  return 1;
}

event adapt(i++;t<He*3600){
  if (j!=3)
    adapt_wavelet({u,v,T1},(double[]){eu,eu,eb},maxlevel,2,{u,v,T1,qt});
}
/**
## Results
The profiles:

~~~gnuplot 
set yr [0:1100]
set xr [286:293]
set ylabel 'Height [m]'
set xlabel 'Temp [K]'
set key box bottom right
plot '0' u 2:1 w l lw 2 t 'Default Lz', \
     '1' u 2:1 w l lw 2 t 'Decreased Lz', \
     '2' u 2:1 w l lw 2 t 'Increased Lz', \
     '3' u 2:1 w l lw 2 t 'Default Lz + equidistant grid'
~~~

The plot looks like a less-cared for version of the one presented in the rebuttle to ref. 1 and 2.
*/
