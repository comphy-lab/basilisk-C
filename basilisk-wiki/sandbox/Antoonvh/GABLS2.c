/**
![The second GABLS intercomparison scenario was based on obervations from the CASES-99 field experiment in Leon, Kansas, USA.](https://www.eol.ucar.edu/rtf/projects/CASES99/stn2.gif)

## The second GABLS intercomparison case
This page discusses the set-up of the GABLS2 case as presented in Van Hooft et al. (under review for GMD). Since there is a extensive discussion on how to set-up an atmopsheric single-column model (SCM) for the GABLS1 case we refer thee reader to the [corresponding page](GABLS1.c) for the general SCM lay-out. Hence, this page focusses on the GABLS2 specific formulations.   
*/
#include "grid/bitree.h"
#include "diffusion.h"
#include "run.h"
/**
After a tedious copy session, we were able to define the expression for the surface temperature. If only there were a more easy way to set a boundary condition for the diurnal cycle... Also a maximum level of refinement is defined using a Macro.
*/
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
/**
Compared to the GABLS1 case, the GABLS2 scenario prescribes the inclusion of 'moist' physics. For this case this means that the buoyuancy of an air parcal is modified by the total water vapor content. This requires us to initialize some number that define some relevant properties of air. 
*/
double Rv= 461.5;     // Gas constant for water vapor
double Rd = 287.;     // Gas constant for dry air
double L= 2500000.;   // Latent heat release for evaporization
double cpd= 1004;     // Density of air.

double exner;         
int nn;
double Up[100];
double Cm,Ch,Cq,qtsat,TK1;

/**
The solver will aim to calculate the evolution of the velocity components ($u,v$), potential temperature $T_1$ and the total water specific humidity $q_t$  
*/
scalar u[],v[],T1[],qt[];

mgstats mgb;
double eu=0.25;
double eb=0.5;

/**
The grid is initialized using 128 cells, and the domain height is set to be 4096 meters. 
*/
int main(){
  init_grid(1<<(maxlevel-2));
  L0=4096;
  X0=0;
  run();
}	
/**
Boundary conditions are set in a similar fashion as was done for the GABLS1 case. 
*/
u[left]=dirichlet(0.);
v[left]=dirichlet(0.);
T1[left]=dirichlet(T1bottom);
T1[right]=neumann((20./3000.));

/**
The simulation is initialized according to the case description. Special care is taken so that the simulation starts with a grid that is consistent with the adaptation algorithm. 
*/
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
  while(adapt_wavelet({u,v,T1},(double[]){eu,eu,eb},maxlevel,4,{u,v,T1,qt}).nf){
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
  /**
  As an intermediate step we calculate the buoyancy field from the potential temperature and the total water specific humidity.  
  */
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
      /**
      The flux of specific humidity is computed based on the humidity of a a water-saturated surface.
      */
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
      kh.x[]=sq(min(0.4*x,70))*(sqrt(sqd.x[]))*fRi.x[];
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
Each hour we perform the output routine ten times. The ouput consists non exclusively of the time series of the 10-meter wund speed $U_10$, profiles of the solution and the evolution of the grid structure. 
*/
event output(t+=360){
  scalar U[];
  int ng=0;
  foreach(){
    U[]=sqrt(sq(u[])+sq(v[]));
    ng++;
  }
  boundary({U});
  static FILE * fp4 = fopen("datogabls.dat","w");
  fprintf(fp4,"%g\t%g\t%d\t%g\t%d\t%g\t%g\n",t,dt,i,interpolate(U,10),ng,qtsat,TK1);

  fflush(fp4);
  fprintf(stdout,"%g\t%g\t%d\t%g\t%d\t%g\t%g\n",t,dt,i,interpolate(U,10),ng,qtsat,TK1);
  static FILE * fp = fopen("prfilesGABLS2.dat","w");
  double yp=0;
  while (yp<L0){
    Point point = locate(yp);
    yp=x;
    fprintf(fp,"%g\t%g\t%g\t%g\t%g\t%g\t%d\n",yp,u[],v[],T1[],qt[],U[],level);
    yp+=Delta/1.5;
  }
  fflush(fp);
  static FILE * fp5 = fopen("gabls2grid.dat","w");
  for (int mm=0.;mm<=4000;mm+=4){
    Point point = locate((double)mm);
    fprintf(fp5,"%d\t",level);
  }
  fprintf(fp5,"\n");
  fflush(fp5);
  
}
/**
Adaptation is based on the waveket-estimated error in the representation of the wind components and temperature. Notice that we do not refine based on the water vapor specific humidity.   
*/
event adapt(i++;t<He*3600){
  adapt_wavelet({u,v,T1},(double[]){eu,eu,eb},maxlevel,2,{u,v,T1,qt});
  
  /** 
  Timestepping is adapted to keep the total number of MG-cycles within reasonable bounds according to a `*Vertrouwen komt te voet en gaat te paard*' scheme. 
  */
  
  if (nn>14)//Quickly reduce timestep if things get rough
    DT=max(DT/(1+((double)nn/10.)),1.);
  if (nn<8)//Slowly increase timestep when timeintegration is easy.
    DT=min(DT*(1+((double)nn/100.)),15.);
}
/**
## Results

In order to check if the simulation has run as expected we plot some output obtained with the code above.

~~~gnuplot Level of refinement
set xr [0:1000]
set ylabel 'time/6min [-]'
set xlabel 'Height/4m [-]'
plot 'gabls2grid.dat' matrix with image
~~~

Looks like a less cared-for and transposed version of fig. 5 in the corresponding manuscript in GMD(D).
*/
