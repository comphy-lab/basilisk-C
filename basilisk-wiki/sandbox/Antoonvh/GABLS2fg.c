/**
# The GABLS2 case using a fixed and equidistant grid approach
comments will follow. Look [here](GABLS2.c).
*/
#include "grid/multigrid1D.h"
#include "diffusion.h"
#include "run.h"


#define T1bottom (((((-25*cos(0.22*((t/3600)+16.) + 0.2)-10.) *(((t/3600.)+16.)<=17.4))+((((t/3600.)+16.)<=30.)*(((t/3600.)+16.)>17.4)*((-0.54*((t/3600.)+16.))+15.2))+((((t/3600.)+16.)<=41.9)*(((t/3600.)+16.)>30)*(-7-(25*cos(0.21*((t/3600.)+16.)+1.8))))+((((t/3600.)+16.)<=53.3)*(((t/3600.)+16.)>41.9)*(-0.37*((t/3600.)+16.)+18.0))+((((t/3600.)+16.)<=65.6)*(((t/3600.)+16.)>53.3)*(-4-25*cos(0.22*((t/3600.)+16.)+2.5)))+((((t/3600.)+16.)>65.6)*(4.4)))+273.15))
#define maxlevel 9
// Holtslag and Boville with Stable F(Ri) according to vdW
#define fris(Ri) (sq((1-(Ri/0.20)))*(Ri<0.20)) //vdW 2017 Critical Ri
//#define fris(x) (1/(1+(10*x*(1+8*x))))
#define friu(Ri) (sqrt(1-(18.*Ri)))            // Holtslag en Boville 1992
#define friubm(Ri,y) ((1-((10*Ri)/(1+75*y*sqrt((x+zo/zo)*fabs(Ri))))))  // Louis 1982
#define friubh(Ri,y) ((1-((15*Ri)/(1+75*y*sqrt((x+zo/zo)*fabs(Ri))))))  // louis 1982

int nn;
double Tref=283.15;
double Ugeo=-8.;
double Vgeo = 3.;
double He= 65.6;
double Rv= 461.5;
double Rd = 287.;
double L= 2500000.;
double cpd= 1004;
double exner;
double Up[100];
double Cm,Ch,Cq,qtsat,TK1;
double zo=0.01;
scalar u[],v[],T1[],qt[];
scalar * tracers = {u,v,T1,qt};
mgstats mgb;
int m = 0;
double eu=0.;
double eb=0.;


int main(){
  init_grid(1<<(maxlevel));
  L0=4096;
  X0=0;
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
}

event Diffusion(i++){
  boundary({u,v,T1,qt});
  nn=0;
  scalar rx[],ry[],rT1[],T1f[],rqt[],b[],w[];
  face vector kh[],sqd[],Ri[],fRi[];
  double CN,U,bbottom,es,qsl;
 
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
      if (b[]>bbottom){
	Cm=CN*fris(((x-zo)*(b[]-(bbottom))/(sq(u[])+sq(v[]))));
	Ch=Cm;
	Cq=Cm*0.025;
      }
      else{
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
  
  // Log total number of Multigrid Cycles
  mgb=diffusion(u,dt,kh,rx);
  nn+=mgb.i;
  mgb=diffusion(v,dt,kh,ry);
  nn+=mgb.i;
  mgb=diffusion(T1,dt,kh,rT1);
  nn+=mgb.i;
  mgb=diffusion(qt,dt,kh,rqt);
  nn+=mgb.i;
}

event output(t+=360){
  scalar U[];
  int ng=0;
  foreach(){
    U[]=sqrt(sq(u[])+sq(v[]));
    ng++;
  }
  boundary({U});
  static FILE * fp4 = fopen("datogablsfg.dat","w");
  fprintf(fp4,"%g\t%g\t%d\t%g\t%d\t%g\t%g\n",t,dt,i,interpolate(U,10),ng,qtsat,TK1);

  fflush(fp4);
  fprintf(stdout,"%g\t%g\t%d\t%g\t%d\t%g\t%g\n",t,dt,i,interpolate(U,10),ng,qtsat,TK1);
  static FILE * fp = fopen("prfilesGABLS2fg.dat","w");
  double yp=0;
  while (yp<L0){
    Point point = locate(yp);
    yp=x;
    fprintf(fp,"%g\t%g\t%g\t%g\t%g\t%g\t%d\n",yp,u[],v[],T1[],qt[],U[],level);
    yp+=Delta/1.5;
  }
  fflush(fp);
  static FILE * fp5 = fopen("gabls2gridfg.dat","w");
  for (int mm=0.;mm<=4000;mm+=4){
    Point point = locate((double)mm);
    fprintf(fp5,"%d\t",level);
  }
  fprintf(fp5,"\n");
  fflush(fp5);
  
}

event adapt(i++;t<He*3600){
  //do not adapt grid but only adapt timestep

  // adapt_wavelet({u,v,T1},(double[]){eu,eu,eb},maxlevel,2,{u,v,T1,qt});
  
  /** Timestepping is adapted to keep the total number of MG-cycles within reasonable bounds
   */
  
  if (nn>14)
    DT=max(DT/(1+((double)nn/10.)),2.);
  if (nn<8)
    DT=min(DT*(1+((double)nn/100.)),15.);
}

