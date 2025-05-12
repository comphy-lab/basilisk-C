/**
# The GABLS1 case with an adaptive-grid Single Column model
On this page the details of the set-up of for the GABLS1 case are presented using a *fixed and equidistant* grid. A large part of the set-up is enherited from the [adaptive-grid set-up](GABLS1.c). Therefore this page will mosty discuss the specific changes made to that set-up to run the case with the equidistant grid. As such, if there is anything unclear about the set-up on this page, you may find more info on the page dedicated to the adaptive-grid set-up, i.e. presented [here](GABLS1.c)     

## General set-up
From the Basilisk source code, we include the one-dimensional multigrid-grid structure for our computations.   
 */
#include "grid/multigrid1D.h"
#include "diffusion.h"
#include "run.h"

#define fris(Ri) (sq((1-(Ri/0.20)))*(Ri<0.20)) //vdW 2017 Critical Ri, Short-tail mixing
//#define fris(x) (1/(1+(10*x*(1+8*x)))) // Long tail mixing
#define friu(Ri) (sqrt(1-(18*Ri)))            // Holtslag en Boville 1992
#define friubm(Ri,y) ((1-((10*Ri)/(1+75*y*sqrt((x+zo/zo)*fabs(Ri))))))  // Louis 1982
#define friubh(Ri,y) ((1-((15*Ri)/(1+75*y*sqrt((x+zo/zo)*fabs(Ri))))))  // louis 1982

#define bbottom (-0.25/26.5*(t/3600))

double zo=0.1;
/**
The maximum level will determined the used number of grid cell troughout the simulation, we set it so that the domain is discretized using ($2^6=)64$ cells.      
*/
int maxlevel = 6;
mgstats mgb;
int nn;
double Up[100],uu[100],vv[100],bb[100];
double Cm,Ch;
int m = 0;
scalar u[],v[],b[];

int main(){
  init_grid(1<<maxlevel);
  L0=400;
  X0=0;
  run();
}	

u[left]=dirichlet(0.);
v[left]=dirichlet(0.);
b[left]=dirichlet(bbottom);
b[right]=neumann(0.01/26.5);

event init(t=0){
  DT=1.;
  foreach(){
    u[]=8;
    v[]=0;
    b[]=(x>100)*(0.01/26.5)*(x-100);
  }
  boundary(all);
}

event Diffusion(i++){
  scalar rx[],ry[],rb[],bf[];
  face vector kh[],sqd[],Ri[],fRi[];
  double CN;
  foreach(){
    rx[]=0.000139*v[];
    rb[]=0;
    ry[]=0.000139*(8-u[]);
    if (x<Delta){
      if (b[]>bbottom){
	Cm=sq(0.4/log((x)/zo))*fris(((x-zo)*(b[]-(bbottom))/(sq(u[])+sq(v[]))));
	Ch=Cm;
      }
      else{
	CN =  sq(0.4/log((x)/zo));
	Cm =CN*friubm((x-zo)*(b[]-(bbottom))/(sq(u[])+sq(v[])),CN);
	Ch =CN*friubh((x-zo)*(b[]-(bbottom))/(sq(u[])+sq(v[])),CN);
      }
      rx[]-=(u[]*Cm*sqrt(sq(u[])+sq(v[])))/Delta;
      ry[]-=(v[]*Cm*sqrt(sq(u[])+sq(v[])))/Delta;
      rb[]-=((b[]-bbottom)*Cm*sqrt(sq(u[])+sq(v[])))/Delta;
    }
  }
  boundary(all);
  foreach_face(){
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
  mgb=diffusion(u,dt,kh,rx);
  nn+=mgb.i;
  mgb=diffusion(v,dt,kh,ry);
  nn+=mgb.i;
  mgb=diffusion(b,dt,kh,rb);
  nn+=mgb.i;
}
/**
## Output
Every ten minutes in physical time we output statistics of our simulation.
*/
event output(t+=360){
  static FILE * fp2 = fopen("GABLScellsMG.dat","w");
  int nnn=0;
  foreach()
    nnn++;
  fprintf(fp2,"%g\t%g\t%d\t%d\n",t,dt,i,nnn);
  fflush(fp2);
  double yp=0;
  static FILE * fp1 = fopen("prfileGABLS10mMG.dat","w");
  while (yp<400){
    Point point = locate(yp);
    yp=x;
    fprintf(fp1,"%g\t%g\t%g\t%g\t%g\t%d\n",yp,u[],v[],b[],sqrt(sq(u[])+sq(v[])),level);
    yp+=Delta/1.5;
  }
  fflush(fp1);
  static FILE * fp5 = fopen("gabls1gridMG.dat","w");
  for (double mm=0.;mm<=400;mm+=3.125){
    Point point = locate((double)mm);
    fprintf(fp5,"%d\t",level);
  }
  fprintf(fp5,"\n");
  fflush(fp5);
}
/**
Furthermore, in the last hour of simulation averaged profiles are calculated, we do this using the solution of each grid cell. 
*/
event avgprof(t=8*3600;i+=20)
{
  scalar U[];
  U[left]=dirichlet(0);
  int ng=0;
  foreach(){
    U[]=sqrt(sq(u[])+sq(v[]));
    ng++;
  }
  boundary(all);
  static FILE * fp = fopen("prfileGABLSfg.dat","w");
  double yp=0.;
  int j=0;
  m++;
  while (yp<400)    {
    Up[j]+=interpolate(U,yp);
    uu[j]+=interpolate(u,yp);
    vv[j]+=interpolate(v,yp);
    bb[j]+=interpolate(b,yp);
    if (t==8*3600)
      fprintf(fp,"%g\t%g\t%g\t%g\t%g\n",yp,uu[j]/m,vv[j]/m,bb[j]/m,(Up[j]/(m)));
    j++;
    yp=yp+400./67.;
  }
  fflush(fp);
}
/**
## Timestep Adaptation
Each timestep the timestep is adapted based on a *Vertrouwen-komt-te-voet-en-gaat-te-paard* strategy. The simulation is stopped when $t = 9 \mathrm{hours}$. Obviously the grid adaptation function has been removed.    
*/
event adapt(i++; t<=9*3600){
  if (nn>14)//Quickly reduce the timestep if things get rough
    DT=max(DT/(1+((double)nn/10.)),2.);
  if (nn<8)//Slowly increase the timestep when time integration is easy.
    DT=min(DT*(1+((double)nn/100.)),15.);
}

/**
## Results

We check if the code above has produced sensible results. Therefore we plot the wind-speed magnitude for everty ten minutes and each model level.

~~~gnuplot
set xr [-1:10]
set xlabel 'Wind-Speed magnitude [m/s]'
set ylabel 'Height [m]'
set key top left box 3
plot 'prfileGABLS10mMG.dat' u 5:1 t 'MG: U'
~~~

That looks familiar...
*/

