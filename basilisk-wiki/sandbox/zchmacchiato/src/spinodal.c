#include "vdw.h"
#include "view.h"
#include "vtk.h"
#include "utils.h"
#include "diverging_cool_warm.h"


#define MU 1e-3
#define LEVEL 8
#define LAMBDA 1e-4
#define rho_c 1.06344
#define R_g 1 
#define theta 0.9
#define p_c 1


double P0(double x)
{
  double rhop;
  rhop=x/rho_c;
  return p_c*rhop*theta*(8/(3-rhop) - 3*rhop/theta);
}

double mynoise(){
  return (1-2.*exp(noise())/exp(1));
}
  
int main()
{
  origin (-0.5, -0.5);
  periodic(right);
  periodic(top);
  init_grid (1 << LEVEL);
  DT = 1e-4;
  run();
}

event init (i = 0)
{

  lambda=LAMBDA;

  foreach()
    {
      foreach_dimension()
      mu[]=MU;
      rho[] = rho_c *( 1.+0.2*noise());
      q.x[] = q.y[] = 0.;
    }
  boundary({mu,rho,q});
}


 event dump_out(i+=1000){
    char ti[100];
    sprintf(ti,"%d",i/1000);  
    dump(ti);
 }
event end (i=10000);
event vtk_out(i+=1000){
   char ti[100];
   sprintf(ti,"%d",i);
   char name[1000]="figure";
   char tale[10]=".vtk";
   strcat(name,ti);
   strcat(name,tale);

   FILE *fp=fopen(name,"w");
   bool linear=true;
   int n=N;
   output_vtk(all,n,fp,linear);
   fclose(fp);
}

event outputfile (t=1.e-3;t+=1e-3) {
  isoline("rho", rho_c);
  squares("rho",  map = mycoolwarm);
  save("rho.mp4");

}


event adapt(i++){
  adapt_wavelet ({rho,q.x,q.y},
    (double[]){5.e-3,3.e-3,3.e-3},LEVEL, 5);
}
