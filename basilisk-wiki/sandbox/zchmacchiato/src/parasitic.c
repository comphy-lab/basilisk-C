#include "grid/multigrid.h"
#include "vdw.h"
#include "view.h"
#include "weno5.h"
#include "utils.h"
#include "vtk.h"
#include "diverging_cool_warm.h"

#define MU 1e-3
#define LEVEL 8
#define LAMBDA 1e-4
#define rho_c 1.
#define theta 0.95
#define p_c 1.
#define D_delta 0.03
#define La 10

double THETAE =  60.*M_PI/180.;

double P0(double x)
{
    double rhop;
    rhop=x/rho_c;

    return p_c*rhop*theta*((8/(3-rhop) - 3*rhop/theta));

}

double myContactAngle(Point point, scalar rho, scalar normGradRho, 
    double THETAE){
    return normGradRho[]/tan(THETAE);
}

int main()
{
  L0 = 2;
  origin (-L0/2.);
  init_grid (1 << LEVEL);
  DT=5e-4;
  lambda=LAMBDA;
  miu= lambda * rho_c *0.15*L0/ La;

  periodic(right);
  periodic(top);
  run();
}

event init (t = 0)
{
    
    foreach()
    {
      double dist = sqrt(sq(x)+sq(y-L0/2))-0.15*L0;
      
      
      rho[] = rho_c *(1. - 0.4*tanh(2.0*dist/D_delta ) );
      mu[] = lambda * rho_c *0.15*L0/ La;
      q.x[] = q.y[] = 0.0;
    }
    boundary(all);


}
 event dump_out(i+=10000){
    char ti[100];
    sprintf(ti,"%d",i/10000);  
    dump(ti);
 }


event movie(t+=6.e-2,last){
  view(ty = -0.5);
  isoline("rho", rho_c);

  stats s = statsf(rho);

  s = statsf(u.x);
  squares("u.x", min = 0.5*s.min, max = 0.5*s.max,  map = mycoolwarm);
  save("ux.mp4");
  s = statsf(rho);
  isoline("rho", rho_c);
  squares("rho", min = 0.6, max = 1.4,  map = mycoolwarm);
  save("rho.mp4");
}


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
event logfile(i+=20){
 scalar e[],ev[],ePF[];
 vector gradrho[];
 foreach()
   foreach_dimension()
     gradrho.x[] = (WENO3_x (point, rho, -1) - WENO3_x (point, rho, 1));
 boundary ((scalar* ){gradrho});

 foreach(){
   e[] = 0.0;//theta*rho[]*log(1/3*rho[]/(1-1/3*rho[])) - 9/8*sq(rho[]);
   ev[] = 0.;
   ePF[] = 0.;
   foreach_dimension(){
     ev[] += 0.5*rho[]*sq(u.x[]);
     ePF[] += 0.5*LAMBDA*sq(gradrho.x[]);
   }
 }
   
 stats s = statsf(e), s1=statsf(ePF), s2 = statsf(ev);
 fprintf(stderr, "%g %g %g %g\n", t, s.sum,s1.sum,s2.sum);
}

event end (t= 50){}

