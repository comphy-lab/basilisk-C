#include "grid/multigrid.h"
#include "vdw.h"
#include "view.h"
#include "diverging_cool_warm.h"
#include "utils.h"
#include "vtk.h"

#define MU 1e-2
#define LEVEL 8
#define LAMBDA 1.e-4
#define rho_c 1 
#define theta 0.95
#define p_c 1
#define D_delta 0.02
#define La 0.01
#define Ca 10
double THETAE =  90.*M_PI/180.;


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
  DT=0.5e-3;
  lambda=LAMBDA;
  miu=lambda * rho_c *L0/ La;
  periodic(right);
  run();
}

event init (t = 0)
{

    foreach()
      {
        double dist;
        if (x>=0){
        dist = x-0.25*L0;
        rho[] = rho_c *(1.0 - 0.4*tanh(2.0*dist/D_delta ) );}
        else
        {dist=x+0.25*L0;
        rho[] = rho_c *(1.0 + 0.4*tanh(2.0*dist/D_delta ) );}

        mu[] = lambda * rho_c *L0/ La;
        q.x[] = q.y[] = 0.;
      }
    boundary(all);
    printf("the velocity = %g\n",Ca*lambda/MU);
    rho[bottom]  = neumann(myContactAngle(point,rho,normGradRho,THETAE));       
    rhop[bottom] = neumann(myContactAngle(point,rhop,normGradRho,THETAE));
    rho[top]  = neumann(myContactAngle(point,rho,normGradRho,THETAE));       
    rhop[top] = neumann(myContactAngle(point,rhop,normGradRho,THETAE));
    u.n[bottom] = dirichlet(0);
    u.t[bottom] = dirichlet(-Ca*lambda/MU);
    u.n[top] = dirichlet(0);
    u.t[top] = dirichlet(Ca*lambda/MU);


}

event dump_out(i+=10000){
   char ti[100];
   sprintf(ti,"%d",i/10000);  
   dump(ti);
}



event movie(t+=1.e-2,last){
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


event vtk_out(i+=10000){
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


event end (t =10 ){}

