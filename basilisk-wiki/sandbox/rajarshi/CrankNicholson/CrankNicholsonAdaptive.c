/**
# Crank-Nicholson diffusion solver on an Adaptive Grid
*/

#include "Poisson-helmholtz4th.h"
#include "CNDLap.h"
#include "utils.h"
#include "run.h"

#define BGHOSTS 2
#define Level 7
#define pi 3.14159265359

scalar T[],Res[];
(const) face vector Diffusion[];
(const) scalar Lambda[];
struct Poisson P;
double diff = 1;
FILE *fp;

int main(){
   L0=1;
   origin(-0.5,-0.5);
   init_grid(1<<Level);
#define EPS 0.01
   
   periodic(left);
   periodic(top);
   run();
}

double GaussianFunction(double x, double x0, double y, double y0, double width, double amplitude){
   return(amplitude*exp(-((x-x0)*(x-x0)/(2*width*width) + (y-y0)*(y-y0)/(2*width*width))));
}

event init (t=0){
   foreach() 
      T[] = GaussianFunction(x,0.25,y,0.25,0.005,100) + GaussianFunction(x,-0.25,y,-0.25,0.005,100);
  
   boundary({T});
   foreach_face()
     Diffusion.x[] = diff;
}

event printdata (t += 0.001; t<0.4) {
    fp = fopen("Temperature-Adaptive.dat","w");
    foreach()
      fprintf(fp, "%g %g %g \n",x,y,T[]);
    fprintf(stdout,"\n\n");
    fclose(fp);
}

event images(t += 1./5000.){

    scalar Lev[];
    foreach()
       Lev[] = level;
    static FILE *fp1 = fopen("Grid-Adaptive.ppm","w");
    output_ppm(Lev,fp1,min=0,max=Level);
    static FILE *fp2 = fopen("Temperature-Adaptive.ppm","w");
    output_ppm(T,fp2,min=0,max=0.1);
}

event end(t=0.01){
    printf("RUN FINISHED (i=%d) & (t=%g) \n",i,t);
}

event adapt(i++) {
    adapt_wavelet ({T}, (double[]){1e-3}, maxlevel = Level);
}

event integration (i++){
    
    dt = 0.00002;
    dt = dtnext(dt);
    printf("\n Time is = %g ",t);
    scalar RHS[];
    DLap({T},{RHS},diff,dt);

    foreach()
      Lambda[] = -(2./(diff*dt));
    
    P.a = T;
    P.b = RHS;
    P.alpha = Diffusion;
    P.lambda = Lambda;
    
    mgstats Sol = poisson(T,RHS,Diffusion,Lambda,1E-03);
    residual ({T},{RHS},{Res},&P);
    output_ppm(T,min=0,max=100,file="Temp.mp4");
}

/**
We produce animations of the vorticity, pressure and kinetic-energy fields 

![Animation of the vorticity field.](TaylorGreen_centered-Euler/Temp.mp4)

*/