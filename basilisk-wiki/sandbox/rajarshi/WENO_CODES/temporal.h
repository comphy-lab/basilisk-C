/**
#Header File : Time-step advancing (Runge Kutta Schemes 3 and 4) !
*/

#include "flux_construction.h"

extern void (* Flux_Computation)  (scalar *, vector *, vector *);

double RungeKutta4(scalar *TL, vector *velL, double Tin, double Tfin, double dt){

 scalar T = TL[0];
 scalar * TList = list_clone ({T,T,T});
 scalar T1 = TList[0], T2 = TList[1], T3 = TList[2], k1[], k2[], k3[], k4[];
 vector vel = velL[0];
 face vector Flux_value[];
 double t = Tin;

 while(t < Tfin){

   t += dt;
   Flux_Computation({T},{vel},{Flux_value});
   foreach(){
     k1[] = 0.;
     foreach_dimension()
        k1[] += (Flux_value.x[]-Flux_value.x[1])/Delta;
     T1[] = T[] + k1[]*dt/2.;  
   }
   boundary({T1});
   
   Flux_Computation({T1},{vel},{Flux_value});
   foreach(){
     k2[] = 0.;
     foreach_dimension()  
        k2[] += (Flux_value.x[]-Flux_value.x[1])/Delta;
     T2[] = T[] + k2[]*dt/2.;     
   }
   boundary({T2});
   
   Flux_Computation({T2},{vel},{Flux_value});
   foreach(){
     k3[] = 0.;
     foreach_dimension()
        k3[] += (Flux_value.x[]-Flux_value.x[1])/Delta;
     T3[] = T[] + k3[]*dt;     
   }
   boundary({T3});
  
   Flux_Computation({T3},{vel},{Flux_value});
   foreach(){
     k4[] = 0.;
     foreach_dimension()
        k4[] += (Flux_value.x[]-Flux_value.x[1])/Delta;
     T[] = T[] + (k1[] + 2.*k2[] + 2.*k3[] + k4[])*dt/6.; 
   }
   boundary({T});
   
 }

 return(t);
}


double RungeKutta3(scalar *TL, vector *velL, double Tin, double Tfin, double dt){
 
 scalar T = TL[0];
 face vector vel = velL[0];
 face vector Flux_value[];
 scalar T1[], T2[];
 double t = Tin; 

 while(t <= Tfin){
   printf("\n\n %g",t);
   t +=dt;
   Flux_Computation({T},{vel},{Flux_value});
   foreach()
     T1[] = T[] - (Flux_value.x[1]-Flux_value.x[])*dt/(3.*Delta);
   T1[left] = periodic();
   boundary({T1});
   
   Flux_Computation({T1},{vel},{Flux_value});
   foreach()  
     T2[] = T[] - (Flux_value.x[1]-Flux_value.x[])*dt/(2.*Delta);
   T2[left] = periodic();
   boundary({T2});
   
   Flux_Computation({T2},{vel},{Flux_value});
   foreach()
     T[] = T[] - (Flux_value.x[1]-Flux_value.x[])*dt/Delta;
   T[left] = periodic();
   boundary({T});
 }
 return(t);
}

double Temporal1(scalar *TL, vector *velL, double Tin, double Tfin, double dt){
 
 scalar T = TL[0];
 face vector vel = velL[0];
 face vector Flux_value[];
 double t = Tin; 

 while(t <= Tfin){
   t+=dt;
   Flux_Computation({T},{vel},{Flux_value});
   foreach()
     T[] = T[] - (Flux_value.x[1]-Flux_value.x[])*dt/(Delta);
   boundary({T});
  }
 return(t);
}
